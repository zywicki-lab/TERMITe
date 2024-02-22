import numpy as np
import pyBigWig
import logging
from scipy.signal import find_peaks
from math import isnan
from subprocess import Popen, PIPE
from os import path
import pybedtools
import pandas as pd
import itertools
from typing import List, Set, Dict, Tuple, Optional, Union, Sequence
from numbers import Number


logging.basicConfig(level = logging.INFO)

class IDR_Exception(Exception):
    pass


class Peak:
	def __init__(self, sample_name : str, chrom : str, prominence : Number,
				left_base : int, right_base : int, height : Number, peak_coverage : np.ndarray, strand : str):
		"""Class representing a single peak as returned by the scipy.signal.find_peaks

		Args:
			sample_name (str): name of the sample
			chrom (str): chromosome
			prominence (float): prominence of the peak (as in scipy.signal.find_peaks)
			left_base (int): the leftmost position of the peak
			right_base (int): the rightmost position of the peak
			height (float): height of the peak (as in scipy.signal.find_peaks)
			peak_coverage (np.ndarray): coverage along the peak
			strand (str): strand
		"""

		self.sample_name = sample_name
		self.chrom = chrom
		self.prominence = prominence
		self.left_base = left_base + 1
		self.right_base = right_base
		self.height = height
		self.peak_coverage = peak_coverage
		self.strand = strand

		self.trimPeak()  # trimms peak ends with insignificant coverage

	def __str__(self):
		return f"{self.chrom}:{self.left_base}-{self.right_base}:prominence {self.prominence}:height {self.height}:values {self.peak_coverage}"

	def __repr__(self):
		return f"{self.chrom}:{self.left_base}-{self.right_base}:prominence {self.prominence}:height {self.height}:values {self.peak_coverage}"

	def trimPeak(self, min_percentage : int = 10):
		"""Trim peak ends with less than [min_percentage] percent of the total peak coverage from the left and right sides.
        Stops on position with >= [min_percentage] percent of the total peak coverage. The minimum length of the peak returned
        is 1 (even if the last position has also < [min_percentage] percent of the full peak coverage).

		Args:
			min_percentage (int, optional): trim peak positions below this value from the left and right sides.
   											Defaults to 10.
		"""

		peak_width = abs(self.right_base - self.left_base)
		if peak_width > 1:

			sum_of_coverage = sum(self.peak_coverage)

			# left trim
			cov = self.peak_coverage
			for i in range(peak_width):
				i_val = self.peak_coverage[i]
				perc = (i_val/sum_of_coverage) * 100
				if perc < min_percentage and len(cov) > 1:
					self.left_base += 1
					cov = np.delete(cov, 0)
				else:
					break
			self.peak_coverage = cov

			peak_width = abs(self.right_base - self.left_base)
			# right trim
			cov = self.peak_coverage
			for i in range(peak_width-1, 0, -1):
				i_val = self.peak_coverage[i]
				perc = (i_val/sum_of_coverage) * 100
				if perc < min_percentage and len(cov) > 1:
					self.right_base -= 1
					cov = np.delete(cov, -1)
				else:
					break
			self.peak_coverage = cov


class Coverage:
	def __init__(self, sample_name : str, chrom : str, chrom_size : int, termseq_coverage : np.ndarray, strand : str):
		"""Store the coverage by the 3' RNA ends

		Args:
			sample_name (str): name of the sample
			chrom (str): chromosome
			chrom_size (int): chromosome length
			termseq_coverage (np.ndarray): coverage by the 3' RNA ends
			strand (str): strand
		"""

		self.sample_name = sample_name
		self.chrom = chrom
		self.chrom_size = chrom_size
		self.termseq_coverage = termseq_coverage
		self.strand = strand


class PeakCaller:
	def __init__(self, bigwig_termseq : str, sample_name : str, strand : str):
		"""Call peaks in the coverage by the 3' RNA ends

		Args:
			bigwig_termseq (str): bigwig file name (coverage by the 3' RNA ends)
			sample_name (str): name of the sample
			strand (str): strand
		"""

		self.bigwig_termseq = bigwig_termseq
		self.sample_name = sample_name
		self.peaks = {}
		self.data = {}
		self.peaks_by_chrom = {}
		self.chrom_sizes = None
		self.strand = strand

		self.readBigWig()  # load data from the bigwig file

	def readBigWig(self):
		"""Load data from the bigwig file
  		"""

		try:
			bw_termseq = pyBigWig.open(self.bigwig_termseq)
		except Exception as e:
			raise FileNotFoundError(f'Failed to load bigwig file: {self.bigwig_termseq}\n{e}')
		
		self.chrom_sizes = bw_termseq.chroms().items()
		for chrom, size in self.chrom_sizes:
			logging.info(f'Reading data for {self.sample_name} and chromosome: {chrom}')
			termseq_coverage = bw_termseq.values(chrom, 0, size, numpy=True)
			self.data[chrom] = Coverage(sample_name=self.sample_name, chrom=chrom, chrom_size=size,
										termseq_coverage=termseq_coverage,
										strand=self.strand)

	def savePeaksToBed(self, file : str):
		"""Save peaks in the standard BED6 file format

		Args:
			file (str): output file name
		"""

		with open(file, "w") as f:
			for chrom, peaks in self.peaks_by_chrom.items():
				for peak in peaks:
					print(f"{chrom}\t{peak.left_base}\t{peak.right_base}\t.\t{peak.prominence}\t{self.strand}", file=f)

	def savePeaksToNarrowPeaks(self, file : str):
		"""Save peaks in the BED narrowPeak file format (ENCODE)

		Args:
			file (str): output file name
		"""

		with open(file, "w") as f:
			for chrom, peaks in self.peaks_by_chrom.items():
				for peak in peaks:
					print(f"{chrom}\t{peak.left_base}\t{peak.right_base}\t.\t0\t{self.strand}\t{peak.prominence}\t-999999\t-1\t-1", file=f)


	def callPeaks(self, min_dist : int, height : Union[Number, np.ndarray, None, Sequence] = None,
               	  threshold : Union[Number, np.ndarray, None, Sequence] = None,
                  distance : Number = None, prominence : Union[Number, np.ndarray, None, Sequence] = None,
                  width : Union[Number, np.ndarray, None, Sequence] = None,
				  wlen : int = None, rel_height : float = 0.5,
      			  plateau_size : Union[Number, np.ndarray, None, Sequence] = None):
		"""Call peaks in the coverage by the 3' RNA ends using scipy.signal.find_peaks

		Args:
			min_dist (int): minimal distance to other peaks. All smaller peaks this far from the higher one will be discarded
			height (Union[Number, np.ndarray, None, Sequence], optional): as in scipy.signal.find_peaks. Defaults to None.
			threshold (Union[Number, np.ndarray, None, Sequence], optional): as in scipy.signal.find_peaks. Defaults to None.
			distance (Number, optional): as in scipy.signal.find_peaks. Defaults to None.
			prominence (Union[Number, np.ndarray, None, Sequence], optional): as in scipy.signal.find_peaks. Defaults to None.
			width (Union[Number, np.ndarray, None, Sequence], optional): as in scipy.signal.find_peaks. Defaults to None.
			wlen (int, optional): as in scipy.signal.find_peaks. Defaults to None.
			rel_height (float, optional): as in scipy.signal.find_peaks. Defaults to 0.5.
			plateau_size (Union[Number, np.ndarray, None, Sequence], optional): as in scipy.signal.find_peaks. Defaults to None.
		"""

		for chrom, coverage in self.data.items():
			logging.info(f'Calling peaks for {self.sample_name} and chromosome: {chrom}')
			peaks, meta = find_peaks(coverage.termseq_coverage, height=height, threshold=threshold,
											distance=distance, prominence=prominence, width=width,
											wlen=wlen, rel_height=rel_height, plateau_size=plateau_size)
			self.peaks[chrom] = peaks, meta

			chrom_peaks = []
			for p in range(len(peaks)):
				prominence = meta["prominences"][p]
				left_base = meta["left_bases"][p]
				right_base = meta["right_bases"][p]
				height = meta["peak_heights"][p]
				peak_coverage = coverage.termseq_coverage[left_base+1:right_base]	
				region_left_base = left_base+1-min_dist if left_base+1-min_dist > 0 else 0
				region_right_base = right_base+min_dist if right_base+min_dist < dict(list(self.chrom_sizes))[chrom] else dict(list(self.chrom_sizes))[chrom]
				max_height_region = max(coverage.termseq_coverage[region_left_base:region_right_base])
				
				if max_height_region > height:
					continue
				peak = Peak(sample_name = self.sample_name, chrom=chrom, prominence=prominence,
							left_base=left_base, right_base=right_base, height=height,
							peak_coverage=peak_coverage, strand=self.strand)
				chrom_peaks.append(peak)
			self.peaks_by_chrom[chrom] = tuple(chrom_peaks)

class IDR:
	def __init__(self, samples : Tuple[str], files : List[str], output_dir : str, idr_threshold : float, strand : str):
		"""Handle the IDR routines

		Args:
			samples (Tuple[str]): names of the samples (biological replicates)
			files (List[str]): file names
			output_dir (str): output directory
			idr_threshold (float): threshold for the IDR value
			strand (str): strand
		"""

		self.samples = samples
		self.files = files
		self.idr_threshold = idr_threshold
		self.strand = strand
		self.output = path.join(output_dir, f"idrs_{samples[0]}_vs_{samples[1]}_{strand}.narrowPeak")
		self.output_bedgraph = path.join(output_dir, f"idrs_{samples[0]}_vs_{samples[1]}.{strand}.bedgraph")
  

	def runIDR(self):
		"""Run the IDR routine

		Raises:
			IDR_Exception: IDR seems not to be available in your $PATH. Aborting.
			IDR_Exception: IDR error.
		"""

		logging.info(f'Running IDR for {" and ".join(self.samples)}')
		try:
			p = Popen((["idr", "--samples"] + [self.files[0], self.files[1]] + ["--output-file", self.output,
						"--plot", "--peak-merge-method", "avg", "--soft-idr-threshold", f"{self.idr_threshold}"]),
				stdout=PIPE, stderr=PIPE)
			stdout, stderr = p.communicate()
			logging.info(stdout.decode(encoding='UTF-8',errors='strict'))
			logging.info(stderr.decode(encoding='UTF-8',errors='strict'))
		except FileNotFoundError as e:
			raise IDR_Exception(f"IDR seems not to be available in your $PATH. Aborting.")
		except OSError:
			raise IDR_Exception(f"IDR exception: {e}.")
			

	@staticmethod
	def combineIDRs(idrs : List, output_file : str, chrom_sizes : List[Tuple]) -> pd.DataFrame:
		"""Combine IDR values from multiple pairwise comparisons between samples
        IDRs are organized per genomic position

		Args:
			idrs (List[IDR]): list of IDR objects
			output_file (str): output file name
			chrom_sizes (List[Tuple]): list of 2-element tuples with chromosome name and its length (chrom, len) or dict_items

		Raises:
			FileNotFoundError: Output of the runIDR method not found

		Returns:
			pd.DataFrame: combined IDRs from pairwise samples comparisons organized by genomic positions
		"""

		logging.info('Combining IDR results from pairwise comparisons')
		for idr in idrs:
			logging.info(f'Creating bedgraph with IDR for {" and ".join(idr.samples)}')
			try:
				f = open(idr.output)
				lines = f.readlines()
				f.close()
			except FileNotFoundError as e:
				raise FileNotFoundError(f"The output of the runIDR method not found ({idr.output})")
			
			bedgraph_content = ""
			for line in lines:
				sl = line.strip().split("\t")
				local_idr = float(sl[10])
				## add local IDR to the output file
				## value in the 11th column (local IDR) is the -log10(Local IDR value)
				back_transformed_local_idr = 1/(10**local_idr)
				bedgraph_content += f"{sl[0]}\t{sl[1]}\t{sl[2]}\t{back_transformed_local_idr}\n"
			
			# add missing genomic positions to the bedgraph and fill them with np.nan
			bedgraph = pybedtools.BedTool(bedgraph_content, from_string=True).sort()	
			genome = {x : (0, y) for x, y in chrom_sizes}
			complement_content = ""
			complement = bedgraph.complement(g=genome)
			for chrom, start, end in complement:
				complement_content += f"{chrom}\t{start}\t{end}\t.\n"
			complement_content += bedgraph_content
			bedgraph = pybedtools.BedTool(complement_content, from_string=True).sort()
			bedgraph.saveas(idr.output_bedgraph)

		bedgraphs = [x.output_bedgraph for x in idrs]
		names = [f"{x.samples[0]}_vs_{x.samples[1]}" for x in idrs]

		if len(bedgraphs) > 1:  # if there are more than 2 replicates (1 pairwise comparison)
			logging.info(f'Unioning IDR bedgraphs')
			bedtool = pybedtools.BedTool()
			unionbedg = bedtool.unionbedg(i=bedgraphs, filler=".", names=names)
			unionbedg_df = unionbedg.to_dataframe(names=(["chromosome", "start", "end"]+names))
			unionbedg_df = unionbedg_df.replace(".", np.nan)
			pybedtools.cleanup(remove_all=True)
			unionbedg_df.to_csv(output_file, sep="\t", index=False)
			return unionbedg_df
		else:
			bedgraph = bedgraph.to_dataframe(names=(["chromosome", "start", "end"]+names))
			bedgraph = bedgraph.replace(".", np.nan)
			bedgraph.to_csv(output_file, sep="\t", index=False)
			pybedtools.cleanup(remove_all=True)
			return bedgraph
		
		
class Result:
	def __init__(self, unionbedg_df : pd.DataFrame, idr_threshold : float, min_no_comp : int, strand : str,
              sample_names : List[str], bigwigs_termseq : List[str], bigwigs_rnaseq : List[str], prefix : str):
		"""Generate the results

		Args:
			unionbedg_df (pd.DataFrame): combined IDRs from pairwise samples comparisons organized by genomic positions - the result of IDR.combineIDRs
			idr_threshold (float): threshold for the IDR value
			min_no_comp (int): minimum required number of pairwise samples comparisons in which IDR must be lower than the specified threshold
			strand (str): strand
			sample_names (List[str]): names of the samples (biological replicates)
			bigwigs_termseq (List[str]): list of the bigwig input files with the coverage by 3' RNA ends
			bigwigs_rnaseq (List[str]): list of the bigwig input files with the RNA-seq read coverage
			prefix (str): prefix for the IDs of the output 3' RNA terminis
		"""

		self.unionbedg_df = unionbedg_df
		self.idr_threshold = idr_threshold
		self.min_no_comp = min_no_comp
		self.idr_filtered_peaks = None
		self.strand = strand
		self.sample_names = sample_names
		self.bigwigs_termseq = bigwigs_termseq
		self.bigwigs_rnaseq = bigwigs_rnaseq
		self.prefix = prefix

		self.selectPeaksByIDR()  # selects peaks with IDR below the specified threshold in at least [--min_no_comp] pairwise comparisons
		self.calculateAverageHeight(self.prefix)  # calculates average peak height between the replicates and finds peak summits
  		# (as a peak summit, we consider positions with the highest average coverage within the peak, if more summits are present
    	# - the one most downstream is selected)

	def selectPeaksByIDR(self):
		"""Selects peaks with IDR below the specified threshold in at least [min_no_comp] pairwise comparisons.
		Results are written to the self.idr_filtered_peaks field with columns ["Count", "Max_IDR", "Min_IDR"].
		Count - the number of pairwise comparisons with IDR below the specified threshold
		Max_IDR - the highest IDR value for all pairwise comparisons
		Min_IDR - the lowest IDR value for all pairwise comparisons
		"""
  
		logging.info(f'Selecting peaks with IDR < {self.idr_threshold}')
		df = self.unionbedg_df
		columns = [x for x in df.columns if x not in ["chromosome", "start", "end"]]
		df = df.replace(np.nan, -99999999)
		df[columns] = df[columns].apply(pd.to_numeric, errors='coerce', axis=1)

		#calculate the number of comparisons with IDR < threshold 
		df["Count"] = df.apply(lambda x: sum( (x[columns] < self.idr_threshold) & (x[columns] >= 0)  ), 1)
		df = df.loc[df['Count'] >= self.min_no_comp]
		df = df.replace(-99999999, np.nan)

		# calculate min and max idr across pairwise comparisons
		df["Max_IDR"] = df[columns].apply(np.nanmax, axis=1)
		df["Min_IDR"] = df[columns].apply(np.nanmin, axis=1)
		self.idr_filtered_peaks = df
		df = self.mergeAdjacentPeaks()
		self.idr_filtered_peaks = df

	def mergeAdjacentPeaks(self) -> pd.DataFrame:
		"""Merge adjacent peaks.

		Rationale: peaks might have different lengths between samples which may lead to the situation that
  		part of the peak is covered by a different number of IDRs < threshold than the rest of the peak.
		In such a case, it will be represented as two adjacent peaks in the IDR.combineIDRs results while it should be a single peak.

		Merged peaks will have:
		minimal value from the ["Max_IDR"] column assigned as a ["Max_IDR"]
		minimal value from the ["Min_IDR"] column assigned as a ["Min_IDR"]
		maximal value from the ["Count"] column assigned as a ["Count"]

		Returns:
			pd.DataFrame: combined IDRs from pairwise samples comparisons organized by genomic positions. Adjacent peaks are merged.
		"""
		
		logging.info(f'Merging adjacent peaks')
		
		rows = [x[1] for x in self.idr_filtered_peaks.iterrows()]
		prv_end = rows[0]["end"]
		prv_start = rows[0]["start"]
		prv_chrom = rows[0]["chromosome"]
		prv_max_idr = rows[0]["Max_IDR"]
		prv_min_idr = rows[0]["Min_IDR"]
		prv_count = rows[0]["Count"]
		d = []
		added = False
		for i in range(1, len(rows)):
			if rows[i]["chromosome"] == prv_chrom and rows[i]["start"] == prv_end:
				prv_end = rows[i]["end"]
				prv_max_idr = min(rows[i]["Max_IDR"], prv_max_idr)
				prv_min_idr = min(rows[i]["Min_IDR"], prv_min_idr)
				prv_count = max(rows[i]["Count"], prv_count)
				added = False
			else:
				modified_row = {"chromosome" : prv_chrom, "start" : prv_start, "end" : prv_end,
								"Max_IDR" : prv_max_idr, "Min_IDR" : prv_min_idr, "Count" : prv_count}
				d.append(modified_row)
				prv_end = rows[i]["end"]
				prv_start = rows[i]["start"]
				prv_chrom = rows[i]["chromosome"]
				prv_max_idr = rows[i]["Max_IDR"]
				prv_min_idr = rows[i]["Min_IDR"]
				prv_count = rows[i]["Count"]
				added = True
		if not added:
			modified_row = {"chromosome" : prv_chrom, "start" : prv_start, "end" : prv_end,
						"Max_IDR" : prv_max_idr, "Min_IDR" : prv_min_idr, "Count" : prv_count}
			d.append(modified_row)
		df = pd.DataFrame(d)
		return df

	def calculateAverageHeight(self, prefix : str):
		"""Calculate average peak height between the replicates and finds peak summits
  		(as a peak summit, we consider positions with the highest average coverage within the peak, if more summits are present
    	- the one most downstream is selected)

		Args:
			prefix (str): prefix for the IDs of the output 3' RNA terminis

		Raises:
			FileNotFoundError: Failed to load the bigwig file
		"""

		logging.info(f'Calculating average peak heights, finding summits')
		d = {}
		for bigwig in self.bigwigs_termseq:
			try:
				bw_termseq = pyBigWig.open(bigwig)
			except Exception as e:
				raise FileNotFoundError(f'Failed to load the bigwig file: {bigwig}\n{e}')
			for chrom, size in bw_termseq.chroms().items():
				coverage = bw_termseq.values(chrom, 0, size, numpy=True)
				if chrom in d:
					d[chrom].append(coverage)
				else:
					d[chrom]=[coverage]
		mean_d = {}
		for chrom, covs in d.items():
			mean_d[chrom] = np.mean(covs, axis=0)
		
		heights, summits, peak_names = [], [], []
		counter = 1
		for index, row in self.idr_filtered_peaks.iterrows():
			c = mean_d[row["chromosome"]][row["start"]: row["end"]]
			c = c if self.strand == "forward" else c[::-1]  # for the reverse strand we need to calculate this 'backwards'
			max_height = max(c)
			ind_max_height = [int(x) for x in np.where(c==max_height)[0]]
			summit = ind_max_height[-1]
			heights.append(max_height)
			summits.append(summit)
			peak_names.append(f"{prefix}_3prime_signal_{self.strand}_{counter}")
			counter += 1
		self.idr_filtered_peaks["Height"] = heights
		self.idr_filtered_peaks["Summit"] = summits
		self.idr_filtered_peaks["Peak_name"] = peak_names

	def calculateTerminationEfficiency(self, window_size : int = 10, min_upstream = 0.25):
		"""Calculate termination efficiency for all termination signals.
		Extracts the average RNA-seq coverage for all replicates starting from [window_size] nt upstream up to
		[window_size] nt upstream, so the total length of the region will be 2*[window_size]+1 with the central position corresponding to the peak summit.
		Termination efficiency is defined as follows:
		T = ( (upstream - downstream)/upstream ) * 100
  
		where upstream is the average value calculated from the coverage of [window_size] nt located upstream from the peak summit and
		downstream is the average value calculated from the coverage of [window_size] nt located downstream from the peak summit

		Args:
			window_size (int, optional): Take this many nucleotides downstream and upstream from the peak summit to calculate the efficiency. Defaults to 10.
			min_upstream (float, optional): if upstream is less than the parameter value, termination efficiency will be set to "."

		Raises:
			FileNotFoundError: Failed to load bigwig file
		"""
  
		logging.info("Calculating termination efficiency")
		if not isinstance(self.bigwigs_rnaseq, list):
			logging.error(f"To calculate termination efficiency, you must use [--rna-seq-coverage] parameter")
			return

		d = {}
		for bigwig in self.bigwigs_rnaseq:
			try:
				bw_rnaseq = pyBigWig.open(bigwig)
			except Exception as e:
				raise FileNotFoundError(f'Failed to load bigwig file: {bigwig}\n{e}')
			for chrom, size in bw_rnaseq.chroms().items():
				coverage = bw_rnaseq.values(chrom, 0, size, numpy=True)
				if chrom in d:
					d[chrom].append(coverage)
				else:
					d[chrom]=[coverage]
		mean_d = {}
		for chrom, covs in d.items():
			mean_d[chrom] = np.mean(covs, axis=0)

		term_efficiencies = []
		strand = "+" if self.strand == "forward" else "reverse"
		for index, row in self.idr_filtered_peaks.iterrows():
			if self.strand == "forward":
				peak_pos = row["start"] + row["Summit"] + 1
			elif self.strand == "reverse":
				peak_pos = row["end"] - row["Summit"]
    
			offset_around_pos = 3
			coverage = mean_d[row["chromosome"]][peak_pos-window_size -1 - offset_around_pos: peak_pos+window_size + offset_around_pos]
			
			coverage = [0 if isnan(x) else x for x in coverage]
			if self.strand == "forward":
				upstream = np.mean(coverage[:window_size])
				downstream = np.mean(coverage[-window_size:])
			elif self.strand == "reverse":
				upstream = np.mean(coverage[-window_size:])
				downstream = np.mean(coverage[:window_size])
    
			if upstream >= min_upstream:
				T = ( (upstream - downstream)/upstream ) * 100
			else:
				T = "."
			term_efficiencies.append(T)
		self.idr_filtered_peaks["Termination_efficiency"] = term_efficiencies

	def extractSequence(self, output_file : str, window_size : int = 10, genome : Union[str, None] = None, min_termination_efficiency : float = 0.0):
		"""Extract the genomic sequence for the region starting [window_size] nt upstream from the peak summit up to [window_size] nt downstream,
		so that the region will have a total length of 2*[window_size]+1 with the central position corresponding to the peak summit.
		If the peak is on the reverse strand, the sequence will be reverse complemented.

		Args:
			output_file (str):output file name
			window_size (int, optional): Windows size describing how many nt will be taken around peak summit to extract the sequence. Defaults to 10.
			genome (Union[str, None], optional): path to the genome fasta file. Defaults to None.
			min_termination_efficiency (float, optional): minimal required termination efficiency to extract the sequence. Defaults to 0.
		"""

		if genome == None:
			logging.error("--genome option not specified. Cannot run extractSequence")
			return
		bed_file = ""
		strand = "+" if self.strand == "forward" else "-"
		for index, row in self.idr_filtered_peaks.iterrows():
			if self.strand == "forward":
				peak_pos = row["start"] + row["Summit"] + 1
			elif self.strand == "reverse":
				peak_pos = row["end"] - row["Summit"]
			seq_start, seq_end = peak_pos - window_size - 1, peak_pos + window_size
			if min_termination_efficiency is not None:
				if "Termination_efficiency" not in self.idr_filtered_peaks.columns:
					logging.error("To filter results based on termination efficiency (--min-term-effi), you must first run calculateTerminationEfficiency. Only supported with --rna-seq-coverage provided")
				elif row['Termination_efficiency'] > min_termination_efficiency:
					bed_file += f"{row['chromosome']}\t{seq_start}\t{seq_end}\t{row['Peak_name']}\t{row['Height']}\t{strand}\n"
			else:
				bed_file += f"{row['chromosome']}\t{seq_start}\t{seq_end}\t{row['Peak_name']}\t{row['Height']}\t{strand}\n"
		bedtool = pybedtools.BedTool(bed_file, from_string=True)
		o = bedtool.sequence(fi=genome, s=True, name=True)
		with open(output_file, "w") as f:
			print(open(o.seqfn).read(), file=f)
			
	def saveResultToBED(self, output_file : str, min_termination_efficiency : float = 0.0):
		"""Save the final results of the pipeline to the BED narrowPeak file (ENCODE)
        Fields:
        1) chromosome name (or another sequence name, e.g., contig, scaffold, etc.)
        2) peak start
        3) peak end
        4) peak name
        5) number of pairwise comparisons with IDR below the specified threshold scaled to the range of <0, 1000> to easily visualize in the IGV (color density)
			e.g. for 3 replicates - 3 pairwise comparisons will be performed (with each possible combination of samples),
   			so that the number of pairwise comparisons with IDR below the specified threshold can be either 1, 2 or 3.
      		After scaling this value to the specified range, 0 will correspond to the value of 1, 500 will correspond to the value of 2
        	and 1000 will correspond to the value of 3.
        6) strand
        7) average peak height
        8) termination efficiency (percentage) or np.nan if the termination efficiency wasn't calculated
        9) max IDR value among the IDRs calculated for each pairwise comparison of samples
        10) peak summit - point source for the peak. 0 corresponds to the most upstream position within the peak, 1 to the next, etc.
        

		Args:
			output_file (str): output file name
			min_termination_efficiency (float, optional): minimal required termination efficiency to save the peak. Defaults to 0.0.
		"""

		no_output_peaks = 0
		if not isinstance(self.idr_filtered_peaks, pd.DataFrame):
			logging.error("selectPeaksByIDR needs to be run before saving results to file")
			return
		strand = "+" if self.strand == "forward" else "-"
		with open(output_file, "w") as f:
			for index, row in self.idr_filtered_peaks.iterrows():
				try:
					if len(self.sample_names) > 2:  #if we have more than 2 replicates 
						normalized_count = int(((row["Count"]-1)/(len(list(itertools.combinations(self.sample_names, 2)))-1)) * 1000)
					else:
						normalized_count = 1000  # if there are only two replicates - there is only one peirwise comparison - nothing to normalize
					if min_termination_efficiency is not None:
						if "Termination_efficiency" not in self.idr_filtered_peaks.columns:
							logging.error("To filter results based on termination efficiency (--min-term-effi), you must first run calculateTerminationEfficiency. Only supported with --rna-seq-coverage provided")
						elif row['Termination_efficiency'] > min_termination_efficiency:
							no_output_peaks += 1
							print(f"{row['chromosome']}\t{row['start']}\t{row['end']}\t{row['Peak_name']}\t{normalized_count}\t{strand}\t{row['Height']:.1f}\t{row['Termination_efficiency']:.1f}\t{row['Max_IDR']:.3f}\t{row['Summit']}", file=f)
					else:
						no_output_peaks += 1
						if "Termination_efficiency" not in self.idr_filtered_peaks.columns:
							print(f"{row['chromosome']}\t{row['start']}\t{row['end']}\t{row['Peak_name']}\t{normalized_count}\t{strand}\t{row['Height']:.1f}\t.\t{row['Max_IDR']:.3f}\t{row['Summit']}", file=f)				
						else:
							print(f"{row['chromosome']}\t{row['start']}\t{row['end']}\t{row['Peak_name']}\t{normalized_count}\t{strand}\t{row['Height']:.1f}\t{row['Termination_efficiency']}\t{row['Max_IDR']:.3f}\t{row['Summit']}", file=f)
				except KeyError as e:
					logging.error(f"Couldn't save the results to a file. Please run result.calculateAverageHeight() before")
					return
		logging.info(f"Found {no_output_peaks} peaks")

	
			



	


		
		
		
		



