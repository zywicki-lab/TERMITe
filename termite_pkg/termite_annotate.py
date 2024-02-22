import pandas as pd
import pybedtools
from typing import List, Set, Dict, Tuple, Optional, Union, Sequence, AnyStr
import re
import argparse
import logging
import pyBigWig
import numpy as np
from math import isnan
import gffutils
from Bio import SeqIO
from Bio.Seq import Seq
from io import StringIO
import subprocess
import os
from pathlib import Path
import ViennaRNA
import string
import random
import pickle


class TransTermHP_Exception(Exception):
    pass

def assign_bigwigs_to_samples(sample_name : str, bigwig_args : List[str]) -> List[str]:
    """Select RNA-seq bigwig files belonging to a given sample.
    Returns None if no RNA-seq bigwigs are provided.

    Args:
        sample_name (str): name of the sample
        bigwig_args (List[str]): list of the bigwig files (from argparse)

    Returns:
        List[str]: bigwigs belonging to a given samples
    """

    if bigwig_args is None:  # when no RNA-seq bigwigs are provided
        return
    for sample in bigwig_args:
        tab = sample.strip().split(":")
        if tab[0] == sample_name:
            return tab[1].strip().split(",")


def downstream_gene(gffutils_database : gffutils.FeatureDB, chromosome : str,
                    strand : str, pos : int, look_on_opposite_strand : bool=False) -> List:
    """Get the nearest downstream gene 5' end coordinate in respect to the pos.

    Args:
        gffutils_database (gffutils.FeatureDB): Connection to the database
        chromosome (str): chromosome/contig ID
        strand (str): strand (+/-) - strand of the terminator
        pos (int): position from which we start looking for the nearest
                       downstream 5' gene end
        look_on_opposite_strand (bool): whether to look on the same or opposite strand

    Returns:
        List[int, str]: (position of the nearest downstream 5' gene end, gene id)
    """

    if strand == "-":
        pos = pos[0] # start
    elif strand == "+":
        pos = pos[1]  # end
    if strand == "+":
        cmd = 'SELECT min(start), id from features where seqid == ' \
            '"{}" AND strand == \"{}\" AND ' \
            'featuretype == \"gene\" AND start > {}'.format(chromosome, "+" if look_on_opposite_strand==False else "-", pos)
        downstream = gffutils_database.execute(cmd)
    else:
        cmd = 'SELECT max(end), id from features where seqid == ' \
            '"{}" AND strand == \"{}\" AND featuretype == \"gene\" ' \
            'AND end < {}'.format(chromosome, "-" if look_on_opposite_strand==False else "+", pos)
        downstream = gffutils_database.execute(cmd)
    return list(downstream.fetchone())
    
def upstream_gene(gffutils_database : gffutils.FeatureDB, chromosome : str,
                    strand : str, pos : int, look_on_opposite_strand : bool=False) -> List:
    """Get the nearest upstream gene 3' end coordinate from the pos.

    Args:
        gffutils_database (gffutils.FeatureDB): Connection to the database
        chromosome (str): chromosome/contig ID
        strand (str): strand (+/-) - strand of the terminator
        pos (int): position from which we start looking for the nearest
                       upstream 5' gene end
        look_on_opposite_strand (bool): whether to look on the same or opposite strand

    Returns:
        List[int, str]: (position of the nearest upstream 3' gene end, gene id)
    """

    if strand == "-":
        pos = pos[1] # end
    elif strand == "+":
        pos = pos[0]  # start
    if strand == "+":
        cmd = 'SELECT max(end), id from features where seqid == ' \
            '"{}" AND strand == \"{}\" AND ' \
            'featuretype == \"gene\" AND end < {}'.format(chromosome, "+" if look_on_opposite_strand==False else "-", pos)
        upstream = gffutils_database.execute(cmd)
    else:
        cmd = 'SELECT min(start), id from features where seqid == ' \
            '"{}" AND strand == \"{}\" AND featuretype == \"gene\" ' \
            'AND start > {}'.format(chromosome, "-" if look_on_opposite_strand==False else "+", pos)
        upstream = gffutils_database.execute(cmd)
    return list(upstream.fetchone())

        
class Sample:
    def __init__(self, sample_name : str, forward_file : str, reverse_file : str, forward_bigwigs : List[str], reverse_bigwigs : List[str],
                forward_rnaseq_bigwigs : List[str], reverse_rnaseq_bigwigs : List[str]) -> None:
        """Store information about the Samples. Multiple replicates might correspond to the sample
        (which should be treated as an experimental condition)

        Args:
            sample_name (str): name of the sample
            forward_file (str): TERMITe output calculated for the forward strand
            reverse_file (str): TERMITe output calculated for the reverse strand
            forward_bigwigs (List[str]): list of the bigwig files (replicates) with the coverage by the 3' RNA ends on the forward strand
            reverse_bigwigs (List[str]): list of the bigwig files (replicates) with the coverage by the 3' RNA ends on the reverse strand
            forward_rnaseq_bigwigs (List[str]): list of the bigwig files (replicates) with the coverage by the whole RNA-seq reads on the forward strand
            reverse_rnaseq_bigwigs (List[str]): list of the bigwig files (replicates) with the coverage by the whole RNA-seq reads on the reverse strand
        """

        self.sample_name = sample_name
        self.forward_file = forward_file
        self.reverse_file = reverse_file
        self.forward_bigwigs = forward_bigwigs
        self.reverse_bigwigs = reverse_bigwigs
        self.termite_results = None
        self.forward_rnaseq_bigwigs = forward_rnaseq_bigwigs
        self.reverse_rnaseq_bigwigs = reverse_rnaseq_bigwigs
        self.mean_fwd_cov = self.extractMeanCoverage(bigwigs=forward_bigwigs)
        self.mean_rev_cov = self.extractMeanCoverage(bigwigs=reverse_bigwigs)
        self.mean_rnaseq_fwd_cov = self.extractMeanCoverage(bigwigs=forward_rnaseq_bigwigs) if forward_rnaseq_bigwigs is not None else None
        self.mean_rnaseq_rev_cov = self.extractMeanCoverage(bigwigs=reverse_rnaseq_bigwigs) if reverse_rnaseq_bigwigs is not None else None
        
    def extractMeanCoverage(self, bigwigs : List[str]) -> Dict[str, np.ndarray]:
        """Calculate the mean coverage from bigwig files

        Args:
            bigwigs (List[str]): list of the bigwig files

        Raises:
            FileNotFoundError: bigwig file not found

        Returns:
            Dict[str, np.ndarray]: mean coverage (value) for each chromosome (key)
        """
        
        d = {}
        for bigwig in bigwigs:
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
        return mean_d
        
    def parseInputData(self):
        """Concatenate TERMITe results for the forward and reverse strands and calculates genomic coordinates of the peak summits.
        Results are saved in the termite_results field.
        """
        
        fwd_data = pd.read_csv(self.forward_file, sep="\t", header=None)
        rev_data = pd.read_csv(self.reverse_file, sep="\t", header=None)
        frames = [fwd_data, rev_data]
        concatenated_frames = pd.concat(frames)
        colnames = ["chromosome", "start", "end", "id", "score", "strand", "average_peak_height", "termination_efficiency", "idr", "summit"]
        concatenated_frames.columns = colnames
        
        # calculate coordinates of the peak summits
        summit_coords = []
        
        for index, row in concatenated_frames.iterrows():
            if row["strand"] == "+":
                peak_summit_coord = row["start"] + row["summit"] + 1
            elif row["strand"] == "-":
                peak_summit_coord = row["end"] - row["summit"]
            summit_coords.append(peak_summit_coord)
            
        concatenated_frames["summit_coordinate"] = summit_coords
        ####
        
        self.termite_results = concatenated_frames
        
    @staticmethod
    def mergeSamplesUsingBedtools(samples : List, output_file : str, window_size : int = 10, min_upstream = 0.25) -> str:
        """Run bedtools merge to merge TERMITe results prepared for different samples (options -c 2,3,4,5,6,7,8,9,10,11, -o collapse, -s).
        Calculate missing average peak heights, peak summits and termination efficiencies when the peak has not been called for a given sample by TERMITe.
         If a peak has not been called, the rest of the combined fields will be filled with nan so that each combined column in the result will have
        the same number of values as the number of samples (values corresponding to the samples will always be in the same order
        - alphabetical based on the sample name)

        Args:
            samples (List): list of Sample objects
            output_file (str): output file name
            window_size (int, optional): take this many nucleotides downstream and upstream from the peak summit to calculate termination efficiency.
                Use the same value as in TERMITe. Defaults to 10.
            min_upstream (float, optional): if upstream is less than the parameter value, termination efficiency will be set to np.nan

        Raises:
            argparse.ArgumentTypeError: Prefixes of the IDs must match sample names [--sample-names]

        Returns:
            str: contents of the output bed file as a string
        """

        with open(output_file, "w") as f:
            header = "chromosome\tstart\tend\tname\tscore\tstrand\tstarts\tends\tid\tscores\taverage_peak_height\ttermination_efficiency\tidr\tsummit_coordinate"
            print(header, file=f)
            bed_out = ""
            data_to_merge = [x.termite_results for x in samples]
            bedtools_input = pd.concat(data_to_merge)
            bedtool = pybedtools.BedTool.from_dataframe(bedtools_input).sort()
            bedtools_merge = bedtool.merge(c="2,3,4,5,6,7,8,9,10,11", o="collapse", s=True)
            
            ### calculate missing average peak heights, peak summits and termination efficiences
            ### if the peak wasn't called in a given sample
            colnames = ["chromosome", "start", "end", "starts", "ends", "id", "score", "strand", "average_peak_height", "termination_efficiency", "idr", "summit", "summit_coordinate"]
            bedtools_merge_df = bedtools_merge.to_dataframe(names=colnames)
            for index, row in bedtools_merge_df.iterrows():
                strand = "forward" if "+" in row["strand"] else "reverse"
                d = {}
                # extract prefixes from ids that should correspond to the sample names
                order = row["id"].strip().split(",")
                order = [re.sub(r'_3prime.+', '', x) for x in order]
                desired_order = sorted([x.sample_name for x in samples])
                
                # check if the prefixes match sample names
                for sample in order:
                    if sample not in desired_order:
                        logging.error(f'Prefixes of the IDs must match sample names [--sample-names] : {", ".join(order)} vs {", ".join(desired_order)} ')
                        raise argparse.ArgumentTypeError(f'Prefixes of the IDs must match sample names [--sample-names] : {", ".join(order)} vs {", ".join(desired_order)}')
                ####
                ### calculate the missing average peak heights, peak summits and termination efficiences
                for i in range(len(order)):
                    for col in ["starts", "ends", "id", "score", "strand", "average_peak_height", "termination_efficiency", "idr", "summit", "summit_coordinate"]:
                        values = str(row[col]).strip().split(",")
                        if not order[i] in d:
                            d[order[i]] = {f"{col}" : values[i]}
                        else:
                            d[order[i]][f"{col}"] = values[i]          
                for sample in desired_order:
                    if sample not in order:  # when the peak hasn't been called for a given sample
                        s = [x for x in samples if x.sample_name == sample][0]  # get the sample object based on sample name
                        if strand == "forward":
                            c = s.mean_fwd_cov[row["chromosome"]][row["start"]: row["end"]]
                        elif strand == "reverse":
                            c = s.mean_rev_cov[row["chromosome"]][row["start"]: row["end"]][::-1]
                        max_height = max(c)
                        ind_max_height = [int(x) for x in np.where(c==max_height)[0]]
                        summit = ind_max_height[-1]
                        peak_summit_coord = row["start"] + summit + 1 if strand == "forward" else row["end"] - summit
                        
                        # termination efficiency
                        if s.mean_rnaseq_fwd_cov is not None and s.mean_rnaseq_rev_cov is not None:
                            offset_around_pos = 3  # ignore this many nucleotides from both sides of POT to calculate the termination efficiency
                            if strand == "forward":
                                coverage = s.mean_rnaseq_fwd_cov[row["chromosome"]][peak_summit_coord-window_size -1 - offset_around_pos: peak_summit_coord+window_size + offset_around_pos]
                            elif strand == "reverse":
                                coverage = s.mean_rnaseq_rev_cov[row["chromosome"]][peak_summit_coord-window_size -1 - offset_around_pos: peak_summit_coord+window_size + offset_around_pos]
                    
                            coverage = [0 if isnan(x) else x for x in coverage]
                            if strand == "forward":
                                upstream = np.mean(coverage[:window_size])
                                downstream = np.mean(coverage[-window_size:])
                            elif strand == "reverse":
                                upstream = np.mean(coverage[-window_size:])
                                downstream = np.mean(coverage[:window_size])
                            # if upstream is equal to 0, set T to "."
                            if upstream >= min_upstream:
                                T = ( (upstream - downstream)/upstream ) * 100
                            else:
                                T = "."
                        else:
                            T = "."
                            
                        # put calculated values in the dict
                        d[sample] = {
                            "starts" : ".",
                            "ends" : ".",
                            "id" : ".",
                            "score" : ".",
                            "strand" : "+" if strand == "forward" else "-",
                            "average_peak_height" : f"{float(max_height):.1f}",
                            "termination_efficiency" : f"{T}",
                            "idr" : ".",
                            "summit" : summit,
                            "summit_coordinate" : peak_summit_coord}
                s = "+" if strand == "forward" else "-"
                output_line = f"{row['chromosome']}\t{row['start']}\t{row['end']}\t.\t.\t{s}\t"
                for column in ["starts", "ends", "id", "score", "average_peak_height", "termination_efficiency", "idr", "summit_coordinate"]:
                    for key in desired_order:
                        if column in ["average_peak_height", "termination_efficiency"]:
                            output_line += f"{float(d[key][column]):.1f}," if d[key][column] != "." else f".,"
                        elif column in ["idr"]:
                            output_line += f"{float(d[key][column]):.3f}," if d[key][column] != "." else f".,"
                        else:
                            output_line += f"{d[key][column]},"
                    output_line = f"{output_line[:-1]}\t"  # remove trailing comma and add tab
                output_line = output_line[:-1]  # remove trailing tab
                print(output_line, file=f)
                bed_out = f"{bed_out}\n{output_line}"
        
        return bed_out
        
        ###
        
class TerminatorHP:
    def __init__(self, id : str, chrom : str, start : int, end : int, strand : str, loc : str, conf : float,
                 hp_score : float, tail_score : float, hairpin : str):
        """Store information about the intrinsic terminator hairpins as returned by the TransTermHP software

        Args:
            id (str): id of the hairpin
            chrom (str): chromosome
            start (int): start coordinate (of the hairpin - tails not included)
            end (int): end coordinate (of the hairpin - tails not included)
            strand (str): strand
            loc (str): localization of the hairpin as described by the TransTermHP
            conf (float): confidence of the intrinsic terminator <0, 100>. Higher the better.
            hp_score (float): hairpin score as calculated by the TransTermHP
            tail_score (float): tail score as calculated by the TransTermHP
            hairpin (str): hairpin structure as str (includes tails)
        """

        self.term_id = id
        self.chromosome = chrom
        self.start = start
        self.end = end
        self.strand = strand 
        self.genomic_loc = loc
        #'Loc' gives the type of region the terminator is in:
        #'G' = in the interior of a gene (at least 50bp from an end),
        #'F' = between two +strand genes,
        #'R' = between two -strand genes,
        #'T' = between the ends of a +strand gene and a -strand gene,
        #'H' = between the starts of a +strand gene and a -strand gene,
        #'N' = none of the above (for the start and end of the DNA)
        #Because of how overlapping genes are handled, these designations are not
        #exclusive. 'G', 'F', or 'R' can also be given in lowercase, indicating that
        #the terminator is on the opposite strand as the region.  Unless the
        #--all-context option is given, only candidate terminators that appear to be in
        #an appropriate genome context (e.g. T, F, R) are output.
        
        self.confidence = conf  # the overall confidence score (which ranges from 0 to 100) is what you probably want to use to assess the quality of a terminator. Higher is better.
        self.hairpin_score = hp_score  # the hairpin score
        self.tail_score = tail_score  # the tail score
        self.hairpin = hairpin 
        repr = self.parse_hairpin_representation()
        self.a_tract, self.hairpin, self.u_tract = repr["A-tract"], repr["Hairpin"], repr["U-tract"]

    def parse_hairpin_representation(self) -> Dict[str, Seq]:
        """Parse TransTermHP hairpin representation and extracts hairpin, possible A-tract and U-tract sequences.
        Returned sequences are always written from 5' to 3' end and complement if on the reverse strand.

        Raises:
            ValueError: wrong hairpin representation

        Returns:
            Dict[str, Seq]: hairpin, possible A-tract and U-tract sequences
        """

        if self.strand == "-":
            self.hairpin = self.hairpin[::-1]
        representation = self.hairpin.strip().split()
        if len(representation) == 4:
            hairpin = self.hairpin[:-15]
            u = self.hairpin[-15:]
            a = "."
        elif len(representation) == 5:
            hairpin = self.hairpin[15:-15]
            u = self.hairpin[-15:]
            a = self.hairpin[:15]
        else:
            raise ValueError(f"Wrong hairpin representation - {representation}")
        if self.strand == "+":
            return {"A-tract" : Seq(a.lower()),
                    "Hairpin" : Seq(hairpin.strip().lower()),
                    "U-tract" : Seq(u.lower())}
        else:
            return {"A-tract" : Seq(a.lower()).complement(),
                    "Hairpin" : Seq(hairpin.strip().lower()).complement(),
                    "U-tract" : Seq(u.lower()).complement()}
     
class AnnotateResults:
    def __init__(self, output_prefix : str, gene_annotations : str, samples : List[Sample], genome_fasta : str, termite_results : str, args : argparse.Namespace):
        """Annotate merged peaks with additional information, such as the closest computationally predicted intrinsic termination hairpins,
        closest genes annotated both upstream and downstream, etc.

        Args:
            output_prefix (str): output prefix name
            gene_annotations (str): GTF file with gene annotations
            genome_fasta (str): FASTA file with the genomic sequence
            samples (List[Sample]): list of Sample objects
            termite_results (str): TERMITe results as returned by Sample.mergeSamplesUsingBedtools
            args (argparse.Namespace): command line arguments parsed using argparse
        """

        self.output_prefix = output_prefix
        self.genome_fasta = genome_fasta
        if args.trans_term_hp:
            logging.info(f"Running TransTermHP")
            self.run_transtermhp()
        self.gene_annotations = gene_annotations
        self.samples = samples
        self.termite_results = termite_results
        self.terminators = self.parseTransTermResults() if args.trans_term_hp == True else None
        self.pyBedTool = self.transTermResultsToBedTool() if args.trans_term_hp == True else None

        if args.trans_term_hp == True:
            self.saveTransTermResultsToTabular()
            df = self.mergeWithTermiteResults()
        else:
            df = StringIO(self.termite_results)
            colnames =  ["chromosome", "start", "end", "name", "score", "strand", "starts", "ends", "termite_id",
                        "termite_score", "average_peak_height", "termination_efficiency","idr", "summit_coordinate"]
        
            df = pd.read_csv(df, sep="\t", header=None, names=colnames)
        
        if args.gene_annotations:
            logging.info(f"Reading gene annotations")
            df = self.add_gene_annotations(df=df)
            df = self.add_overlapping_features(df=df)
        
        df = self.assignRefPOT(df, samples=samples, transtermhp=args.trans_term_hp, max_distance_to_hairpin=args.max_distance_to_hairpin,
                               min_distance_to_hairpin=args.min_distance_to_hairpin)
        if args.rna_seq_coverage_forward is not None:
            logging.info(f"Calculating termination efficiencies")
            df = self.calculate_termination_efficiency_around_POT(df=df, samples=samples, window_size=10, min_upstream=args.min_upstream)
        if args.genome is not None:
            logging.info(f"Extracting sequences around POTs")
            df = self.extractSequence(df=df, upstream=args.upstream_nt, downstream=args.downstream_nt, output_file_prefix=output_prefix)

        
        if args.bed is not None:
            df = self.intersectBEDs(df=df, beds=args.bed)
        
        self.df = df
        
                


    def run_transtermhp(self) -> None:
        """Run TransTermHP and saves results in the file"

        Raises:
            TransTermHP_Exception: TransTermHP seems not to be available in your $PATH. Aborting.
            TransTermHP_Exception: TransTermHP exception.
            TransTermHP_Exception: TransTermHP finished with the return code other than 0.
            TransTermHP_Exception: Couldn't save TransTermHP results to the temporary file.

        Returns:
            None: None
        """

        # create file with the coordinates of 'fakegenes' like suggested in the TransTermHP
        n = 1
        with open(f"{self.output_prefix}_fake_genes.coords", "w") as f:
            fasta = SeqIO.parse(open(self.genome_fasta),'fasta')
            for fasta_record in fasta:
                name, sequence = fasta_record.id, fasta_record.seq
                chrom = name.split()[0]
                print(f"fakegene{n}\t1\t2\t{chrom}", file = f)
                print(f"fakegene{n+1}\t{len(sequence)-1}\t{len(sequence)}\t{chrom}", file = f)
                n += 2
        ##

        # run transtermhp
        try:
            script_path = os.path.dirname(os.path.realpath(__file__))
            p = subprocess.Popen(["transterm", "--all-context", "-p", f"{script_path}/expterm.dat", f"{self.genome_fasta}",
                                f"{self.output_prefix}_fake_genes.coords"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            out, err = p.communicate()
            out, err = out.decode("utf-8"), err.decode("utf-8")
        except FileNotFoundError as e:
            raise TransTermHP_Exception(f"TransTermHP seems not to be available in your $PATH. Aborting.")
        except OSError as e:
            raise TransTermHP_Exception(f"TransTermHP exception: {e}.")
        
        if p.returncode != 0:
            logging.error(err)
            raise TransTermHP_Exception(f"TransTermHP finished with the return code: {p.returncode}.")
        else:
            logging.info(f"TransTermHP finished successfully and returned code: {p.returncode}")
            
        try:
            with open(f"{self.output_prefix}_transterm_predictions.tt", "w") as tt_out:
                print(out, file = tt_out)
        except Exception as e:
            raise TransTermHP_Exception(f"Couldn't save TransTermHP results to the temporary file: {e}.")

        return
        
        
    def saveTransTermResultsToTabular(self) -> str:
        """Save TransTermHP results in the tabular file format (tsv) in the file "{output_prefix}_transterm.tab"

        Returns:
            str: text representation of the output file contents
        """

        with open(f"{self.output_prefix}_transterm.tab", "w") as f:
            output = ""
            for terminator in self.terminators:
                output_line = f"{terminator.chromosome}\t{terminator.start}\t{terminator.end}\t{terminator.term_id}\t{terminator.confidence}\t{terminator.strand}\t{terminator.genomic_loc}\t{terminator.hairpin_score}\t{terminator.tail_score}\t{terminator.a_tract}\t{terminator.hairpin}\t{terminator.u_tract}"
                print(output_line, file=f)
                output = f"{output}\n{output_line}"
        return output
    
    def transTermResultsToBedTool(self) -> pybedtools.BedTool:
        """Save TransTermHP results in the tabular file format (tsv) in the file "{output_prefix}_transterm.tab" and returns pyBEDTool object

        Returns:
            pybedtools.BedTool: pyBEDTool object with the TransTermHP results
        """

        output = ""
        for terminator in self.terminators:
            output_line = f"{terminator.chromosome}\t{terminator.start}\t{terminator.end}\t{terminator.term_id}\t{terminator.confidence}\t{terminator.strand}\t{terminator.genomic_loc}\t{terminator.hairpin_score}\t{terminator.tail_score}\t{terminator.a_tract}\t{terminator.hairpin}\t{terminator.u_tract}"
            output = f"{output}\n{output_line}"
        return pybedtools.BedTool(output, from_string=True)
    
    def saveTransTermResultsToBed(self) -> str:
        """Save TransTermHP results in the BED file format in the file "{output_prefix}_transterm.tab"

        Returns:
            str: text representation of the output file contents
        """

        with open(f"{self.output_prefix}_transterm.tab", "w") as f:
            output = ""
            for terminator in self.terminators:
                output_line = f"{terminator.chromosome}\t{terminator.start}\t{terminator.end}\t{terminator.term_id}\t{terminator.confidence}\t{terminator.strand}"
                output = f"{output}\n{output_line}"
        return output
        
    def parseTransTermResults(self) -> List[TerminatorHP]:
        """Parse the output of the TransTermHP software with computational genome-wide predictions of the intrinsic terminators (*.tt).

        Raises:
            argparse.FileNotFoundError: Couldn't find the output file of the TransTermHP (*.tt)

        Returns:
            List[TerminatorHP]: list of TerminatorHP objects
        """

        loc_desc = {'G' : 'in the interior of a gene (at least 50bp from an end)',
                    'F' : 'between two +strand genes',
                    'R' : 'between two -strand genes',
                    'T' : 'between the ends of a +strand gene and a -strand gene',
                    'H' : 'between the starts of a +strand gene and a -strand gene',
                    'N' : 'not appicable (start and end of DNA)',
                    'g' : 'on the opposite strand to the interior of a gene (at least 50bp from an end)',
                    'f' : 'between two +strand genes on the opposite strand',
                    'r' : 'between two -strand genes on the opposite strand'}
        result = []
        try:
            with open(f"{self.output_prefix}_transterm_predictions.tt") as f:
                prv_line_was_term = False
                for line in f:
                    sl = line.strip()
                    if sl.startswith("SEQUENCE"):
                        chromosome = sl.split()[1]  # chromosome name is read until the first whitespace (the rest is ignored)
                    if prv_line_was_term == True:
                        termHP = sl
                        terminator = TerminatorHP(id=term_id, chrom=chromosome, start=start, end=end, strand=strand, loc=genomic_loc, conf=conf,
                                                  hp_score=hp_score, tail_score=tail_score, hairpin=termHP)
                        result.append(terminator)
                    if sl.startswith("TERM"):
                        sl = sl.split()
                        term_id = "_".join(sl[:2])
                        start = min(int(sl[2]), int(sl[4]))
                        end = max(int(sl[2]), int(sl[4]))
                        strand = sl[5]
                        genomic_loc = ""
                        for loc in sl[6]:
                            genomic_loc += f'{loc_desc[loc]},'
                        genomic_loc = genomic_loc[:-1]
                        conf = sl[7]
                        hp_score = f"{float(sl[8]):.1f}"
                        tail_score = f"{float(sl[9]):.1f}"
                        prv_line_was_term = True
                    else:
                        prv_line_was_term = False
        except FileNotFoundError:
            raise argparse.FileNotFoundError(f"No such file: {self.input_file}")
        
        return result
            
    def mergeWithTermiteResults(self) -> pd.DataFrame:
        """Annotate merged peaks with the closest TransTermHP hairpins (bedtools closest -s -D a -t all)

        Args:
            termite_results (str): TERMITe results as returned by Sample.mergeSamplesUsingBedtools

        Returns:
            pd.DataFrame: results of the bedtools closest in the form of a DataFrame
        """

        termite_bedtool = pybedtools.BedTool(self.termite_results, from_string=True).sort()

        termite_bedtool.saveas('termite_bedtool.bed')
        p = subprocess.Popen(["sort", "-V", "-k1,1", "-k2,2", f"{self.output_prefix}_transterm.tab"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = p.communicate()
        out, err = out.decode("utf-8"), err.decode("utf-8")
        with open(f"{self.output_prefix}_transterm.sorted.tab", "w") as f:
            print(out, file=f)
        result = termite_bedtool.closest(f"{self.output_prefix}_transterm.sorted.tab", s=True, D="a", t="all", id=True)
        result.saveas("closest.bed")
        colnames =  ["chromosome", "start", "end", "name", "score", "strand", "starts", "ends", "termite_id", "termite_score", "average_peak_height", "termination_efficiency",
                     "idr", "summit_coordinate", "transtermhp_hairpin_chromosome", "transtermhp_closest_hairpin_start", "transtermhp_closest_hairpin_end",
                     "transtermhp_id", "transtermhp_confidence", "transtermhp_strand", "transtermhp_location", "transtermhp_hairpin_score", "transtermhp_tail_score", "transtermhp_a_tract", "transtermhp_hairpin", "transtermhp_u_tract", "distance_to_hairpin"]
        results_df = result.to_dataframe(names=colnames)

        del results_df["transtermhp_hairpin_chromosome"]
        del results_df["transtermhp_strand"]
        del results_df["distance_to_hairpin"]
        
        return results_df
    
    def add_gene_annotations(self, df : pd.DataFrame) -> pd.DataFrame:
        """Annotate the results with information about the closest genes both upstream and downstream from the merged peak on the same and opposite strand.
        Annotation is done using the gffutils package.
        
        Args:
            df (pd.DataFrame): input dataframe to be annotated.
            
        Returns:
            pd.DataFrame: output data frame with the annotation."""

        db = gffutils.create_db(self.gene_annotations, dbfn=f'{self.gene_annotations}.db', force=True, keep_order=True, merge_strategy='merge', sort_attribute_values=True, disable_infer_genes=True, disable_infer_transcripts=True)
        distances_to_upstream_gene_same_strand, distances_to_downstream_gene_same_strand = [], []
        distances_to_upstream_gene_opposite_strand, distances_to_downstream_gene_opposite_strand = [], []
        overlapping_genes_same_strand, overlapping_genes_opposite_strand = [], []
        for index, row in df.iterrows():
            chromosome, start, end, strand = row["chromosome"], row["start"], row["end"], row["strand"]
            overlapping_features = [x.id for x in db.region((chromosome, start, end), strand=strand, featuretype="gene")]
            overlapping_gene_same_strand = "." if len(overlapping_features) == 0 else ",".join(overlapping_features)
            overlapping_features_opposite_strand = [x.id for x in db.region((chromosome, start, end), strand="+" if strand == "-" else "-", featuretype="gene")]
            overlapping_gene_opposite_strand = "." if len(overlapping_features_opposite_strand) == 0 else ",".join(overlapping_features_opposite_strand)
            upstream_gene_same_strand = upstream_gene(gffutils_database=db, chromosome=chromosome, strand=strand, pos=(start, end), look_on_opposite_strand=False)
            upstream_gene_opposite_strand = upstream_gene(gffutils_database=db, chromosome=chromosome, strand=strand, pos=(start, end), look_on_opposite_strand=True)
            downstream_gene_same_strand = downstream_gene(gffutils_database=db, chromosome=chromosome, strand=strand, pos=(start,end), look_on_opposite_strand=False)
            downstream_gene_opposite_strand = downstream_gene(gffutils_database=db, chromosome=chromosome, strand=strand, pos=(start, end), look_on_opposite_strand=True)
            
            try:
                distance_to_upstream_gene_same_strand = start - upstream_gene_same_strand[0] if strand == "+" else upstream_gene_same_strand[0] - end
            except TypeError: ## upstream_gene_same_strand not found (None) - beginning of the chromosome
                distance_to_upstream_gene_same_strand = "."
            try:
                distance_to_downstream_gene_same_strand = downstream_gene_same_strand[0] - end if strand == "+" else start - downstream_gene_same_strand[0]
            except TypeError:
                distance_to_downstream_gene_same_strand = "."
            try:
                distance_to_upstream_gene_opposite_strand = start - upstream_gene_opposite_strand[0] if strand == "+" else upstream_gene_opposite_strand[0] - end
            except TypeError:
                distance_to_upstream_gene_opposite_strand = "."
            try:
                distance_to_downstream_gene_opposite_strand = downstream_gene_opposite_strand[0] - end if strand == "+" else start - downstream_gene_opposite_strand[0]
            except TypeError:
                distance_to_downstream_gene_opposite_strand = "."
            

            distances_to_upstream_gene_same_strand.append(f"{upstream_gene_same_strand[1]} - {distance_to_upstream_gene_same_strand}bp")
            distances_to_downstream_gene_same_strand.append(f"{downstream_gene_same_strand[1]} - {distance_to_downstream_gene_same_strand}bp")
            distances_to_upstream_gene_opposite_strand.append(f"{upstream_gene_opposite_strand[1]} - {distance_to_upstream_gene_opposite_strand}bp")
            distances_to_downstream_gene_opposite_strand.append(f"{downstream_gene_opposite_strand[1]} - {distance_to_downstream_gene_opposite_strand}bp")
            overlapping_genes_same_strand.append(overlapping_gene_same_strand)
            overlapping_genes_opposite_strand.append(overlapping_gene_opposite_strand)
        
        df["overlapping_gene"] = overlapping_genes_same_strand
        #df["overlapping_gene_on_the_opposite_strand"] = overlapping_genes_opposite_strand
        df["upstream_gene"] = [x if x != "None - nanbp" else "." for x in distances_to_upstream_gene_same_strand]
        df["downstream_gene"] = [x if x != "None - nanbp" else "." for x in distances_to_downstream_gene_same_strand]
        #df["distance_to_upstream_gene_on_the_opposite_strand"] = distances_to_upstream_gene_opposite_strand
        #df["distance_to_downstream_gene_on_the_opposite_strand"] = distances_to_downstream_gene_opposite_strand
        
        return df
    
    def add_overlapping_features(self, df : pd.DataFrame) -> pd.DataFrame:
        """Annotate the results with the overlapping feature types (like CDS) from the GTF input file. Annotation is done using the gffutils package.
        
        Args:
            df (pd.DataFrame): input dataframe to be annotated.
            
        Returns:
            pd.DataFrame: output data frame with the annotation."""

        db = gffutils.create_db(self.gene_annotations, dbfn=f'{self.gene_annotations}.db', force=True, keep_order=True, merge_strategy='merge', sort_attribute_values=True, disable_infer_genes=True, disable_infer_transcripts=True)
        overlapping_features_per_3prime_end = []
        for index, row in df.iterrows():
            chromosome, start, end, strand = row["chromosome"], row["start"], row["end"], row["strand"]
            overlapping_features = ",".join(sorted([x.featuretype for x in db.region((chromosome, start, end), strand=strand)]))
            overlapping_features = "." if overlapping_features == '' else overlapping_features
            overlapping_features_per_3prime_end.append(overlapping_features)
        df["overlapping_feature_types"] = overlapping_features_per_3prime_end
        return df
    
    def assignRefPOT(self, df : pd.DataFrame, samples: List[Sample], transtermhp : bool, max_distance_to_hairpin : int = 10, min_distance_to_hairpin : int = 0) -> pd.DataFrame:
        """Assign the reference POT (Point of Observed [RNA] Termini) by an averaging signal from all samples and all replicates
        (in the genomic coordinates of the merged peak), then signal summit within this region is treated as the reference POT.
        If more than one position with maximal signal within the region has the same average signal value, the one most downstream is selected.

        Args:
            df (pd.DataFrame): output of the 
             method. DataFrame to annotate.
            samples (List): list of Sample objects
            transtermhp (bool): True if transtermhp option has been used
            max_distance_to_hairpin (int, optional): maximum distance between the POT and the intrinsic terminator hairpin located upstream. Defaults to 10.
            min_distance_to_hairpin (int, optional): minimum distance between the POT and the intrinsic terminator hairpin located upstream. Defaults to 0. 


        Raises:
            FileNotFoundError: failed to load one of the bigwig files with the coverage by the RNA 3' ends

        Returns:
            pd.DataFrame: annotated data frame (with POT column)
        """

        pots = []
        for index, row in df.iterrows():
            chromosome = row["chromosome"]
            start = row["start"]
            end = row["end"]
            strand = row["strand"]
            coverages = []
            for sample in self.samples:
                bigwigs = sample.forward_bigwigs if strand == "+" else sample.reverse_bigwigs
                for bigwig in bigwigs:
                    try:
                        bw_termseq = pyBigWig.open(bigwig)
                    except Exception as e:
                        raise FileNotFoundError(f'Failed to load the bigwig file: {bigwig}\n{e}')
                    coverage = bw_termseq.values(chromosome, start, end, numpy=True)
                    coverages.append(coverage)
            c = np.mean(coverages, axis=0)
            max_height = max(c)
            ind_max_height = [int(x) for x in np.where(c==max_height)[0]]
            summit = ind_max_height[-1]
            peak_summit_coord = start + summit + 1
            pots.append(peak_summit_coord)
        df["POT"] = pots

        if transtermhp == True:
            distances = [int(row["transtermhp_closest_hairpin_start"]) - int(row["POT"]) -1 if row["strand"] == "-" else int(row["POT"]) - int(row["transtermhp_closest_hairpin_end"]) - 1 for _, row in df.iterrows()]
            df["transtermhp_distance_to_hairpin"] = distances
            df["transtermhp"] = ["+" if row["transtermhp_distance_to_hairpin"] <= max_distance_to_hairpin and row["transtermhp_distance_to_hairpin"] >= min_distance_to_hairpin else "-" for ind, row in df.iterrows()]

            
            modified_u_tract = []
            for _, row in df.iterrows():
                try:
                    dist = int(row["transtermhp_distance_to_hairpin"])
                    if dist < 0:  # possible when the POT overlaps the hairpin
                        modified_u_tract.append(row["transtermhp_u_tract"])
                        continue
                    u = row["transtermhp_u_tract"][:dist]
                    u = f"{u}{row['transtermhp_u_tract'][dist].upper()}{row['transtermhp_u_tract'][dist+1:]}"
                    modified_u_tract.append(u)
                except IndexError:
                    modified_u_tract.append(row["transtermhp_u_tract"])
            df["transtermhp_u_tract"] = modified_u_tract

            modified_hairpins = []   # mark with the capital letter if the 3' end lies within the hairpin
            for _, row in df.iterrows():
                try:
                    dist = int(row["transtermhp_distance_to_hairpin"])
                    if dist >= 0:  # possible when the POT overlaps the hairpin
                        modified_hairpins.append(row["transtermhp_hairpin"])
                        continue
                    else:
                        dist = abs(dist)
                    modified_hairpin = ""
                    ind = 0
                    hrp = row["transtermhp_hairpin"][::-1]
                    for nucl in hrp:
                        if nucl != ' ' and nucl != "-":
                            ind += 1
                            if ind == dist:
                                modified_hairpin = f"{modified_hairpin}{nucl.upper()}"
                            else:
                                modified_hairpin = f"{modified_hairpin}{nucl}"
                        else:
                            modified_hairpin = f"{modified_hairpin}{nucl}"
                    modified_hairpins.append(modified_hairpin[::-1])
                except IndexError:
                    modified_hairpins.append(row["transtermhp_hairpin"])
            df["transtermhp_hairpin"] = modified_hairpins
            
            ### if the hairpin is to far from the POT, set all the fields to "."
            df["transtermhp_closest_hairpin_start"] = ["." if row["transtermhp"] == "-" else row["transtermhp_closest_hairpin_start"] for ind, row in df.iterrows()]
            df["transtermhp_closest_hairpin_end"] = ["." if row["transtermhp"] == "-" else row["transtermhp_closest_hairpin_end"] for ind, row in df.iterrows()]
            df["transtermhp_id"] = ["." if row["transtermhp"] == "-" else row["transtermhp_id"] for ind, row in df.iterrows()]
            df["transtermhp_confidence"] = ["." if row["transtermhp"] == "-" else row["transtermhp_confidence"] for ind, row in df.iterrows()]
            df["transtermhp_location"] = ["." if row["transtermhp"] == "-" else row["transtermhp_location"] for ind, row in df.iterrows()]
            df["transtermhp_hairpin_score"] = ["." if row["transtermhp"] == "-" else row["transtermhp_hairpin_score"] for ind, row in df.iterrows()]
            df["transtermhp_tail_score"] = ["." if row["transtermhp"] == "-" else row["transtermhp_tail_score"] for ind, row in df.iterrows()]
            df["transtermhp_hairpin"] = ["." if row["transtermhp"] == "-" else row["transtermhp_hairpin"] for ind, row in df.iterrows()]
            df["transtermhp_closest_hairpin_start"] = ["." if row["transtermhp"] == "-" else row["transtermhp_closest_hairpin_start"] for ind, row in df.iterrows()]
            df["transtermhp_u_tract"] = ["." if row["transtermhp"] == "-" else row["transtermhp_u_tract"] for ind, row in df.iterrows()]
            df["transtermhp_a_tract"] = ["." if row["transtermhp"] == "-" else row["transtermhp_a_tract"] for ind, row in df.iterrows()]
            df["transtermhp_distance_to_hairpin"] = ["." if row["transtermhp"] == "-" else row["transtermhp_distance_to_hairpin"] for ind, row in df.iterrows()]
            ##
        
        return df
    
    def calculate_termination_efficiency_around_POT(self, df: pd.DataFrame, samples: List[Sample], window_size : int = 10, min_upstream = 0.25) -> pd.DataFrame:
        """Calculate termination efficiency for each replicate using the reference Point of Observed [RNA] Termini (POT) calculated with the assignRefPOT method

        Args:
            df (pd.DataFrame): output of the assignRefPOT method. DataFrame to annotate.
            samples (List[Sample]): list of Sample objects
            window_size (int, optional): Take this many nucleotides downstream and upstream from the calculated Reference POT to calculate termination efficiency. Defaults to 10.
            min_upstream (float, optional): if upstream is less than the parameter value, termination efficiency will be set to np.nan
        
        Raises:
            FileNotFoundError: failed to load one of the bigwig files with the coverage by the RNA 3' ends

        Returns:
            pd.DataFrame: annotated data frame (with termination efficiencies for each replicate)
        """
        upstreams, downstreams = {}, {}
        result = {s.sample_name : [] for s in samples}
        for index, row in df.iterrows():
            chromosome = row["chromosome"]
            start = row["start"]
            end = row["end"]
            pot = row["POT"]
            strand = "forward" if "+" in row["strand"] else "reverse"
            
            for s in samples:
                term_efficiencies = []
                for i in range(len(s.forward_rnaseq_bigwigs)):
                    # termination efficiency
                    d = {}
                    if strand == "forward":
                        bigwig = s.forward_rnaseq_bigwigs[i]
                    elif strand == "reverse":
                        bigwig = s.reverse_rnaseq_bigwigs[i]
                    try:
                        bw_rnaseq = pyBigWig.open(bigwig)
                    except Exception as e:
                        raise FileNotFoundError(f'Failed to load the bigwig file: {bigwig}\n{e}')
                    coverage = bw_rnaseq.values(chromosome, 0, bw_rnaseq.chroms()[chromosome], numpy=True)
                    
                    offset_around_pos = 3
                    coverage = coverage[pot-window_size -1 - offset_around_pos : pot+window_size + offset_around_pos]
            
                    coverage = [0 if isnan(x) else x for x in coverage]
                    if strand == "forward":
                        upstream = np.mean(coverage[:window_size])
                        downstream = np.mean(coverage[-window_size:])
                    elif strand == "reverse":
                        upstream = np.mean(coverage[-window_size:])
                        downstream = np.mean(coverage[:window_size])
                    # if upstream is equal to 0, set T to 0
                    upstreams[f"{chromosome}_{start}_{end}_{pot}_{strand}_{s}"] = upstream
                    downstreams[f"{chromosome}_{start}_{end}_{pot}_{strand}_{s}"] = downstream
                    if upstream >= min_upstream:
                        T = ( (upstream - downstream)/upstream ) * 100
                    else:
                        T = "."
                    term_efficiencies.append( f"{T}")
                result[s.sample_name].append(",".join(term_efficiencies))
        for sample_name in result:
            df[f"termination_efficiencies_{sample_name}"] = result[sample_name]
        with open('upstreams.pickle', 'wb') as handle:
            pickle.dump(upstreams, handle, protocol=pickle.HIGHEST_PROTOCOL)
        with open('downstreams.pickle', 'wb') as handle:
            pickle.dump(downstreams, handle, protocol=pickle.HIGHEST_PROTOCOL)
        return df
            

    def extractSequence(self, df : pd.DataFrame, upstream : int, downstream : int, output_file_prefix : str) -> pd.DataFrame:
        """Extract the genomic sequence for the region starting [upstream] nt upstream from the reference POT up to [downstream] nt downstream,
		so that the region will have a total length of upstream+downstream+1 with the central position corresponding to the reference POT (Point of Observed [RNA] Termini).
		If the peak is on the reverse strand, the sequence will be reverse complemented.

		Args:
            df (pd.DataFrame): DataFrame as returned by the assignRefPOT method (with POT) column
			output_file_prefix (str): output file name prefix (.fasta will be added as an extension)
			downstream (int): how many nt will be taken downstream from the POT to extract the sequence.
            upstream (int): how many nt will be taken upstream from the POT to extract the sequence.

        Returns:
            pd.DataFrame: annotated data frame (with Sequence column - POT is marked in uppercase)
		"""

        bed_file = ""
        chrom_sizes = {}
        with open(f"{self.genome_fasta}.fai") as f:
            for line in f:
                sl = line.strip().split("\t")
                chrom_sizes[sl[0]] = int(sl[1])
                
        for index, row in df.iterrows():
            chromosome = row["chromosome"]
            start = row["start"]
            end = row["end"]
            strand = row["strand"]
            pot = row["POT"]
            seq_start = pot - upstream - 1 if strand == "+" else pot - downstream - 1
            seq_start = 0 if seq_start < 0 else seq_start
            seq_end = pot + downstream if strand == "+" else pot + upstream
            seq_end = chrom_sizes[chromosome] if seq_end > chrom_sizes[chromosome] else seq_end
            bed_file += f"{chromosome}\t{seq_start}\t{seq_end}\t{row['termite_id']}\t.\t{strand}\n"
        bedtool = pybedtools.BedTool(bed_file, from_string=True)
        
        o = bedtool.sequence(fi=self.genome_fasta, s=True, name=True)
        with open(f"{output_file_prefix}.fasta", "w") as f:
            print(open(o.seqfn).read(), file=f)
            
        fasta_sequences = SeqIO.parse(open(f"{output_file_prefix}.fasta"),'fasta')
        seqs = {}
        for fasta in fasta_sequences:
            name, sequence = fasta.id.strip().split("::")[0], str(fasta.seq).lower()
            seqs[name] = name
            if fasta.id.strip().endswith("(+)"):
                strand = "+"
            elif fasta.id.strip().endswith("(-)"):
                strand == "-"
            else:
                raise ValueError(f'Invalid strand - {fasta.id.strip()[-2]}')
            pot = upstream if strand == "+" else downstream
            s = sequence[:pot] + sequence[pot].upper() + sequence[pot+1:]
            seqs[name] = s
        
        sequences = []
        for index, row in df.iterrows():
            sequences.append(seqs[row['termite_id']])
            
        df["sequence"] = sequences
        
        return df
    
    def intersectBEDs(self, df : pd.DataFrame, beds : List[str], strand_specific : bool = True) -> pd.DataFrame:
        """Check if the TERMITe's peak overlaps with the features from BED(s).
        Adds column marked with the name assigned to the BED file(s) and fills the column with "+" and "-" values,
        when the overlap has been or hasn't been found respectively.

        Args:
            df (pd.DataFrame): DataFrame with annotated TERMITe's results
            beds (List[str]): list of the bed files
            strand_specific (bool, optional): True if the search for the overlaps is to be strand-specific; False otherwise. Defaults to True.

        Returns:
            pd.DataFrame: data frame with additional columns describing overlaps
        """

        d = {}
        for bed_file in beds:
            sl = bed_file.strip().split(":")
            header, f = sl
            d[header] = []
            bedtool_a = pybedtools.BedTool.from_dataframe(df).sort()
            bedtool_b = pybedtools.BedTool(f)
            a_and_b = bedtool_a.intersect(bedtool_b, s=True) if strand_specific else bedtool_a.intersect(bedtool_b)
            for interval in a_and_b:
                d[header].append(interval[8])
        
            df[header] = [ "+" if row["termite_id"] in d[header] else "-" for x, row in df.iterrows() ]
            
        return df


class Terminator():
    """Describe Terminator structure and properties as returned by the TransTermHP
    """
    def __init__(self, _id, start, end, strand, confidence, hairpin_score, tail_score, representation):
        self.id = _id
        self.start = start
        self.end = end
        self.strand = strand # we should keep only those on positive strand as they correspond to the same strand on which we observe termseq signal
        self.confidence = confidence
        self.hairpin_score = hairpin_score
        self.tail_score = tail_score
        self.representation = representation.lower()
 
    
class RNAfold:
    """Predict the intrinsic terminator hairpin structures in the vicinity of the TERMITe signal.
    Prediction is made by running TransTermHP with very permissive options increasing the sensitivity of the method.
    TransTermHP scans sequences corresponding to the region starting 80nt upstream and ending 30nt downstream from the POT,
    so the sequence is 111nt long. The best one is selected from the predicted hairpins on the same strand as the POT.  
    To select the best hairpin, the following criteria are analyzed (checked in the specified order until only one hairpin remains):
    - hairpin confidence
    - hairpin score
    - tail score.
    Then RNAfold from the ViennaRNA package is run to re-fold the heuristically obtained hairpin structure (from TransTermHP).
    Class checks if the MFE structure is indeed the hairpin structure with energy <= -3 and
    if the POT lies at most 10nt downstream from the hairpin.
    Predicted intrinsic terminators meeting these criteria are reported in the final DataFrame.
    """

    def __init__(self, df : pd.DataFrame,  genome_fasta : str, max_distance_to_hairpin : int = 10, min_distance_to_hairpin : int = 0):
        """Run the class logic.

        Args:
            df (pd.DataFrame): pandas DataFrame to annotate
            genome_fasta (str): path to the fasta file with the genomic sequence
            max_distance_to_hairpin (int, optional): maximum distance between the POT and the intrinsic terminator hairpin located upstream. Defaults to 10.
            min_distance_to_hairpin (int, optional): minimum distance between the POT and the intrinsic terminator hairpin located upstream. Defaults to 0. 
        """

        self.df = df
        self.genome_fasta = genome_fasta
        self.random_file_name = "".join(random.choices(string.ascii_uppercase + string.digits + string.ascii_lowercase, k=20))
        logging.info("Running the RNAfold pipeline")
        self.extractSequence()
        self.calculate_MFE_RNAfold(max_distance_to_hairpin, min_distance_to_hairpin)
        
        fl = Path(f"{self.random_file_name}.fasta")
        if fl.is_file():
            os.system(f"rm {self.random_file_name}.fasta")
            
 
    def extractSequence(self, upstream : int = 80, downstream : int = 30) -> pd.DataFrame:
        """Extract the genomic sequence for the region starting [upstream] nt upstream from the reference POT up to [downstream] nt downstream,
		so that the region will have a total length of upstream+downstream+1 with the central position corresponding to the reference POT (Point of Observed [RNA] Termini).
		If the peak is on the reverse strand, the sequence will be reverse complemented. 

		Args:
			downstream (int): how many nt will be taken downstream from the POT to extract the sequence.
            upstream (int): how many nt will be taken upstream from the POT to extract the sequence.

        Returns:
            pd.DataFrame: annotated dataframe (with Sequence-TransTermHP column - POT is marked in uppercase)
		"""

        bed_file = ""
        for index, row in self.df.iterrows():
            chromosome = row["chromosome"]
            start = row["start"]
            end = row["end"]
            strand = row["strand"]
            pot = row["POT"]
            seq_start = pot - upstream - 1 if strand == "+" else pot - downstream - 1
            seq_start = max(seq_start, 0)
            seq_end = pot + downstream if strand == "+" else pot + upstream
            seq_end = max(seq_end, 0)
            bed_file += f"{chromosome}\t{seq_start}\t{seq_end}\t{row['termite_id']}\t.\t{strand}\n"
        bedtool = pybedtools.BedTool(bed_file, from_string=True)
        o = bedtool.sequence(fi=self.genome_fasta, s=True, name=True)
        
        
        with open(f"{self.random_file_name}.fasta", "w") as f:
            print(open(o.seqfn).read(), file=f)
            
        fasta_sequences = SeqIO.parse(open(f"{self.random_file_name}.fasta"),'fasta')
        seqs = {}
        for fasta in fasta_sequences:
            name, sequence = fasta.id.strip().split("::")[0], str(fasta.seq).lower()
            seqs[name] = name
            if fasta.id.strip().endswith("(+)"):
                strand = "+"
            elif fasta.id.strip().endswith("(-)"):
                strand == "-"
            else:
                raise ValueError(f'Invalid strand - {fasta.id.strip()[-2]}')
            pot = upstream if strand == "+" else downstream
            s = sequence[:pot] + sequence[pot].upper() + sequence[pot+1:]
            seqs[name] = s
        
        sequences = []
        for index, row in self.df.iterrows():
            sequences.append(seqs[row['termite_id']])
        
        self.df["temp_sequence"] = sequences
    
    def parse_transterm_output(self, output : str) -> List[Terminator]:
        """Parse the TransTermHP output file and returns the list of Terminator class objects.

        Args:
            output (str): file path to the TransTermHP results

        Returns:
           List[Terminator]: list of Terminator class objects representing terminators predicted by the TransTermHP
        """
        
        terminators = []
        lines = output.split("\n")
        for line_ind in range(0, len(lines)):
            if lines[line_ind].strip().startswith("TERM"):
                info = lines[line_ind].strip().split()
                representation = lines[line_ind+1].strip()
                if len(representation.strip().split()) == 5 and len(representation.strip().split()[0]) != 15:
                    representation = "   ".join(representation.strip().split()[-4:])
                term = Terminator(_id = f"{info[0]}_{info[1]}",
                                start = int(info[2]),
                                end = int(info[4]),
                                strand = info[5],
                                confidence = int(info[7]),
                                hairpin_score = float(info[8]),
                                tail_score = float(info[9]),
                                representation = representation)
                terminators.append(term)
        return terminators
    
    def select_hairpin(self, terminators : List[Terminator]) -> Terminator:
        """Select the best terminator predicted for a given POT with permissive criteria.
        To select the best hairpin, the following criteria are analyzed (checked in the specified order until only one hairpin remains):
        - hairpin confidence
        - hairpin score
        - tail score.

        Args:
            terminators (List[Terminator]): list of Terminator class objects representing terminators predicted by the TransTermHP

        Returns:
            Terminator: the best terminator predicted for a given POT
        """
        
        s = False
        hairpins_pos_strand = [x for x in terminators if x.strand == "+" or x.strand == None]
        if len(hairpins_pos_strand) == 0: return None
        max_confidence = max([x.confidence for x in hairpins_pos_strand])
        candidates = [x for x in hairpins_pos_strand if x.confidence == max_confidence]
        if len(candidates) == 0: return None
        min_hp_score = min([x.hairpin_score for x in candidates])
        candidates = [x for x in candidates if x.hairpin_score == min_hp_score]
        if len(candidates) == 0: return None
        min_tail_score = min([x.tail_score for x in candidates])
        candidates = [x for x in candidates if x.tail_score == min_tail_score]
        if len(candidates) == 0: return None
        if len(candidates) != 1:
            logging.warn("More than 1 hairpin was found with the exact confidence and hairpin score.")
        return candidates[0]
    
    
    def evaluate_hairpin(self, ss : str, hairpin : str, hairpin_end : int, row : pd.Series, energy : float, upstream_nt : int = 80, max_distance_to_hairpin : int = 10,
                         min_distance_to_hairpin : int = 0) -> Tuple[bool, int]:
        """Evaluate the hairpin predicted by the TransTermHP and RNAfold and check if the predicted MFE structure is indeed the hairpin structure with energy <= -3 and
        the POT lying at most [max_distance_to_hairpin] nt downstream from the hairpin.

        Args:
            ss (str): the MFE secondary structure in a dot-bracket format as folded by the RNAfold
            hairpin (str): hairpin structure as returned by the TransTermHP
            hairpin_end (int): end coordinate of the hairpin
            row (pd.Series): pandas DataFrame row with the POT to be annotated
            energy (float): min energy of the hairpin structure to be considered
            upstream_nt (int, optional): length of the sequence that is scanned for the terminators by the TransTermHP upstream from the POT. Defaults to 80.
            max_distance_to_hairpin (int, optional): maximum distance between the POT and the intrinsic terminator hairpin located upstream. Defaults to 10.
            min_distance_to_hairpin (int, optional): minimum distance between the POT and the intrinsic terminator hairpin located upstream. Defaults to 0. 

        Raises:
            ValueError: wrong hairpin representation

        Returns:
            Tuple[bool, int]: the first element of the tuple is True if the hairpin is meeting the specified criteria; False otherwise,
                              the second element is the distance from the POT to the hairpin located upstream if >= [max_distance_to_hairpin] nt; None otherwise
        """
        
        hp = hairpin.split()
        if len(hp) != 3:
            raise ValueError("Wrong hairpin representation: ", hp)
        stem1, loop, stem2 = hp[0], hp[1], hp[2]    
        if energy <= -3:
            if not ")" in ss[:len(stem1)] and not "(" in ss[-len(stem2):]:
                if re.search("\)\.+\(", ss) == None:
                    if (upstream_nt - int(hairpin_end) <= max_distance_to_hairpin and upstream_nt - int(hairpin_end) >= min_distance_to_hairpin):
                        return True, upstream_nt - int(hairpin_end)
        return False, None
    
    
    def calculate_MFE_RNAfold(self, max_distance_to_hairpin : int = 10, min_distance_to_hairpin : int = 0) -> None:
        """RNAfold from the ViennaRNA package is run to re-fold the heuristically obtained hairpin structure (from TransTermHP).
        Method checks if the MFE structure is indeed the hairpin structure with energy <= -3 and
        if the POT lies at [--max-distance-to-hairpin] nt downstream from the hairpin.
        Predicted intrinsic terminators meeting these criteria are reported in the final DataFrame (self.df).
        
        Args:
            max_distance_to_hairpin (int, optional): maximum distance between the POT and the intrinsic terminator hairpin located upstream. Defaults to 10.
            min_distance_to_hairpin (int, optional): minimum distance between the POT and the intrinsic terminator hairpin located upstream. Defaults to 0. 

        Raises:
            TransTermHP_Exception: TransTermHP seems not to be available in your $PATH. Aborting.
            TransTermHP_Exception: TransTermHP exception.
            TransTermHP_Exception: TransTermHP finished with the return code other than 0.
            ValueError: wrong hairpin representation.
            AttributeError: no hairpin was found using the specified criteria.

        Returns:
            None: None
        """

        d = {}
        mfes, sss, U_tract, A_tract, hp_seqs, starts, ends, res, distances = [], [], [], [], [], [], [], [], []

        for index, row in self.df.iterrows(): 
            with open("sequence.fasta", "w") as f:
                print(f">sequence\n{row['temp_sequence']}", file=f)
            
                
            try:
                script_path = os.path.dirname(os.path.realpath(__file__))
                #p = subprocess.Popen(["transterm", "--all-context", "--min-conf=0", "--min-stem=4", "--min-loop=4",
                #                        "--max-len=30", "--max-loop=20", "-p", f"{script_path}/newexpterm.dat", "sequence.fasta",
                #                    f"{script_path}/fakegenes.coords"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                p = subprocess.Popen(["transterm", "--all-context", "--min-conf=0", "--min-stem=4", "--min-loop=4",
                                        "--max-len=40", "--uwin-require=1", "--uwin-size=6", "--max-loop=20", "-p", f"{script_path}/newexpterm_L_and_I_shaped.dat", "sequence.fasta",
                                    f"{script_path}/fakegenes.coords"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                out, err = p.communicate()
                out, err = out.decode("latin-1"), err.decode("latin-1")
            except FileNotFoundError as e:
                raise TransTermHP_Exception(f"TransTermHP seems not to be available in your $PATH. Aborting.")
            except OSError as e:
                raise TransTermHP_Exception(f"TransTermHP exception: {e}.")
            
            if p.returncode != 0:
                if "Bad gene coordinates" in err:
                    mfes.append(".")
                    sss.append(".")
                    U_tract.append(".")
                    A_tract.append(".")
                    hp_seqs.append(".")
                    starts.append(".")
                    ends.append(".")
                    res.append("-")
                    distances.append(".")
                    continue
                logging.error(err)
                raise TransTermHP_Exception(f"TransTermHP finished with the return code: {p.returncode}")

            terminators = self.parse_transterm_output(out)
            terminator = self.select_hairpin(terminators)
            start = None if terminator is None else str(terminator.start)
            end = None if terminator is None else str(terminator.end)
            
            try:
                representation = terminator.representation.strip().split()
                if len(representation) == 4:
                    hairpin = terminator.representation[:-15]
                    u = terminator.representation[-15:]
                    a = None
                elif len(representation) == 5:
                    hairpin = terminator.representation[15:-15]
                    u = terminator.representation[-15:]
                    a = terminator.representation[:15]
                else:
                    raise ValueError(f"Wrong hairpin representation - {terminator.representation}.")

                hp_seq = hairpin.replace(" ", "").replace("-", "")
                ss, mfe = ViennaRNA.RNA.fold(hp_seq)
                evl, distance_to_hairpin = self.evaluate_hairpin(ss, hairpin, end, row, mfe, 80, max_distance_to_hairpin, min_distance_to_hairpin)
                if evl is False:
                    raise AttributeError()
                mfes.append(mfe)
                sss.append(ss)
                hp_seqs.append(hp_seq)
                U_tract.append(u)
                A_tract.append(a)
                starts.append(start)
                ends.append(end)
                distances.append(distance_to_hairpin)
                if evl:
                    res.append("+")
                else:
                    res.append("-")

            except AttributeError:
                mfes.append(".")
                sss.append(".")
                U_tract.append(".")
                A_tract.append(".")
                hp_seqs.append(".")
                starts.append(".")
                ends.append(".")
                res.append("-")
                distances.append(".")
            
        self.df["rnafold_a_tract"] = A_tract
        self.df["rnafold_hairpin"] = hp_seqs
        self.df["rnafold_u_tract"] = U_tract
        self.df["rnafold_hairpin_structure"] = sss
        self.df["rnafold_distance_to_hairpin"] = distances
        self.df["rnafold_energy"] = mfes
        self.df["rnafold"] = res
        
        del self.df["temp_sequence"]
        


              
        
        
            
            
            

            
        
        
    
           
            
