import argparse
from pathlib import Path
import logging
import itertools
from typing import List, Set, Dict, Tuple, Optional, Union, Sequence
import sys
import subprocess
import os

class RNAfold_Exception(Exception):
    pass

class BEDTools_Exception(Exception):
    pass

class IDR_Exception(Exception):
    pass

def parse_arguments_termite():

	parser = argparse.ArgumentParser(description="""TERMITe: a tool for identifying and analyzing TERMInaTion signals in the term-seq data.""", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	
	subprasers = parser.add_subparsers(dest='command')
	stable_rna_ends_subparser = subprasers.add_parser('find_stable_rna_ends', help="Identification of stable 3'RNA ends based on the term-seq data.",  epilog="Thank you for using TERMITe.")
	stable_rna_ends_subparser.add_argument('--rna-3prime-ends', nargs="+", type=str, required=True,
						help="bigWig files with normalized coverage of RNA 3' ends", dest='rna_3prime_ends', metavar=("TERMSEQ-BIGWIG"))
	stable_rna_ends_subparser.add_argument('--rna-seq-coverage', nargs="+", type=str, required=False, default=None,
						help="bigWig files with normalized RNA-seq coverage. Files corresponding to the same samples should be placed in the same order as for --rna-3prime-ends",
      					dest='rnaseq', metavar=("RNASEQ-BIGWIG"))
	stable_rna_ends_subparser.add_argument('--sample-names', nargs="+", type=str, required=True,
						help="sample names. Should be placed in the same order as corresponding files in --rna-3prime-ends and --rna-seq-coverage",
      					dest='sample_names', metavar=("SAMPLE-NAME"))
	stable_rna_ends_subparser.add_argument('--out-dir', type=str, required=True,
						help="output directory", dest='out_dir')
	stable_rna_ends_subparser.add_argument('--strand', type=str, required=True, choices=["forward", "reverse"],
						help="strand", dest='strand')
	stable_rna_ends_subparser.add_argument('--idr-threshold', type=float, required=False, default=0.05,
						help="IDR threshold. Default is 0.05", dest='idr')
	stable_rna_ends_subparser.add_argument('--min-no-comp', type=int, required=False, default=1,
						help="IDR must be smaller than the threshold in a minimum [--min_no_comp] pairwise comparisons between samples. The total number of comparisons equals the number of possible two-element combinations of samples provided by --sample_names. Default is 1", dest='min_no_comp')
	stable_rna_ends_subparser.add_argument('--min-distance', type=int, required=False, default=20,
						help="the minimum distance between peaks. The peak will not be reported if there is a higher signal +/- [--min-distance] nt from the peak. Default is 20", dest='min_dist')
	stable_rna_ends_subparser.add_argument('--min-term-effi', type=float, required=False, default=None,
						help="filter out signals with termination efficiency below this threshold (percentage). Must be used together with [--rna-seq-coverage]. Default is None", dest='min_effi')
	stable_rna_ends_subparser.add_argument('--min-upstream', type=float, required=False, default=0.25,
						help="do not calculate termination efficiency if the average upstream signal is less than [--min-upstream]. Default is 0.25", dest='min_upstream')
	stable_rna_ends_subparser.add_argument('--genome', type=str, required=False, default=None,
						help="genome FASTA file to extract the sequences from (+/- [--window] nt from the peak summit). Sequences around the peaks on the reverse strand will be reverse complemented. Currently, this tool does not support gzipped fasta files. Default is None", dest='genome')
	stable_rna_ends_subparser.add_argument('--window-size', type=int, required=False, default=10,
						help="windows size describing how many nt will be taken around peak summit to extract the sequence. Must be used together with --genome to have an effect. By default is set to 10, so that 10nt upstream and 10nt downstream from the summit will be extracted (21nt in total). Default is 10", dest='window_size')
	stable_rna_ends_subparser.add_argument('--name-prefix', type=str, required=True, help="add this prefix to the ID of each output 3' RNA termini", dest='name_prefix')
	
	annotate_subparser = subprasers.add_parser('annotate', help="Annotation of stable 3'RNA ends identified by the TERMITe.", epilog="Thank you for using TERMITe.")
	annotate_subparser.add_argument('--termite-out-forward', nargs="+", type=str, required=True,
						help="list of the TERMITe output files (narrowPeak) calculated for the forward strand", dest='termite_out_forward', metavar=("TERMITe-OUTPUT-FWD"))
	annotate_subparser.add_argument('--termite-out-reverse', nargs="+", type=str, required=True,
						help="list of the TERMITe output files (narrowPeak) calculated for the reverse strand. Files belonging to the same samples should be specified in the same order as in the --termite-out-forward",
	  					dest='termite_out_reverse', metavar=("TERMITe-OUTPUT_REV"))
	annotate_subparser.add_argument('--sample-names', nargs="+", type=str, required=True,
						help="sample names. Should be placed in the same order as corresponding files in --termite-out-forward and --termite-out-reverse. Should be the same as --name-prefix(es) provided in the TERMITe module",
	  					dest='sample_names', metavar=("SAMPLE-NAME"))
	annotate_subparser.add_argument('--rna-3prime-ends-forward', nargs="+", type=str, required=True,
						help="bigWig files with the genomic coverage by the 3'RNA ends on the forward strand. Should be provided as follows: 'sample:path_to_file1,path_to_file2,...,path_to_fileN'",
	  					dest='rna_3prime_ends_forward', metavar=("TERMSEQ-BIGWIG-FWD"))
	annotate_subparser.add_argument('--rna-3prime-ends-reverse', nargs="+", type=str, required=True,
						help="bigWig files with the genomic coverage by the 3'RNA ends on the reverse strand. Should be provided as follows: 'sample:path_to_file1,path_to_file2,...,path_to_fileN'",
	  					dest='rna_3prime_ends_reverse', metavar=("TERMSEQ-BIGWIG-REV"))
	annotate_subparser.add_argument('--rna-seq-coverage-forward', nargs="+", type=str, required=False,
						help="bigWig files with the genomic coverage by the whole RNA-seq reads on the forward strand. Should be provided as follows: 'sample:path_to_file1,path_to_file2,...,path_to_fileN'",
	  					dest='rna_seq_coverage_forward', metavar=("RNASEQ-BIGWIG-FWD"))
	annotate_subparser.add_argument('--rna-seq-coverage-reverse', nargs="+", type=str, required=False,
						help="bigWig files with the genomic coverage by the whole RNA-seq reads on the reverse strand. Should be provided as follows: 'sample:path_to_file1,path_to_file2,...,path_to_fileN'",
	  					dest='rna_seq_coverage_reverse', metavar=("RNASEQ-BIGWIG-REV"))
	annotate_subparser.add_argument('--intersect', nargs="+", type=str, required=False,
						help="report 3' RNA end overlaps with the features in the [--intersect] BED file(s). Multiple bed files can be provided using this option. All files should be specified as 'NAME:PATH_TO_THE_BED_FILE', where NAME will be used as the column name in the TERMITe report and therefore must be unique. BED files should be in a standard BED6 format",
	  					dest='bed', metavar=("BED"))
	annotate_subparser.add_argument('--min-upstream', type=float, required=False, default=0.25,
						help="do not calculate termination efficiency if the average upstream signal is less than [--min-upstream]. Default is 0.25", dest='min_upstream')
	annotate_subparser.add_argument('--max-distance-to-hairpin', type=int, required=False, default=10,
						help="the maximum distance between the POT and the upstream terminator hairpin", dest='max_distance_to_hairpin')
	annotate_subparser.add_argument('--min-distance-to-hairpin', type=int, required=False, default=0,
						help="the minimum distance between the POT and the intrinsic end of the upstream terminator hairpin. Negative values are possible when the POT overlaps the hairpin", dest='min_distance_to_hairpin')
	annotate_subparser.add_argument('--trans-term-hp', action="store_true", help="annotate results with the intrinsic terminator hairpins predicted by the TransTermHP software. To run the TransTermHP, you need to specify the --genome option", dest='trans_term_hp')
	annotate_subparser.add_argument('--rnafold', action="store_true", help="annotate results with the intrinsic terminator hairpins predicted by the pipeline based on TransTermHP and RNAfold. To use this option, you need to specify the --genome option", dest='rnafold')
	annotate_subparser.add_argument('--gene-annotations', type=str, required=False, help="GTF/GFF file with gene annotations", default=None, dest='gene_annotations')
	annotate_subparser.add_argument('--genome', type=str, required=False, help="reference genome in the FASTA file format", default=None, dest='genome')
	annotate_subparser.add_argument('--upstream-nt', type=int, required=False, default=10, help="extract the sequence consisting of this many nt upstream from the POT, by default is set to 10, so that [--upstream-nt] nt upstream and [--downstream-nt] nt downstream from the POT will be extracted ([--upstream-nt] + 1 + [--downstream-nt] nt in total). Default is 10",
					 	dest='upstream_nt')
	annotate_subparser.add_argument('--downstream-nt', type=int, required=False, default=10, help="extract the sequence consisting of this many nt upstream from the POT, by default is set to 10, so that [--upstream-nt] nt upstream and [--downstream-nt] nt downstream from the POT will be extracted ([--upstream-nt] + 1 + [--downstream-nt] nt in total). Default is 10",
					 	dest='downstream_nt')
	annotate_subparser.add_argument('--output', type=str, required=True, help="a prefix of the output file name", dest='output')

	if len(sys.argv)==1:
		parser.print_help(sys.stderr)
		sys.exit(1)
	args = parser.parse_args()

	if args.command == "find_stable_rna_ends":
		_ = validate_args_termite_main_module(args)
	elif args.command == "annotate":
		_ = validate_args_termite_annotate_module(args)

	return args

def validate_args_termite_main_module(args : argparse.Namespace):
	"""Validates arguments passed to the main TERMITe's module.

	Args:
		args (argparse.Namespace): argparse object (argument parser)

	Raises:
		argparse.ArgumentTypeError: --min_upstream must be greater or equal than 0.
		argparse.ArgumentTypeError: The number of elements provided in the --rna-3prime-ends, --rna-seq-coverage, and --sample-names must be equal.
		argparse.ArgumentTypeError: TERMITe needs to be run with at least two replicates [--sample-names].
		argparse.ArgumentTypeError: Sample names provided with --sample-names must be unique.
		argparse.ArgumentTypeError: The number of elements provided in the --rna-3prime-ends, --rna-seq-coverage, and --sample-names must be equal.
		argparse.ArgumentTypeError: No such file.
		argparse.ArgumentTypeError: --min-term-effi must be used together with --rna-seq-coverage.
		argparse.ArgumentTypeError: --min-term-effi must be in range <0, 100> (percentage).
		argparse.ArgumentTypeError: No such directory.
		argparse.ArgumentTypeError: str '_vs_' cannot be included in the sample name.
		argparse.ArgumentTypeError: IDR threshold (--idr-threshold) must be in the range (0, 1].
		argparse.ArgumentTypeError: --min-no-comp must be in range <1, number of pairwise comparisons>.
		argparse.ArgumentTypeError: --min-distance must be >= 0.
	"""

	try:
		p = subprocess.Popen(["idr", "--help"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		out, err = p.communicate()
		out, err = out.decode("utf-8"), err.decode("utf-8")
	except FileNotFoundError as e:
		raise IDR_Exception(f"idr seems not to be available in your $PATH. Aborting.")
	except OSError as e:
		raise IDR_Exception(f"idr exception: {e}.")	
 
	if args.min_upstream < 0:
		logging.error('--min_upstream must be greater or equal to 0.')
		raise argparse.ArgumentTypeError("--min_upstream must be greater or equal to 0.")

	if not (len(args.sample_names) == len(args.rna_3prime_ends)):
		logging.error('The number of elements provided in the --rna-3prime-ends, --rna-seq-coverage, and --sample-names must be equal.')
		raise argparse.ArgumentTypeError("The number of elements provided in the --rna-3prime-ends, --rna-seq-coverage, and --sample-names must be equal.")

	if len(args.sample_names) < 2:
		logging.error('TERMITe needs to be run with at least two replicates [--sample-names].')
		raise argparse.ArgumentTypeError("TERMITe needs to be run with at least two replicates [--sample-names].")

	if len(args.sample_names) != len(set(args.sample_names)):
		logging.error('Sample names provided with --sample-names must be unique.')
		raise argparse.ArgumentTypeError('Sample names provided with --sample-names must be unique.')
		
	if isinstance(args.rnaseq, list):
		if len(args.rnaseq) != len(args.sample_names):
			logging.error('The number of elements provided in the --rna-3prime-ends, --rna-seq-coverage, and --sample-names must be equal.')
			raise argparse.ArgumentTypeError("The number of elements provided in the --rna-3prime-ends, --rna-seq-coverage, and --sample-names must be equal.")
		for file in args.rnaseq:
			my_file = Path(file)
			if not my_file.is_file():
				raise argparse.ArgumentTypeError(f"No such file: {file}.")

	if args.min_effi is not None and args.rnaseq is None:
		logging.error('--min-term-effi must be used together with --rna-seq-coverage.')
		raise argparse.ArgumentTypeError("--min-term-effi must be used together with --rna-seq-coverage.")

	if args.min_effi is not None and (args.min_effi < 0 or args.min_effi > 100):
		logging.error('--min-term-effi must be in range <0, 100> (percentage).')
		raise argparse.ArgumentTypeError("--min-term-effi must be in range <0, 100> (percentage).")

	for file in args.rna_3prime_ends:
		my_file = Path(file)
		if not my_file.is_file():
			raise argparse.ArgumentTypeError(f"No such file: {file}.")

	if args.genome is not None:
		my_file = Path(args.genome)
		if not my_file.is_file():
			raise argparse.ArgumentTypeError(f"No such file: {args.genome}.")

	out_dir = Path(args.out_dir)
	if not out_dir.is_dir():
		os.system(f"mkdir -p {out_dir}")

	for name in args.sample_names:
		if "_vs_" in name:
			raise argparse.ArgumentTypeError(f"str '_vs_' cannot be included in the sample name.")

	if args.idr <= 0 or args.idr > 1:
		raise argparse.ArgumentTypeError(f"IDR threshold (--idr-threshold) must be in the range (0, 1].")

	if args.min_no_comp < 1 or args.min_no_comp > len(list(itertools.combinations(args.sample_names, 2))):
		raise argparse.ArgumentTypeError(f"--min-no-comp must be in range <1, number of pairwise comparisons>.")

	if args.min_dist < 0:
		raise argparse.ArgumentTypeError(f"--min-distance must be >= 0.")
	return

def validate_bigwigs(forward_bigwigs : list, reverse_bigwigs : list, sample_names : list) -> None:
	"""Validates the input bigwig files.

	Args:
		forward_bigwigs (list): the list of the bigwig files for the forward strand
		reverse_bigwigs (list): the list of the bigwig files for the reverse strand
		sample_names (list): the list of sample names

	Raises:
		argparse.ArgumentTypeError: The number of elements provided in the --sample-names, --rna-3prime-ends-forward, --rna-3prime-ends-reverse, --rna-seq-coverage-forward and --rna-seq-coverage-reverse must be equal.
		argparse.ArgumentTypeError: At least one of --rna-3prime-ends-forward or --rna-seq-coverage-forward has a wrong format.
		argparse.ArgumentTypeError: Sample name provided in --rna-3prime-ends-forward or --rna-seq-coverage-forward not found in --sample-names.
		argparse.ArgumentTypeError: At least one of --rna-3prime-ends-reverse or --rna-seq-coverage-reverse has a wrong format.
		argparse.ArgumentTypeError: Sample name provided in --rna-3prime-ends-reverse or --rna-seq-coverage-reverse not found in --sample-names.
		argparse.ArgumentTypeError: -rna-3prime-ends-forward, -rna-3prime-ends-reverse, --rna-seq-coverage-forward and --rna-seq-coverage-reverse must specify the same number of input files.
	"""

	d_fwd, d_rev = {}, {}
	if not(len(forward_bigwigs) == len(reverse_bigwigs) and len(reverse_bigwigs) == len(sample_names)):
		logging.error('The number of elements provided in the --sample-names, --rna-3prime-ends-forward, --rna-3prime-ends-reverse, --rna-seq-coverage-forward and --rna-seq-coverage-reverse must be equal.')
		raise argparse.ArgumentTypeError('The number of elements provided in the --sample-names, --rna-3prime-ends-forward, --rna-3prime-ends-reverse, --rna-seq-coverage-forward and --rna-seq-coverage-reverse must be equal.')
	
	for sample in forward_bigwigs:
		splitted_arg = sample.strip().split(":")
		if len(splitted_arg) != 2:
			logging.error(f'At least one of --rna-3prime-ends-forward or --rna-seq-coverage-forward has a wrong format - {sample}.')
			raise argparse.ArgumentTypeError(f'At least one of --rna-3prime-ends-forward or --rna-seq-coverage-forward has a wrong format - {sample}.')
		s_name = splitted_arg[0]
		if s_name not in sample_names:
			logging.error(f'Sample name - {s_name} - provided in --rna-3prime-ends-forward or --rna-seq-coverage-forward not found in --sample-names.')
			raise argparse.ArgumentTypeError(f'Sample name - {s_name} - provided in --rna-3prime-ends-forward or --rna-seq-coverage-forward not found in --sample-names.')
		d_fwd[s_name] = splitted_arg[1].strip().split(",")
	
	for sample in reverse_bigwigs:
		splitted_arg = sample.strip().split(":")
		if len(splitted_arg) != 2:
			logging.error(f'At least one of --rna-3prime-ends-reverse or --rna-seq-coverage-reverse has a wrong format - {sample}.')
			raise argparse.ArgumentTypeError(f'At least one of --rna-3prime-ends-reverse or --rna-seq-coverage-reverse has a wrong format - {sample}.')
		s_name = splitted_arg[0]
		if s_name not in sample_names:
			logging.error(f'Sample name - {s_name} - provided in --rna-3prime-ends-reverse or --rna-seq-coverage-reverse not found in --sample-names.')
			raise argparse.ArgumentTypeError(f'Sample name - {s_name} - provided in --rna-3prime-ends-reverse or --rna-seq-coverage-reverse not found in --sample-names.')
		d_rev[s_name] = splitted_arg[1].strip().split(",")
	
	# check if the same number of input bigwigs was provided for forward and reverse strands
	for sample in sample_names:
		if len(d_fwd[sample]) != len(d_rev[sample]):
			logging.error(f'-rna-3prime-ends-forward, -rna-3prime-ends-reverse, --rna-seq-coverage-forward and --rna-seq-coverage-reverse must specify the same number of input files.')
			raise argparse.ArgumentTypeError(f'-rna-3prime-ends-forward, -rna-3prime-ends-reverse, --rna-seq-coverage-forward and --rna-seq-coverage-reverse must specify the same number of input files.')
	
	# check if the input files exist
	bigwigs = list(d_rev.values()) + list(d_fwd.values())
	bigwigs = list(itertools.chain(*bigwigs))
	for file in bigwigs:
		my_file = Path(file)
		if not my_file.is_file():
			raise argparse.ArgumentTypeError(f"No such file: {file}")
			
def validate_BED6_file(file : str):
	"""Validates the input bed6 files.

	Args:
		file (str): path to the input BED6 file.

	Raises:
		argparse.ArgumentTypeError: The file is not in the standard BED6 file format.
		argparse.ArgumentTypeError: The start coordinate in the file is not a number.
		argparse.ArgumentTypeError: The end coordinate in the file is not a number.
		argparse.ArgumentTypeError: The score value in the file is not a number nor a ".".
		argparse.ArgumentTypeError: The strand in the file should be one of ".", "+", "-".
	"""

	with open(file) as f:
		for line in f:
			sl = line.strip().split()
			if len(sl) != 6:
				logging.error(f'The {file} is not in the standard BED6 file format. Line: {line}.')
				raise argparse.ArgumentTypeError(f'The {file} is not in the standard BED6 file format. Line: {line}.')
			try:
				start = int(sl[1])
			except ValueError:
				logging.error(f'The start coordinate in the {file} file is not a number. Line: {line}.')
				raise argparse.ArgumentTypeError(f'The start coordinate in the {file} file is not a number. Line: {line}.')
			try:
				end = int(sl[2])
			except ValueError:
				logging.error(f'The end coordinate in the {file} file is not a number. Line: {line}.')
				raise argparse.ArgumentTypeError(f'The end coordinate in the {file} file is not a number. Line: {line}.')
			try:
				if sl[4] != ".":
					score = int(sl[4])
			except ValueError:
				logging.error(f'The score value in the {file} file is not a number, neither a ".". Line: {line}.')
				raise argparse.ArgumentTypeError(f'The score value in the {file} file is not a number, neither a ".". Line: {line}.')
			if sl[5] != "." and sl[5] != "+" and sl[5] != "-":
				logging.error(f'The strand in the {file} file should be one of ".", "+", "-". Line: {line}.')
				raise argparse.ArgumentTypeError(f'The strand in the {file} file should be one of ".", "+", "-". Line: {line}.')
                
                
			
def validate_args_termite_annotate_module(args : argparse.Namespace) -> None:
	"""Validates arguments passed to the TERMITe's annotate module.

	Args:
		args (argparse.Namespace): argparse object (argument parser)

	Raises:
		argparse.ArgumentTypeError: --min_upstream must be greater or equal than 0.
		argparse.ArgumentTypeError: The number of elements provided in the --termite-out-forward, --termite-out-reverse, and --sample-names must be equal.
		argparse.ArgumentTypeError: --trans-term-hp must be used together with --genome.
  		argparse.ArgumentTypeError: "--max_distance_to_hairpin must be greater or equal than 1.
  		argparse.ArgumentTypeError: --rnafold must be used together with --genome.
		argparse.ArgumentTypeError: Sample names provided with --sample-names must be unique.
		argparse.ArgumentTypeError: Wrong format of the --intersect argument.
		argparse.ArgumentTypeError: BED file NAME's length is 0.
		argparse.ArgumentTypeError: BED file PATH's length is 0.
		argparse.ArgumentTypeError: NAMES of the BED files must be unique.
		argparse.ArgumentTypeError: No such file.
	"""

	if args.rnafold:
		try:
			p = subprocess.Popen(["RNAfold", "--help"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
			out, err = p.communicate()
			out, err = out.decode("utf-8"), err.decode("utf-8")
		except FileNotFoundError as e:
			raise RNAfold_Exception(f"RNAfold seems not to be available in your $PATH. Aborting.")
		except OSError as e:
			raise RNAfold_Exception(f"RNAfold exception: {e}.")

	try:
		p = subprocess.Popen(["bedtools", "--help"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		out, err = p.communicate()
		out, err = out.decode("utf-8"), err.decode("utf-8")
	except FileNotFoundError as e:
		raise BEDTools_Exception(f"BEDTools seems not to be available in your $PATH. Aborting.")
	except OSError as e:
		raise BEDTools_Exception(f"BEDTools exception: {e}.")

	if args.max_distance_to_hairpin < 1:
		logging.error('--max_distance_to_hairpin must be greater or equal to 1.')
		raise argparse.ArgumentTypeError("--max_distance_to_hairpin must be greater or equal to 1.")
	
	if args.min_upstream < 0:
		logging.error('--min-upstream must be greater or equal to 0.')
		raise argparse.ArgumentTypeError("--min_upstream must be greater or equal to 0.")

	if not (len(args.sample_names) == len(args.termite_out_reverse) and len(args.sample_names) == len(args.termite_out_forward)):
		logging.error('The number of elements provided in the --termite-out-forward, --termite-out-reverse, and --sample-names must be equal.')
		raise argparse.ArgumentTypeError("The number of elements provided in the --termite-out-forward, --termite-out-reverse, and --sample-names must be equal.")

	if args.trans_term_hp and args.genome == None:
		logging.error('--trans-term-hp must be used together with --genome.')
		raise argparse.ArgumentTypeError('--trans-term-hp must be used together with --genome.')

	if args.rnafold and args.genome == None:
		logging.error('--rnafold must be used together with --genome.')
		raise argparse.ArgumentTypeError('--rnafold must be used together with --genome.')

	if len(args.sample_names) != len(set(args.sample_names)):
		logging.error('Sample names provided with --sample-names must be unique.')
		raise argparse.ArgumentTypeError('Sample names provided with --sample-names must be unique.')

	bed_names, bed_files = [], []
	if args.bed is not None:
		for bed_file in args.bed:
			sl = bed_file.strip().split(":")
			if len(sl) != 2:
				logging.error(f'Wrong format of the --intersect argument - {bed_file}.')
				raise argparse.ArgumentTypeError(f'Wrong format of the --intersect argument - {bed_file}.')
			if len(sl[0]) == 0:
				logging.error(f'BED file NAME\'s length is 0.')
				raise argparse.ArgumentTypeError(f'BED file NAME\'s length is 0.')
			if len(sl[1]) == 0:
				logging.error(f'BED file PATH\'s length is 0.')
				raise argparse.ArgumentTypeError(f'BED file PATH\'s length is 0.')
			bed_names.append(sl[0])
			bed_files.append(sl[1])
			validate_BED6_file(sl[1])
   
	if len(bed_names) != len(set(bed_names)):
		logging.error(f'NAMES of the BED files must be unique: {bed_names}.')
		raise argparse.ArgumentTypeError(f'NAMES of the BED files must be unique: {bed_names}.')

	files_to_test =  args.termite_out_reverse + args.termite_out_forward + [args.gene_annotations, args.genome] + bed_files
	files_to_test = [x for x in files_to_test if x is not None]
	for file in files_to_test:
		my_file = Path(file)
		if not my_file.is_file():
			raise argparse.ArgumentTypeError(f"No such file: {file}.")

	validate_bigwigs(forward_bigwigs=args.rna_3prime_ends_forward, reverse_bigwigs=args.rna_3prime_ends_reverse, sample_names=args.sample_names)
	if args.rna_seq_coverage_forward is not None or args.rna_seq_coverage_reverse is not None:
		validate_bigwigs(forward_bigwigs=args.rna_seq_coverage_forward, reverse_bigwigs=args.rna_seq_coverage_reverse, sample_names=args.sample_names)

	return

