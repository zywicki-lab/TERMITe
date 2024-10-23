# TERMITe

<b>a tool for identifying and analyzing transcription termination signals based on bacterial term-seq data</b>

<hr>
<br>

TERMITe identifies stable 3'RNA ends in data from replicated Term-seq experiments, annotates them to understand their origin and extract intrinsic transcription termination signals.

Part of the TERMITe software is based on the algorithm implemented in the `Term-seq peak-caller` [Adams PP, et al. 2021; https://github.com/NICHD-BSPC/termseq-peaks] but with several essential adjustments, improvements, and a module dedicated to the annotation of the identified stable 3'RNA ends.

`Philip P Adams Gabriele Baniulyte Caroline Esnault Kavya Chegireddy Navjot Singh Molly Monge Ryan K Dale Gisela Storz Joseph T Wade (2021) Regulatory roles of Escherichia coli 5' UTR and ORF-internal RNAs detected by 3' end mapping eLife 10:e62438. https://doi.org/10.7554/eLife.62438`

<br>
<br>

## Introduction

Term-seq is a high-throughput sequencing method allowing researchers to sequence only 3' ends of the RNA molecules in bacteria. TERMITe was developed to help researchers identify stable 3'RNA ends in somewhat noisy datasets. The noise in the term-seq data can represent intermediate RNA decay or processing products, unproperly mapped sequencing reads, ends of RNAs that have been fragmented during the library preparation, and other sequencing or library preparation-related artifacts. Furthermore, TERMITe contains a dedicated module for annotating the identified 3'RNA ends, providing extensive information about their genomic and transcriptomic context and the presence of predicted intrinsic terminators in the vicinity of the term-seq signal. In addition, the User can quickly obtain information about the termination efficiency at a given termination site, surrounding sequence, and more.

Two variations of the Term-seq protocol are described in the following papers:

`Mondal S, Yakhnin AV, Sebastian A, Albert I, Babitzke P. NusA-dependent transcription termination prevents misregulation of global gene expression. Nat Microbiol. 2016 Jan 11;1:15007. doi: 10.1038/nmicrobiol.2015.7. PMID: 27571753; PMCID: PMC5358096.`

`Dar D, Shamir M, Mellin JR, Koutero M, Stern-Ginossar N, Cossart P, Sorek R. Term-seq reveals abundant ribo-regulation of antibiotics resistance in bacteria. Science. 2016 Apr 8;352(6282):aad9822. doi: 10.1126/science.aad9822. PMID: 27120414; PMCID: PMC5756622.`

<br>
<br>

## Installation

Due to a significant number of external dependencies, we highly recommend using Docker to install and use TERMITe. Docker can be downloaded from [this link](https://www.docker.com).

You can directly download the docker image (AMD64) from docker hub

```{bash}
docker pull jankdocker1994/termite
docker tag jankdocker1994/termite termite
```

or use the following commands to clone the TERMITe repository and build a docker container

```{bash}
git clone https://github.com/zywicki-lab/TERMITe.git
cd TERMITe
docker build . -t termite
```

<hr>

## Help

To see the manual and available options, run

	docker run --rm -v $(pwd):/data --platform=linux/amd64 termite --help

where $(pwd) is a location of a working directory, usually containing data for analysis.

Running TERMITe is a two-step process. First, you need to invoke the module dedicated to the identification of the stable 3'RNA ends (`find_stable_rna_ends`) and then optionally annotate the results (`annotate`).

	docker run --rm -v $(pwd):/data --platform=linux/amd64 termite find_stable_rna_ends --help

	docker run --rm -v $(pwd):/data --platform=linux/amd64 termite annotate --help

<br>
<br>

## Input

TERMITe requires at least two biological replicates to identify statistically reproducible peaks in the Term-seq data.

Input files include:
- (required) the normalized genomic coverage by the 3' RNA ends calculated based on the Term-seq reads mapped to the reference genome (bigwig)
- (optionally) the normalized genomic coverage by the RNA-seq reads mapped to the reference genome (bigwig). RNA-seq samples must correspond to the Term-seq samples and be prepared from the same biological material.
- (optionally) the fasta file with the reference genome
- (optionally) the GTF/GFF3 file with genome annotations (must be compatible with the fasta file)

Both RNA-seq read coverage and coverage by 3' RNA ends should be computed for each strand separately.

Coverage by 3' RNA ends can be easily calculated using `bamCoverage` from the deeptools package as follows:

	bamCoverage --bam input_file.bam -o output_file.bw --normalizeUsing CPM --outFileFormat bigwig --binSize 1 --Offset 1 --samFlagExclude 256 --filterRNAstrand forward

	bamCoverage --bam input_file.bam -o output_file.bw --normalizeUsing CPM --outFileFormat bigwig --binSize 1 --Offset 1 --samFlagExclude 256 --filterRNAstrand reverse

`--Offset 1` will guarantee that to compute the coverage, only the 5' ends of the reads will be used (corresponding to the 3' ends of RNA molecules)

RNA-seq coverage can be computed on BAM files containing mapped RNA-seq reads almost the same way. The only suggested change is to remove the `--Offset` flag.

<hr>
<br>
<br>
<br>


## Running TERMITe find_stable_rna_ends

<b>a module dedicated for identification of stable 3'RNA ends in the data from replicated Term-seq experiments</b>

<br>
<br>

#### Help

To see the manual and available options, run

	docker run --rm -v $(pwd):/data --platform=linux/amd64 termite find_stable_rna_ends --help

<hr><br>

	usage: termite find_stable_rna_ends [-h] --rna-3prime-ends TERMSEQ-BIGWIG [TERMSEQ-BIGWIG ...] [--rna-seq-coverage RNASEQ-BIGWIG [RNASEQ-BIGWIG ...]] --sample-names SAMPLE-NAME [SAMPLE-NAME ...] --out-dir OUT_DIR --strand {forward,reverse} [--idr-threshold IDR] [--min-no-comp MIN_NO_COMP] [--min-distance MIN_DIST] [--min-term-effi MIN_EFFI] [--min-upstream-cpm MIN_UPSTREAM] [--genome GENOME] [--window-size WINDOW_SIZE] --name-prefix NAME_PREFIX


##### Required options

<b>--rna-3prime-ends</b> : bigWig files with normalized coverage of RNA 3' ends

<b>--sample-names</b> : sample names. Should be placed in the same order as corresponding files in --rna-3prime-ends and --rna-seq-coverage

<b>--out-dir OUT_DIR</b> : output directory

<b>--strand {forward,reverse}</b> : strand

<b>--name-prefix NAME_PREFIX</b> : add this prefix to the ID of each output 3' RNA termini

##### Additional options

<b>-h</b>, <b>--help</b> : show this help message and exit
  
<b>--rna-seq-coverage</b> : bigWig files with normalized RNA-seq coverage. Files corresponding to the same samples should be placed in the same order as for --rna-3prime-ends
			
<b>--idr-threshold</b> : IDR threshold. Default is 0.05  

<b>--min-no-comp</b> : IDR must be smaller than the threshold in a minimum [--min_no_comp] pairwise comparisons between samples. The total number of comparisons equals the number of possible two-element combinations of samples provided by --sample_names. Default is 1 

<b>--min-distance</b> : the minimum distance between peaks. The peak will not be reported if there is a higher signal +/- [--min-distance] nt from the peak. Default is 20
  
<b>--min-term-effi</b> : filter out signals with termination efficiency below this threshold (percentage). Must be used together with [--rna-seq-coverage]. Default is None

<b>--min-upstream</b> : do not calculate termination efficiency if the average upstream signal is less than [--min-upstream]. Default is 0.25
  
<b>--genome</b> : genome FASTA file to extract the sequences from (+/- [--window] nt from the peak summit). Sequences around the peaks on the reverse strand will be reverse complemented.Currently, this tool does not support gzipped fasta files. Default is None

<b>--window-size</b> : windows size describing how many nt will be taken around peak summit to extract the sequence. Must be used together with --genome to have an effect. By default is set to 10, so that 10nt upstream and 10nt downstream from the summit will be extracted (21nt in total). Default is 10

<br>
<br>

#### Algorithm

1) Find peaks using default parameters using `scipy.signal.find_peaks` (for each of the replicates separately) - peaks will represent both actual termination signal and noise
3) Trim peak ends with less than 10 percent of the total peak coverage from the left and right sides. Stop on position with >= 10 percent of the complete peak coverage. The minimum length of the peak returned is 1 (even if the last position also has < 10 percent of the total peak coverage). This results in much shorter peaks which much better correspond to the termination signals in the data (mainly from the intrinsic terminators -> sharp peaks)
4) Run IDR routines for each possible combination of two replicates (`https://github.com/nboley/idr`)
5) Extract the peaks with the IDR below the specified threshold [--idr-threshold] in at least [--min-no-comp]
6) (optional) Calculate termination efficiency for all termination signals. Extracts the average RNA-seq coverage for all replicates starting from 13 nt upstream to 13 nt downstream. The total length of the region will be 27, with the central position corresponding to the peak summit. Peak summit and three consecutive positions, both upstream and downstream, are not considered for the calculations. Termination efficiency is defined as follows:

`T = ( (upstream - downstream)/upstream ) * 100`
  
where

*upstream is the average value calculated from the coverage of 10 nt located upstream from the peak summit
*downstream is the average value calculated from the coverage of 10 nt located downstream from the peak summit
If _upstream_ is smaller than `--min-upstream`, the termination efficiency will be set to ".".
6) (optional) Extract the genomic sequence for the region starting [--window_size] nt upstream from the peak summit up to [--window_size] nt downstream so that the region will have a total length of 2*[--window_size]+1 with the central position corresponding to the peak summit. If the peak is on the reverse strand, the sequence will be reverse complemented

<br>
<br>
                     
#### Output

Please note that all outputs of the TERMITe software use a 0-based coordinate system (as in the BED file format).

The main output file is in the BED narrowPeak format and is named `3prime_RNA_termini_rep1-rep2[...-repN]_[strand].narrowPeak`.
There is an additional FASTA file name with sequences extracted around peak summits (see *Algorithm*).

Fields:

1) chromosome name (or another sequence name, e.g., contig, scaffold, etc.)
2) peak start
3) peak end
4) peak name
5) number of pairwise comparisons with IDR below the specified threshold scaled to the range of <0, 1000> to quickly visualize in the IGV (color density)
			e.g. for three replicates - three pairwise comparisons will be performed (with each possible combination of samples),
   		so that the number of pairwise comparisons with IDR below the specified threshold can be either 1, 2 or 3.
      After scaling this value to the specified range, 0 will correspond to the value of 1, 500 will correspond to the value of 2
      and 1000 will correspond to the value of 3.
6) strand
7) average peak height
8) termination efficiency (percentage) or "." if the termination efficiency wasn't calculated. To calculate termination efficiency <b>--rna-seq-coverage</b> needs to be c.
9) max IDR value among the IDRs calculated for each pairwise comparison of samples
10) peak summit - point source for the peak. 0 corresponds to the most upstream position within the peak, 1 to the next, etc.

<hr>
<br>
<br>
<br>

## Running TERMITe annotate

The goal of the TERMITe annotate is to:
1) annotate identified 3'RNA ends, providing extensive information about their genomic and transcriptomic context, presence of predicted intrinsic terminators in the vicinity of the term-seq signal, termination efficiency at a given termination site, surrounding sequences, and more.
2) create a summarized report for TERMITe runs dedicated to multiple experimental conditions.

<br>
<br>

#### Help

To see the manual and available options, run

	docker run --rm -v $(pwd):/data --platform=linux/amd64 termite annotate --help

<hr><br>

	usage: termite annotate [-h] --termite-out-forward TERMITe-OUTPUT-FWD [TERMITe-OUTPUT-FWD ...] --termite-out-reverse TERMITe-OUTPUT_REV [TERMITe-OUTPUT_REV ...] --sample-names SAMPLE-NAME [SAMPLE-NAME ...] --rna-3prime-ends-forward TERMSEQ-BIGWIG-FWD [TERMSEQ-BIGWIG-FWD ...] --rna-3prime-ends-reverse TERMSEQ-BIGWIG-REV [TERMSEQ-BIGWIG-REV ...] [--rna-seq-coverage-forward RNASEQ-BIGWIG-FWD [RNASEQ-BIGWIG-FWD ...]] [--rna-seq-coverage-reverse RNASEQ-BIGWIG-REV [RNASEQ-BIGWIG-REV ...]] [--intersect BED [BED ...]] [--min-upstream MIN_UPSTREAM] [--max-distance-to-hairpin MAX_DISTANCE_TO_HAIRPIN] [--min-distance-to-hairpin MIN_DISTANCE_TO_HAIRPIN] [--trans-term-hp] [--rnafold] [--gene-annotations GENE_ANNOTATIONS] [--genome GENOME] [--upstream-nt UPSTREAM_NT] [--downstream-nt DOWNSTREAM_NT] --output OUTPUT
			 
##### Required options

<b>--termite-out-forward</b> : list of the TERMITe output files (narrowPeak) calculated for the forward strand

<b>--termite-out-reverse</b> : list of the TERMITe output files (narrowPeak) calculated for the reverse strand. Files belonging to the same samples should be specified in the same order as in the --termite-out-forward

<b>--sample-names</b> : sample names. Should be placed in the same order as corresponding files in --termite-out-forward and --termite-out-reverse. Should be the same as --name-prefix(es) provided in the TERMITe module

<b>--rna-3prime-ends-forward</b> : bigWig files with the genomic coverage by the 3'RNA ends on the forward strand. Should be provided as follows: 'sample:path_to_file1,path_to_file2,...,path_to_fileN'

<b>--rna-3prime-ends-reverse</b> : bigWig files with the genomic coverage by the 3'RNA ends on the reverse strand. Should be provided as follows: 'sample:path_to_file1,path_to_file2,...,path_to_fileN'

<b>--output</b> : a prefix of the output file name

##### Additional options

<b>--rna-seq-coverage-forward</b> : bigWig files with the genomic coverage by the whole RNA-seq reads on the forward strand. Should be provided as follows: 'sample:path_to_file1,path_to_file2,...,path_to_fileN'

<b>--rna-seq-coverage-reverse</b> : bigWig files with the genomic coverage by the whole RNA-seq reads on the reverse strand. Should be provided as follows: 'sample:path_to_file1,path_to_file2,...,path_to_fileN'

<b>--intersect</b> : report 3' RNA end overlaps with the features in the [--intersect] BED file(s). Multiple bed files can be provided using this option. All files should be specified as 'NAME:PATH_TO_THE_BED_FILE', where NAME will be used as the column name in the TERMITe report and therefore must be unique. BED files should be in a standard BED6 format

<b>--min-upstream</b> : do not calculate termination efficiency if the average upstream signal is less than [--min-upstream]. Default is 0.25

<b>--max-distance-to-hairpin</b> : the maximum distance between the POT (Point of Observed RNA Termini) and the upstream terminator hairpin

<b>--min-distance-to-hairpin</b> : the minimum distance between the POT and the intrinsic end of the upstream terminator hairpin. Negative values are possible when the POT overlaps the hairpin

<b>--trans-term-hp</b> : annotate results with the intrinsic terminator hairpins predicted by the TransTermHP software. To run the TransTermHP, you need to specify the --genome option

<b>--rnafold</b> : annotate results with the intrinsic terminator hairpins predicted by the pipeline based on TransTermHP and RNAfold. To use this option, you need to specify the --genome option

<b>--gene-annotations</b> : GTF/GFF file with gene annotations to annotate genomic context of the stable 3' RNA end

<b>--genome</b> : reference genome in the FASTA file format

<b>--upstream-nt</b> : extract the sequence consisting of this many nt upstream from the POT, by default is set to 10, so that [--upstream-nt] nt upstream and [--downstream-nt] nt downstream from the POT will be extracted ([--upstream-nt] + 1 + [--downstream-nt] nt in total). Default is 10

<b>--downstream-nt</b> : extract the sequence consisting of this many nt upstream from the POT, by default is set to 10, so that [--upstream-nt] nt upstream and [--downstream-nt] nt downstream from the POT will be extracted ([--upstream-nt] + 1 + [--downstream-nt] nt in total). Default is 10

<br>
<br>

#### Algorithm

1) Merging peaks from multiple TERMITe's runs representing various experimental conditions using bedtools merge (options -s and -o collapse).
2) TERMITe calculates the reference Point of Observed [RNA] Termini (POT) by averaging signals from all samples and all replicates (in the genomic coordinates of the merged peak and extracting the signal summit (the highest point) within this region. If more than one position with maximal signal within the region has the same average signal value, the one most downstream is selected.
3) If the [--trans-term-hp] option is used, TERMITe runs TransTermHP to predict intrinsic terminators genome-wide [--genome] (present in any genomic context) using default parameters and scoring scheme (see the output fields 17-27 and 34). TERMITe checks:
    - which terminators might be assigned to an identified stable RNA 3'
    - what is the distance between the terminator hairpin and the POT
    - what is the sequence of the predicted U-tract, hairpin, and the possible A-tract (marking the position of the POT with an uppercase letter)
    - what is the confidence assigned to the terminator by the TransTermHP
4) If the [--gene-annotations] option is employed, with help from the `gffutils` package, TERMITe looks for:
    - the potential gene overlapping the peak (on the same strand)
    - the closest upstream gene (on the same strand; reporting the gene id and the distance between the gene and the peak)
    - the closest downstream gene (on the same strand; reporting the gene id and the distance between the gene and the peak)
    - overlapping feature types from the GTF/GFF3 file, like CDS, UTRs, etc. (on the same strand; at least 1bp of the overlap is required)
5) If the [--genome] option is invoked, TERMITe extracts the genomic sequence surrounding the POT (the User can specify the number of nucleotides both downstream and upstream from the POT using `--upstream-nt` and `--downstream-nt` ), marking the POT with an uppercase letter
6) TERMITe calculates the termination efficiency at a given POT if both the `--rna-seq-coverage-forward` and `--rna-seq-coverage-reverse` options are used (see the output section for a detailed description)
7) Optionally, TERMITe looks for the overlaps between the peak and features supplied in the form of one or more BED files. This option can be used to annotate the results, e.g., with known terminators downloaded from the reference databases or established as a result of other experiments (on the same strand; at least 1bp of the overlap is required)
8) To further examine the possibility of forming the intrinsic terminator hairpin upstream to the POT, TERMITe runs TransTermHP with very permissive options (--all-context, --min-conf=0, --min-stem=4, --min-loop=4, --max-len=40, --max-loop=20, --uwin-size=6, --uwin-require=1) to maximally increase the sensitivity of the search. Furthermore, only sequences starting 80 nt upstream and ending 30 nt downstream from the POT are scanned (111 nt in total). The best terminator predicted for a given POT is selected. To select the best hairpin, the following criteria are analyzed (checked in the specified order until only one hairpin remains):
        - hairpin confidence
        - hairpin score
        - tail score.
   Then RNAfold from the ViennaRNA package is run to re-fold the heuristically obtained hairpin structure (from TransTermHP). TERMITe checks if the obtained MFE structure is indeed the hairpin structure with energy <= -3 and if the POT lies at most [--max-distance-to_hairpin] nt downstream from the hairpin. Intrinsic terminators found with these parameters are reported similarly to those from 3) (see the output section 27-33). This option can be used to also examine the possibility of forming I-shaped non-canonical terminator hairpins (e.g. with poor or no U-tract's).

<br>
<br>

#### Output

Please note that all outputs of the TERMITe software use a 0-based coordinate system (as in the BED file format).

The TERMITe annotate module results in a tabular file with multiple columns (tab-separated; .tsv) and a name prefix defined in the command line ([--output]).

Fields:
1) chromosome - chromosome name (or another sequence name, e.g., contig, scaffold, etc.)
2) start - peak start
3) end - peak end
4) POT - the reference Point of Observed [RNA] Termini - calculated by averaging signals from all samples and all replicates (in the genomic coordinates of the merged peak - columns 1,2,3 and 6) and extracting the signal summit (the highest point) within this region. If more than one position with maximal signal within the region has the same average signal value, the one most downstream is selected.
5) sequence - the genomic sequence for the region starting [--upstream-nt] nt upstream from the POT up to [--downstream-nt] nt downstream. The sequence is in lowercase except for the reference POT, which is in uppercase. If the peak is on the reverse strand, the sequence will be reverse complemented - always written from 5' to 3'. The same sequence is saved in the fasta file (all nucleotides in uppercase) - see the explanation below.
6) strand - strand on which the peak is identified
7) termite_id - identifiers of the peaks called by the TERMITe in each run
8) termite_score - number of pairwise comparisons with IDR below the specified threshold; scaled to the range of <0, 1000> to easily visualize in the IGV (color density), e.g., for three replicates - three pairwise comparisons will be performed (with each possible combination of samples) so that the number of pairwise comparisons with IDR below the specified threshold can be either 1, 2 or 3. After scaling this value to the specified range, 0 will correspond to the value of 1, 500 will correspond to the value of 2 and 1000 will correspond to the value of 3. (order as in the termite_id field (7))
9) average_peak_height - average peak heights for each evaluated sample. Even if the peak has not been called for a particular sample - the average peak height is calculated within the coordinates of the merged peak - fields 1,2,3 (order as in the termite_id field (7))
10) termination_efficiency - termination_efficiency (percentage). Even if a peak has not been called for a particular sample - termination efficiency is calculated by establishing the peak summit (as in TERMITe find_stable_rna_ends) for the coordinates of the merged peak - fields 1,2,3 (order as in 7-termite_id). See the TERMITe find_stable_rna_ends documentation on how the termination efficiency is calculated
11) idr - Irreproducible Discovery Rate. Max IDR value among the IDRs calculated for each pairwise comparison of replicates (see TERMITe find_stable_rna_ends). (order as in 7-termite_id)
12) summit_coordinate - genomic coordinate of the peak summit. Even the if a peak has not been called for a particular sample - the peak summit is calculated within the coordinates of the merged peak. (order as in 7-termite_id)
13) overlapping_gene - id of the gene overlapping (at least 1nt-overlap required) the merged peak. "." if none. Gene annotations are taken from the input GTF file 
14) upstream_gene - id and the minimal possible distance between either start or end coordinates of the merged peak and the closest upstream gene on the same strand
15) downstream_gene - id and minimal possible distance between either start or end coordinates of the merged peak and the closest downstream gene on the same strand
16) overlapping_feature_types - list of feature types overlapping the merged peak (as in the GTF file); at least 1-nt overlap is required
17) transtermhp_closest_hairpin_start - the start of the closest intrinsic terminator hairpin predicted by the TransTermHP software
18) transtermhp_closest_hairpin_end - the end of the closest intrinsic terminator hairpin predicted by the TransTermHP software
19) transtermhp_id - id of the closest intrinsic terminator hairpin predicted by the TransTermHP software
20) transtermhp_confidence - confidence of the closest intrinsic terminator hairpin predicted by the TransTermHP software
21) transtermhp_hairpin_score - hairpin score of the closest intrinsic terminator hairpin predicted by the TransTermHP software
22) transtermhp_tail_score - tail score of the closest intrinsic terminator predicted by the TransTermHP software
23) transtermhp_a_tract - sequence of the possible A-tract predicted by the TransTermHP software (always written from 5' to 3' end, complement if on the reverse strand)
24) transtermhp_hairpin - sequence and structure of the hairpin predicted by the TransTermHP software (always written from 5' to 3' end, complement if on the reverse strand)
25) transtermhp_u_tract - sequence of the U-tract predicted by the TransTermHP software (always written from 5' to 3' end, complement if on the reverse strand)
26) transtermhp_distance_to_hairpin -  the minimal possible distance between the merged peak's start or end coordinates and the hairpin predicted by the TransTermHP. If negative - 3' RNA end is upstream from the hairpin end; if positive - downstream
27) rnafold_a_tract - sequence of the possible A-tract predicted by the TransTermHP + RNAfold pipeline (always written from 5' to 3' end, complement if on the reverse strand)
28) rnafold_hairpin - sequence and structure of the hairpin predicted by the TransTermHP + RNAfold pipeline (always written from 5' to 3' end, complement if on the reverse strand)
29) rnafold_u_tract - sequence of the U-tract predicted by the TransTermHP + RNAfold pipeline (always written from 5' to 3' end, complement if on the reverse strand)
transtermhp_closest_hairpin_start - closest intrinsic terminator hairpin start (as in results of the TransTermHP) annotated on the opposite strand
30) rnafold_hairpin_structure - MFE structure of the hairpin presented in 28), folded using the RNAfold software from the ViennaRNA package
transtermhp_closest_hairpin_start - closest intrinsic terminator hairpin start (as in results of the TransTermHP) annotated on the opposite strand
31) rnafold_distance_to_hairpin -  the minimal possible distance between the merged peak's start or end coordinates and the hairpin predicted by the TransTermHP + RNAfold pipeline. If negative - 3' RNA end is upstream from the hairpin end; if positive - downstream
32) rnafold_energy - the energy of the secondary structure in 31)
33) rnafold - "+" if the RNAfold software confirmed the possible intrinsic termination hairpin in the POT's genomic neighborhood defined by both the [--min-distance-to-hairpin] and the [--max-distance-to-hairpin] options; "-" otherwise
34) transtermhp - "+" if the TransTermHP software confirmed the possible intrinsic termination hairpin in the POT's genomic neighborhood defined by both the [--min-distance-to-hairpin] and the [--max-distance-to-hairpin] options; "-" otherwise
35-...) Additional columns describing existing overlaps ("+") with the features from the supplied BED files [--intersect] (at least 1-bp overlap is required, strand-specific search) 

PLEASE NOTE:
- field 5 is available only if the  `--genome` option is used (column will be empty otherwise)
- field 10 is available only if both the  `--rna-seq-coverage-forward` and `--rna-seq-coverage-reverse` options are used (column will be empty otherwise)
- fields 13-16 are available only if the  `--gene-annotations` option is used (columns will be empty otherwise)
- fields 17-26 are available only if the  `--trans-term-hp` option is used and are set to "." if no intrinsic terminator hairpin has been predicted in the POT's genomic neighborhood defined by both the [--min-distance-to-hairpin] and the [--max-distance-to-hairpin] options (columns will be empty otherwise)
- fields 27-32 are available only if the  `--rnafold` option is used and are set to "." if no intrinsic terminator hairpin has been predicted in the POT's genomic neighborhood defined by both the [--min-distance-to-hairpin] and the [--max-distance-to-hairpin] options (columns will be empty otherwise)
- fields 33 and 34 are available only if the  `--rnafold` or `--trans-term-hp` options are used (respectively) (columns will be empty otherwise)


Additionally, TERMITe extracts the genomic sequence for the region starting [--upstream-nt] nt upstream from the POT up to [--downstream-nt] nt downstream,
so that the region will have a total length of upstream+downstream+1 with the central position corresponding to the reference POT (Point of Observed [RNA] Termini).
If the peak is on the reverse strand, the sequence will be reverse complemented - always written from 5' to 3'. Sequences are saved in the fasta file with the name prefix defined from the command line ([--output]).

<hr>
<br>
<br>
<br>

## Test case

All files used in this test case are available at https://github.com/zywicki-lab/TERMITe-dev (`sample_data.tar.gz`).

Extract files and swith working directory
```{bash}
tar -xzvf sample_data.tar.gz
cd sample_data
```

Term-seq, altogether with the paired RNA-seq data, were downloaded from the NCBI SRA database (PRJEB12568).

<b>Organism:</b> Bacillus subtilis subsp. subtilis str. 168
<br>
<b>Description:</b> B. subtilis "grown in LB media at 37c to OD600 = 0.1-0.2; culture supplemented with Lincomycin (Lm) 0.5 g/ml for 15min before harvesting and RNA extraction"
<br>
<b>Samples:</b> ERR1248394 (RNA-seq); ERR1248372 (Term-seq, replicate 1), ERR1248373 (Term-seq, replicate 2), ERR1248374 (Term-seq, replicate 3)
<br>
<b>Publication:</b> `Dar D, Shamir M, Mellin JR, Koutero M, Stern-Ginossar N, Cossart P, Sorek R. Term-seq reveals abundant ribo-regulation of antibiotics resistance in bacteria. Science. 2016 Apr 8;352(6282):aad9822. doi: 10.1126/science.aad9822. PMID: 27120414; PMCID: PMC5756622.`
<br>
<br>



Two bigWig files are present in the sample data for each replicate, representing normalized coverage by either 3' RNA ends or coverage by RNA-seq reads, for Term-seq and RNA-seq samples, respectively. One of them (*.f.bw) stores the coverage on the forward strand, the other (*.r.bw) on the reverse strand. All bigWig files were prepared using 'bamCoverage' tool as described in the `Input` section above.

<br>
<br>
<b>Additional files:</b><br>
`NC_000964.3.fa` - reference genome, Bacillus subtilis subsp. subtilis str. 168 (GCA_000009045, genome assembly ASM904v1) downloaded from the Ensembl database (bacteria.ensembl.org)

`NC_000964.3.gff` - gene annotations, Bacillus subtilis subsp. subtilis str. 168 (GCA_000009045, genome assembly ASM904v1) downloaded from the Ensembl database (bacteria.ensembl.org)

`mandell_et_al_2021_suppl_file_2.bed` - genomic coordinates of intrinsic termination sites identified by Term-seq in the wild-type Bacillus subtilis strain converted to the BED6 file format. Data from `Mandell ZF, Oshiro RT, Yakhnin AV, Vishwakarma R, Kashlev M, Kearns DB, Babitzke P. NusG is an intrinsic transcription termination factor that stimulates motility and coordinates gene expression with NusA. Elife. 2021 Apr 9;10:e61880. doi: 10.7554/eLife.61880. PMID: 33835023; PMCID: PMC8060035`.

Term-seq reads in the FASTQ file format were downloaded from the NCBI SRA database, preprocessed using `trimmomatic` (with the _MINLEN:25_ option, removing TruSeq3-SE adapters; Bolger et al. 2014), and mapped with `bowtie` (-l 20 --sam --best --strata -m 1; Langmead et al. 2009). SAM files were converted to BAMs and sorted using `samtools` (Li et al. 2009; Danecek et al. 2021).

<hr>
<br>
<br>
<br>

To identify stable 3' RNA ends (represented by reproducible bigWig signal across the replicates) using TERMITe, run:

```{bash}
docker run --rm -v $(pwd):/data --platform=linux/amd64 termite find_stable_rna_ends \
	--rna-3prime-ends /data/ERR1248372.f.bw /data/ERR1248373.f.bw /data/ERR1248374.f.bw \
	--rna-seq-coverage /data/ERR1248394.f.bw /data/ERR1248394.f.bw /data/ERR1248394.f.bw \
	--sample-names Lm1 Lm2 Lm3 \
	--out-dir /data/termite_results \
	--strand forward \
	--name-prefix BS_Lm \
	--min-no-comp 3
```
<hr>

where $(pwd) is a location of a working directory, usually containing data for analysis.

```{bash}
docker run --rm -v $(pwd):/data --platform=linux/amd64 termite find_stable_rna_ends \
	--rna-3prime-ends /data/ERR1248372.r.bw /data/ERR1248373.r.bw /data/ERR1248374.r.bw \
	--rna-seq-coverage /data/ERR1248394.r.bw /data/ERR1248394.r.bw /data/ERR1248394.r.bw \
	--sample-names Lm1 Lm2 Lm3 \
	--out-dir /data/termite_results \
	--strand reverse \
	--name-prefix BS_Lm \
	--min-no-comp 3
```
`--min-no-comp 3` will guarantee that the identified 3' RNA ends are reproducible across all the possible pairwise replicate comparisons, resulting in highly-confident peaks.

The first step of the TERMITe pipeline will result in the two main output files - `3prime_RNA_termini_Lm1-Lm2-Lm3_forward.narrowPeak` and `3prime_RNA_termini_Lm1-Lm2-Lm3_reverse.narrowPeak` - both in the BED NarrowPeak file format (see the `Output` section for a detailed description). Those files describe genomic positions of the peaks called on Term-seq data, representing stable 3' RNA ends identified by the TERMITe software.

<hr>
<br>

To annotate the results with the `annotate` module, run:
```{bash}
docker run --rm -v $(pwd):/data --platform=linux/amd64 termite annotate \
	--termite-out-forward /data/termite_results/3prime_RNA_termini_Lm1-Lm2-Lm3_forward.narrowPeak \
	--termite-out-reverse /data/termite_results/3prime_RNA_termini_Lm1-Lm2-Lm3_reverse.narrowPeak \
	--sample-names BS_Lm \
	--rna-3prime-ends-forward BS_Lm:/data/ERR1248372.f.bw,/data/ERR1248373.f.bw,/data/ERR1248374.f.bw \
	--rna-3prime-ends-reverse BS_Lm:/data/ERR1248372.r.bw,/data/ERR1248373.r.bw,/data/ERR1248374.r.bw \
	--rna-seq-coverage-forward BS_Lm:/data/ERR1248394.f.bw,/data/ERR1248394.f.bw,/data/ERR1248394.f.bw \
	--rna-seq-coverage-reverse BS_Lm:/data/ERR1248394.r.bw,/data/ERR1248394.r.bw,/data/ERR1248394.r.bw \
	--gene-annotations /data/NC_000964.3.gff \
	--genome /data/NC_000964.3.fa \
	--output /data/termite_results/Lm \
	--trans-term-hp \
	--rnafold \
	--intersect Mandell_et_al_2021:/data/mandell_et_al_2021_suppl_file_2.bed \
	--min-upstream 5
```
To find out more about those options, please consult the `Running TERMITe annotate` section above. The main output will be the `termite_results/Lm.tsv` file (see the `Output` section for a detailed description). An additional FASTA file containing sequences around the identified Points of Termination will be saved as `termite_results/Lm.fasta`. 
	
<hr>
<br>
<br>
<br>

<b>How do we compare those results to other samples?</b><br>

Let's run TERMITe again on data from the control experiment.

<b>Description:</b> B. subtilis "grown in LB media at 37c to OD600 = 0.1-0.2; culture supplemented with control media for 15min before harvesting and RNA extraction."
<br>
<b>Samples:</b> ERR1248393 (RNA-seq); ERR1248370 (Term-seq, replicate 1), ERR1248371 (Term-seq, replicate 2)
<br>
```{bash}
docker run --rm -v $(pwd):/data --platform=linux/amd64 termite find_stable_rna_ends \
	--rna-3prime-ends /data/ERR1248370.r.bw /data/ERR1248371.r.bw \
	--rna-seq-coverage /data/ERR1248393.r.bw /data/ERR1248393.r.bw \
	--sample-names Ctrl1 Ctrl2 \
	--out-dir /data/termite_results \
	--strand reverse \
	--name-prefix BS_Ctrl \
	--min-no-comp 1
```
<hr>
```{bash}
docker run --rm -v $(pwd):/data --platform=linux/amd64 termite find_stable_rna_ends \
	--rna-3prime-ends /data/ERR1248370.f.bw /data/ERR1248371.f.bw \
	--rna-seq-coverage /data/ERR1248393.f.bw /data/ERR1248393.f.bw \
	--sample-names Ctrl1 Ctrl2 \
	--out-dir /data/termite_results \
	--strand forward \
	--name-prefix BS_Ctrl \
	--min-no-comp 1
```
	
Please note that `--min-no-comp` has been set to 1, as there are only two replicates from the control experiment and, therefore, only one possible pairwise comparison between replicates.

Now, it is possible to summarize the results in a single report using the `annotate` module.

```{bash}
docker run --rm -v $(pwd):/data --platform=linux/amd64 termite annotate \
--termite-out-forward \
	/data/termite_results/3prime_RNA_termini_Lm1-Lm2-Lm3_forward.narrowPeak \
	/data/termite_results/3prime_RNA_termini_Ctrl1-Ctrl2_forward.narrowPeak \
--termite-out-reverse \
	/data/termite_results/3prime_RNA_termini_Lm1-Lm2-Lm3_reverse.narrowPeak \
	/data/termite_results/3prime_RNA_termini_Ctrl1-Ctrl2_reverse.narrowPeak \
--sample-names \
	BS_Lm \
	BS_Ctrl \
--rna-3prime-ends-forward \
	BS_Lm:/data/ERR1248372.f.bw,/data/ERR1248373.f.bw,/data/ERR1248374.f.bw \
	BS_Ctrl:/data/ERR1248370.f.bw,/data/ERR1248371.f.bw \
--rna-3prime-ends-reverse \
	BS_Lm:/data/ERR1248372.r.bw,/data/ERR1248373.r.bw,/data/ERR1248374.r.bw \
	BS_Ctrl:/data/ERR1248370.r.bw,/data/ERR1248371.r.bw \
--rna-seq-coverage-forward \
	BS_Lm:/data/ERR1248394.f.bw,/data/ERR1248394.f.bw,/data/ERR1248394.f.bw \
	BS_Ctrl:/data/ERR1248393.f.bw,/data/ERR1248393.f.bw \
--rna-seq-coverage-reverse \
	BS_Lm:/data/ERR1248394.r.bw,/data/ERR1248394.r.bw,/data/ERR1248394.r.bw \
	BS_Ctrl:/data/ERR1248393.r.bw,/data/ERR1248393.r.bw \
--gene-annotations /data/NC_000964.3.gff \
--genome /data/NC_000964.3.fa \
--output /data/termite_results/summary \
--trans-term-hp --rnafold \
--intersect Mandell_et_al_2021:/data/mandell_et_al_2021_suppl_file_2.bed \
--min-upstream 5
```

The main results are stored in `termite_results/summary.tsv` and `termite_results/summary.fasta` files (see the `Output` section for a detailed description).
	



	

<hr>
<br>
<br>
<br>

## Literature

To check the details about the results of the TransTermHP, please read:

`Carl Kingsford, Kunmi Ayanbule, and Steven Salzberg. Rapid, accurate, computational discovery of Rho-independent transcription terminators illuminates their relationship to DNA uptake. Genome Biology 8:R22 (2007)`; http://transterm.ccb.jhu.edu

With particular emphasis on the sections: `Algorithm to search for candidate terminators` and `Function to evaluate the quality of terminators`.

<hr>

To get acquainted with other software employed in the TERMITe's pipeline, please check:

`Cock PJ, Antao T, Chang JT, Chapman BA, Cox CJ, Dalke A, Friedberg I, Hamelryck T, Kauff F, Wilczynski B, de Hoon MJ. Biopython: freely available Python tools for computational molecular biology and bioinformatics. Bioinformatics. 2009 Jun 1;25(11):1422-3. doi: 10.1093/bioinformatics/btp163. Epub 2009 Mar 20. PMID: 19304878; PMCID: PMC2682512.`

`Dale RK, Pedersen BS, Quinlan AR. Pybedtools: a flexible Python library for manipulating genomic datasets and annotations. Bioinformatics. 2011 Dec 15;27(24):3423-4. doi: 10.1093/bioinformatics/btr539. Epub 2011 Sep 23. PMID: 21949271; PMCID: PMC3232365.`

`Li Q, Brown JB, Huang H, Bickel PJ. Measuring reproducibility of high-throughput experiments, Annals of Applied Statistics. 2011 Vol. 5, No. 3, 1752-1779. doi: 10.1214/11-AOAS466`; https://github.com/nboley/idr

`Lorenz R, Bernhart SH, Höner Zu Siederdissen C, Tafer H, Flamm C, Stadler PF, Hofacker IL. ViennaRNA Package 2.0. Algorithms Mol Biol. 2011 Nov 24;6:26. doi: 10.1186/1748-7188-6-26. PMID: 22115189; PMCID: PMC3319429.`

`Quinlan AR, Hall IM. BEDTools: a flexible suite of utilities for comparing genomic features. Bioinformatics. 2010 Mar 15;26(6):841-2. doi: 10.1093/bioinformatics/btq033. Epub 2010 Jan 28. PMID: 20110278; PMCID: PMC2832824.`

`Quinlan AR. BEDTools: The Swiss-Army Tool for Genome Feature Analysis. Curr Protoc Bioinformatics. 2014 Sep 8;47:11.12.1-34. doi: 10.1002/0471250953.bi1112s47. PMID: 25199790; PMCID: PMC4213956.`

http://daler.github.io/gffutils/

https://github.com/deeptools/pyBigWig

<hr>

Software used in the TERMITe's example use-case:

`Bolger AM, Lohse M, Usadel B. Trimmomatic: a flexible trimmer for Illumina sequence data. Bioinformatics. 2014 Aug 1;30(15):2114-20. doi: 10.1093/bioinformatics/btu170. Epub 2014 Apr 1. PMID: 24695404; PMCID: PMC4103590.`

`Danecek P, Bonfield JK, Liddle J, Marshall J, Ohan V, Pollard MO, Whitwham A, Keane T, McCarthy SA, Davies RM, Li H. Twelve years of SAMtools and BCFtools. Gigascience. 2021 Feb 16;10(2):giab008. doi: 10.1093/gigascience/giab008. PMID: 33590861; PMCID: PMC7931819.`

`Langmead B, Trapnell C, Pop M, Salzberg SL. Ultrafast and memory-efficient alignment of short DNA sequences to the human genome. Genome Biol. 2009;10(3):R25. doi: 10.1186/gb-2009-10-3-r25. Epub 2009 Mar 4. PMID: 19261174; PMCID: PMC2690996.`

`Li H, Handsaker B, Wysoker A, Fennell T, Ruan J, Homer N, Marth G, Abecasis G, Durbin R; 1000 Genome Project Data Processing Subgroup. The Sequence Alignment/Map format and SAMtools. Bioinformatics. 2009 Aug 15;25(16):2078-9. doi: 10.1093/bioinformatics/btp352. Epub 2009.`



<hr>
<br>
<br>
<br>

## Contribute

If you notice any errors and mistakes or would like to suggest some new features, please use Github's issue tracking system to report them. You are also welcome to send a pull request with your corrections and suggestions.

## License

This project is licensed under the GNU General Public License v3.0 license terms.

## Funding

This work was supported by the Faculty of Biology Dean’s Grant (Adam Mickiewicz University in Poznan, Poland) `[GDWB-11/2020 to J.K.]`

