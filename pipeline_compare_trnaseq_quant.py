"""
====================
Compare tRNA Seq quantification methods pipeline
====================


Overview
========
This pipeline compares alternative approaches for quantification from tRNA-Seq
reads using simulated data


Input
-----
Reads are imported by placing files or linking to files in the :term:
`working directory`.
The following suffixes/file types are possible:

fastq.gz
   Single-end reads in fastq format.
fastq.1.gz, fastq.2.gzf
   Paired-end reads in fastq format.
   The two fastq files must be sorted by read-pair.
Code
====
"""

###################################################
# load modules
###################################################

# import ruffus
from ruffus import transform, suffix, regex, merge, \
    follows, mkdir, originate, add_inputs, jobs_limit, split, \
    subdivide, formatter, collate


# import ruffus
from ruffus import *

# import useful standard python modules
import sys
import os

import pandas as pd

import pickle
# import re
# import shutil
# import sqlite3
# import glob

# import cgatcore modules
#import cgatcore.experiment as E
from cgatcore import pipeline as P
import cgatcore.iotools as iotools

# import module for RNASeq functions. Allows absraction of the pipeline task
# and re-use of core functions
# import ModuleRNASeq as RNASeq

import simulatetrna.alignmentSummary as alignmentSummary
import simulatetrna.simulateReads as simulateReads
import simulatetrna.fasta as fasta
import simulatetrna.bam as bam

import PipelineCompareTrnaSeq as CompareTrnaSeq
# load options from the config file
PARAMS = P.get_parameters(
    ["%s/pipeline.yml" % os.path.splitext(__file__)[0],
    "pipeline.yml"])

conda_base_env = PARAMS['conda_base_env']

# Helper functions mapping tracks to conditions, etc
# determine the location of the input files (reads).
try:
    PARAMS["input"]
except KeyError:
    PARAMS["input"] = "."

# define input files. Here we allow single or paired end
SEQUENCESUFFIXES = ("*.fastq.1.gz", "*.fastq.gz", "*.fastq")

SEQUENCEFILES = tuple([os.path.join(PARAMS["input"], suffix_name)
                       for suffix_name in SEQUENCESUFFIXES])

SEQUENCEFILES_REGEX = regex(r"(.*\/)*(\S+).(fastq.1.gz|fastq.gz|fastq)")


###################################################
# Prepare inputs
###################################################


@mkdir('trna_sequences.dir')
@originate('trna_sequences.dir/hg38-tRNAs.fa')
def downloadtRNAsequences(outfile):

    url = 'http://gtrnadb.ucsc.edu/genomes/eukaryota/Hsapi38/hg38-tRNAs.fa'

    statement = '''
    wget %(url)s
    -O %(outfile)s
    ''' % locals()
    # Very small task so running locally is OK, even on HPC!
    P.run(statement, submit=False)


@transform(downloadtRNAsequences,
           regex('trna_sequences.dir/(\S+).fa'),
           r'trna_sequences.dir/\1.sanitised.fa')
def santitisetRNAsequences(infile, outfile):
    fasta.santitiseGtRNAdbFasta(infile, outfile)


@transform(santitisetRNAsequences,
           regex('trna_sequences.dir/(\S+).sanitised.fa'),
           r'trna_sequences.dir/\1.sanitised.fa.bwt')
def buildBWAIndex(infile, outfile):
    'Index the tRNA sequences for BWA'

    statement = 'bwa index %(infile)s' % locals()

    # Very small task so running locally is OK, even on HPC!
    P.run(statement, submit=False)


@transform(santitisetRNAsequences,
           regex('trna_sequences.dir/(\S+).sanitised.fa'),
           r'trna_sequences.dir/\1.sanitised.fa.1.bt2')
def buildBowtie2Index(infile, outfile):
    'Index the tRNA sequences for Bowtie2'

    outfile_base = iotools.snip(outfile, '.1.bt2')
    statement = 'bowtie2-build %(infile)s %(outfile_base)s' % locals()

    # Very small task so running locally is OK, even on HPC!
    P.run(statement, submit=False)


@transform(santitisetRNAsequences,
           regex('trna_sequences.dir/(\S+).sanitised.fa'),
           r'trna_sequences.dir/\1.sanitised.fa.idx')
def buildSegemehlIndex(infile, outfile):
    'Index the tRNA sequences for Segemehl'

    statement = 'segemehl.x -x %(outfile)s -d %(infile)s' % locals()

    # Very small task so running locally is OK, even on HPC!
    P.run(statement, submit=False)
###################################################
# Learn mutation and truncation signatures
###################################################

print(SEQUENCEFILES)
@mkdir('mut_trunc_sig.dir')
@transform(SEQUENCEFILES,
           SEQUENCEFILES_REGEX,
           add_inputs(buildBWAIndex),
           r"mut_trunc_sig.dir/\2_bwaMemReportSingle.bam")
def mapWithBWAMEMSingleReport(infiles, outfile):
    '''
    Map reads to the tRNA sequences, reporting only one alignment per read
    '''
    infile, bwa_index = infiles

    bwa_index = iotools.snip(bwa_index, '.bwt')

    tmp_file = P.get_temp_filename()

    statement = '''
    bwa mem -k 10 -T 15
    %(bwa_index)s %(infile)s
    > %(tmp_file)s 2> %(outfile)s.log;
    samtools flagstat %(tmp_file)s > %(outfile)s.flagstat;
    samtools sort -o %(outfile)s -O BAM %(tmp_file)s;
    samtools index %(outfile)s;
    rm -f %(tmp_file)s
    ''' % locals()

    P.run(statement)


@mkdir('mut_trunc_sig.dir')
@transform(SEQUENCEFILES,
           SEQUENCEFILES_REGEX,
           add_inputs(buildBowtie2Index),
           r"mut_trunc_sig.dir/\2_bowtie2ReportSingle.bam")
def mapWithBowtie2SingleReport(infiles, outfile):
    '''
    Map reads to the tRNA sequences, reporting only one alignment per read
    '''
    infile, bowtie2_index = infiles

    bowtie2_index = iotools.snip(bowtie2_index, '.1.bt2')

    tmp_file = P.get_temp_filename()

    statement = '''
    bowtie2
    --min-score G,1,8
    --local
    -D 20 -R 3 -N 1 -L 10
    -i S,1,0.5
    -x %(bowtie2_index)s
    -U %(infile)s
    -S %(tmp_file)s
    2> %(outfile)s.log;
    samtools flagstat %(tmp_file)s > %(outfile)s.flagstat;
    samtools sort -o %(outfile)s -O BAM %(tmp_file)s;
    samtools index %(outfile)s;
    rm -f %(tmp_file)s
    ''' % locals()

    P.run(statement)

@mkdir('mut_trunc_sig.dir')
@transform(SEQUENCEFILES,
           SEQUENCEFILES_REGEX,
           add_inputs(buildSegemehlIndex, santitisetRNAsequences),
           r"mut_trunc_sig.dir/\2_SegemehlReportSingle.bam")
def mapWithSegemehlSingleReport(infiles, outfile):
    '''
    Map reads to the tRNA sequences, reporting only one alignment per read
    '''

    infile, segemehl_index, trna_fasta = infiles

    tmp_file = P.get_temp_filename()

    statement = '''
    segemehl.x
    -A 85 -E 500 -r 1
    -i %(segemehl_index)s
    -d %(trna_fasta)s
    -q %(infile)s
    > %(tmp_file)s
    2> %(outfile)s.log;
    samtools flagstat %(tmp_file)s > %(outfile)s.flagstat;
    samtools sort -o %(outfile)s -O BAM %(tmp_file)s;
    samtools index %(outfile)s;
    rm -f %(tmp_file)s
    ''' % locals()

    P.run(statement)

# @collate((mapWithBWAMEMSingleReport,
#           mapWithBowtie2SingleReport,
#           mapWithSegemehlSingleReport),
@collate((mapWithBWAMEMSingleReport,
          mapWithBowtie2SingleReport),
       regex('mut_trunc_sig.dir/(\S+)_(bowtie2ReportSingle|bwaMemReportSingle|SegemehlReportSingle).bam'),
       r'mut_trunc_sig.dir/\1.merged.bam')
def mergeSingleReports(infiles, outfile):

    infiles = ' '.join(infiles)

    statement = '''
    samtools merge -f -o %(outfile)s %(infiles)s;
    samtools index %(outfile)s
    ''' % locals()

    P.run(statement)

# @transform((mapWithBWAMEMSingleReport,
#             mapWithBowtie2SingleReport,
#             mapWithSegemehlSingleReport),
@transform((mapWithBWAMEMSingleReport,
            mapWithBowtie2SingleReport),
           suffix('.bam'),
           add_inputs(santitisetRNAsequences),
           '.summariseAlignments.pickle')
def summariseIndividualAlignments(infiles, outfile):
    '''
    Use alignmentSummary.clustalwtrnaAlignmentSummary class to summarise the
    alignments and learn the mutational & truncation signatures
    '''

    infile, trna_fasta = infiles

    CompareTrnaSeq.summariseAlignments(infile, trna_fasta, outfile, submit=True)

@transform(mergeSingleReports,
           suffix('.bam'),
           add_inputs(santitisetRNAsequences),
           '.summariseAlignments.pickle')
def summariseMergedAlignments(infiles, outfile):
    '''
    Use alignmentSummary.clustalwtrnaAlignmentSummary class to summarise the
    merged alignments and learn the mutational & truncation signatures
    '''

    infile, trna_fasta = infiles

    CompareTrnaSeq.summariseAlignments(infile, trna_fasta, outfile, submit=True)


###################################################
# simulate reads
###################################################
@mkdir('simulation_dummy_files')
@originate(['simulation_dummy_files/%s' % x for x in range(0, PARAMS['simulation_n'])])
def create_files(output_file):
    with open(output_file, "w"):
        pass

@mkdir('simulations.dir')
@transform(santitisetRNAsequences,
           regex('trna_sequences.dir/(\S+).sanitised.fa'),
           r'simulations.dir/\1.simulation_uniform.fastq.gz')
def simulation_uniform_no_errors(infile, outfile):

    n_reads = int(PARAMS['simulation_uniform_reads'])

    gt = simulateReads.make_gt(infile, 1, 0)

    adjustment_factor = n_reads / sum(gt.values())

    gt = {x:round(gt[x] * adjustment_factor) for x in gt.keys()}

    with open(outfile + '.gt', 'w') as out:
        for k,v in gt.items():
            out.write('%s\t%s\n' % (k, v))

    simulateReads.simulate_reads(
        infile,
        outfile,
        gt,
        error_rate=0,
        mutation_threshold=None,
        truncate=False,
        alignment_summary=None,
        summary_level='anticodon')


#### HERE!! Needs update so simulate from each input file, e.g not transform but combination of dummy file and summarisedMergedAlignments ###
@mkdir('simulations.dir')
@transform(create_files,
           regex('simulation_dummy_files/(\d+)'),
           add_inputs(summariseMergedAlignments, santitisetRNAsequences),
           r'simulations.dir/\1.simulation_withseqer.fastq.gz')
def simulation_with_sequencing_errors(infiles, outfile):

    simulation_n, alignment_summary_picklefile, infile = infiles

    n_reads = int(PARAMS['simulation_withmut_reads'])

    gt = simulateReads.make_gt(infile, 1000, 500)

    adjustment_factor = n_reads / sum(gt.values())

    gt = {x:round(gt[x] * adjustment_factor) for x in gt.keys()}

    with open(outfile + '.gt', 'w') as out:
        for k,v in gt.items():
            out.write('%s\t%s\n' % (k, v))

    alignment_summary=pickle.load(open(alignment_summary_picklefile, 'rb'))

    simulateReads.simulate_reads(
        infile,
        outfile,
        gt,
        error_rate=1/100,
        mutation_threshold=None,
        truncate=False,
        alignment_summary=None,
        summary_level='anticodon')


@mkdir('simulations.dir')
@transform(create_files,
           regex('simulation_dummy_files/(\d+)'),
           add_inputs(summariseMergedAlignments, santitisetRNAsequences),
           r'simulations.dir/\1.simulation_withmut.fastq.gz')
def simulation_with_mutation_errors(infiles, outfile):

    simulation_n, alignment_summary_picklefile, infile = infiles

    n_reads = int(PARAMS['simulation_withmut_reads'])

    gt = simulateReads.make_gt(infile, 1000, 500)

    adjustment_factor = n_reads / sum(gt.values())

    gt = {x:round(gt[x] * adjustment_factor) for x in gt.keys()}

    with open(outfile + '.gt', 'w') as out:
        for k,v in gt.items():
            out.write('%s\t%s\n' % (k, v))

    alignment_summary=pickle.load(open(alignment_summary_picklefile, 'rb'))

    simulateReads.simulate_reads(
        infile,
        outfile,
        gt,
        error_rate=1/100,
        mutation_threshold=0.1,
        truncate=False,
        alignment_summary=alignment_summary,
        summary_level='anticodon')

@mkdir('simulations.dir')
@transform(create_files,
           regex('simulation_dummy_files/(\d+)'),
           add_inputs(summariseMergedAlignments, santitisetRNAsequences),
           r'simulations.dir/\1.simulation_withtrunc.fastq.gz')
def simulation_with_truncations(infiles, outfile):

    simulation_n, alignment_summary_picklefile, infile = infiles

    n_reads = int(PARAMS['simulation_withmut_reads'])

    gt = simulateReads.make_gt(infile, 1000, 500)

    adjustment_factor = n_reads / sum(gt.values())

    gt = {x:round(gt[x] * adjustment_factor) for x in gt.keys()}

    with open(outfile + '.gt', 'w') as out:
        for k,v in gt.items():
            out.write('%s\t%s\n' % (k, v))

    alignment_summary=pickle.load(open(alignment_summary_picklefile, 'rb'))

    simulateReads.simulate_reads(
        infile,
        outfile,
        gt,
        error_rate=1/100,
        mutation_threshold=0.1,
        truncate=True,
        alignment_summary=alignment_summary,
        summary_level='anticodon')


###################################################
# quantify
###################################################
@mkdir('quant.dir')
@transform((simulation_uniform_no_errors,
            simulation_with_sequencing_errors,
            simulation_with_mutation_errors,
            simulation_with_truncations),
# @transform(simulation_uniform_no_errors,
           regex('simulations.dir/(\S+).(simulation_\S+).fastq.gz'),
           add_inputs(buildBWAIndex),
           r'quant.dir/\1.\2.bwa.bam')
def alignWithBWA(infiles, outfile):

    infile, bwa_index= infiles

    bwa_index = iotools.snip(bwa_index, '.bwt')

    tmp_file = P.get_temp_filename()

    statement = '''
    bwa mem
    -k 10
    -T 15
    -a
    %(bwa_index)s
    %(infile)s
    > %(tmp_file)s
    2> %(outfile)s.log;
    samtools flagstat %(tmp_file)s > %(outfile)s.flagstat;
    '''

    P.run(statement)

    bam.filter_sam(tmp_file, outfile)

    os.unlink(tmp_file)


@mkdir('quant.dir')
@transform((simulation_uniform_no_errors,
            simulation_with_sequencing_errors,
            simulation_with_mutation_errors,
            simulation_with_truncations),
# @transform(simulation_uniform_no_errors,
           regex('simulations.dir/(\S+).(simulation_\S+).fastq.gz'),
           add_inputs(buildBowtie2Index),
           r'quant.dir/\1.\2.bowtie2.bam')
def alignWithBowtie2(infiles, outfile):

    infile, bowtie2_index = infiles

    bowtie2_index = iotools.snip(bowtie2_index, '.1.bt2')

    tmp_file = P.get_temp_filename()

    # -a = report all alignments
    statement = '''
    bowtie2
    --min-score G,1,8
    --local
    -a
    -D 20 -R 3 -N 1 -L 10
    -i S,1,0.5
    -x %(bowtie2_index)s
    -U %(infile)s
    -S %(tmp_file)s
    2> %(outfile)s.log;
    samtools flagstat %(tmp_file)s > %(outfile)s.flagstat;
    ''' % locals()

    P.run(statement)

    bam.filter_sam(tmp_file, outfile)

    os.unlink(tmp_file)


@transform((alignWithBWA,
            alignWithBowtie2),
           regex('quant.dir/(\S+).(simulation_\S+)\.(\S+).bam'),
           [r'quant.dir/\1.\2.\3.gene_count_individual',
           r'quant.dir/\1.\2.\3.gene_count_isodecoder',
           r'quant.dir/\1.\2.\3.gene_count_anticodon'])
def quantDiscreteCounts(infile, outfiles):

    outfile_individual, outfile_isodecoder, outfile_anticodon = outfiles

    CompareTrnaSeq.tally_read_counts(
    infile,
    outfile_individual,
    outfile_isodecoder,
    outfile_anticodon,
    submit=True)

@transform((alignWithBWA,
            alignWithBowtie2),
           regex('quant.dir/(\S+).(simulation_\S+)\.(\S+).bam'),
           add_inputs(santitisetRNAsequences),
           r'quant.dir/\1.\2.\3/quant.sf')
def quantWithSalmon(infiles, outfile):

    infile, trna_fasta = infiles

    outfile_dir = os.path.dirname(outfile)

    if not os.path.exists(outfile_dir):
        os.makedirs(outfile_dir)

    statement = '''
    salmon quant -t
    %(trna_fasta)s
    -l A
    -a %(infile)s
    -o %(outfile_dir)s;
    '''

    P.run(statement)


###################################################
# compare
###################################################
@collate(quantWithSalmon,
       regex(r'quant.dir/(\S+).(simulation_\S+)\.(\S+)/quant.sf'),
       add_inputs(r'simulations.dir/\1.\2.fastq.gz.gt'),
       [r'quant.dir/\2.SalmonCompareTruthEstimate.tsv',
        r'quant.dir/\2.SalmonCompareTruthEstimateIsodecoder.tsv',
        r'quant.dir/\2.SalmonCompareTruthEstimateAnticodon.tsv'])
def compareTruthEstimateSalmon(infiles, outfiles):

    outfile_individual, outfile_isodecoder, outfile_anticodon = outfiles

    all_counts_vs_truth = []

    for infile_pair in infiles:

        estimate, truth = infile_pair

        estimate_counts = pd.read_csv(estimate, sep='\t')[['Name', 'NumReads']]
        truth_counts = pd.read_csv(truth, sep='\t', header=None,
                                   names=['Name', 'truth'])

        counts_vs_truth = pd.merge(estimate_counts, truth_counts, on='Name')

        simulation_n = os.path.basename(truth).split('.')[0]
        counts_vs_truth['simulation_n']=simulation_n

        quant_method = os.path.dirname(estimate).split('.')[-1] + '_salmon'
        counts_vs_truth['quant_method']=quant_method

        all_counts_vs_truth.append(counts_vs_truth)

    all_counts_vs_truth = pd.concat(all_counts_vs_truth)
    all_counts_vs_truth.to_csv(outfile_individual, sep='\t', index=False)

    all_counts_vs_truth['Name'] = ['-'.join(x.split('-')[0:4]) for x in all_counts_vs_truth.Name]
    all_counts_vs_truth_isodecoder = all_counts_vs_truth.groupby(
        ['Name', 'simulation_n', 'quant_method']).agg(
            {'NumReads':'sum', 'truth':'sum'}).reset_index()
    all_counts_vs_truth_isodecoder.to_csv(outfile_isodecoder, sep='\t', index=False)

    all_counts_vs_truth['Name'] = ['-'.join(x.split('-')[0:3]) for x in all_counts_vs_truth.Name]
    all_counts_vs_truth_ac = all_counts_vs_truth.groupby(
        ['Name', 'simulation_n', 'quant_method']).agg(
            {'NumReads':'sum', 'truth':'sum'}).reset_index()
    all_counts_vs_truth_ac.to_csv(outfile_anticodon, sep='\t', index=False)





@collate(quantDiscreteCounts,
         regex('quant.dir/(\S+).(simulation_\S+)\.(\S+).gene_count_.*'),
         add_inputs(r'simulations.dir/\1.\2.fastq.gz.gt'),
         [r'quant.dir/\2.DiscreteCompareTruthEstimate.tsv',
          r'quant.dir/\2.DiscreteCompareTruthEstimateIsodecoder.tsv',
          r'quant.dir/\2.DiscreteCompareTruthEstimateAnticodon.tsv'])
def compareTruthEstimateDiscreteCounts(infiles, outfiles):

    outfile_individual, outfile_isodecoder, outfile_anticodon = outfiles

    all_counts_vs_truth = {'individual':[], 'isodecoder':[], 'anticodon':[]}

    for infile_set in infiles:

        estimates, truth = infile_set

        estimate_individual, estimate_isodecoder, estimate_anticodon = estimates

        simulation_n = os.path.basename(truth).split('.')[0]
        quant_method = os.path.basename(estimate_individual).split('.')[-2] + '_discrete'

        truth_counts = pd.read_csv(truth, sep='\t', header=None,
                                   names=['Name', 'truth'])

        counts_vs_truth = CompareTrnaSeq.mergeDiscreteEstimateWithTruth(
            truth_counts, estimate_individual, simulation_n, quant_method)
        all_counts_vs_truth['individual'].append(counts_vs_truth)

        truth_counts['Name'] = ['-'.join(x.split('-')[0:4]) for x in truth_counts.Name]
        truth_counts_isodecoder = truth_counts.groupby(
            ['Name']).agg(
                {'truth':'sum'}).reset_index()
        counts_vs_truth_isodecoder = CompareTrnaSeq.mergeDiscreteEstimateWithTruth(
            truth_counts_isodecoder, estimate_isodecoder, simulation_n, quant_method)
        all_counts_vs_truth['isodecoder'].append(counts_vs_truth_isodecoder)

        truth_counts['Name'] = ['-'.join(x.split('-')[0:3]) for x in truth_counts.Name]
        truth_counts_anticodon = truth_counts.groupby(
            ['Name']).agg(
                {'truth':'sum'}).reset_index()
        counts_vs_truth_anticodon = CompareTrnaSeq.mergeDiscreteEstimateWithTruth(
            truth_counts_anticodon, estimate_anticodon, simulation_n, quant_method)
        all_counts_vs_truth['anticodon'].append(counts_vs_truth_anticodon)


    all_counts_vs_truth_individual = pd.concat(all_counts_vs_truth['individual'])
    all_counts_vs_truth_individual.to_csv(outfile_individual, sep='\t', index=False)

    all_counts_vs_truth_isodecoder = pd.concat(all_counts_vs_truth['isodecoder'])
    all_counts_vs_truth_isodecoder.to_csv(outfile_isodecoder, sep='\t', index=False)

    all_counts_vs_truth_anticodon = pd.concat(all_counts_vs_truth['anticodon'])
    all_counts_vs_truth_anticodon.to_csv(outfile_anticodon, sep='\t', index=False)

# Function to merge comparisons for discrete and salmon-based quant

# Plotting functions
###################################################
# Targets
###################################################

@follows(buildBWAIndex)
def build_input():
    'prepare inputs'
    pass


@follows(summariseMergedAlignments)
def summarise_alignments():
    'learn the mutation & truncation signatures'
    pass


@follows(simulation_uniform_no_errors,
         simulation_with_sequencing_errors,
         simulation_with_mutation_errors,
         simulation_with_truncations)
def simulate():
    'simulate tRNA reads from basic to more realistic'
    pass


@follows(summariseIndividualAlignments,
         summariseMergedAlignments)
def summariseAlignments():
    'summarise the alignments with real reads'
    pass


@follows(quantWithSalmon)
def quant():
    'Quantify from the simulated reads'
    pass

@follows(compareTruthEstimateSalmon)
def compare():
    'compare observed and ground truth'
    pass

@follows(summariseAlignments,
         quant,
         compare)
def full():
    'run it all!'
    pass


###################################################
# Making pipline command-line friendly
###################################################

# Facilitate command line parsing
def main(argv=None):
    if argv is None:
        argv = sys.argv
    P.main(argv)


if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
