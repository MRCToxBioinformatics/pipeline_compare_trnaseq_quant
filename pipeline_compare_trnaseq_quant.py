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




###################################################
# Learn mutation and truncation signatures
###################################################

print(SEQUENCEFILES)
@mkdir('mut_trunc_sig.dir')
@transform(SEQUENCEFILES,
           SEQUENCEFILES_REGEX,
           add_inputs(buildBWAIndex),
           r"mut_trunc_sig.dir/\2.bam")
def mapWithBWAMEMSingleReport(infiles, outfile):
    '''
    Map reads to the tRNA sequences, reporting only one alignment per read
    '''
    infile, bwa_index = infiles

    bwa_index = iotools.snip(bwa_index, '.bwt')

    tmp_file = P.get_temp_filename()

    statement = '''
    bwa mem %(bwa_index)s %(infile)s > %(tmp_file)s 2> %(outfile)s.log;
    samtools sort -o %(outfile)s -O BAM %(tmp_file)s;
    samtools index %(outfile)s;
    rm -f %(tmp_file)s
    ''' % locals()

    P.run(statement)


@transform(mapWithBWAMEMSingleReport,
           suffix('.bam'),
           add_inputs(santitisetRNAsequences),
           '.summariseAlignments.pickle')
def summariseAlignments(infiles, outfile):
    '''
    Use alignmentSummary.clustalwtrnaAlignmentSummary class to summarise the
    alignments and learn the mutational & truncation signatures
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


@mkdir('simulations.dir')
@transform(create_files,
           regex('simulation_dummy_files/(\d+)'),
           add_inputs(summariseAlignments, santitisetRNAsequences),
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
           add_inputs(summariseAlignments, santitisetRNAsequences),
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
           add_inputs(summariseAlignments, santitisetRNAsequences),
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
           regex('simulations.dir/(\S+).(simulation_\S+).fastq.gz'),
           add_inputs(buildBWAIndex),
           r'quant.dir/\1.\2.')
def quantWithBWASalmon(infiles, outfile):

    infile, bwa_index = infiles

    bwa_index = iotools.snip(bwa_index, '.bwt')


    statement ='bwa mem -a %(bwa_index)s %(infile)s > %(tmp_file)s 2> %(outfile)s.log;'

    P.run(statement)

    

@mkdir('mut_trunc_sig.dir')
@transform(SEQUENCEFILES,
           SEQUENCEFILES_REGEX,
           add_inputs(buildBWAIndex),
           r"mut_trunc_sig.dir/\2.bam")
def mapWithBWAMEMSingleReport(infiles, outfile):
    '''
    Map reads to the tRNA sequences, reporting only one alignment per read
    '''
    infile, bwa_index = infiles

    bwa_index = iotools.snip(bwa_index, '.bwt')

    tmp_file = P.get_temp_filename()

    statement = '''
    bwa mem %(bwa_index)s %(infile)s > %(tmp_file)s 2> %(outfile)s.log;
    samtools sort -o %(outfile)s -O BAM %(tmp_file)s;
    samtools index %(outfile)s;
    rm -f %(tmp_file)s
    ''' % locals()

    P.run(statement)

###################################################
# Meta-task targets
###################################################

# build_input =  prepare inputs
@follows(buildBWAIndex)
def build_input():
    pass

# summarise_alignments = learn the mutation & truncation signatures
@follows(summariseAlignments)
def summarise_alignments():
    pass

# simulate = simulate tRNA reads from basic to more realistic
@follows(simulation_uniform_no_errors,
         simulation_with_sequencing_errors,
         simulation_with_mutation_errors,
         simulation_with_truncations)
def simulate():
    pass

simulation_uniform_no_errors

# full = run it all!
@follows(summariseAlignments)
def full():
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
