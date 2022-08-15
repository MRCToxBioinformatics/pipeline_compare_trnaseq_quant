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

from ruffus.combinatorics import product

# import ruffus
from ruffus import *

# import useful standard python modules
import sys
import os
import re
import pickle
import shutil
import glob

import pandas as pd
import pysam

from cgatcore import pipeline as P
import cgatcore.iotools as iotools

import simulatetrna.alignmentSummary as alignmentSummary
import simulatetrna.simulateReads as simulateReads
import simulatetrna.fasta as fasta
import simulatetrna.bam as bam

import PipelineCompareTrnaSeq as CompareTrnaSeq
# load options from the config file
PARAMS = P.get_parameters(
    ["%s/pipeline.yml" % os.path.splitext(__file__)[0],
    "pipeline.yml"])

# Helper functions mapping tracks to conditions, etc
# determine the location of the input files (reads).
try:
    PARAMS["input"]
except KeyError:
    PARAMS["input"] = "."

print(PARAMS)
for k,v in PARAMS.items():
    print(k,v)

# define input files. Here we allow single or paired end
SEQUENCESUFFIXES = ("*.fastq.1.gz", "*.fastq.gz", "*.fastq")

SEQUENCEFILES = tuple([os.path.join(PARAMS["input"], suffix_name)
                       for suffix_name in SEQUENCESUFFIXES])

SEQUENCEFILES_REGEX = regex(r"(.*\/)*(\S+).(fastq.1.gz|fastq.gz|fastq)")

###################################################
# Prepare inputs
###################################################

@mkdir('trna_sequences.dir')
@originate('trna_sequences.dir/mt_sequences.fa')
def updateMtFastaNaming(outfile):
    CompareTrnaSeq.updateMtFastaNaming(PARAMS['trna_mt_sequences_infile'], outfile)

@merge((PARAMS['trna_sequences_infile'], updateMtFastaNaming),
       'trna_sequences.dir/trna_sequences_all.fa')
def mergeNucMtSequences(infiles, outfile):

    with open(outfile, 'w') as outf:
        for infile in infiles:
            with open(infile, 'r') as inf:
                for line in inf:
                    outf.write(line)

@transform(mergeNucMtSequences,
           regex('.*/(\S+).fa'),
           r'trna_sequences.dir/\1.fa.bwt')
def buildBWAIndex(infile, outfile):
    'Index the tRNA sequences for BWA'

    prefix = iotools.snip(outfile, '.bwt')
    statement = 'bwa index -p %(prefix)s %(infile)s' % locals()

    job_condaenv=PARAMS['conda_base_env']

    P.run(statement)


@transform(mergeNucMtSequences,
       regex('.*/(\S+).fa'),
       r'trna_sequences.dir/\1.fa.1.bt2')
def buildBowtie2Index(infile, outfile):
    'Index the tRNA sequences for Bowtie2'

    outfile_base = iotools.snip(outfile, '.1.bt2')

    statement = 'bowtie2-build %(infile)s %(outfile_base)s' % locals()

    job_condaenv=PARAMS['conda_base_env']

    P.run(statement)

@transform(mergeNucMtSequences,
           regex('.*/(\S+).fa'),
           r'trna_sequences.dir/\1_gsnap_index')
def buildGSNAPIndex(infile, outfile):
    'Index the tRNA sequences for GSNAP'

    outdir = os.path.dirname(outfile)
    outfile_base = os.path.basename(outfile)

    statement = 'gmap_build -q 1 -D %(outdir)s -d %(outfile_base)s %(infile)s' % locals()

    job_condaenv=PARAMS['conda_base_env']

    P.run(statement)


@transform(mergeNucMtSequences,
           regex('.*/(\S+).fa'),
           r'trna_sequences.dir/\1.fa.idx')
def buildSegemehlIndex(infile, outfile):
    'Index the tRNA sequences for Segemehl'

    statement = 'segemehl.x -x %(outfile)s -d %(infile)s' % locals()

    job_condaenv=PARAMS['conda_base_env']

    P.run(statement)

###################################################
# Learn mutation and truncation signatures
###################################################

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

    job_options = PARAMS['cluster_options'] + " -t 1:00:00"
    job_condaenv=PARAMS['conda_base_env']

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

    job_options = PARAMS['cluster_options'] + " -t 1:00:00"
    job_condaenv=PARAMS['conda_base_env']

    P.run(statement)

SEQUENCEFILES_REGEX2 = regex(r"(.*\/)*(YAMATseq_BT20_A).(fastq.1.gz|fastq.gz|fastq)")

@mkdir('mut_trunc_sig.dir')
@transform(SEQUENCEFILES,
           SEQUENCEFILES_REGEX,
           add_inputs(buildGSNAPIndex),
           r"mut_trunc_sig.dir/\2_GSNAPReportSingle.bam")
def mapWithGSNAPSingleReport(infiles, outfile):
    '''
    Map reads to the tRNA sequences, reporting only one alignment per read
    '''
    infile, gsnap_index = infiles

    gsnap_index_dir = os.path.dirname(gsnap_index)
    gsnap_index_base = os.path.basename(gsnap_index)

    tmp_file_sam = P.get_temp_filename()
    tmp_file_single_sam = P.get_temp_filename()

    statement = '''
    gsnap
    --gunzip
    -D %(gsnap_index_dir)s
    -d %(gsnap_index_base)s
    --format sam
    --genome-unk-mismatch 0
    --ignore-trim-in-filtering 1
    %(infile)s
    -o %(tmp_file_sam)s
    2> %(outfile)s.log;
    samtools flagstat %(tmp_file_sam)s > %(outfile)s.flagstat;
    ''' % locals()

    job_options = PARAMS['cluster_options'] + " -t 1:00:00"
    job_condaenv=PARAMS['conda_base_env']

    P.run(statement)

    CompareTrnaSeq.keep_random_alignment(tmp_file_sam, tmp_file_single_sam, outfile)

    os.unlink(tmp_file_sam)
    os.unlink(tmp_file_single_sam)


@mkdir('mut_trunc_sig.dir')
@transform(SEQUENCEFILES,
           SEQUENCEFILES_REGEX,
           add_inputs(mergeNucMtSequences),
           r"mut_trunc_sig.dir/\2_SHRiMPReportSingle.bam")
def mapWithSHRiMP2SingleReport(infiles, outfile):
    '''
    Map reads to the tRNA sequences, reporting only one alignment per read
    '''
    infile, trna_sequences = infiles

    tmp_file_fasta = P.get_temp_filename()
    tmp_file_sam = P.get_temp_filename()
    tmp_file_single_sam = P.get_temp_filename()

    job_threads = PARAMS['shrimp_threads']

    # Need to convert fastq to fasta for SHRiMP
    statement = '''
    zcat < %(infile)s |
    sed 's/^@//g'|
    paste - - - - |
    awk 'BEGIN { FS="\\t" } {print ">"$1"\\n"$2}' >
    %(tmp_file_fasta)s;

    gmapper
    %(tmp_file_fasta)s
    %(trna_sequences)s
    -N %(job_threads)s
    --strata  --report 1000
    --sam-unaligned
    --mode mirna
    > %(tmp_file_sam)s
    2> %(outfile)s.log ;
    samtools flagstat %(tmp_file_sam)s > %(outfile)s.flagstat;
    ''' % locals()

    job_options = PARAMS['cluster_options'] + " -t 3:00:00"
    job_condaenv=PARAMS['conda_base_env']

    P.run(statement)

    CompareTrnaSeq.keep_random_alignment(tmp_file_sam, tmp_file_single_sam, outfile)

    os.unlink(tmp_file_fasta)
    os.unlink(tmp_file_sam)
    os.unlink(tmp_file_single_sam)


# Not run
@mkdir('mut_trunc_sig.dir')
@transform(SEQUENCEFILES,
           SEQUENCEFILES_REGEX,
           add_inputs(buildSegemehlIndex, mergeNucMtSequences),
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

    job_options = PARAMS['cluster_options'] + " -t 5:00:00"
    job_condaenv=PARAMS['conda_base_env']

    P.run(statement)


@collate((mapWithBWAMEMSingleReport,
          mapWithBowtie2SingleReport,
          mapWithSHRiMP2SingleReport),
       regex('mut_trunc_sig.dir/(\S+)_(bwaMemReportSingle|bowtie2ReportSingle|SHRiMPReportSingle).bam'),
       r'mut_trunc_sig.dir/\1.merged.bam')
def mergeSingleReports(infiles, outfile):

    infiles = ' '.join(infiles)

    tmp_file = P.get_temp_filename()

    statement = '''
    samtools merge -f -o %(tmp_file)s %(infiles)s;
    samtools sort %(tmp_file)s > %(outfile)s;
    samtools index %(outfile)s;
    rm -f %(tmp_file)s
    ''' % locals()

    job_options = PARAMS['cluster_options'] + " -t 1:00:00"
    job_condaenv=PARAMS['conda_base_env']

    P.run(statement)


@transform((mapWithBWAMEMSingleReport,
            mapWithBowtie2SingleReport,
            mapWithSHRiMP2SingleReport,),
           suffix('.bam'),
           add_inputs(mergeNucMtSequences),
           '.summariseAlignments.pickle')
def summariseIndividualAlignments(infiles, outfile):
    '''
    Use alignmentSummary.clustalwtrnaAlignmentSummary class to summarise the
    alignments and learn the mutational & truncation signatures
    '''

    infile, trna_fasta = infiles

    job_options = PARAMS['cluster_options'] + " -t 1:00:00"
    job_condaenv=PARAMS['conda_base_env']

    CompareTrnaSeq.summariseAlignments(infile, trna_fasta, outfile, submit=True)

@transform(mergeSingleReports,
           suffix('.bam'),
           add_inputs(mergeNucMtSequences),
           '.summariseAlignments.pickle')
def summariseMergedAlignments(infiles, outfile):
    '''
    Use alignmentSummary.clustalwtrnaAlignmentSummary class to summarise the
    merged alignments and learn the mutational & truncation signatures
    '''

    infile, trna_fasta = infiles

    job_options = PARAMS['cluster_options'] + " -t 1:00:00"
    job_condaenv=PARAMS['conda_base_env']

    CompareTrnaSeq.summariseAlignments(infile, trna_fasta, outfile, submit=True)


###################################################
# simulate reads
###################################################

@mkdir('simulations.dir')
@transform(summariseMergedAlignments,
           formatter('mut_trunc_sig.dir/(?P<input_file>\S+).merged.summariseAlignments.pickle'),
           add_inputs(mergeNucMtSequences),
           'simulations.dir/{input_file[0]}.0.simulation_uniform.fastq.gz')
def simulation_uniform(infiles, outfile):

    alignment_summary_picklefile, infile = infiles

    n_reads = int(PARAMS['simulation_uniform_reads'])

    gt = simulateReads.make_gt(infile, n_reads, 1, 0, filter='\-Und\-')

    alignment_summary=pickle.load(open(alignment_summary_picklefile, 'rb'))

    job_options = PARAMS['cluster_options'] + " -t 1:00:00"
    job_condaenv=PARAMS['conda_base_env']


    CompareTrnaSeq.simulate_reads(
        infile=infile,
        outfile=outfile,
        ground_truth=gt,
        error_rate=1/100,
        mutation_threshold=0.1,
        truncate=True,
        alignment_summary=alignment_summary,
        summary_level='anticodon',
        submit=True,
        outfile_gt=outfile + '.gt')


@mkdir('simulation_dummy_files')
@originate(['simulation_dummy_files/%s' % x for x in range(0, PARAMS['simulation_n'])])
def create_files(output_file):
    with open(output_file, "w"):
        pass

@mkdir('simulations.dir')
@product(summariseMergedAlignments,
         formatter('mut_trunc_sig.dir/(?P<input_file>\S+).merged.summariseAlignments.pickle'),
         create_files,
         formatter('simulation_dummy_files/(?P<simulation_n>\d+)'),
         add_inputs(mergeNucMtSequences),
         'simulations.dir/{input_file[0][0]}.{simulation_n[1][0]}.simulation_withtrunc.fastq.gz')
def simulation_with_truncations(infiles, outfile):

    alignment_summary_picklefile = infiles[0][0]
    infile = infiles[1]

    n_reads = int(PARAMS['simulation_withmut_reads'])

    gt = simulateReads.make_gt(infile, n_reads, 10, 5, filter='\-Und\-')

    alignment_summary=pickle.load(open(alignment_summary_picklefile, 'rb'))

    job_options = PARAMS['cluster_options'] + " -t 1:00:00"
    job_condaenv=PARAMS['conda_base_env']

    CompareTrnaSeq.simulate_reads(
        infile=infile,
        outfile=outfile,
        ground_truth=gt,
        error_rate=1/100,
        mutation_threshold=0.1,
        truncate=True,
        alignment_summary=alignment_summary,
        summary_level='anticodon',
        submit=True,
        outfile_gt=outfile +'.gt')



###################################################
# align
###################################################

# not run
@mkdir('quant.dir')
@transform((simulation_uniform,
            simulation_with_truncations),
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

    job_options = PARAMS['cluster_options'] + " -t 1:00:00"
    job_condaenv=PARAMS['conda_base_env']

    P.run(statement)

    CompareTrnaSeq.filter_sam(tmp_file, outfile, submit=True)

    os.unlink(tmp_file)


@mkdir('quant.dir')
@transform((simulation_uniform,
            simulation_with_truncations),
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

    job_options = PARAMS['cluster_options'] + " -t 1:00:00"
    job_condaenv=PARAMS['conda_base_env']

    P.run(statement)

    CompareTrnaSeq.filter_sam(tmp_file, outfile, submit=True)

    os.unlink(tmp_file)


@mkdir('quant.dir')
@transform((simulation_uniform,
            simulation_with_truncations),
           regex('simulations.dir/(\S+).(simulation_\S+).fastq.gz'),
           add_inputs(mergeNucMtSequences),
           r'quant.dir/\1.\2.shrimp.bam')
def alignWithSHRiMP(infiles, outfile):

    infile, trna_sequences = infiles

    tmp_file_fasta = P.get_temp_filename()
    tmp_file_sam = P.get_temp_filename()

    job_threads = PARAMS['shrimp_threads']

    # Need to convert fastq to fasta for SHRiMP
    statement = '''
    zcat < %(infile)s |
    sed 's/^@//g'|
    paste - - - - |
    awk 'BEGIN { FS="\\t" } {print ">"$1"\\n"$2}' >
    %(tmp_file_fasta)s;

    gmapper
    %(tmp_file_fasta)s
    %(trna_sequences)s
    -N %(job_threads)s
    --strata  --report 1000
    --sam-unaligned
    --mode mirna
    > %(tmp_file_sam)s
    2> %(outfile)s.log ;
    samtools flagstat %(tmp_file_sam)s > %(outfile)s.flagstat;
    ''' % locals()

    job_options = PARAMS['cluster_options'] + " -t 1:00:00"
    job_condaenv=PARAMS['conda_base_env']

    P.run(statement)

    CompareTrnaSeq.filter_sam(tmp_file_sam, outfile, submit=True)
    os.unlink(tmp_file_fasta)
    os.unlink(tmp_file_sam)


@mkdir('quant.dir')
@transform((simulation_uniform,
            simulation_with_truncations),
           regex('simulations.dir/(\S+).(simulation_\S+).fastq.gz'),
           add_inputs(buildGSNAPIndex),
           r'quant.dir/\1.\2.gsnap.bam')
def alignWithGSNAP(infiles, outfile):

    infile, gsnap_index = infiles

    gsnap_index_dir = os.path.dirname(gsnap_index)
    gsnap_index_base = os.path.basename(gsnap_index)

    tmp_file = P.get_temp_filename()

    statement = '''
    gsnap
    --gunzip
    -D %(gsnap_index_dir)s
    -d %(gsnap_index_base)s
    --format sam
    --genome-unk-mismatch 0
    --ignore-trim-in-filtering 1
    %(infile)s
    -o %(tmp_file_sam)s
    2> %(outfile)s.log;
    samtools flagstat %(tmp_file_sam)s > %(outfile)s.flagstat;
    ''' % locals()

    job_options = PARAMS['cluster_options'] + " -t 1:00:00"
    job_condaenv=PARAMS['conda_base_env']

    P.run(statement)

    CompareTrnaSeq.filter_sam(tmp_file, outfile, submit=True)

    os.unlink(tmp_file)

###################################################
# quantify
###################################################

@transform((alignWithBowtie2,
            alignWithSHRiMP),
           regex('quant.dir/(\S+).(simulation_\S+)\.(\S+).bam'),
           [r'quant.dir/\1.\2.\3.gene_count_individual',
           r'quant.dir/\1.\2.\3.gene_count_isodecoder',
           r'quant.dir/\1.\2.\3.gene_count_anticodon'])
def quantDiscreteCounts(infile, outfiles):

    outfile_individual, outfile_isodecoder, outfile_anticodon = outfiles

    job_options = PARAMS['cluster_options'] + " -t 1:00:00"
    job_condaenv=PARAMS['conda_base_env']

    CompareTrnaSeq.tally_read_counts(
    infile,
    outfile_individual,
    outfile_isodecoder,
    outfile_anticodon,
    submit=True)

@transform((alignWithBowtie2,
            alignWithSHRiMP),
           regex('quant.dir/(\S+).(simulation_\S+)\.(\S+).bam'),
           [r'quant.dir/\1.\2.\3.gene_count_fractional_individual',
           r'quant.dir/\1.\2.\3.gene_count_fractional_isodecoder',
           r'quant.dir/\1.\2.\3.gene_count_fractional_anticodon'])
def quantFractionalCounts(infile, outfiles):

    outfile_individual, outfile_isodecoder, outfile_anticodon = outfiles

    job_options = PARAMS['cluster_options'] + " -t 1:00:00"
    job_condaenv=PARAMS['conda_base_env']

    CompareTrnaSeq.tally_read_counts(
    infile,
    outfile_individual,
    outfile_isodecoder,
    outfile_anticodon,
    discrete=False,
    submit=True)

@transform((alignWithBowtie2,
            alignWithSHRiMP),
           regex('quant.dir/(\S+).(simulation_\S+)\.(\S+).bam'),
           [r'quant.dir/\1.\2.\3.gene_count_mapq10_individual',
           r'quant.dir/\1.\2.\3.gene_count_mapq10_isodecoder',
           r'quant.dir/\1.\2.\3.gene_count_mapq10_anticodon'])
def quantDiscreteCountsMAPQ10(infile, outfiles):

    outfile_individual, outfile_isodecoder, outfile_anticodon = outfiles

    job_options = PARAMS['cluster_options'] + " -t 1:00:00"
    job_condaenv=PARAMS['conda_base_env']

    CompareTrnaSeq.tally_read_counts(
    infile,
    outfile_individual,
    outfile_isodecoder,
    outfile_anticodon,
    min_mapq=10,
    submit=True)


@transform((alignWithBowtie2,
            alignWithSHRiMP),
           regex('quant.dir/(\S+).(simulation_\S+)\.(\S+).bam'),
           [r'quant.dir/\1.\2.\3.gene_count_no_multi_individual',
           r'quant.dir/\1.\2.\3.gene_count_no_multi_isodecoder',
           r'quant.dir/\1.\2.\3.gene_count_no_multi_anticodon'])
def quantDiscreteCountsNoMultimapping(infile, outfiles):

    outfile_individual, outfile_isodecoder, outfile_anticodon = outfiles

    job_options = PARAMS['cluster_options'] + " -t 1:00:00"
    job_condaenv=PARAMS['conda_base_env']

    CompareTrnaSeq.tally_read_counts(
    infile,
    outfile_individual,
    outfile_isodecoder,
    outfile_anticodon,
    allow_multimapping=False,
    submit=True)

@transform((alignWithBowtie2,
            alignWithSHRiMP),
           regex('quant.dir/(\S+).(simulation_\S+)\.(\S+).bam'),
           [r'quant.dir/\1.\2.\3.gene_count_random_single_individual',
           r'quant.dir/\1.\2.\3.gene_count_random_single_isodecoder',
           r'quant.dir/\1.\2.\3.gene_count_random_single_anticodon'])
def quantDiscreteCountsRandomSingle(infile, outfiles):

    outfile_individual, outfile_isodecoder, outfile_anticodon = outfiles

    job_options = PARAMS['cluster_options'] + " -t 1:00:00"
    job_condaenv=PARAMS['conda_base_env']

    CompareTrnaSeq.tally_read_counts(
    infile,
    outfile_individual,
    outfile_isodecoder,
    outfile_anticodon,
    random_single=True,
    submit=True)


@transform((alignWithBowtie2,
            alignWithSHRiMP),
           regex('quant.dir/(\S+).(simulation_\S+)\.(\S+).bam'),
           add_inputs(mergeNucMtSequences),
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
    job_options = PARAMS['cluster_options'] + " -t 1:00:00"
    job_condaenv=PARAMS['conda_base_env']

    P.run(statement)


@follows(quantWithSalmon) # avoid mimseq running alongside salmon!
@jobs_limit(1)
@mkdir('mimseq.dir')
@collate((simulation_uniform, simulation_with_truncations),
         regex('simulations.dir/(\S+?)_(\S+).(simulation_\S+).fastq.gz'),
         r'mimseq.dir/\1_\3/counts/Anticodon_counts_raw.txt')
def quantWithMimSeq(infiles, outfile):

    nuc_trnas = PARAMS['trna_sequences_infile']
    mt_trnas = PARAMS['trna_mt_sequences_infile']
    trna_scan = PARAMS['trna_scan_infile']
    job_threads = PARAMS['mimseq_threads']

    tmp_outdir = P.get_temp_dir(clear=True)
    tmp_stdouterr = P.get_temp_filename()

    # create a dummy sample file so that all samples can be
    # quantified in one mimseq run
    tmp_sample_data = P.get_temp_filename(dir='.')
    with open(tmp_sample_data, 'w') as outf:
        for ix, infile in enumerate(infiles):
            if (ix % 2) == 0:
                condition = 'condition1'
            else:
                condition = 'condition2'
            outf.write('./%s\t%s\n' % (infile, condition))


    outdir = os.path.dirname(os.path.dirname(outfile))

    statement = '''
    mimseq
    --search vsearch
    --trnas %(nuc_trnas)s
    --mito-trnas %(mt_trnas)s
    --trnaout %(trna_scan)s
    --cluster-id 0.97
    --threads %(job_threads)s
    --min-cov 0.0005
    --max-mismatches 0.075
    --control-condition condition1
    -n results
    --out-dir %(tmp_outdir)s
    --max-multi 4
    --remap
    --remap-mismatches 0.05
    %(tmp_sample_data)s
    > %(tmp_stdouterr)s
    2>&1;
    rm -rf %(outdir)s;
    mkdir %(outdir)s;
    mv %(tmp_outdir)s/* %(outdir)s;
    mv %(tmp_stdouterr)s %(outdir)s/stdout_stderr
    ''' % locals()

    job_options = PARAMS['cluster_options'] + " -t 3:00:00"
    job_condaenv=PARAMS['mimseq_conda_env_name']

    P.run(statement)


###################################################
# compare ground truth and quantification estimates
###################################################
@transform(quantWithMimSeq,
           regex('mimseq.dir/(\S+?)_(simulation_\S+)/counts/Anticodon_counts_raw.txt'),
           add_inputs(r'simulations.dir/\1*.simulation_withtrunc.fastq.gz.gt'),
           [r'mimseq.dir/\1.\2.MimseqCompareTruthEstimateMimseqIsodecoder.tsv',
            r'mimseq.dir/\1.\2.MimseqCompareTruthEstimateAnticodon.tsv'])
def compareTruthEstimateMimseq(infiles, outfiles):
    '''mimseq reports isodecoder-level quantification with some isodecoders still merged together,
    hence isodecoder file is suffixed with MimSeqIsodecoder.tsv to distinguish if from
    the other Isodecoder.tsv outputs
    '''

    isodecoder_out, anticodon_out = outfiles

    infile = infiles[0]
    truths = infiles[1:]

    job_options = PARAMS['cluster_options'] + " -t 1:00:00"
    job_condaenv=PARAMS['conda_base_env']

    CompareTrnaSeq.compareMimSeq(infile, truths, isodecoder_out, anticodon_out, submit=True)



@collate(compareTruthEstimateMimseq,
         regex('mimseq.dir/(\S+)\.(simulation_\S+)\.MimseqCompareTruthEstimateMimseqIsodecoder.tsv'),
       [r'quant.dir/\2.Mimseq.CompareTruthEstimateMimseqIsodecoder.tsv',
        r'quant.dir/\2.Mimseq.CompareTruthEstimateAnticodon.tsv'])
def mergeCompareTruthEstimateMimseq(infiles, outfiles):
    outfile_isodecoder, outfile_anticodon = outfiles
    infiles_isodecoder = [x[0] for x in infiles]
    infiles_anticodon = [x[1] for x in infiles]

    job_options = PARAMS['cluster_options'] + " -t 1:00:00"
    job_condaenv=PARAMS['conda_base_env']

    CompareTrnaSeq.mergeCompareTruthEstimateMimseq(
        infiles_isodecoder, outfile_isodecoder, infiles_anticodon, outfile_anticodon, submit=True)


@merge(compareTruthEstimateMimseq,
       'quant.dir/mimseq_isodecoder_maps.tsv')
def mergeMimSeqIsodecoderMaps(infiles, outfile):

    infiles = [x[0] for x in infiles]

    job_options = PARAMS['cluster_options'] + " -t 1:00:00"
    job_condaenv=PARAMS['conda_base_env']

    CompareTrnaSeq.mergeMimSeqIsodecoderMaps(infiles, outfile, submit=True)

@collate(quantWithSalmon,
       regex(r'quant.dir/(\S+)\.(\S+).(simulation_\S+)\.(\S+)/quant.sf'),
       add_inputs(r'simulations.dir/\1.\2.\3.fastq.gz.gt'),
       [r'quant.dir/\3.Salmon.CompareTruthEstimate.tsv',
        r'quant.dir/\3.Salmon.CompareTruthEstimateIsodecoder.tsv',
        r'quant.dir/\3.Salmon.CompareTruthEstimateAnticodon.tsv'])
def compareTruthEstimateSalmon(infiles, outfiles):

    outfile_individual, outfile_isodecoder, outfile_anticodon = outfiles

    job_options = PARAMS['cluster_options'] + " -t 1:00:00"
    job_condaenv=PARAMS['conda_base_env']

    CompareTrnaSeq.compareTruthEstimateSalmon(
        infiles, outfile_individual, outfile_isodecoder, outfile_anticodon, submit=True)


@collate((quantFractionalCounts, quantDiscreteCounts, quantDiscreteCountsMAPQ10,
          quantDiscreteCountsRandomSingle, quantDiscreteCountsNoMultimapping),
         regex('quant.dir/(\S+).(simulation_\S+)\.(\S+).(gene_count.*)_.*'),
         add_inputs(r'simulations.dir/\1.\2.fastq.gz.gt'),
         [r'quant.dir/\2.\4_Decision.CompareTruthEstimate.tsv',
          r'quant.dir/\2.\4_Decision.CompareTruthEstimateIsodecoder.tsv',
          r'quant.dir/\2.\4_Decision.CompareTruthEstimateAnticodon.tsv'])
def compareTruthEstimateDecisionCounts(infiles, outfiles):

    outfile_individual, outfile_isodecoder, outfile_anticodon = outfiles

    job_options = PARAMS['cluster_options'] + " -t 1:00:00"
    job_condaenv=PARAMS['conda_base_env']

    CompareTrnaSeq.compareTruthEstimateDecisionCounts(
        infiles, outfile_individual, outfile_isodecoder, outfile_anticodon, submit=True)


@transform((compareTruthEstimateDecisionCounts, compareTruthEstimateSalmon),
           regex('quant.dir/(\S+).CompareTruthEstimate.tsv'),
           add_inputs(mergeMimSeqIsodecoderMaps),
           r'quant.dir/\1.CompareTruthEstimateMimseqIsodecoder.tsv')
def makeMimseqIsodecoderQuant(infiles, outfile):
    infile = infiles[0][1]
    mapping_file = infiles[1]

    job_options = PARAMS['cluster_options'] + " -t 1:00:00"
    job_condaenv=PARAMS['conda_base_env']

    CompareTrnaSeq.makeMimseqIsodecoderQuant(infile, mapping_file, outfile, submit=True)
# Function to merge comparisons for decision, salmon-based and mimseq quant

###################################################
# Alignment confusion matrices
###################################################
@follows(compareTruthEstimateMimseq)
@mkdir('truth2assignment.dir')
@transform((alignWithBowtie2,
            alignWithSHRiMP),
           regex('quant.dir/(\S+?)_(\S+).(simulation_\S+)\.(\S+).bam'),
           add_inputs(r'mimseq.dir/\1.\3.MimseqCompareTruthEstimateMimseqIsodecoder.tsv.mapping'),
               r'truth2assignment.dir/\1_\2.\3.\4.truth2assignment.tsv')
def getTruth2Assignment(infiles, outfile):
    'Get the tally of ground truths to assignments'

    infile, isodecoder_mapping = infiles

    job_options = PARAMS['cluster_options'] + " -t 1:00:00"
    job_condaenv=PARAMS['conda_base_env']

    CompareTrnaSeq.getTruth2Assignment(infile, isodecoder_mapping, outfile, submit=True)


@follows(compareTruthEstimateMimseq)
@mkdir('truth2assignment.dir')
@transform((simulation_uniform, simulation_with_truncations),
           regex('simulations.dir/(\S+?)_(\S+).(simulation_\S+).fastq.gz'),
           add_inputs(r'mimseq.dir/\1_\3/counts/Isodecoder_counts_raw.txt',
                      r'mimseq.dir/\1.\3.MimseqCompareTruthEstimateMimseqIsodecoder.tsv.mapping',
                      r'mimseq.dir/\1_\3/align/\1_\2.\3.unpaired_uniq.bam'),
           r'truth2assignment.dir/\1_\2.\3.mimseq.truth2assignment.tsv')
def getTruth2AssignmentMimSeq(infiles, outfile):
    'Get the tally of ground truths to assignments'

    infile_sim, mimseq_isodecoder_counts, isodecoder_mapping, infile = infiles

    job_options = PARAMS['cluster_options'] + " -t 1:00:00"
    job_condaenv=PARAMS['conda_base_env']

    CompareTrnaSeq.getTruth2AssignmentMimSeq(
        infile, mimseq_isodecoder_counts, isodecoder_mapping, outfile, submit=True)


@collate((getTruth2Assignment, getTruth2AssignmentMimSeq),
         regex('truth2assignment.dir/(\S+)\.(\d+)\.(simulation_\S+)\.(\S+)\.truth2assignment.tsv'),
         r'truth2assignment.dir/truth2assignment.\3.tsv')
def mergeTruth2Assignment(infiles, outfile):

    job_options = PARAMS['cluster_options'] + " -t 1:00:00"
    job_condaenv=PARAMS['conda_base_env']

    CompareTrnaSeq.mergeTruth2Assignment(infiles, outfile, submit=True)

    infiles_anticodon  = [re.sub('.tsv$', '_anticodon.tsv', x) for x in infiles]
    outfile_anticodon = re.sub('.tsv$', '_anticodon.tsv', outfile)
    CompareTrnaSeq.mergeTruth2Assignment(infiles_anticodon, outfile_anticodon, submit=True)

    infiles_isodecoder  = [re.sub('.tsv$', '_isodecoder.tsv', x) for x in infiles]
    # no isodecoder level output from mimseq
    infiles_isodecoder = [x for x in infiles_isodecoder if os.path.exists(x)]
    outfile_isodecoder = re.sub('.tsv$', '_isodecoder.tsv', outfile)
    CompareTrnaSeq.mergeTruth2Assignment(infiles_isodecoder, outfile_isodecoder, submit=True)

    infiles_mimseq_isodecoder  = [re.sub('.tsv$', '_mimseq_isodecoder.tsv', x) for x in infiles]
    outfile_mimseq_isodecoder = re.sub('.tsv$', '_mimseq_isodecoder.tsv', outfile)
    CompareTrnaSeq.mergeTruth2Assignment(infiles_mimseq_isodecoder, outfile_mimseq_isodecoder, submit=True)
####################################################




###################################################
# Targets
###################################################


@follows(simulation_uniform,
          simulation_with_truncations)
def simulate():
    'simulate tRNA reads from basic to more realistic'
    pass


@follows(summariseIndividualAlignments, summariseMergedAlignments)
def summariseAlignments():
    'summarise the alignments with real reads'
    pass


@follows(quantWithSalmon, quantWithMimSeq)
def quant():
    'Quantify from the simulated reads'
    pass

@follows(compareTruthEstimateSalmon,
         compareTruthEstimateDecisionCounts,
         compareTruthEstimateMimseq,
         mergeCompareTruthEstimateMimseq,
         makeMimseqIsodecoderQuant,
         mergeTruth2Assignment)
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
