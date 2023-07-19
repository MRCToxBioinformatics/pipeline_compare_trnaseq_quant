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
import shutil
import re
import pickle
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
@originate('trna_sequences.dir/nuc_sequences.fa')
def filterNucFasta(outfile):
    CompareTrnaSeq.filterFasta(PARAMS['trna_sequences_infile'], outfile)

@mkdir('trna_sequences.dir')
@originate('trna_sequences.dir/mt_sequences.fa')
def filterMtFasta(outfile):
    CompareTrnaSeq.filterFasta(PARAMS['trna_mt_sequences_infile'], outfile)

@mkdir('trna_sequences.dir')
@transform(filterMtFasta,
           suffix('.fa'),
           '.filtered.fa')
def updateMtFastaNaming(infile, outfile):
    CompareTrnaSeq.updateMtFastaNaming(infile, outfile)

@merge((filterNucFasta, updateMtFastaNaming),
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
# Simulation full length no mutations and compare alignment strategies

@transform(mergeNucMtSequences,
           suffix('.fa'),
           '_simulation_no_error.fastq.gz')
def simulation_no_errors(infile, outfile):

    n_reads = int(PARAMS['simulation_uniform_reads'])

    gt = simulateReads.make_gt(infile, n_reads, 1, 0, filter='\-Und\-')

    job_options = PARAMS['cluster_options'] + " -t 2:00:00"
    job_condaenv=PARAMS['conda_base_env']


    CompareTrnaSeq.simulate_reads(
        infile=infile,
        outfile=outfile,
        ground_truth=gt,
        mutation_threshold=None,
        error_rate=0.001,
        truncate=False,
        outfile_gt=outfile + '.gt',
        submit=True,
        job_options=job_options)


no_errors_alignment_params_files = []
for d_value in P.as_list(PARAMS['bowtie2_d_values']):
    for l_value in P.as_list(PARAMS['bowtie2_l_values']):
        for n_value in P.as_list(PARAMS['bowtie2_n_values']):
            no_errors_alignment_params_files.append(
                'no_errors_aligment_params.dir/d%s_l%s_n%s_bowtie2.tsv' % (d_value, l_value, n_value))


for k_value in P.as_list(PARAMS['bwamem_k_values']):
    for r_value in P.as_list(PARAMS['bwamem_r_values']):
        no_errors_alignment_params_files.append(
        'no_errors_aligment_params.dir/k%s_r%s_bwamem.tsv' % (k_value, r_value))


@mkdir('no_errors_aligment_params.dir')
@originate(no_errors_alignment_params_files)
def create_no_errors_sim_alignment_param_files(output_file):
    with open(output_file, "w"):
        pass



@mkdir('simulation_no_errors.dir')
@transform(create_no_errors_sim_alignment_param_files,
           regex('no_errors_aligment_params.dir/(\S+)_bwamem.tsv'),
           add_inputs(simulation_no_errors, buildBWAIndex),
           r"simulation_no_errors.dir/bwamem_\1_ReportSingle.bam")
def mapSimulationNoErrorsBwamemSingleReport(infiles, outfile):
    '''
    Map reads to the tRNA sequences, reporting only one alignment per read
    '''

    params_dummy, infile, bwa_index = infiles
    
    bwa_index = iotools.snip(bwa_index, '.bwt')

    tmp_file = P.get_temp_filename()

    k_value, r_value = [x[1:] for x in os.path.basename(params_dummy).strip('_bwamem.tsv').split('_')]

    statement = '''
    bwa mem
    -k %(k_value)s
    -r %(r_value)s
    -T 15
    %(bwa_index)s %(infile)s
    > %(tmp_file)s 2> %(outfile)s.log;
    samtools flagstat %(tmp_file)s > %(outfile)s.flagstat;
    samtools sort -o %(outfile)s -O BAM %(tmp_file)s;
    samtools index %(outfile)s;
    rm -f %(tmp_file)s
    ''' % locals()

    job_options = PARAMS['cluster_options'] + " -t 2:00:00"
    job_condaenv=PARAMS['conda_base_env']

    P.run(statement)


@mkdir('simulation_no_errors.dir')
@transform(create_no_errors_sim_alignment_param_files,
           regex('no_errors_aligment_params.dir/(\S+)_bowtie2.tsv'),
           add_inputs(simulation_no_errors, buildBowtie2Index),
           r"simulation_no_errors.dir/bowtie2_\1_ReportSingle.bam")
def mapSimulationNoErrorsBowtie2SingleReport(infiles, outfile):
    '''
    Map reads to the tRNA sequences, reporting only one alignment per read
    '''

    params_dummy, infile, bowtie2_index = infiles

    d_value, l_value, n_value = [x[1:] for x in os.path.basename(params_dummy).strip('_bowtie2.tsv').split('_')]

    bowtie2_index = iotools.snip(bowtie2_index, '.1.bt2')

    tmp_file = P.get_temp_filename()

    statement = '''
    bowtie2
    --min-score G,1,8
    --local
    -D %(d_value)s
    -N %(n_value)s
    -L %(l_value)s
    -R 3
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

    job_options = PARAMS['cluster_options'] + " -t 8:00:00"
    job_condaenv=PARAMS['conda_base_env']

    P.run(statement)


@mkdir('truth2assignment_no_errors.dir')
@transform((mapSimulationNoErrorsBowtie2SingleReport, mapSimulationNoErrorsBwamemSingleReport),
           regex('simulation_no_errors.dir/(\S+)_ReportSingle.bam'),
           r'truth2assignment_no_errors.dir/\1_truth2assignment.tsv')
def getTruth2AssignmentNoErrors(infile, outfile):
    'Get the tally of ground truths to assignments'

    job_options = PARAMS['cluster_options'] + " -t 2:00:00"
    job_condaenv=PARAMS['conda_base_env']

    CompareTrnaSeq.getTruth2Assignment(infile, outfile, submit=True, job_options=job_options)

    
    
@mkdir('final_results.dir')
@collate(getTruth2AssignmentNoErrors,
         regex('truth2assignment_no_errors.dir/(bowtie2|bwamem)_(\S+)_truth2assignment.tsv'),
         r'final_results.dir/\1_truth2assignment_no_errors.tsv', r'\1')
def mergeTruth2AssignmentNoErrors(infiles, outfile, aligner):

    job_options = PARAMS['cluster_options'] + " -t 1:00:00"
    job_condaenv=PARAMS['conda_base_env']

    CompareTrnaSeq.mergeTruth2AssignmentNoErrors(infiles, outfile, aligner, submit=True, job_options=job_options)

    infiles_anticodon  = [re.sub('.tsv$', '_anticodon.tsv', x) for x in infiles]
    outfile_anticodon = re.sub('.tsv$', '_anticodon.tsv', outfile)
    CompareTrnaSeq.mergeTruth2AssignmentNoErrors(infiles_anticodon, outfile_anticodon, aligner, submit=True, job_options=job_options)

    infiles_isodecoder  = [re.sub('.tsv$', '_isodecoder.tsv', x) for x in infiles]
    outfile_isodecoder = re.sub('.tsv$', '_isodecoder.tsv', outfile)
    CompareTrnaSeq.mergeTruth2AssignmentNoErrors(infiles_isodecoder, outfile_isodecoder, aligner, submit=True, job_options=job_options)


####
#

@mkdir('multimapping.dir')
@transform(SEQUENCEFILES,
           SEQUENCEFILES_REGEX,
           add_inputs(buildBowtie2Index),
           r"multimapping.dir/\2_bowtie2ReportAll.bam")
def alignWithBowtie2Multimapping(infiles, outfile):

    infile, bowtie2_index = infiles

    bowtie2_index = iotools.snip(bowtie2_index, '.1.bt2')

    tmp_file = P.get_temp_filename()

    job_threads = PARAMS['bowtie2_threads']

    statement = '''
    bowtie2
    --min-score G,1,8
    --local
    -a
    -D 100 -R 3 -N 1 -L 10
    -i S,1,0.5
    -p %(job_threads)s
    -x %(bowtie2_index)s
    -U %(infile)s
    -S %(tmp_file)s
    2> %(outfile)s.log;
    samtools flagstat %(tmp_file)s > %(outfile)s.flagstat;
    ''' % locals()

    job_options = PARAMS['cluster_options'] + " -t 8:00:00"
    job_condaenv=PARAMS['conda_base_env']

    P.run(statement)

    CompareTrnaSeq.filter_sam(tmp_file, outfile, submit=True, job_options=job_options)

    os.unlink(tmp_file)




###################################################
# Learn mutation and truncation signatures


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

    job_options = PARAMS['cluster_options'] + " -t 2:00:00"
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
    -D 100 -R 3 -N 1 -L 10
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

    job_options = PARAMS['cluster_options'] + " -t 8:00:00"
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

    job_options = PARAMS['cluster_options'] + " -t 2:00:00"
    job_condaenv=PARAMS['conda_base_env']

    P.run(statement)

    CompareTrnaSeq.keep_random_alignment(tmp_file_sam, tmp_file_single_sam, outfile, submit=True, job_options=job_options)

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


#@collate((mapWithBWAMEMSingleReport,
#          mapWithBowtie2SingleReport,
#          mapWithSHRiMP2SingleReport),
@collate((mapWithBWAMEMSingleReport,
          mapWithBowtie2SingleReport),
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

    job_options = PARAMS['cluster_options'] + " -t 2:00:00"
    job_condaenv=PARAMS['conda_base_env']

    P.run(statement)


#@transform((mapWithBWAMEMSingleReport,
#            mapWithBowtie2SingleReport,
#            mapWithSHRiMP2SingleReport,),
@transform((mapWithBWAMEMSingleReport,
            mapWithBowtie2SingleReport),
           suffix('.bam'),
           add_inputs(mergeNucMtSequences),
           '.summariseAlignments.pickle')
def summariseIndividualAlignments(infiles, outfile):
    '''
    Use alignmentSummary.clustalwtrnaAlignmentSummary class to summarise the
    alignments and learn the mutational & truncation signatures
    '''

    infile, trna_fasta = infiles

    job_options = PARAMS['cluster_options'] + " -t 2:00:00"
    job_condaenv=PARAMS['conda_base_env']

    CompareTrnaSeq.summariseAlignments(infile, trna_fasta, outfile, submit=True, job_options=job_options)


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

    job_options = PARAMS['cluster_options'] + " -t 2:00:00"
    job_condaenv = PARAMS['conda_base_env']

    CompareTrnaSeq.summariseAlignments(infile, trna_fasta, outfile, submit=True, job_options=job_options)

@collate(summariseMergedAlignments,
         regex('mut_trunc_sig.dir/(\S+?)_.*.merged.summariseAlignments.pickle'),
         r'mut_trunc_sig.dir/\1_common_genes.tsv')
def defineCommonGenesPerMethod(infiles, outfile):
    mut_dict = pickle.load(open(infiles[0], 'rb'))
    common_genes = set(mut_dict.alignment_coordinates.keys())

    for infile in infiles[1:]:
        mut_dict = pickle.load(open(infile, 'rb'))
        common_genes = common_genes.intersection(set(mut_dict.alignment_coordinates.keys()))

    with open(outfile, 'w') as outf:
        for gene in common_genes:
            outf.write('%s\n' % gene)


###################################################
# compare errors to known modifications
###################################################

@mkdir('final_results.dir')
@transform(mergeNucMtSequences,
           regex('trna_sequences.dir/trna_sequences_all.fa'),
           r'final_results.dir/modomics_mapped_to_fasta.tsv')
def mapModomics2fasta(infile, outfile):
    '''
    Take the Modomics modifications and map them to the tRNA sequences
    '''

    modification_index = PARAMS['modification_index']
    modomics_json = PARAMS['modomics_json']

    job_options = PARAMS['cluster_options'] + " -t 2:00:00"
    job_condaenv=PARAMS['conda_base_env']

    CompareTrnaSeq.mapModomics2fasta(
        infile, modomics_json, modification_index, outfile,
        submit=False, job_options=job_options)

@follows(mapModomics2fasta)
@merge(summariseMergedAlignments,
       ['final_results.dir/mutations_vs_modomics.tsv',
        'final_results.dir/truncations_summary.tsv',
        'final_results.dir/truncations_vs_modomics.tsv'])
def mergeMutationProfileModomics(infiles, outfiles):

    outfile_mutations, outfile_read_end, outfile_truncations = outfiles
    modomics_positions = 'final_results.dir/modomics_mapped_to_fasta.tsv'

    job_options = PARAMS['cluster_options'] + " -t 2:00:00"
    job_condaenv=PARAMS['conda_base_env']

    CompareTrnaSeq.mergeMutationProfileModomics(
        infiles, modomics_positions,
        outfile_mutations, outfile_read_end, outfile_truncations,
        submit=False, job_options=job_options)


###################################################
# simulate reads
###################################################

@follows(defineCommonGenesPerMethod)
@mkdir('simulations.dir')
@transform(SEQUENCEFILES,
           SEQUENCEFILES_REGEX,
           r'simulations.dir/\2.0.simulation_null.fastq.gz')
def simulation_null(infile, outfile):
    '''create softlinks to real data in same filestructure format as simulations so
    they can be processed by the same quantification tasks
    '''

    shutil.copyfile(os.path.abspath(infile), outfile)

@follows(defineCommonGenesPerMethod)
@mkdir('simulations.dir')
@transform(summariseMergedAlignments,
           regex('mut_trunc_sig.dir/(\S+?)_(\S+).merged.summariseAlignments.pickle'),
           add_inputs(mergeNucMtSequences, r'mut_trunc_sig.dir/\1_common_genes.tsv'),
           r'simulations.dir/\1_\2.0.simulation_uniform.fastq.gz')
def simulation_uniform(infiles, outfile):

    alignment_summary_picklefile, infile, common_genes_inf = infiles


    common_genes = set()
    with open(common_genes_inf, 'r') as inf:
        for line in inf:
            common_genes.add(line.strip(),)

    n_reads = int(PARAMS['simulation_uniform_reads'])

    gt = simulateReads.make_gt(infile, n_reads, 1, 0, filter='\-Und\-', genes=common_genes)

    alignment_summary=pickle.load(open(alignment_summary_picklefile, 'rb'))

    job_options = PARAMS['cluster_options'] + " -t 2:00:00"
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
        outfile_gt=outfile + '.gt',
        submit=True,
        job_options=job_options)


@mkdir('simulation_dummy_files')
@originate(['simulation_dummy_files/%s' % x for x in range(0, PARAMS['simulation_n'])])
def create_files(output_file):
    with open(output_file, "w"):
        pass


@follows(defineCommonGenesPerMethod)
@mkdir('simulations.dir')
@product(summariseMergedAlignments,
         formatter('mut_trunc_sig.dir/(?P<trna_seq_method>\S+?)_(?P<sample>\S+).merged.summariseAlignments.pickle'),
         create_files,
         formatter('simulation_dummy_files/(?P<simulation_n>\d+)'),
         add_inputs(mergeNucMtSequences, 'mut_trunc_sig.dir/{trna_seq_method[0][0]}_common_genes.tsv'),
         'simulations.dir/{trna_seq_method[0][0]}_{sample[0][0]}.{simulation_n[1][0]}.simulation_realistic.fastq.gz')
def simulation_realistic(infiles, outfile):

    alignment_summary_picklefile = infiles[0][0]
    common_genes_inf = infiles[2]

    common_genes = set()
    with open(common_genes_inf, 'r') as inf:
        for line in inf:
            common_genes.add(line.strip(),)

    n_reads = int(PARAMS['simulation_withmut_reads'])

    infile = infiles[1]

    gt = simulateReads.make_gt(infile, n_reads, 10, 5, filter='\-Und\-', genes=common_genes)

    alignment_summary=pickle.load(open(alignment_summary_picklefile, 'rb'))

    job_options = PARAMS['cluster_options'] + " -t 2:00:00"
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
        outfile_gt=outfile +'.gt',
        submit=True,
        job_options=job_options)



###################################################
# align
###################################################

# not run
@mkdir('quant.dir')
@transform((simulation_null,
            simulation_uniform,
            simulation_realistic),
           regex('simulations.dir/(\S+).(simulation_\S+).fastq.gz'),
           add_inputs(buildBWAIndex),
           r'quant.dir/\1.\2.bwa.bam')
def alignWithBWA(infiles, outfile):

    infile, bwa_index= infiles

    bwa_index = iotools.snip(bwa_index, '.bwt')

    tmp_file = P.get_temp_filename()

    job_threads = PARAMS['bwa_threads']

    statement = '''
    bwa mem
    -k 10
    -T 15
    -a
    -t %(job_threads)s
    %(bwamem_index)s
    %(infile)s
    > %(tmp_file)s
    2> %(outfile)s.log;
    samtools flagstat %(tmp_file)s > %(outfile)s.flagstat;
    '''

    job_options = PARAMS['cluster_options'] + " -t 8:00:00"
    job_condaenv=PARAMS['conda_base_env']

    P.run(statement)

    CompareTrnaSeq.filter_sam(tmp_file, outfile, submit=True, job_options=job_options)

    os.unlink(tmp_file)


@mkdir('quant.dir')
@transform((simulation_null,
            simulation_uniform,
            simulation_realistic),
           regex('simulations.dir/(\S+).(simulation_\S+).fastq.gz'),
           add_inputs(buildBowtie2Index),
           r'quant.dir/\1.\2.bowtie2.bam')
def alignWithBowtie2(infiles, outfile):

    infile, bowtie2_index = infiles

    bowtie2_index = iotools.snip(bowtie2_index, '.1.bt2')

    tmp_file = P.get_temp_filename()

    job_threads = PARAMS['bowtie2_threads']

    # -a = report all alignments
    statement = '''
    bowtie2
    --min-score G,1,8
    --local
    -a
    -D 20 -R 3 -N 1 -L 10
    -i S,1,0.5
    -p %(job_threads)s
    -x %(bowtie2_index)s
    -U %(infile)s
    -S %(tmp_file)s
    2> %(outfile)s.log;
    samtools flagstat %(tmp_file)s > %(outfile)s.flagstat;
    ''' % locals()

    job_options = PARAMS['cluster_options'] + " -t 8:00:00"
    job_condaenv=PARAMS['conda_base_env']

    P.run(statement)

    CompareTrnaSeq.filter_sam(tmp_file, outfile, submit=True, job_options=job_options)

    os.unlink(tmp_file)


@mkdir('quant.dir')
@transform((simulation_null,
            simulation_uniform,
            simulation_realistic),
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

    job_options = PARAMS['cluster_options'] + " -t 4:00:00"
    job_condaenv=PARAMS['conda_base_env']

    P.run(statement)

    CompareTrnaSeq.filter_sam(tmp_file_sam, outfile, submit=True, job_options=job_options)
    os.unlink(tmp_file_fasta)
    os.unlink(tmp_file_sam)


@mkdir('quant.dir')
@transform((simulation_null,
            simulation_uniform,
            simulation_realistic),
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

    job_options = PARAMS['cluster_options'] + " -t 4:00:00"
    job_condaenv=PARAMS['conda_base_env']

    P.run(statement)

    CompareTrnaSeq.filter_sam(tmp_file, outfile, submit=True, job_options=job_options)

    os.unlink(tmp_file)

###################################################
# quantify
###################################################

#@transform((alignWithBowtie2,
#            alignWithSHRiMP),
@transform(alignWithBowtie2,
           regex('quant.dir/(\S+).(simulation_\S+)\.(\S+).bam'),
           [r'quant.dir/\1.\2.\3.gene_count_decision.individual',
           r'quant.dir/\1.\2.\3.gene_count_decision.isodecoder',
           r'quant.dir/\1.\2.\3.gene_count_decision.anticodon'])
def quantDiscreteCounts(infile, outfiles):

    outfile_individual, outfile_isodecoder, outfile_anticodon = outfiles

    job_options = PARAMS['cluster_options'] + " -t 2:00:00"
    job_condaenv=PARAMS['conda_base_env']

    CompareTrnaSeq.tally_read_counts(
        infile,
        outfile_individual,
        outfile_isodecoder,
        outfile_anticodon,
        submit=True,
        job_options=job_options)

@transform(alignWithBowtie2,
           regex('quant.dir/(\S+).(simulation_\S+)\.(\S+).bam'),
           [r'quant.dir/\1.\2.\3.gene_count_fractional.individual',
           r'quant.dir/\1.\2.\3.gene_count_fractional.isodecoder',
           r'quant.dir/\1.\2.\3.gene_count_fractional.anticodon'])
def quantFractionalCounts(infile, outfiles):

    outfile_individual, outfile_isodecoder, outfile_anticodon = outfiles

    job_options = PARAMS['cluster_options'] + " -t 2:00:00"
    job_condaenv=PARAMS['conda_base_env']

    CompareTrnaSeq.tally_read_counts(
    infile,
    outfile_individual,
    outfile_isodecoder,
    outfile_anticodon,
    discrete=False,
        submit=True,
     job_options=job_options)


# not shrimp as it doesn't report MAPQ as parameterised here
@transform(alignWithBowtie2,
           regex('quant.dir/(\S+).(simulation_\S+)\.(\S+).bam'),
           [r'quant.dir/\1.\2.\3.gene_count_mapq10.individual',
           r'quant.dir/\1.\2.\3.gene_count_mapq10.isodecoder',
           r'quant.dir/\1.\2.\3.gene_count_mapq10.anticodon'])
def quantDiscreteCountsMAPQ10(infile, outfiles):

    outfile_individual, outfile_isodecoder, outfile_anticodon = outfiles

    job_options = PARAMS['cluster_options'] + " -t 2:00:00"
    job_condaenv=PARAMS['conda_base_env']

    CompareTrnaSeq.tally_read_counts(
    infile,
    outfile_individual,
    outfile_isodecoder,
    outfile_anticodon,
    min_mapq=10,
        submit=True,
         job_options=job_options)


@transform(alignWithBowtie2,
           regex('quant.dir/(\S+).(simulation_\S+)\.(\S+).bam'),
           [r'quant.dir/\1.\2.\3.gene_count_no_multi.individual',
           r'quant.dir/\1.\2.\3.gene_count_no_multi.isodecoder',
           r'quant.dir/\1.\2.\3.gene_count_no_multi.anticodon'])
def quantDiscreteCountsNoMultimapping(infile, outfiles):

    outfile_individual, outfile_isodecoder, outfile_anticodon = outfiles

    job_options = PARAMS['cluster_options'] + " -t 2:00:00"
    job_condaenv=PARAMS['conda_base_env']

    CompareTrnaSeq.tally_read_counts(
    infile,
    outfile_individual,
    outfile_isodecoder,
    outfile_anticodon,
    allow_multimapping=False,
        submit=True,
     job_options=job_options)

@transform(alignWithBowtie2,
           regex('quant.dir/(\S+).(simulation_\S+)\.(\S+).bam'),
           [r'quant.dir/\1.\2.\3.gene_count_random_single.individual',
           r'quant.dir/\1.\2.\3.gene_count_random_single.isodecoder',
           r'quant.dir/\1.\2.\3.gene_count_random_single.anticodon'])
def quantDiscreteCountsRandomSingle(infile, outfiles):

    outfile_individual, outfile_isodecoder, outfile_anticodon = outfiles

    job_options = PARAMS['cluster_options'] + " -t 2:00:00"
    job_condaenv=PARAMS['conda_base_env']

    CompareTrnaSeq.tally_read_counts(
    infile,
    outfile_individual,
    outfile_isodecoder,
    outfile_anticodon,
    random_single=True,
        submit=True,
     job_options=job_options)


@transform(alignWithBowtie2,
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
    job_options = PARAMS['cluster_options'] + " -t 2:00:00"
    job_condaenv=PARAMS['conda_base_env']

    P.run(statement)


@follows(quantWithSalmon) # avoid mimseq running alongside salmon! (only needed for non-HPC runs)
@jobs_limit(PARAMS['mimseq_max_sim_tasks'])
@mkdir('mimseq.dir')
@collate((simulation_uniform, simulation_realistic, simulation_null),
         regex('simulations.dir/(\S+?)_(\S+).(simulation_\S+).fastq.gz'),
         add_inputs(filterNucFasta, filterMtFasta),
         r'mimseq.dir/\1_\3/counts/Anticodon_counts_raw.txt')
def quantWithMimSeq(infiles, outfile):

    nuc_trnas, mt_trnas  = infiles[0][1:3]
    infiles = [x[0] for x in infiles]

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


    #--search vsearch

    statement = '''
    mimseq
    -t %(nuc_trnas)s
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

    job_options = PARAMS['cluster_options'] + " -t 15:00:00"
    job_condaenv=PARAMS['mimseq_conda_env_name']

    P.run(statement)


#####################################################
# concatenate quantification estiamtes from real data
#####################################################

@mkdir('final_results.dir')
@collate((quantFractionalCounts, quantDiscreteCounts, quantDiscreteCountsMAPQ10,
          quantDiscreteCountsRandomSingle, quantDiscreteCountsNoMultimapping),
         regex('quant.dir/(\S+).(simulation_null)\.(\S+).(gene_count_.*).*'),
         [r'final_results.dir/\2.\4_Decision.ConcatenateEstimate.tsv',
          r'final_results.dir/\2.\4_Decision.ConcatenateEstimateIsodecoder.tsv',
          r'final_results.dir/\2.\4_Decision.ConcatenateEstimateAnticodon.tsv'])
def concatenateEstimateDecisionCounts(infiles, outfiles):

    job_options = PARAMS['cluster_options'] + " -t 1:00:00"
    job_condaenv=PARAMS['conda_base_env']

    CompareTrnaSeq.concatenateEstimateDecisionCounts(
        infiles,
        outfiles,
        submit=True,  job_options=job_options)


@mkdir('final_results.dir')
@collate(quantWithSalmon,
       regex(r'quant.dir/(\S+)\.(\S+).(simulation_null)\.(\S+)/quant.sf'),
       [r'final_results.dir/\3.Salmon.ConcatenateEstimate.tsv',
        r'final_results.dir/\3.Salmon.ConcatenateEstimateIsodecoder.tsv',
        r'final_results.dir/\3.Salmon.ConcatenateEstimateAnticodon.tsv'])
def concatenateEstimateSalmon(infiles, outfiles):

    outfile_individual, outfile_isodecoder, outfile_anticodon = outfiles

    job_options = PARAMS['cluster_options'] + " -t 1:00:00"
    job_condaenv=PARAMS['conda_base_env']

    CompareTrnaSeq.concatenateEstimateSalmon(
        infiles,
        outfile_individual, outfile_isodecoder, outfile_anticodon,
        submit=True, job_options=job_options)

###################################################
# compare ground truth and quantification estimates
###################################################
@transform(quantWithMimSeq,
           regex('mimseq.dir/(\S+?)_(simulation_uniform|simulation_realistic)/counts/Anticodon_counts_raw.txt'),
           add_inputs(r'simulations.dir/\1*.\2.fastq.gz.gt'),
           [r'mimseq.dir/\1.\2.MimseqCompareTruthEstimateMimseqIsodecoder.tsv',
            r'mimseq.dir/\1.\2.MimseqCompareTruthEstimateAnticodon.tsv'])
def compareTruthEstimateMimseq(infiles, outfiles):
    '''mimseq reports isodecoder-level quantification with some isodecoders still merged together,
    hence isodecoder file is suffixed with MimSeqIsodecoder.tsv to distinguish it from
    the other Isodecoder.tsv outputs
    '''

    isodecoder_out, anticodon_out = outfiles

    infile = infiles[0]
    truths = infiles[1:]

    job_options = PARAMS['cluster_options'] + " -t 2:00:00"
    job_condaenv=PARAMS['conda_base_env']

    CompareTrnaSeq.compareMimSeq(infile, truths, isodecoder_out, anticodon_out, submit=True, job_options=job_options)


@mkdir('final_results.dir')
@collate(compareTruthEstimateMimseq,
         regex('mimseq.dir/(\S+)\.(simulation_\S+)\.MimseqCompareTruthEstimateMimseqIsodecoder.tsv'),
       [r'final_results.dir/\2.Mimseq.CompareTruthEstimateMimseqIsodecoder.tsv',
        r'final_results.dir/\2.Mimseq.CompareTruthEstimateAnticodon.tsv'])
def mergeCompareTruthEstimateMimseq(infiles, outfiles):
    outfile_isodecoder, outfile_anticodon = outfiles
    infiles_isodecoder = [x[0] for x in infiles]
    infiles_anticodon = [x[1] for x in infiles]

    job_options = PARAMS['cluster_options'] + " -t 1:00:00"
    job_condaenv=PARAMS['conda_base_env']

    CompareTrnaSeq.mergeCompareTruthEstimateMimseq(
        infiles_isodecoder, outfile_isodecoder, infiles_anticodon, outfile_anticodon, submit=True, job_options=job_options)


@collate(compareTruthEstimateMimseq,
         regex('mimseq.dir/(\S+?).(simulation_\S+).MimseqCompareTruthEstimateMimseqIsodecoder.tsv'),
         r'quant.dir/\2_mimseq_isodecoder_maps.tsv')
def mergeMimSeqIsodecoderMaps(infiles, outfile):

    infiles = [x[0] for x in infiles]

    job_options = PARAMS['cluster_options'] + " -t 1:00:00"
    job_condaenv=PARAMS['conda_base_env']

    CompareTrnaSeq.mergeMimSeqIsodecoderMaps(infiles, outfile, submit=True, job_options=job_options)

@mkdir('final_results.dir')
@collate(quantWithSalmon,
       regex(r'quant.dir/(\S+)\.(\S+).(simulation_uniform|simulation_realistic)\.(\S+)/quant.sf'),
       add_inputs(r'simulations.dir/\1.\2.\3.fastq.gz.gt'),
       [r'final_results.dir/\3.Salmon.CompareTruthEstimate.tsv',
        r'final_results.dir/\3.Salmon.CompareTruthEstimateIsodecoder.tsv',
        r'final_results.dir/\3.Salmon.CompareTruthEstimateAnticodon.tsv'])
def compareTruthEstimateSalmon(infiles, outfiles):

    outfile_individual, outfile_isodecoder, outfile_anticodon = outfiles

    job_options = PARAMS['cluster_options'] + " -t 1:00:00"
    job_condaenv=PARAMS['conda_base_env']

    CompareTrnaSeq.compareTruthEstimateSalmon(
        infiles, outfile_individual, outfile_isodecoder, outfile_anticodon, submit=True, job_options=job_options)

@mkdir('final_results.dir')
@collate((quantFractionalCounts, quantDiscreteCounts, quantDiscreteCountsMAPQ10,
          quantDiscreteCountsRandomSingle, quantDiscreteCountsNoMultimapping),
         regex('quant.dir/(\S+).(simulation_uniform|simulation_realistic)\.(\S+).(gene_count_.*).*'),
         add_inputs(r'simulations.dir/\1.\2.fastq.gz.gt'),
         [r'final_results.dir/\2.\4_Decision.CompareTruthEstimate.tsv',
          r'final_results.dir/\2.\4_Decision.CompareTruthEstimateIsodecoder.tsv',
          r'final_results.dir/\2.\4_Decision.CompareTruthEstimateAnticodon.tsv'])
def compareTruthEstimateDecisionCounts(infiles, outfiles):

    outfile_individual, outfile_isodecoder, outfile_anticodon = outfiles

    job_options = PARAMS['cluster_options'] + " -t 1:00:00"
    job_condaenv=PARAMS['conda_base_env']

    CompareTrnaSeq.compareTruthEstimateDecisionCounts(
        infiles, outfile_individual, outfile_isodecoder, outfile_anticodon, submit=True,  job_options=job_options)

@mkdir('final_results.dir')
@follows(mergeMimSeqIsodecoderMaps)
@transform((compareTruthEstimateDecisionCounts, compareTruthEstimateSalmon),
           regex('final_results.dir/(simulation_\S+)\.(Salmon|gene_count_\S+).CompareTruthEstimate.tsv'),
           add_inputs(r'quant.dir/\1_mimseq_isodecoder_maps.tsv'),
           r'final_results.dir/\1.\2.CompareTruthEstimateMimseqIsodecoder.tsv')
def makeMimseqIsodecoderQuant(infiles, outfile):
    infile = infiles[0][1]
    mapping_file = infiles[1]

    job_options = PARAMS['cluster_options'] + " -t 1:00:00"
    job_condaenv=PARAMS['conda_base_env']

    CompareTrnaSeq.makeMimseqIsodecoderQuant(infile, mapping_file, outfile, submit=True,  job_options=job_options)
# Function to merge comparisons for decision, salmon-based and mimseq quant

#####################################################
# Summary of multiple mapping and agreement wit truth
#####################################################
@mkdir('multiple_mapped_summary.dir')
@transform(alignWithBowtie2,
           regex('quant.dir/(\S+?)_(\S+).(simulation_uniform|simulation_realistic)\.(\S+).bam'),
           r'multiple_mapped_summary.dir/\1_\2.\3.\4.tsv')
def summariseMultimappedTruth2Assignment(infile, outfile):
    'Summarise reads as single/multi mapped and correct/incorrect at the anticodon level'

    job_options = PARAMS['cluster_options'] + " -t 2:00:00"
    job_condaenv=PARAMS['conda_base_env']

    CompareTrnaSeq.summariseMultimappedTruth2Assignment(infile, outfile, submit=True, job_options=job_options)


@mkdir('final_results.dir')
@collate(summariseMultimappedTruth2Assignment,
         regex('multiple_mapped_summary.dir/(\S+)\.(\d+)\.(simulation_\S+)\.(\S+)\.tsv'),
         r'final_results.dir/multiple_mapped_summary.\3.tsv')
def mergeSummariseMultimapped(infiles, outfile):

    job_options = PARAMS['cluster_options'] + " -t 1:00:00"
    job_condaenv=PARAMS['conda_base_env']

    CompareTrnaSeq.mergesummariseMultimapped(infiles, outfile, submit=True, job_options=job_options)

###################################################
# Alignment confusion matrices
###################################################



@follows(compareTruthEstimateMimseq)
@mkdir('truth2assignment.dir')
@transform(alignWithBowtie2,
           regex('quant.dir/(\S+?)_(\S+).(simulation_uniform|simulation_realistic)\.(\S+).bam'),
           add_inputs(r'mimseq.dir/\1.\3.MimseqCompareTruthEstimateMimseqIsodecoder.tsv.mapping'),
               r'truth2assignment.dir/\1_\2.\3.\4.truth2assignment.tsv')
def getTruth2Assignment(infiles, outfile):
    'Get the tally of ground truths to assignments'

    infile, isodecoder_mapping = infiles

    job_options = PARAMS['cluster_options'] + " -t 2:00:00"
    job_condaenv=PARAMS['conda_base_env']

    CompareTrnaSeq.getTruth2Assignment(infile, outfile, isodecoder_mapping, submit=True, job_options=job_options)


@follows(compareTruthEstimateMimseq)
@mkdir('truth2assignment.dir')
@transform((simulation_uniform, simulation_realistic),
           regex('simulations.dir/(\S+?)_(\S+).(simulation_uniform|simulation_realistic).fastq.gz'),
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
        infile, mimseq_isodecoder_counts, outfile, isodecoder_mapping, submit=True, job_options=job_options)

@mkdir('final_results.dir')
@collate((getTruth2Assignment, getTruth2AssignmentMimSeq),
         regex('truth2assignment.dir/(\S+)\.(\d+)\.(simulation_\S+)\.(\S+)\.truth2assignment.tsv'),
         r'final_results.dir/truth2assignment.\3.tsv')
def mergeTruth2Assignment(infiles, outfile):

    job_options = PARAMS['cluster_options'] + " -t 1:00:00"
    job_condaenv=PARAMS['conda_base_env']

    CompareTrnaSeq.mergeTruth2Assignment(infiles, outfile, submit=True, job_options=job_options)

    infiles_anticodon  = [re.sub('.tsv$', '_anticodon.tsv', x) for x in infiles]
    outfile_anticodon = re.sub('.tsv$', '_anticodon.tsv', outfile)
    CompareTrnaSeq.mergeTruth2Assignment(infiles_anticodon, outfile_anticodon, submit=True, job_options=job_options)

    infiles_isodecoder  = [re.sub('.tsv$', '_isodecoder.tsv', x) for x in infiles]
    # no isodecoder level output from mimseq
    infiles_isodecoder = [x for x in infiles_isodecoder if os.path.exists(x)]
    outfile_isodecoder = re.sub('.tsv$', '_isodecoder.tsv', outfile)
    CompareTrnaSeq.mergeTruth2Assignment(infiles_isodecoder, outfile_isodecoder, submit=True, job_options=job_options)

    infiles_mimseq_isodecoder  = [re.sub('.tsv$', '_mimseq_isodecoder.tsv', x) for x in infiles]
    outfile_mimseq_isodecoder = re.sub('.tsv$', '_mimseq_isodecoder.tsv', outfile)
    CompareTrnaSeq.mergeTruth2Assignment(infiles_mimseq_isodecoder, outfile_mimseq_isodecoder, submit=True, job_options=job_options)
####################################################




###################################################
# Targets
###################################################


@follows(simulation_uniform,
         simulation_realistic)
def simulate():
    'simulate tRNA reads from basic to more realistic'
    pass


@follows(summariseIndividualAlignments,
         summariseMergedAlignments,
         summariseMultimappedTruth2Assignment,
         mergeMutationProfileModomics)
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
         mergeTruth2AssignmentNoErrors,
         makeMimseqIsodecoderQuant,
         mergeTruth2Assignment)
def compare():
    'compare observed and ground truth'
    pass


@follows(concatenateEstimateDecisionCounts,
         concatenateEstimateSalmon)
def quantifyReal():
    '''quantify from the real raw data. This task by itself will not include
    mimseq quantification'''
    pass


@follows(quantifyReal,
         summariseAlignments,
         quant,
         compare)
def full():
    'run it all'
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
