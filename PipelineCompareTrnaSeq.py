import simulatetrna.simulateReads as simulateReads
import simulatetrna.alignmentSummary as alignmentSummary
import simulatetrna.bam as bam
from cgatcore.pipeline import cluster_runnable
import pysam
import pandas as pd
from collections import Counter, defaultdict
import re
import os
from Bio import SeqIO

def updateMtFastaNaming(infile, outfile):
    with open(outfile, 'w') as outf:
        for trna_record in SeqIO.parse(infile, "fasta"):
            _, species, __, aa, ac = trna_record.id.split('|')
            trna_record.id = '-'.join([species + '_MTtRNA', aa, ac, '1', '1'])
            SeqIO.write(trna_record, outf, "fasta")

@cluster_runnable
def mergeMimSeqIsodecoderMaps(infiles, outfile):

    isodecoder_maps = []
    for mapping_file in infiles:
        isodecoder_mapping = pd.read_csv(
            mapping_file + '.mapping', sep='\t', header=None, names=['isodecoder','mimseq_isodecoder'])
        isodecoder_mapping['trna_method'] = os.path.basename(mapping_file).split('.')[0]
        isodecoder_maps.append(isodecoder_mapping)

    isodecoder_maps = pd.concat(isodecoder_maps)

    isodecoder_maps.to_csv(outfile, sep='\t', index=False)

def cat_mimseq_compare_files(infiles, outfile):
    with open(outfile, 'w') as outf:
        out_header = True
        for fn in infiles:
            with open(fn, 'r') as inf:
                header = next(inf)
                if out_header:
                    outf.write(header)
                    out_header = False
                for line in inf:
                    outf.write(line)

@cluster_runnable
def mergeCompareTruthEstimateMimseq(infiles_isodecoder, outfile_isodecoder,
                                    infiles_anticodon, outfile_anticodon):

    cat_mimseq_compare_files(infiles_isodecoder, outfile_isodecoder)
    cat_mimseq_compare_files(infiles_anticodon, outfile_anticodon)

@cluster_runnable
def filter_sam(tmp_file, outfile):
    bam.filter_sam(tmp_file, outfile)

@cluster_runnable
def simulate_reads(infile,
                   outfile,
                   ground_truth,
                   error_rate,
                   mutation_threshold,
                   truncate,
                   alignment_summary,
                   summary_level,
                   outfile_gt=None):

    ground_truth = simulateReads.simulate_reads(
        infile, outfile, ground_truth, error_rate,
        mutation_threshold, truncate, alignment_summary, summary_level)

    if outfile_gt is not None:
        with open(outfile_gt, 'w') as out:
            for k,v in ground_truth.items():
                out.write('%s\t%s\n' % (k, v))

@cluster_runnable
def keep_random_alignment(tmp_file_sam, tmp_file_single_sam, outfile):
    bam.keep_random_alignment(tmp_file_sam, tmp_file_single_sam)
    pysam.sort(tmp_file_single_sam, "-o", outfile)
    pysam.index(outfile)

@cluster_runnable
def summariseAlignments(infile, trna_fasta, outfile):
    '''
    Use alignmentSummary.clustalwtrnaAlignmentSummary class to summarise the
    alignments and learn the mutational & truncation signatures
    '''

    alignment_summary = alignmentSummary.clustalwtrnaAlignmentSummary()

    alignment_summary.build(infile, trna_fasta)

    alignment_summary.pickle(outfile)


@cluster_runnable
def tally_read_counts(bam_infile, outfile_individual, outfile_isodecoder, outfile_anticodon,
                      random_single=False,
                      allow_multimapping=True,
                      min_mapq=0,
                      discrete=True):

    inbam = pysam.Samfile(bam_infile, 'r')

    read_counts = {'individual_sequence':Counter(),
                   'isodecoder':Counter(),
                   'anticodon':Counter()}

    for read_group in bam.iterate_reads(inbam, allow_multimapping):

        if min_mapq > 0:
            # MAPQ 255 is unassigned
            assignments = [read.reference_name for read in read_group if read.mapq>=min_mapq and read.mapq!=255]

        assignments = [read.reference_name for read in read_group if read.reference_name is not None]

        assignments = set(assignments)

        # No assignments for read after filtering
        if len(assignments) == 0:
            continue

        # Just one assigned gene or just taking one random assignment
        if len(assignments) == 1 or random_single:
            assignment = assignments.pop()
            read_counts['individual_sequence'][assignment] += 1
            read_counts['anticodon'][re.sub('\d$', '', '-'.join(assignment.split('-')[0:3]).replace('tRX', 'tRNA'))] +=1
            read_counts['isodecoder']['-'.join(assignment.split('-')[0:4])] +=1

        # multiple gene assignments
        else:
            anticodons = set([re.sub('\d$', '', '-'.join(x.split('-')[0:3]).replace('tRX', 'tRNA')) for x in assignments])
            isodecoders = set(['-'.join(x.split('-')[0:4]) for x in assignments])

            if discrete:
                if len(anticodons)==1:
                    read_counts['anticodon'][anticodons.pop()] +=1

                if len(isodecoders)==1:
                    read_counts['isodecoder'][isodecoders.pop()] +=1
            else:
                n_isodecoders = len(isodecoders)
                fractional_assignment_isodecoders = 1/n_isodecoders

                for isodecoder in isodecoders:
                    read_counts['isodecoder'][isodecoder] += fractional_assignment_isodecoders

                n_anticodons = len(anticodons)
                fractional_assignment_anticodons = 1/n_anticodons

                for anticodon in anticodons:
                    read_counts['anticodon'][anticodon] += fractional_assignment_anticodons

    with open(outfile_individual, 'w') as outf:
        for k,v in read_counts['individual_sequence'].items():
            outf.write(k + '\t' + str(v) + '\n')

    with open(outfile_isodecoder, 'w') as outf:
        for k,v in read_counts['isodecoder'].items():
            outf.write(k + '\t' + str(v) + '\n')

    with open(outfile_anticodon, 'w') as outf:
        for k,v in read_counts['anticodon'].items():
            outf.write(k + '\t' + str(v) + '\n')

def readMimSeqQuant(infile, drop_cols, melt_cols):

    'read a quant file from mimseq, drop unwanted columns, melt and annotate'

    mimseq_quant = pd.read_csv(infile, sep='\t')

    mimseq_quant = mimseq_quant.drop(
        drop_cols, axis=1).melt(
            id_vars=melt_cols, value_name='NumReads', var_name='filepath')

    mimseq_quant['input_file'] = [os.path.basename(x).split('.')[0] for x in mimseq_quant['filepath']]

    mimseq_quant['simulation_n'] = [os.path.basename(x).split('.')[1] for x in mimseq_quant['filepath']]

    mimseq_quant['quant_method']='mimseq'
    return(mimseq_quant)

@cluster_runnable
def compareMimSeq(infile, truths, isodecoder_out, anticodon_out):
    '''
    Compare mimseq quant with ground truth at isodecoder and anticodon level
    In both case, we need to update the ground truths to account for the internal working
    of mimseq
    '''
    all_truth_counts = []

    for truth_infile in truths:
        truth = pd.read_csv(truth_infile, sep='\t', header=None, names=['Name', 'truth'])
        truth['input_file'] = os.path.basename(truth_infile).split('.')[0]
        truth['simulation_n']= os.path.basename(truth_infile).split('.')[1]
        all_truth_counts.append(truth)

    truth_counts = pd.concat(all_truth_counts)

    # Isodecoder-level comparison
    mimseq_iso_quant = readMimSeqQuant(infile.replace('Anticodon','isodecoder'),
                                       drop_cols=['Single_isodecoder', 'parent', 'size'],
                                       melt_cols=['isodecoder'])

    # Mimseq reports some isodecoders as sum of multiple isodecoders
    # Below, we identify the mapping from isodecoder to mimseq isodecoder
    isodecoders2mimseqisodecoder = {}

    for isodecoder in mimseq_iso_quant['isodecoder']:
        isodecoder_ixs = isodecoder.split('-')[-1].split('/')

        if len(isodecoder_ixs)>1:
            isodecoder_prefix = '-'.join(isodecoder.split('-')[0:-1])
            for isodecoder_ix in isodecoder_ixs:
                if 'tRX' in isodecoder_ix:
                    isodecoder_ix = isodecoder_ix.replace('tRX', '')
                    key = isodecoder_prefix.replace('tRNA', 'tRX') + '-' + isodecoder_ix
                else:
                    key = isodecoder_prefix + '-' + isodecoder_ix

                isodecoders2mimseqisodecoder[key] = isodecoder
        else:
            isodecoders2mimseqisodecoder[isodecoder] = isodecoder

    # And then update the truth counts for isodecoders to use the mimseq isodecoder name
    truth_counts_isodecoder = truth_counts.copy()

    truth_counts_isodecoder['Name'] = ['-'.join(x.split('-')[0:4]) for x in truth_counts_isodecoder.Name]

    mimseq_isodecoder = []
    for x in truth_counts_isodecoder['Name']:
        if 'Und' in x:
            mimseq_isodecoder.append('Und')
            continue
        x = re.sub('\d-([ATCGN]{3})', r'-\1', x.replace('_MT', '_mito_').replace('tRX', 'tRNA'))
        if x in isodecoders2mimseqisodecoder:
            mimseq_isodecoder.append(isodecoders2mimseqisodecoder[x])
        else:
            mimseq_isodecoder.append(isodecoders2mimseqisodecoder[x.replace('tRNA', 'tRX')])


    truth_counts_isodecoder['mimseq_isodecoder'] = mimseq_isodecoder

    truth_counts_isodecoder = truth_counts_isodecoder.groupby(
        ['mimseq_isodecoder', 'input_file', 'simulation_n']).agg(
            {'truth':'sum'}).reset_index()

    mimseq_vs_truth_iso = mimseq_iso_quant.merge(
        truth_counts_isodecoder,
        right_on=['mimseq_isodecoder', 'input_file', 'simulation_n'],
        left_on=['isodecoder', 'input_file', 'simulation_n'])

    mimseq_vs_truth_iso = mimseq_vs_truth_iso.drop(
        ['filepath', 'mimseq_isodecoder'], axis=1).rename(
            columns={'isodecoder':'Name'})

    mimseq_vs_truth_iso.to_csv(isodecoder_out, sep='\t', index=False)

    # write out the isodecoder mapping to a file so we can use this to also
    # aggregate the other quantifications by the same isodecoder ids
    with open(isodecoder_out + '.mapping', 'w') as outf:
        truth_counts_isodecoders = set(['-'.join(x.split('-')[0:4]) for x in truth_counts.Name])
        for k, v in isodecoders2mimseqisodecoder.items():
            if k not in truth_counts_isodecoders:
                k = k.replace('tRX', 'tRNA').replace('_mito_', '_MT')
                if k not in truth_counts_isodecoders:
                    r = re.compile(re.sub('-([ATCGN]{3})', r'\\d-\1', k))
                    matches = list(filter(r.match, truth_counts_isodecoders))
                    if len(matches)==1:
                        k = matches[0]
                    else:
                        raise ValueError('More than one possible mimseq isodecoder to input sequence name!')

            outf.write('%s\t%s\n' % (k,v))

    # Anticodon-level comparison
    truth_counts_anticodon =  truth_counts.copy()

    truth_counts_anticodon['Name'] = ['-'.join(x.split('-')[0:3]) for x in truth_counts_anticodon.Name]


    # Mimseq internally replaces:
    # _MT -> _mito_
    # tRX -> tRNA
    # Below we update anticodon name in truth so we can merge
    truth_counts_anticodon['Name_mimseq'] = [
        re.sub('\d-([ATCGN]{3})', r'-\1', x.replace('_MT', '_mito_')) if '_MT' in x else x for x in truth_counts_anticodon['Name']]

    truth_counts_anticodon['Name_mimseq'] = [
        x.replace('tRX', 'tRNA') for x in truth_counts_anticodon['Name_mimseq']]

    truth_counts_anticodon = truth_counts_anticodon.groupby(
        ['Name_mimseq', 'input_file', 'simulation_n']).agg(
            {'truth':'sum'}).reset_index()

    mimseq_ac_quant = readMimSeqQuant(infile,
                                      drop_cols=['size'],
                                      melt_cols=['Anticodon'])

    mimseq_vs_truth_ac = mimseq_ac_quant.merge(truth_counts_anticodon,
                                               left_on=['Anticodon', 'input_file', 'simulation_n'],
                                               right_on=['Name_mimseq', 'input_file', 'simulation_n'])

    mimseq_vs_truth_ac = mimseq_vs_truth_ac.drop(
        ['filepath'], axis=1).rename(
            columns={'Anticodon':'Name'})

    mimseq_vs_truth_ac.to_csv(anticodon_out, sep='\t', index=False)



def mergeDecisionEstimateWithTruth(truth, estimate_infile, input_file, simulation_n, quant_method):
    estimate_counts = pd.read_csv(estimate_infile , sep='\t', header=None,
                                  names=['Name', 'NumReads'])
    counts_vs_truth = pd.merge(estimate_counts, truth, on='Name', how='right')
    counts_vs_truth['NumReads'] = counts_vs_truth['NumReads'].fillna(0)

    counts_vs_truth['input_file']=input_file
    counts_vs_truth['simulation_n']=simulation_n
    counts_vs_truth['quant_method']=quant_method
    return(counts_vs_truth)



@cluster_runnable
def getTruth2Assignment(infile, isodecoder_mapping, outfile):

    truth2assignment = defaultdict(Counter)

    for read_group in bam.iterate_reads(pysam.Samfile(infile)):

        assignments = set([read.reference_name for read in read_group if read.reference_name is not None])

        # all the read group queries are the same,
        # so just take one and strip the trailing '_x',
        # where x can be any length number
        truth = re.sub('_\d+$', '', read_group.pop().query_name)

        # If multiple assignments, split read equally
        n_assignments = len(assignments)
        fractional_assignment = 1/n_assignments

        for assignment in assignments:
            truth2assignment[truth][assignment]+=fractional_assignment

    with open(outfile, 'w') as outf:
        outf.write('\t'.join(('truth', 'assignment', 'count')) + '\n')
        for k,v in truth2assignment.items():
            for k2, c in v.items():
                outf.write('\t'.join(map(str, (k, k2, c))) + '\n')

    t2a_df = t2a_df = pd.read_csv(outfile, sep='\t')
    t2a_df['truth_anticodon'] = [re.sub('\d$', '', '-'.join(x.split('-')[0:3]).replace('tRX', 'tRNA')) for x in t2a_df['truth']]
    t2a_df['assignment_anticodon'] = [re.sub('\d$', '', '-'.join(x.split('-')[0:3]).replace('tRX', 'tRNA')) for x in t2a_df['assignment']]

    t2a_df_anticodon = t2a_df.groupby(['truth_anticodon', 'assignment_anticodon']).agg({'count':'sum'}).reset_index()
    t2a_df_anticodon.to_csv(re.sub('.tsv$', '_anticodon.tsv', outfile), sep='\t')

    t2a_df['truth_isodecoder'] = ['-'.join(x.split('-')[0:4]) for x in t2a_df['truth']]
    t2a_df['assignment_isodecoder'] = ['-'.join(x.split('-')[0:4]) for x in t2a_df['assignment']]

    t2a_df_isodecoder = t2a_df.groupby(['truth_isodecoder', 'assignment_isodecoder']).agg({'count':'sum'}).reset_index()
    t2a_df_isodecoder.to_csv(re.sub('.tsv$', '_isodecoder.tsv', outfile), sep='\t')

    # read in the mimseq isodecoder mapping file to map from each tRNA to the mimseq isodecoder
    trna2mimseqcluster = pd.read_csv(isodecoder_mapping, sep='\t', header=None, names=['trna', 'mimseq_cluster'])
    trna2isodecoder = defaultdict(str, {t:i for t,i in zip(trna2mimseqcluster['trna'], trna2mimseqcluster['mimseq_cluster'])})

    with open(re.sub('.tsv$', '_mimseq_isodecoder.tsv', outfile), 'w') as outf:
        for x,y in trna2isodecoder.items():
            outf.write("%s\t%s\n" % (x, y))

    t2a_df_isodecoder['truth_isodecoder'] = [trna2isodecoder[x] for x in t2a_df_isodecoder['truth_isodecoder']]
    t2a_df_isodecoder['assignment_isodecoder'] = [trna2isodecoder[x] for x in t2a_df_isodecoder['assignment_isodecoder']]

    t2a_df_mimseq_isodecoder = t2a_df_isodecoder.groupby(['truth_isodecoder', 'assignment_isodecoder']).agg({'count':'sum'}).reset_index()
    t2a_df_mimseq_isodecoder.to_csv(re.sub('.tsv$', '_mimseq_isodecoder.tsv', outfile), sep='\t', index=False)



@cluster_runnable
def getTruth2AssignmentMimSeq(infile, mimseq_isodecoder_counts, isodecoder_mapping, outfile):

    truth2assignment = defaultdict(Counter)

    for read in pysam.Samfile(infile):
        assignment = read.reference_name

        if read.reference_name is None:
            continue

        # all the read group queries are the same,
        # so just take one and strip the trailing '_x',
        # where x can be any length number
        truth = re.sub('-\d+_\d+$', '', read.query_name)

        truth2assignment[truth][assignment] += 1

    with open(outfile, 'w') as outf:
        outf.write('\t'.join(('truth', 'assignment', 'count')) + '\n')
        for k,v in truth2assignment.items():
            for k2, c in v.items():
                outf.write('\t'.join(map(str, (k, k2, c))) + '\n')

    t2a_df = pd.read_csv(outfile, sep='\t')

    t2a_df['truth_anticodon'] = [
        re.sub('\d$', '', '-'.join(x.split('-')[0:3]).replace('tRX', 'tRNA')) for x in t2a_df['truth']]

    t2a_df['assignment_anticodon'] = [
        re.sub('\d$', '', '-'.join(x.split('-')[0:3]).replace('tRX', 'tRNA').replace('mito_', 'MT')) for x in t2a_df['assignment']]

    t2a_df_anticodon = t2a_df.groupby(['truth_anticodon', 'assignment_anticodon']).agg({'count':'sum'}).reset_index()
    t2a_df_anticodon.to_csv(re.sub('.tsv$', '_anticodon.tsv', outfile), sep='\t', index=False)

    mimseq_isodecoder_table = pd.read_csv(mimseq_isodecoder_counts, sep='\t')
    parent2isodecoder = {p:i for i,p in zip(mimseq_isodecoder_table['isodecoder'], mimseq_isodecoder_table['parent'])}

    # read in the mimseq isodecoder mapping file to map from each tRNA to the mimseq isodecoder
    trna2mimseqcluster = pd.read_csv(isodecoder_mapping, sep='\t', header=None, names=['trna', 'mimseq_cluster'])
    trna2isodecoder = defaultdict(str, {t:i for t,i in zip(trna2mimseqcluster['trna'], trna2mimseqcluster['mimseq_cluster'])})

    with open(re.sub('.tsv$', '_mimseq_isodecoder.tsv', outfile), 'w') as outf:
        for x,y in trna2isodecoder.items():
            outf.write("%s\t%s\n" % (x,y))

    t2a_df['truth_isodecoder'] = [trna2isodecoder[x] for x in t2a_df['truth']]
    t2a_df['assignment_isodecoder'] = [parent2isodecoder[re.sub('-\d$', '', x)] for x in t2a_df['assignment']]

    t2a_df_isodecoder = t2a_df.groupby(['truth_isodecoder', 'assignment_isodecoder']).agg({'count':'sum'}).reset_index()
    t2a_df_isodecoder.to_csv(re.sub('.tsv$', '_mimseq_isodecoder.tsv', outfile), sep='\t', index=False)


@cluster_runnable
def mergeTruth2Assignment(infiles, outfile):
    dfs = []
    for infile in infiles:
        input_file = os.path.basename(infile).split('.')[0]
        simulation_n = os.path.basename(infile).split('.')[1]
        quant_method = os.path.basename(infile).split('.')[3]

        df = pd.read_csv(infile, sep='\t')
        df['input_file']=input_file
        df['simulation_n']=simulation_n
        df['quant_method']=quant_method

        dfs.append(df)

    df = pd.concat(dfs)
    df.to_csv(outfile, sep='\t')

@cluster_runnable
def compareTruthEstimateSalmon(
    infiles, outfile_individual, outfile_isodecoder, outfile_anticodon):

    all_counts_vs_truth = []

    for infile_pair in infiles:

        estimate, truth = infile_pair

        estimate_counts = pd.read_csv(estimate, sep='\t')[['Name', 'NumReads']]
        truth_counts = pd.read_csv(truth, sep='\t', header=None,
                                   names=['Name', 'truth'])

        counts_vs_truth = pd.merge(estimate_counts, truth_counts, on='Name', how='right')
        counts_vs_truth['NumReads'] = counts_vs_truth['NumReads'].fillna(0)

        input_file = os.path.basename(truth).split('.')[0]
        counts_vs_truth['input_file']=input_file

        simulation_n = os.path.basename(truth).split('.')[1]
        counts_vs_truth['simulation_n']=simulation_n

        quant_method = os.path.dirname(estimate).split('.')[-1] + '_salmon'
        counts_vs_truth['quant_method']=quant_method

        all_counts_vs_truth.append(counts_vs_truth)

    all_counts_vs_truth = pd.concat(all_counts_vs_truth)
    all_counts_vs_truth.to_csv(outfile_individual, sep='\t', index=False)

    all_counts_vs_truth['Name'] = ['-'.join(x.split('-')[0:4]) for x in all_counts_vs_truth.Name]
    all_counts_vs_truth_isodecoder = all_counts_vs_truth.groupby(
        ['Name', 'simulation_n', 'quant_method', 'input_file']).agg(
            {'NumReads':'sum', 'truth':'sum'}).reset_index()
    all_counts_vs_truth_isodecoder.to_csv(outfile_isodecoder, sep='\t', index=False)

    all_counts_vs_truth['Name'] = [re.sub('\d$', '', '-'.join(x.split('-')[0:3]).replace('tRX', 'tRNA')) for x in all_counts_vs_truth.Name]
    all_counts_vs_truth_ac = all_counts_vs_truth.groupby(
        ['Name', 'simulation_n', 'quant_method', 'input_file']).agg(
            {'NumReads':'sum', 'truth':'sum'}).reset_index()

    all_counts_vs_truth_ac.to_csv(outfile_anticodon, sep='\t', index=False)


@cluster_runnable
def compareTruthEstimateDecisionCounts(infiles, outfile_individual, outfile_isodecoder, outfile_anticodon):
    all_counts_vs_truth = {'individual':[], 'isodecoder':[], 'anticodon':[]}

    for infile_set in infiles:

        estimates, truth = infile_set

        estimate_individual, estimate_isodecoder, estimate_anticodon = estimates

        input_file = os.path.basename(truth).split('.')[0]
        simulation_n = os.path.basename(truth).split('.')[1]
        quant_method = os.path.basename(estimate_individual).split('.')[-2] + '_decision'

        if 'mapq10' in estimate_individual:
            quant_method += '_mapq10'

        truth_counts = pd.read_csv(truth, sep='\t', header=None,
                                   names=['Name', 'truth'])

        counts_vs_truth = mergeDecisionEstimateWithTruth(
            truth_counts, estimate_individual, input_file, simulation_n, quant_method)

        all_counts_vs_truth['individual'].append(counts_vs_truth)

        truth_counts['Name'] = ['-'.join(x.split('-')[0:4]) for x in truth_counts.Name]
        truth_counts_isodecoder = truth_counts.groupby(
            ['Name']).agg(
                {'truth':'sum'}).reset_index()
        counts_vs_truth_isodecoder = mergeDecisionEstimateWithTruth(
            truth_counts_isodecoder, estimate_isodecoder, input_file, simulation_n, quant_method)
        all_counts_vs_truth['isodecoder'].append(counts_vs_truth_isodecoder)

        truth_counts['Name'] = [re.sub('\d$', '', '-'.join(x.split('-')[0:3]).replace('tRX', 'tRNA')) for x in truth_counts.Name]
        truth_counts_anticodon = truth_counts.groupby(
            ['Name']).agg(
                {'truth':'sum'}).reset_index()
        counts_vs_truth_anticodon = mergeDecisionEstimateWithTruth(
            truth_counts_anticodon, estimate_anticodon, input_file, simulation_n, quant_method)
        all_counts_vs_truth['anticodon'].append(counts_vs_truth_anticodon)


    all_counts_vs_truth_individual = pd.concat(all_counts_vs_truth['individual'])
    all_counts_vs_truth_individual.to_csv(outfile_individual, sep='\t', index=False)

    all_counts_vs_truth_isodecoder = pd.concat(all_counts_vs_truth['isodecoder'])
    all_counts_vs_truth_isodecoder.to_csv(outfile_isodecoder, sep='\t', index=False)

    all_counts_vs_truth_anticodon = pd.concat(all_counts_vs_truth['anticodon'])
    all_counts_vs_truth_anticodon.to_csv(outfile_anticodon, sep='\t', index=False)


@cluster_runnable
def makeMimseqIsodecoderQuant(infile, mapping_file, outfile):

    isodecoder_mapping = pd.read_csv(mapping_file, sep='\t')

    isodecoder_quant = pd.read_csv(infile, sep='\t')

    isodecoder_quant['trna_method'] = [x.split('_')[0] for x in isodecoder_quant['input_file']]

    isodecoder_quant = isodecoder_quant.merge(isodecoder_mapping,
                                              left_on=['Name', 'trna_method'],
                                                  right_on=['isodecoder', 'trna_method'])

    mimseq_isodecoder_quant = isodecoder_quant.groupby(
            ['mimseq_isodecoder', 'input_file', 'simulation_n', 'quant_method']).agg(
                {'NumReads':'sum', 'truth':'sum'}).reset_index().rename(
                    columns={'mimseq_isodecoder':'Name'})

    mimseq_isodecoder_quant.to_csv(outfile, sep='\t', index=False)
