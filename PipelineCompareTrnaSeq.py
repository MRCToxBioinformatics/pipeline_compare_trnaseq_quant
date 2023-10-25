import simulatetrna.simulateReads as simulateReads
import simulatetrna.alignmentSummary as alignmentSummary
import simulatetrna.bam as bam

from cgatcore.pipeline import cluster_runnable

import pysam

import numpy as np
import pandas as pd

from collections import Counter, defaultdict
from itertools import permutations

import re
import os
import pickle

from Bio import SeqIO
from Bio import pairwise2
from Bio.Seq import Seq

import networkx as nx
from networkx.algorithms import community

import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.colors import LinearSegmentedColormap


@cluster_runnable
def getGraph(bamfile, nx_graph_outfile, edge_weights_outfile, edge_weights_table_outfile):
    '''
    bamfile = BAM with all top scoring alignments for each read, sorted in read order. Can
    be obtained by e.g using -a with bowtie2 and then using simulatetrna.bam.filtersam()

    nx_graph_outfile = Where to save the nx.graph

    edge_weights_outfile = where to save the dictionary of weighted edges
    '''

    # Create dictionaries to hold the nx.Graphs and node and egde counts
    G = {'gid':nx.Graph(),
         'tid':nx.Graph(),
         'ac':nx.Graph(),
         'aa':nx.Graph()}

    # Calculate the edge weights
    node_counts = node_counts = {
        'gid':Counter(),
        'tid':Counter(),
        'ac':Counter(),
        'aa':Counter()}

    edge_weights = {
        'gid':Counter(),
        'tid':Counter(),
        'ac':Counter(),
        'aa':Counter()}

    inreads = pysam.Samfile(bamfile, 'r')

    # Iterate over alignments
    for alignments in bam.iterate_reads(inreads, allow_multimapping=True, remove_trx=True, remove_mt=True):

        alignments_reference = {
            'gid': set([read.reference_name.split('_')[-1] for read in alignments]),
            'tid': set(['-'.join(read.reference_name.split('_')[-1].split('-')[0:4]) for read in alignments]),
            'ac': set(['-'.join(read.reference_name.split('_')[-1].split('-')[0:3]) for read in alignments]),
            'aa': set(['-'.join(read.reference_name.split('_')[-1].split('-')[0:2]) for read in alignments])
         }
        
        
        for level in alignments_reference.keys():
            
            alignments = alignments_reference[level]
        
            if len(alignments) == 1:
                # Single alignment, add the read to the corresponding node
                # read = next(iter(alignments))

                node = list(alignments)[0]

                node_counts[level][node] += 1

                G[level].add_node(node)

            else:
                # Multiple alignments, add each read to each edge for all pairs of nodes

                for ref1, ref2 in permutations(alignments, 2):

                    edge_key = frozenset((ref1, ref2))

                    # Check if the edge key exists in the edge_weights dictionary

                    edge_weights[level][edge_key] += 1

                    G[level].add_edge(ref1, ref2)
                
    # recalculate the weights using the node counts. Weights are now bounded 0-1
    recalculated_weights = {}
    for level in edge_weights.keys():
        recalculated_weights[level] = {}
        for edge, weight in edge_weights[level].items():
            node1, node2 = edge
            recalculated_weights[level][edge] = weight / (
                weight + node_counts[level][node1] + node_counts[level][node2])


    pickle.dump(recalculated_weights, open(edge_weights_outfile, 'wb'))
    pickle.dump(G, open(nx_graph_outfile, 'wb'))
    
    with open(edge_weights_table_outfile, 'w') as outf:
        outf.write('bam\tlevel\tnode1\tnode2\tnode1_count\tnode2_count\tedge_weight\tedge_reweighted\n')
        for level in edge_weights.keys():
            for edge, weight in edge_weights[level].items():
                node1, node2 = edge
                outf.write('\t'.join(map(str, (bamfile,
                                               level,
                                               node1,
                                               node2,
                                               node_counts[level][node1],
                                               node_counts[level][node2],
                                               weight,
                                               recalculated_weights[level][edge]))) + '\n')


def custom_cmap():
    # Create a custom colormap that transitions from white (0.0) to blue (1.0)
    cmap_colors = [(1, 1, 1), (0, 0, 1)]
    return LinearSegmentedColormap.from_list("custom_cmap", cmap_colors)


@cluster_runnable
def plot_heatmap(edge_weights_infile,
                 nxGraph_infile,
                 output_filename,
                 level='ac',
                 dividing_lines=True,
                 cmap_col=custom_cmap()):

    '''
    Plot a heatmap for the multimapping between tRNAs using the output from getGraph

    edge_weights = normalised edge weights
    output_filename = filepath for outfile
    level = multimapping level to plot (aa=aminoacid, ac=anticodon, tid=transcriptID, gid=genomeLocusID)
    '''

    edge_weights = pickle.load(open(edge_weights_infile, 'rb'))
    G = pickle.load(open(nxGraph_infile, 'rb'))

    allowed_levels = ['aa', 'ac', 'gid', 'tid'] 

    if level not in allowed_levels:
        raise ValueError('level must be one of %s' % ', '.join(allowed_levels))

    # Get the list of nodes and their order in alphabetical order
    nodes = sorted(list(G[level].nodes()))
    num_nodes = len(nodes)

    # Create an empty matrix to store the edge weights
    edge_matrix = np.zeros((num_nodes, num_nodes))

    # Populate the edge matrix with the logarithmic transformed edge weights
    for edge in G[level].edges():
        node1 = nodes.index(edge[0])
        node2 = nodes.index(edge[1])
        recalc_edge_weight = edge_weights[level].get(frozenset(edge), 0)
        edge_matrix[node1, node2] = np.log1p(recalc_edge_weight)  # Apply logarithmic transformation

    # Get the list of amino acid labels from the nodes
    amino_acids = [node.split('-')[1] for node in nodes]

    # Remove 'eColiLys' amino acid if present
    if 'eColiLys' in amino_acids:
        idx = amino_acids.index('eColiLys')
        amino_acids.pop(idx)
        edge_matrix = np.delete(edge_matrix, idx, axis=0)
        edge_matrix = np.delete(edge_matrix, idx, axis=1)

    
    # Get unique amino acids and their corresponding indices
    unique_amino_acids, amino_acid_indices = np.unique(amino_acids, return_inverse=True)

    
    num_unique_amino_acids = len(unique_amino_acids)

    # Create the heatmap plot with the custom colormap
    plt.figure(figsize=(7, 7))
    heatmap = plt.imshow(edge_matrix, cmap=cmap_col, interpolation='nearest', vmin=0.0, vmax=1.0)

    # Add color scale legend and adjust its size using the 'shrink' parameter
    plt.colorbar(heatmap, shrink=0.5)

    # Find the indices where the amino acid groups change
    group_change_indices = np.where(np.diff(amino_acid_indices) != 0)[0]

    if dividing_lines:
    
        # Loop through the group change indices and draw dividing lines
        for idx in group_change_indices:
            plt.axhline(y=idx+0.5, color='grey', linewidth=0.1, alpha=0.2)
            plt.axvline(x=idx+0.5, color='grey', linewidth=0.1, alpha=0.2)

        if level != 'aa':
            # Get the list of amino acid labels from the nodes
            ac = [node.split('-')[2] for node in nodes]

            # Get unique amino acids and their corresponding indices
            unique_ac, ac_indices = np.unique(ac, return_inverse=True)

            # Find the indices where the amino acid groups change
            ac_change_indices = np.where(np.diff(ac_indices) != 0)[0]

            # Loop through the group change indices and draw dividing lines
            for idx in ac_change_indices:
                plt.axhline(y=idx+0.5, color='lightgrey', linewidth=0.05, alpha=0.2)
                plt.axvline(x=idx+0.5, color='lightgrey', linewidth=0.05, alpha=0.2)

        group_change_indices_inc_start_end = np.append(np.append(-1, group_change_indices), len(amino_acid_indices))

        start_end = [x for x in zip(group_change_indices_inc_start_end[0::1], group_change_indices_inc_start_end[1::1])]

        for start, end in start_end:
            start = start+0.5
            end = end+0.5
            #Create a Rectangle patch
            rect = patches.Rectangle((start, end), end-start, start-end,
                                     edgecolor='darkgrey', linewidth=1, alpha=1, facecolor='none')

            # Add the patch to the Axes
            plt.gca().add_patch(rect)

    # Add amino acid labels to the x and y axes
    plt.xticks(np.arange(num_unique_amino_acids), [])
    plt.yticks(np.arange(num_unique_amino_acids), [])

    previous = None
    previous_ix = None

    # Loop through unique amino acids and add labels at their midpoints on x and y axes
    previous = amino_acid_indices[0]
    previous_ix = [0]
    i = 0
    for ix, aai in enumerate(amino_acid_indices[1:], start=1):

        if aai == previous:
            previous_ix.append(ix)
        else:
            current_index = np.mean(previous_ix)  # Midpoint of assigned indexes
            amino_acid = unique_amino_acids[i]
            plt.text(current_index, len(edge_matrix)+1, amino_acid, ha='center', va='top', fontsize=8, color='black', rotation=90)
            plt.text(-2, current_index, amino_acid, ha='right', va='center', fontsize=7, color='black')

            i += 1
            
            previous = aai
            previous_ix = [ix]

    current_index = np.mean(previous_ix)  # Midpoint of assigned indexes
    amino_acid = unique_amino_acids[i]
    plt.text(current_index, len(edge_matrix)+1, amino_acid, ha='center', va='top', fontsize=8, color='black', rotation=90)
    plt.text(-2, current_index, amino_acid, ha='right', va='center', fontsize=7, color='black')            

    # Remove ticks from both the x and y axes
    plt.tick_params(axis='both', which='both', length=0)

    # Save the figure with higher resolution
    plt.savefig(output_filename, format='pdf', dpi=500)

    plt.close()


def filterFasta(infile, outfile):
    with open(outfile, 'w') as outf:
        for trna_record in SeqIO.parse(infile, "fasta"):
            if len(trna_record.seq) < 100:
                SeqIO.write(trna_record, outf, "fasta")


def updateMtFastaNaming(infile, outfile):
    with open(outfile, 'w') as outf:
        for trna_record in SeqIO.parse(infile, "fasta"):
            _, species, __, aa, ac = trna_record.id.split('|')
            trna_record.id = '-'.join([species + '_MTtRNA', aa, ac, '1', '1'])
            SeqIO.write(trna_record, outf, "fasta")


@cluster_runnable
def mapModomics2fasta(fasta_infile,
                      modomics_json,
                      modification_index,
                      outfile,
                      alignment_score_threshold = 65):
    '''
    Take the Modomics modifications and map them to the tRNA sequences

    For this, we need the MODOMICS modifications (in json format) and
    a file mapping the MODOMICS code to nucleotide (modification_index) and,
    finally, the fasta file of tRNA sequences
    '''

    mod_code_to_base = {}
    mod_code_to_mod_name = {}

    with open(modification_index, 'r') as inf:
        header = inf.readline()
        for line in inf:
            mod_code, name, short_name, code, moiety_type = [x for x in line.strip('"').split('","')]
            if code == '""':
                code = '"'
            if code != '':
                mod_code_to_base[code] = mod_code[-1:]
                mod_code_to_mod_name[code] = short_name

    sequences_df = pd.read_json(modomics_json).transpose()

    fasta_sequences = SeqIO.parse(open(fasta_infile),'fasta')

    anticodon_sequences = defaultdict(
        lambda: defaultdict(
            lambda: defaultdict(list)))

    for entry in fasta_sequences:
        entry.seq = entry.seq.upper()
        species, amino_acid, anticodon, transcript_id, genome_id = entry.name.split('-')
        species = species.replace('_tRNA', '').replace('_', ' ')
        anticodon_sequences[species][amino_acid][anticodon].append(entry)

    outf = open(outfile, 'w')
    outf.write('\t'.join(('modomics_position', 'fasta_position', 'modification', 'fasta_entry')) + '\n')

    for row in sequences_df.itertuples():
        if re.search('mitochondrion', row.organellum):
            continue

        #print(row.subtype, row.anticodon)
        sequence_no_mod = ''
        mod_positions = {}
        for pos, base in enumerate(row.seq):
            if base not in ['C', 'A', 'G', 'U']:
                if base not in mod_code_to_base:
                    raise ValueError(sprintf(
                        'No mapping available from character to modification: %s' % base))
                mod_positions[pos] = mod_code_to_mod_name[base]
                base = mod_code_to_base[base]
            sequence_no_mod+=base

        keep_alignments = []

        for anticodon in anticodon_sequences[row.organism][row.subtype]:
            for entry in anticodon_sequences[row.organism][row.subtype][anticodon]:
                seq1 = Seq(sequence_no_mod)
                seq2 = Seq(str(entry.seq).replace('T', 'U'))

                # Finding similarities
                try:
                    alignment = pairwise2.align.localms(seq1, seq2, 1, -1, -1, -.1, one_alignment_only=True)[0]

                except:
                    print(seq1)
                    print(seq2)
                    raise ValueError()

                if alignment.score > alignment_score_threshold:
                    keep_alignments.append((entry, alignment))


        for keep_alignment in keep_alignments:

            keep_entry, keep_alignment = keep_alignment

            mod_pos = 0
            fasta_pos = 0
            mod_to_fasta_pos = {}
            for mod_seq_base, fasta_seq_base in zip(keep_alignment.seqA,
                                                    keep_alignment.seqB):
                if(mod_seq_base==fasta_seq_base):
                    mod_to_fasta_pos[mod_pos] = fasta_pos
                    mod_pos +=1
                    fasta_pos +=1
                elif(fasta_seq_base=='-'):
                    mod_pos +=1
                elif(mod_seq_base=='-'):
                    fasta_pos +=1
                else:
                    mod_pos +=1
                    fasta_pos +=1

            for mod_position in mod_positions:
                if(mod_position in mod_to_fasta_pos):
                    outf.write('\t'.join(map(str,
                        [mod_position,
                         mod_to_fasta_pos[mod_position],
                         mod_positions[mod_position],
                         keep_entry.name])) + '\n')

    outf.close()

@cluster_runnable
def mergeMutationProfileModomics(pickle_infiles, modomics_positions, outfile_mutations, outfile_read_end, outfile_truncations):

    mutation_df = []
    start_df = []

    for infile in pickle_infiles:

        alignment_summary = pickle.load(open(infile, 'rb'))

        sample_name = os.path.basename(infile).replace('.merged.summariseAlignments.pickle', '')
        quant_method = sample_name.split('_')[0]
        sample = '_'.join(sample_name.split('_')[1:]).replace('Hsap_', '')

        mutations = []
        starts = []

        for trna in alignment_summary.contig_base_frequencies:
            #if trna != toi: # save runtime & memory
            #    continue
            for position, all_base_counts in alignment_summary.contig_base_frequencies[trna].items():
                for ref_base, base_counts in all_base_counts.items():
                    total_counts = sum(base_counts.values())
                    if total_counts>0:
                        for obs_base, count in base_counts.items():
                            if(obs_base != ref_base):
                                mutations.append([sample, quant_method, trna,
                                                  position, ref_base, obs_base,
                                                  count, total_counts,
                                                  count/total_counts])
                for start_counts in alignment_summary.alignment_coordinates[trna].values():
                    for start, count in start_counts.items():
                        starts.append([sample, quant_method, trna, start, count])

        mutation_df.append(pd.DataFrame.from_records(
            mutations, columns=['sample', 'quant_method', 'trna', 'position', 'ref_base',
                                'obs_base', 'count', 'total_counts', 'mutation_rate']))

        start_df.append(pd.DataFrame.from_records(
            starts, columns=['sample', 'quant_method', 'trna', 'start', 'frequency']))

    ###
    # merge the mutations and modomics positions
    mutation_df = pd.concat(mutation_df)
    
    # can be multiple rows with the same modification information
    modomics_modifications = pd.read_table(modomics_positions).drop_duplicates()

    keep_trnas = set(modomics_modifications['fasta_entry'])
    mutation_df_with_modomics = mutation_df[[x in keep_trnas for x in mutation_df['trna']]]

    mutation_df_with_modomics_modifications = pd.merge(mutation_df_with_modomics,
                 modomics_modifications,
                 left_on=['position', 'trna'],
                 right_on=['fasta_position', 'fasta_entry'],
                 how='left')

    mutation_df_with_modomics_modifications.to_csv(outfile_mutations, sep='\t', index=False)
    ###

    ###
    # merge the start positions (end of read) and modomics positions
    start_df = pd.concat(start_df)

    start_df_summary = start_df.groupby(['quant_method', 'start']).agg(
        total_frequency=pd.NamedAgg(column='frequency', aggfunc=sum)).reset_index()

    start_df_summary.to_csv(outfile_read_end, sep='\t', index=False)


    modomics_modifications['0'] = modomics_modifications['fasta_position']
    modomics_modifications['-1'] = [x+1 for x in modomics_modifications['0']]
    modomics_modifications['+1'] = [x-1 for x in modomics_modifications['0']]
    
    modomics_modifications_plus_minus_one = pd.melt(
        modomics_modifications[['fasta_entry', 'modification', '0', '-1', '+1']],
        id_vars=['fasta_entry', 'modification'],
        value_vars=['0', '-1', '+1'],
        var_name='relative_to_mod_position',
        value_name='fasta_position')

    start_df = start_df.groupby(['sample', 'quant_method', 'trna', 'start']).agg({'frequency' : 'sum'}).reset_index()

    start_df_with_modomics = start_df[[x in keep_trnas for x in start_df['trna']]]

    start_df_with_modomics_modifications = pd.merge(
        start_df_with_modomics,
        modomics_modifications_plus_minus_one,
        left_on=['start', 'trna'],
        right_on=['fasta_position', 'fasta_entry'],
        how='left')

    start_df_with_modomics_modifications.to_csv(outfile_truncations, sep='\t', index=False)
    ###
    
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
                   alignment_summary=None,
                   summary_level='anticodon',
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

        assignments = [read for read in read_group if read.reference_name is not None]

        if min_mapq > 0:
            # MAPQ 255 is unassigned
            assignments = [read.reference_name for read in read_group if read.mapq>=min_mapq and read.mapq!=255]
        else:
            assignments = [read.reference_name for read in read_group]

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
    mimseq_quant['tally_method']='mimseq'
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
    mimseq_iso_quant = readMimSeqQuant(infile.replace('Anticodon','Isodecoder'),
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
        x = re.sub('\d-([ATCGN]{3})', r'-\1', x.replace('_MT', '_mito_'))#.replace('tRX', 'tRNA'))
        if x in isodecoders2mimseqisodecoder:
            mimseq_isodecoder.append(isodecoders2mimseqisodecoder[x])
        else:
            try:
                mimseq_isodecoder.append(isodecoders2mimseqisodecoder[x.replace('tRNA', 'tRX')])
            except:
                mimseq_isodecoder.append(x)


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
                k = k.replace('tRX', 'tRNA').replace('-tRNA', '-').replace('_mito_', '_MT')
                if k not in truth_counts_isodecoders:
                    r = re.compile(re.sub('-([ATCGN]{3})', r'\\d-\1', k))
                    matches = list(filter(r.match, truth_counts_isodecoders))
                    if len(matches)==1:
                        k = matches[0]
                    else:
                        raise ValueError('Not just one possible mimseq isodecoder to input sequence name! %s : %s' %(k, ','.join(matches)))

            outf.write('%s\t%s\n' % (k,v))

    # Anticodon-level comparison
    truth_counts_anticodon =  truth_counts.copy()

    truth_counts_anticodon['Name'] = ['-'.join(x.split('-')[0:3]) for x in truth_counts_anticodon['Name']]

    truth_counts_anticodon['Name'] = [
        re.sub('\d-([ATCGN]{3})', r'-\1', x.replace('tRX', 'tRNA')) for x in truth_counts_anticodon['Name']]

    truth_counts_anticodon = truth_counts_anticodon.groupby(
        ['Name', 'input_file', 'simulation_n']).agg(
            {'truth':'sum'}).reset_index()

    mimseq_ac_quant = readMimSeqQuant(infile,
                                      drop_cols=['size'],
                                      melt_cols=['Anticodon'])
    # Mimseq internally replaces:
    # _MT -> _mito_
    # tRX -> tRNA
    # Below we update anticodon name in mimseq output so we can merge
    mimseq_ac_quant['Anticodon'] = [x.replace('_mito_', '_MT') for x in mimseq_ac_quant['Anticodon']]

    mimseq_vs_truth_ac = truth_counts_anticodon.merge(mimseq_ac_quant,
                                               right_on=['Anticodon', 'input_file', 'simulation_n'],
                                               left_on=['Name', 'input_file', 'simulation_n'])


    mimseq_vs_truth_ac = mimseq_vs_truth_ac.drop(['filepath', 'Anticodon'], axis=1)

    mimseq_vs_truth_ac.to_csv(anticodon_out, sep='\t', index=False)



def mergeDecisionEstimateWithTruth(truth, estimate_infile, input_file, simulation_n, quant_method, tally_method):
    estimate_counts = pd.read_csv(estimate_infile , sep='\t', header=None,
                                  names=['Name', 'NumReads'])
    counts_vs_truth = pd.merge(estimate_counts, truth, on='Name', how='right')
    counts_vs_truth['NumReads'] = counts_vs_truth['NumReads'].fillna(0)

    counts_vs_truth['input_file']=input_file
    counts_vs_truth['simulation_n']=simulation_n
    counts_vs_truth['quant_method']=quant_method
    counts_vs_truth['tally_method']=tally_method
    return(counts_vs_truth)



@cluster_runnable
def summariseMultimappedTruth2Assignment(infile, outfile):

    # events keys:values are assignments:anticodons:correct:count
    events = defaultdict(lambda: defaultdict(Counter))

    for read_group in bam.iterate_reads(pysam.Samfile(infile)):

        assignments = set([read.reference_name for read in read_group if read.reference_name is not None])

        # all the read group queries are the same,
        # so just take one and strip the trailing '_x',
        # where x can be any length number
        truth = re.sub('_\d+$', '', read_group.pop().query_name)
        truth_anticodon = re.sub('\d$', '', '-'.join(truth.split('-')[0:3]).replace('tRX', 'tRNA'))

        assignment_anticodons = [re.sub('\d$', '', '-'.join(x.split('-')[0:3]).replace('tRX', 'tRNA')) for x in assignments]

        if len(assignments) == 1:
            assignment_anticodon = assignment_anticodons.pop()
            if assignment_anticodon == truth_anticodon:
                events['single']['single']['correct'] += 1
            else:
                events['single']['single']['incorrect'] += 1
        else:
            if len(set(assignment_anticodons)) == 1:
                assignment_anticodon = assignment_anticodons.pop()
                if assignment_anticodon == truth_anticodon:
                    events['multiple']['single']['correct'] += 1
                else:
                    events['multiple']['single']['incorrect'] += 1
            else:
                events['multiple']['multiple']['incorrect'] += 1

    outf = open(outfile, 'w')
    outf.write('\t'.join(('alignments', 'anticodons', 'agreement', 'count')) + '\n')
    for k,v in events.items():
        for k2, v2 in v.items():
            for k3, c in v2.items():
                outf.write('\t'.join(map(str, (k, k2, k3, c))) + '\n')
    outf.close()

@cluster_runnable
def getTruth2Assignment(infile, outfile, isodecoder_mapping=None):

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

    if isodecoder_mapping is not None:
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

#mimtRNAseq_Hsap_iPSC_rep2.0.simulation_uniform.bowtie2.tsv
@cluster_runnable
def mergesummariseMultimapped(infiles, outfile):
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
    df.to_csv(outfile, sep='\t', index=False)


@cluster_runnable
def mergeTruth2AssignmentNoErrors(infiles, outfile, aligner):
    dfs = []
    for infile in infiles:
        
        df = pd.read_csv(infile, sep='\t')

        if aligner == 'bowtie2':
            d_value = os.path.basename(infile).split('_')[1][1:]
            l_value = os.path.basename(infile).split('_')[2][1:]
            n_value = os.path.basename(infile).split('_')[3][1:]

            df['D']=d_value
            df['L']=l_value
            df['N']=n_value
        elif aligner == 'bwamem':
            k_value = os.path.basename(infile).split('_')[1][1:]
            r_value = os.path.basename(infile).split('_')[2][1:]

            df['K']=k_value
            df['R']=r_value
        dfs.append(df)

    df = pd.concat(dfs)
    df.to_csv(outfile, sep='\t', index=False)


@cluster_runnable
def mergeTruth2Assignment(infiles, outfile):
    dfs = []
    for infile in infiles:
        input_file = os.path.basename(infile).split('.')[0]
        simulation_n = os.path.basename(infile).split('.')[1]
        quant_method = os.path.basename(infile).split('.')[3]

        if quant_method == 'mimseq':
            tally_method == 'mimseq'
        else:
            tally_method = os.path.basename(infile).split('.')[4].replace('gene_count_', '')

        df = pd.read_csv(infile, sep='\t')
        df['input_file']=input_file
        df['simulation_n']=simulation_n
        df['quant_method']=quant_method
        df['tally_method']=tally_method

        dfs.append(df)

    df = pd.concat(dfs)
    df.to_csv(outfile, sep='\t', index=False)


@cluster_runnable
def concatenateEstimateSalmon(infiles, outfile_individual, outfile_isodecoder, outfile_anticodon):

    all_counts = []

    for estimate in infiles:

        estimate_counts = pd.read_csv(estimate, sep='\t')[['Name', 'NumReads']]

        input_file = os.path.basename(os.path.dirname(estimate)).split('.')[0]
        quant_method = os.path.dirname(estimate).split('.')[-1]

        estimate_counts['input_file']=input_file
        estimate_counts['quant_method']=quant_method
        estimate_counts['tally_method']='salmon'

        all_counts.append(estimate_counts)

    all_counts_df = pd.concat(all_counts)
    all_counts_df.to_csv(outfile_individual, sep='\t', index=False)

    all_counts_df['Name'] = ['-'.join(x.split('-')[0:4]) for x in all_counts_df.Name]
    all_counts_isodecoder = all_counts_df.groupby(
        ['Name', 'quant_method', 'tally_method', 'input_file']).agg(
            {'NumReads':'sum'}).reset_index()
    all_counts_isodecoder.to_csv(outfile_isodecoder, sep='\t', index=False)

    all_counts_df['Name'] = [re.sub('\d$', '', '-'.join(x.split('-')[0:3]).replace('tRX', 'tRNA')) for x in all_counts_df.Name]
    all_counts_ac = all_counts_df.groupby(
        ['Name', 'quant_method', 'tally_method', 'input_file']).agg(
            {'NumReads':'sum'}).reset_index()
    all_counts_ac.to_csv(outfile_anticodon, sep='\t', index=False)


@cluster_runnable
def concatenateEstimateDecisionCounts(infiles, outfiles):
    all_counts = {'individual':[], 'isodecoder':[], 'anticodon':[]}

    for infile in infiles:

        estimates = {'individual': infile[0],
                     'isodecoder': infile[1],
                     'anticodon': infile[2]}

        input_file = os.path.basename(infile[0]).split('.')[0]
        quant_method = os.path.basename(infile[0]).split('.')[-3]
        tally_method = os.path.basename(infile[0]).split('.')[-2].replace('gene_count_', '')

        for estimate_level in estimates.keys():

            estimate_counts = pd.read_csv(estimates[estimate_level], sep='\t', header=None,
                                          names=['Name', 'NumReads'])

            estimate_counts['input_file']=input_file
            estimate_counts['quant_method']=quant_method
            estimate_counts['tally_method']=tally_method

            all_counts[estimate_level].append(estimate_counts)

    outfiles = {'individual': outfiles[0],
                 'isodecoder': outfiles[1],
                 'anticodon': outfiles[2]}

    # we can re-use the keys from the last dictionary of filepaths
    for estimate_level in estimates.keys():
        all_counts_df = pd.concat(all_counts[estimate_level])
        all_counts_df.to_csv(outfiles[estimate_level], sep='\t', index=False)


@cluster_runnable
def compareTruthEstimateSalmon(infiles, outfile_individual, outfile_isodecoder, outfile_anticodon):

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

        quant_method = os.path.dirname(estimate).split('.')[-1]
        counts_vs_truth['quant_method']=quant_method

        counts_vs_truth['tally_method']='salmon'

        all_counts_vs_truth.append(counts_vs_truth)

    all_counts_vs_truth = pd.concat(all_counts_vs_truth)
    all_counts_vs_truth.to_csv(outfile_individual, sep='\t', index=False)

    all_counts_vs_truth['Name'] = ['-'.join(x.split('-')[0:4]) for x in all_counts_vs_truth.Name]
    all_counts_vs_truth_isodecoder = all_counts_vs_truth.groupby(
        ['Name', 'simulation_n', 'quant_method', 'tally_method', 'input_file']).agg(
            {'NumReads':'sum', 'truth':'sum'}).reset_index()
    all_counts_vs_truth_isodecoder.to_csv(outfile_isodecoder, sep='\t', index=False)

    all_counts_vs_truth['Name'] = [re.sub('\d$', '', '-'.join(x.split('-')[0:3]).replace('tRX', 'tRNA')) for x in all_counts_vs_truth.Name]
    all_counts_vs_truth_ac = all_counts_vs_truth.groupby(
        ['Name', 'simulation_n', 'quant_method', 'tally_method', 'input_file']).agg(
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
        quant_method = os.path.basename(estimate_individual).split('.')[-3]
        tally_method = (os.path.basename(estimate_individual).split('.')[-2]).replace('gene_count_', '')

        truth_counts = pd.read_csv(truth, sep='\t', header=None,
                                   names=['Name', 'truth'])

        counts_vs_truth = mergeDecisionEstimateWithTruth(
            truth_counts, estimate_individual, input_file, simulation_n, quant_method, tally_method)

        all_counts_vs_truth['individual'].append(counts_vs_truth)

        truth_counts['Name'] = ['-'.join(x.split('-')[0:4]) for x in truth_counts.Name]
        truth_counts_isodecoder = truth_counts.groupby(
            ['Name']).agg(
                {'truth':'sum'}).reset_index()
        counts_vs_truth_isodecoder = mergeDecisionEstimateWithTruth(
            truth_counts_isodecoder, estimate_isodecoder, input_file, simulation_n, quant_method, tally_method)
        all_counts_vs_truth['isodecoder'].append(counts_vs_truth_isodecoder)

        truth_counts['Name'] = [re.sub('\d$', '', '-'.join(x.split('-')[0:3]).replace('tRX', 'tRNA')) for x in truth_counts.Name]
        truth_counts_anticodon = truth_counts.groupby(
            ['Name']).agg(
                {'truth':'sum'}).reset_index()
        counts_vs_truth_anticodon = mergeDecisionEstimateWithTruth(
            truth_counts_anticodon, estimate_anticodon, input_file, simulation_n, quant_method, tally_method)
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
    isodecoder_quant.to_csv(outfile, sep='\t', index=False)

    isodecoder_quant['trna_method'] = [x.split('_')[0] for x in isodecoder_quant['input_file']]

    isodecoder_quant = isodecoder_quant.merge(isodecoder_mapping,
                                              left_on=['Name', 'trna_method'],
                                                  right_on=['isodecoder', 'trna_method'])


    mimseq_isodecoder_quant = isodecoder_quant.groupby(
            ['mimseq_isodecoder', 'input_file', 'simulation_n', 'quant_method', 'tally_method']).agg(
                {'NumReads':'sum', 'truth':'sum'}).reset_index().rename(
                    columns={'mimseq_isodecoder':'Name'})

    mimseq_isodecoder_quant.to_csv(outfile, sep='\t', index=False)
