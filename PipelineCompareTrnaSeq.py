import simulatetrna.alignmentSummary as alignmentSummary
import simulatetrna.bam as bam
from cgatcore.pipeline import cluster_runnable

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
def tally_read_counts(bam_infile, outfile_individual, outfile_isodecoder, outfile_anticodon):

    inbam = pysam.Samfile(bam_infile, 'r')

    read_counts = {'individual_sequence':Counter(),
                   'isodecoder':Counter(),
                   'anticodon':Counter()}

    for read_group in bam.iterate_reads(inbam):

        assignments = set([read.reference_name for read in read_group if read.reference_name is not None])

        # No assignments for read
        if len(assignments) == 0:
            continue

        # Just one assigned gene
        if len(assignments)==0:
            assignment = assignments.pop()
            read_counts['individual_sequence'][assignment] +=1
            read_counts['anticodon']['-'.join(assignment.split('-')[0:3])] +=1
            read_counts['isodecoder']['-'.join(assignment.split('-')[0:4])] +=1

        # multiple gene assignments
        else:
            anticodons = set(['-'.join(x.split('-')[0:3]) for x in assignments])
            isodecoders = set(['-'.join(x.split('-')[0:4]) for x in assignments])

            if len(anticodons)==1:
                read_counts['anticodon'][anticodons.pop()] +=1

            if len(isodecoders)==1:
                read_counts['isodecoder'][isodecoders.pop()] +=1

    with open(outfile_individual, 'w') as outf:
        for k,v in read_counts['individual_sequence'].items():
            outf.write(k + '\t' + str(v) + '\n')

    with open(outfile_isodecoder, 'w') as outf:
        for k,v in read_counts['isodecoder'].items():
            outf.write(k + '\t' + str(v) + '\n')

    with open(outfile_anticodon, 'w') as outf:
        for k,v in read_counts['anticodon'].items():
            outf.write(k + '\t' + str(v) + '\n')
