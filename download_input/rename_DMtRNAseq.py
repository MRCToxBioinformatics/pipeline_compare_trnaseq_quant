import os
infile = './filereport_read_run_PRJNA593498_tsv.txt'
with open(infile, 'r') as inf:

    header = inf.readline()
    for line in inf:
        line = line.strip().split('\t')
        fastq_filename = os.path.join(line[2], line[2] + '.fastq.gz')
        updated_filename = 'DMtRNAseq_' + line[6] + '.fastq.gz'

        if os.path.exists(fastq_filename):
            os.rename(fastq_filename, updated_filename)
