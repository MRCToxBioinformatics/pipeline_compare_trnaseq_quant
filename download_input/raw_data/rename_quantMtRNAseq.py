import os

infile = './filereport_read_run_PRJNA593498_quantMtRNAseq_10090_Cortex_Cerebellum_tsv.txt'

with open(infile, 'r') as inf:

    header = inf.readline()
    for line in inf:
        line = line.strip().split('\t')
        fastq_filename = line[0] + '.fastq.gz'

        if os.path.exists(fastq_filename):
            print(fastq_filename)
            print(line)
            updated_filename = 'quantMtRNAseq_' + line[4].split(': ')[2].replace(' Mus musculus OTHER', '') + '.fastq.gz'

            os.rename(fastq_filename, updated_filename)
