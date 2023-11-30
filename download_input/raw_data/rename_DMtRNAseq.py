import os
infile = 'filereport_read_run_PRJNA277309_tsv_tRNA.txt'
with open(infile, 'r') as inf:

    header = inf.readline()
    for line in inf:
        line = line.strip().split('\t')
        condition = line[9].split(' ')[1]
        replicate_number = line[9][-1]
        updated_filename = 'DMtRNAseq_%s_%s.fastq.gz' % (condition, replicate_number)
        fastq_filename = line[3] + '.fastq.gz'

        print(fastq_filename, updated_filename)
        
        if os.path.exists(fastq_filename):
            os.rename(fastq_filename, updated_filename)

