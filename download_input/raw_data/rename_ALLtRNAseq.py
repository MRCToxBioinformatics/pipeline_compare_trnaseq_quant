import os
infile = './filereport_read_run_PRJNA775872_tsv_hESC.txt'
with open(infile, 'r') as inf:

    for line in inf:
        line = line.strip().split('\t')
        fastq_filename = line[3]  + '.fastq.gz'

        sample_name = line[12]

        cell = sample_name.split(' ')[0]
        timepoint = sample_name.split(' ')[1][-1]
        replicate = sample_name[-1]
 
        lane = line[8][-1]

        updated_filename = 'ALLtRNAseq_%s_%s_%s-l%s.fastq.gz' % (cell, timepoint, replicate, lane) 

        print(fastq_filename, updated_filename)

        if os.path.exists(fastq_filename):
            os.rename(fastq_filename, updated_filename)
