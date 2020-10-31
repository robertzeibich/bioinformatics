#!/usr/bin/python

import sys
import getopt
import gzip

def usage(): print("""usage: fastq_paired_end_quality_filter [-h] [-q N] [-p N] [-i INFILE] [-o OUTFILE] [-o OUTFILE]
by R. Zeibich (robertzeibich@yahoo.de)

   [-h]          = This helpful help screen.
   [-q N]        = Minimum quality score to keep.
   [-p N]        = Minimum percent of bases that must have [-q] quality.
   [-i INFILE]   = Compressed FASTQ input file (e.g. -i fastq.gz).
   [-o OUTFILES] = FASTQs output files (e.g. -o out1.fastq -o out2.fastq).""")
   
opts,args = getopt.getopt(sys.argv[1:],'hi:o:q:p:')

opts_dict = {}
for o,a in opts:
    opts_dict[o] = []
for o,a in opts:
    opts_dict[o].append(a) 

if '-h' in opts_dict.keys():
    usage();sys.exit()
elif '-q' not in opts_dict.keys():
    usage();sys.exit("Minimum quality score missing")
elif '-p' not in opts_dict.keys():
    usage();sys.exit("Minimum percent of bases that must have [-q] quality missing")
elif '-i' not in opts_dict.keys():
    usage();sys.exit("Compressed input fastq.gz file missing")


fastq_in = opts_dict['-i'][0]
q = int(opts_dict['-q'][0])
p = int(opts_dict['-p'][0])

fastq1_out = 'out1.fastq'
fastq2_out = 'out2.fastq'

if '-o' in opts_dict.keys():
    if len(opts_dict['-o']) == 2:
        fastq1_out = opts_dict['-o'][0]
        fastq2_out = opts_dict['-o'][1]

def phred33ToQ(qual):
    """Turn Phred+33 ASCII-encoded quality into Q"""
    return ord(qual)-33
    
def fastq_paired_end_quality_filter(fastq_in,q,p,fastq1_out,fastq2_out):
          
    f = open(fastq1_out, "w")
    f.close()
    
    f = open(fastq2_out, "w")
    f.close()
    
    with gzip.open(fastq_in) as fq:
        while True:
            
            head1 = fq.readline().rstrip().decode("utf-8")
            seq1 = fq.readline().rstrip().decode("utf-8")
            ph1 = fq.readline().rstrip().decode("utf-8")
            qual1 = fq.readline().rstrip().decode("utf-8")
               
            head2 = fq.readline().rstrip().decode("utf-8")
            seq2 = fq.readline().rstrip().decode("utf-8")
            ph2 = fq.readline().rstrip().decode("utf-8")
            qual2 = fq.readline().rstrip().decode("utf-8")
            
            qual_len = len(qual1)
            allowance = qual_len - int(qual_len*p/100)
            
            count1 = 0
            for val in qual1:
                qual_val = phred33ToQ(val)
                if qual_val < q:
                    count1+=1
                if count1 > allowance:
                    break
            
            count2 = 0  
            if count1 <= allowance:
                for val in qual2:
                    qual_val = phred33ToQ(val)
                    if qual_val < q:
                        count2+=1
                    if count2 > allowance:
                        break
                    
            if len(seq1) == 0 or len(seq2) == 0:
                break
                
            if count1 <= allowance and count2 <= allowance:  
                with open(fastq1_out,'a') as fq1:
                    fq1.write(head1+'\n')
                    fq1.write(seq1+'\n')
                    fq1.write(ph1+'\n')
                    fq1.write(qual1+'\n')
                    fq1.close()
                with open(fastq2_out,'a') as fq2:
                    fq2.write(head2+'\n')
                    fq2.write(seq2+'\n')
                    fq2.write(ph2+'\n')
                    fq2.write(qual2+'\n')
                    fq2.close()
        fq.close()

fastq_paired_end_quality_filter(fastq_in,q,p,fastq1_out,fastq2_out)
