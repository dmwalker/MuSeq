# Written by David Walker 1502013

import os
import numpy as np
from os.path import isfile, join
from optparse import OptionParser


def bitflag(n):
	bit=bin(n)[2:]
	while len(bit)< 12:
		bit = '0'+bit
	return bit[::-1]

parser = OptionParser()

# Lets give it the location of the files we will be running
parser.add_option("-d", "--dir", dest="work_dir", action="store", type="string", default=None, help="Working directory to put the output")
# General settings
parser.add_option("-g", "--g", dest='ref_genome', action="store", type="string", default=None, help="Pathway to the host (MG1655) reference genome for alignment.")
parser.add_option("-m", "--mu", dest='mu_genome', action="store", type="string", default=None, help="Pathway to the Mu reference genome for alignment.")
parser.add_option("-r1", dest="read_1", action="store", default=None, help="Pathway to the first read file in the paired end reads.")
parser.add_option("-r2", dest="read_2", action="store", default=None, help="Pathway to the second read file in the paired end reads.")
parser.add_option("-o", dest="output", action="store", default=None, help="Output naming scheme for the results.")
parser.add_option("-v", "--verbose", dest="verbose", action="store_true", default=False, help="Run the program verbosely")


# parse command line arguments
(options, args)=parser.parse_args()
w_dir = options.work_dir
ref_genome = options.ref_genome
mu_genome = options.mu_genome
r1 = options.read_1
r2 = options.read_2
output=options.output
verbose=options.verbose


# Combine ref and Mu genomes
if verbose:
    print('Catenating reference and Mu genome into one file.')
os.system('cat %s %s >%s/ref-mu_genome.fasta' % (ref_genome, mu_genome, w_dir))

# Index genome for alignment
if verbose:
    print('Indexing reference and mu genome for alignment.')
os.system("bwa index %s/ref-mu_genome.fasta")


# Align sequences with ref and Mu genomes
align = "bwa mem -t 20 %s/ref-mu_genome.fasta %s %s" % (w_dir, r1, r2)
view = "samtools view -buS"
sort = "samtools sort -O BAM -o %s/%s" % (w_dir, output + 'all_reads_alignment.bam')) 
os.system('%s | %s | %s' % (align, view, sort))

# Index alignments ln 152
os.system('samtools index %s/%s' % (w_dir, r1.replace('.fastq', '.all_reads_alignment.bam')))

# extract the reads aligning with Mu
align = 'bwa mem -T 20 -o %s/tmp.mu_align.sam %s %s %s ' % (w_dir, mu_genome, r1, r2)
os.system(align)

# Clean up samfiles
last=''
reads = {}
R1 = open("%s/%s.fq1" % (w_dir, output), "w")
R2 = open("%s/%s.fq2" % (w_dir, output), "w")
while open('%s/tmp.mu_align.sam' % w_dir) as f:
    for line in f:
        try:
            qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual = line.rsplit('\t')
            if last and qname != last:
                try:
                    R1.write(reads[1])
                    R2.write(reads[2])
                    del reads[1], reads[2]
                except KeyError:
                    pass

            last = qname
            flag = bitflag(flag)
            if int(flag[2]) and int(flag[3]):
                pass
            if int(flag[8]) or int(flag[11]):
                pass
            elif int(flag[4]):
                print('reverse')
            r = 1 if int(flag)[6] else r=2
            info = "\@%s\n%s\n+\n%s\n" % (qname, seq, qual)
            reads[r] = info
            
        except IndexError:
            pass
            
R1.close()
R2.close()

rds = "%s/%s.fq1 %s/%s.fq2" % (w_dir, output, w_dir, output)

align = "bwa mem -t 20 %s/ref-mu_genome.fasta %s" % (w_dir, rds)
view = "samtools view -buS"
sort = "samtools sort -O BAM -o %s/%s" % (w_dir, output + '.ref-mu.sorted.bam')) 
os.system('%s | %s | %s' % (align, view, sort))
