# Written by David Walker 1502013
# Version 1.2

import os
import sys
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
parser.add_option("--r1", dest="read_1", action="store", default=None, help="Pathway to the first read file in the paired end reads.")
parser.add_option("--r2", dest="read_2", action="store", default=None, help="Pathway to the second read file in the paired end reads.")
parser.add_option("-o", dest="output", action="store", default=None, help="Output naming scheme for the results.")
parser.add_option("-e", dest="enhanced", action="store_true", default=False, help="Flag for target enriched samples")
parser.add_option("-e", dest="self_insertion", action="store_true", default=False, help="Flag for monitoring Mu self insertions.")
parser.add_option("-v", "--verbose", dest="verbose", action="store_true", default=False, help="Run the program verbosely")
parser.add_option("-H", "--long-help", dest="longhelp", action="store_true", default=False, help="Display long help")


# parse command line arguments
(options, args)=parser.parse_args()
w_dir = options.work_dir
ref_genome = options.ref_genome
mu_genome = options.mu_genome
r1 = options.read_1
r2 = options.read_2
output=options.output
verbose=options.verbose
enhanced=options.enhanced
self_insertion=option.self_insertions
longhelp=options.longhelp

if longhelp or len(sys.argv) == 1 or w_dir is None:
    scriptname = sys.argv[0]
    print("Long help")
    parser.print_help()
    print('Version 1.2')
    exit()
##    print("You need to supply a directory with fasta formated DNA sequences like the following:")
##    print("\n genome-parse.py -d \"\\C:\\Research\\14101611seq\\ \n")
##    print("Note this wil only work with fasta formating at the moment.")

with open(ref_genome) as ref:
    ref_info = ref.readline().split()
    host = ref_info[0][1:]

with open(mu_genome) as mu:
    mu_info = mu.readline().split()
    mu = mu_info[0][1:]
'''
# Combine ref and Mu genomes
if verbose:
    print('Catenating reference and Mu genome into one file.')
os.system('cat %s %s >%s/ref-mu_genome.fasta' % (ref_genome, mu_genome, w_dir))
# Index genome for alignment
if verbose:
    print('Indexing reference and mu genome for alignment.')
os.system("bwa index %s/ref-mu_genome.fasta" % w_dir)
# Align sequences with ref and Mu genomes
align = "bwa mem -t 20 %s/ref-mu_genome.fasta %s %s" % (w_dir, r1, r2)
view = "samtools view -buS"
sort = "samtools sort -O BAM -o %s/%s" % (w_dir, output + '.all_reads_alignment.bam')
os.system('%s | %s | %s' % (align, view, sort))
# Index alignments
os.system('samtools index %s/%s' % (w_dir, output + '.all_reads_alignment.bam'))
# extract the reads aligning with Mu
align = 'bwa mem -T 20 -o %s/tmp.mu_align.sam %s %s %s ' % (w_dir, mu_genome, r1, r2)
os.system(align)
# Clean up samfiles
last=''
reads = {}
R1 = open("%s/%s.fq1" % (w_dir, output), "w")
R2 = open("%s/%s.fq2" % (w_dir, output), "w")
with open('%s/tmp.mu_align.sam' % w_dir) as f:
    for line in f:
        try:
            qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual = line.rsplit('\t')[:11]
            if last and qname != last:
                try:
                    R1.write(reads[1])
                    R2.write(reads[2])
                    del reads[1], reads[2]
                except KeyError:
                    pass
            last = qname
            flag = bitflag(int(flag))
            if int(flag[2]) and int(flag[3]):
                pass
            if int(flag[8]) or int(flag[11]):
                pass
            if int(flag[6]):
                r=1
            else:
                r=2
            info = "\@%s\n%s\n+\n%s\n" % (qname, seq, qual)
            reads[r] = info
            
        except IndexError as ie:
            print(ie)
            pass
        except ValueError as ve:
            print(ve)
            pass
            
R1.close()
R2.close()
rds = "%s/%s.fq1 %s/%s.fq2" % (w_dir, output, w_dir, output)
# sorting out the reads
align = "bwa mem -t 20 %s/ref-mu_genome.fasta %s" % (w_dir, rds)
view = "samtools view -buS"
sort = "samtools sort -O BAM -o %s/%s" % (w_dir, output + '.ref-mu.sorted.bam')
os.system('%s | tee %s/%s.ref-mu.sam | %s | %s' % (align, w_dir, output, view, sort))
# index the reads
index = "samtools index %s/%s" % (w_dir, output + '.ref-mu.sorted.bam')
os.system(index)
# extract identify reads that have both Mu and genome alignment
junctions = open('%s/%s.mu-junctions.sam' % (w_dir, output), "w")
print('Writing Mu junctions...')
with open('%s/%s.ref-mu.sam' % (w_dir, output)) as f:
    for line in f:
        try:
            qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual = line.rsplit('\t')[:11]
            if rname != rnext and rnext != '=':
                junctions.write(line)
                
        except ValueError as ve:
             print(ve)
             if len(line.rsplit('\t')) < 6:
                junctions.write(line)
             pass
        except IndexError as ie:
             print(ie)
             pass
'''
# modifying the reads
view = "samtools view -h -S %s/%s -o %s/%s.informative.sam" % (w_dir,output+'.mu-junctions.sam', w_dir, output)
os.system(view)

view = "samtools view -buS %s/%s.informative.sam" % (w_dir, output)
sort = "samtools sort -O BAM -o %s/%s.informative.sorted.bam" % (w_dir, output)
index = "samtools index %s/%s.informative.sorted.bam" % (w_dir, output)
if verbose:
    print('%s | %s' % (view,sort))
os.system('%s | %s' % (view, sort))
os.system(index)

os.system('samtools view -h -o %s/%s.informative.sorted.sam %s/%s.informative.sorted.bam' % (w_dir, output, w_dir, output))
if verbose:
    print('Sorted sam file is now ready for identifying inserts...')

mu_insert_array = np.zeros(4641625)
print('Writing Mu Insertion Locations')
if enhanced:
        mu_positions = open('%s/%s.enriched_sites.sam' % (w_dir, output), "w")
        with open('%s/%s.informative.sorted.sam' % (w_dir, output)) as f:
            info=[]
            old_seq = 0
            for line in f:
                try:
                    qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual = line.rsplit('\t')[:11]
                    if 'TTCAGAGCTTC' in seq and rname == 'U00096.3':
                        if seq != old_seq:
                            cig = cigar.split('S')# The CIGAR string will start with a soft mismatch that is the y-linker adapter and then follow up with E. coli DNA
                            cig1 = cig[1].split('M')
                            pos = int(pos) + int(cig1[0])
                            mu_positions.write(line)
                            info.append([qname, pos, rname])
                        old_seq = seq
                        
                except ValueError as ve:
                     print(ve)
                     if len(line.rsplit('\t')) < 6:
                        enriched.write(line)
                     pass
                except IndexError as ie:
                     print(ie)
                     pass
            np.savetxt('%s/%s.nucleotide_positions.txt' % (w_dir, output), np.array(info), fmt='%s')


elif not enhanced:        
        mu_positions = open('%s/%s.mu_inserts.txt' % (w_dir, output), 'w')
        with open('%s/%s.informative.sorted.sam' % (w_dir, output)) as f:
            info=[]
            old_seq=0
            for line in f:
                try:
                    qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual = line.rsplit('\t')[:11]
                    flag = bitflag(int(flag))
                    if rname == host:
                        if int(flag[4]):
                            mu_insert = int(pos)
                        else:
                            match = re.match(r'[^M]*', cigar)
                            match = match.group()
                            if 'S' or 'H' in match:
                                match = re.split(r'[S,H]', match)[-1]
                            mu_insert = int(pos) + int(match)
                        mu_insert_array[mu_insert] += 1
                        mu_positions.write(line)
                        info.append([qname, mu_insert, rname])

                except ValueError as ve:
                    if verbose:
                       print(ve)
                    pass
                except IndexError as ie:
                    if verbose:
                       print(mu_insert)
                       print(ie)
                    pass

            np.savetxt('%s/%s.nucleotide_positions.txt' % (w_dir, output), np.array(info), fmt='%s')

mu_positions.close()
mu_inserts = np.transpose([np.arange(4641625), mu_insert_array])
np.savetxt('%s/%s.transposition.counts.txt' % (w_dir, output), mu_inserts, fmt='%d', delimiter='\t')

if self_insertion:
        junctions = open('%s/%s.mu-self_inserts.sam' % (w_dir, output), "w")
        print('Writing Mu junctions...')
        with open('%s/%s.ref-mu.sam' % (w_dir, output)) as f:
            info =[]
            for line in f:
                try:
                    qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual = line.rsplit('\t')[:11]
                    if 'TTCAGAGCTTC' in seq and rname == 'AF083977.1':
                        if rnext == '=' and np.abs(int(tlen))>500:
                            if pos == '1':
                                info.append([pnext, rname])
                            junctions.write(line)
                except ValueError as ve:
                     print(ve)
                     if len(line.rsplit('\t')) < 6:
                        junctions.write(line)
                     pass
                except IndexError as ie:
                     print(ie)
                     pass
            np.savetxt('%s/%s.self_insertions.txt' % (w_dir, output), np.array(info), fmt='%s')
        junctions.close()

print('Finished tabulating Mu insertions...')


