import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('name',help='sample name')
parser.add_argument('ref_fasta', type = argparse.FileType('r'),help ='path to human genome reference file in FASTA format')
parser.add_argument('sample_bam',type = argparse.FileType('r'),help='path to sample file in BAM format ')
parser.add_argument('union',type = argparse.FileType('r'),help='path to list of SNVs in tsv format')
parser.add_argument('--baseQ',nargs = '?',default=0,help= 'integer for baseQ threshold')
parser.add_argument('--mapQ',nargs = '?',default=0,help= 'integer for mapQ threshold')
parser.add_argument('--baseQ_file',nargs = '?',type = argparse.FileType('r'),help= 'path to baseQ thresholds file in tsv format')
parser.add_argument('--mapQ_file',nargs = '?',type = argparse.FileType('r'),help= 'path to mapQ.tsv threshold file in tsv format')
args = parser.parse_args()

baseQ = False
mapQ = False
name = args.name
ref = args.ref_fasta
bam = args.sample_bam
union = args.union
q = args.baseQ
Q = args.mapQ
if args.baseQ_file:
    with open(os.path.abspath(args.baseQ_file)) as f:
        baseQ_list= f.read().splitlines()
        baseQ = True

if args.mapQ_file:
    with open(os.path.abspath(args.mapQ_file)) as f:
        mapQ_list= f.read().splitlines()
        mapQ = True
#case 1
if baseQ and mapQ:
    i = 0
    with open(union ,'r') as u:
        for line in inf:
            snv =line.split(".")
            chrom = snv[0]
            pos = snv[1]
            os.system('samtools mpileup -r' + chrom + ':' + pos + '-' + pos + '-q '+ baseQ_list[i] + '-Q ' + mapQ_list[i] + ' -x -f ' + ref + ' ' + bam + ' -s -O  >> ' + name + '.txt')
            i += 1
# case 2
if not baseQ and  mapQ:
    i = 0
    with open(union ,'r') as u:
        for line in inf:
            snv =line.split(".")
            chrom = snv[0]
            pos = snv[1]
            os.system('samtools mpileup -r' + chrom + ':' + pos + '-' + pos + '-q '+ q + '-Q ' + mapQ_list[i] + ' -x -f ' + ref + ' ' + bam + ' -s -O  >> ' + name + '.txt')
            i += 1

# case 3
if baseQ and not mapQ:
    i = 0
    with open(union ,'r') as u:
        for line in inf:
            snv =line.split(".")
            chrom = snv[0]
            pos = snv[1]
            os.system('samtools mpileup -r' + chrom + ':' + pos + '-' + pos + '-q '+ baseQ_list[i] + '-Q ' + Q + ' -x -f ' + ref + ' ' + bam + ' -s -O  >> ' + name + '.txt')
            i += 1
# case 4

if not baseQ and not mapQ:
    with open(union ,'r') as u:
        for line in inf:
            snv =line.split(".")
            chrom = snv[0]
            pos = snv[1]
            os.system('samtools mpileup -r' + chrom + ':' + pos + '-' + pos + '-q '+ q + '-Q ' + Q + ' -x -f ' + ref + ' ' + bam + ' -s -O  >> ' + name + '.txt')
