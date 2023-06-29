# -*- coding: utf-8 -*-

"""
Parse a bam file for barcodes
bam file index sorted and indexed
The program is written with illumina pair reads in mind.
It will extract the read (and its mate) containing a 
barcode sequence. Shuld work for single end reads as well

usage
python extract_barcodes_def.py bamFile barcode_A barcode_B 1 False
python extract_barcodes_def.py bamFile barcode_A barcode_B 1 True
python extract_barcodes_def.py bamFile barcode_A barcode_B 2


--python extract_barcodes_def.py bamFile barcode_A barcode_B 1 False--
--python extract_barcodes_def.py bamFile barcode_A barcode_B 1 True--
Search for a barcode in forward and reverse complement
This is the case of phenotyping experiment with iRNA
The first argument is the barcode to search in forward (TCGCGAGGC)
and the second argument is the reverse complement of the barcode (GCCTCGCGA).

If the fourth argument != True:
3 outputs file
_F.bam  -> Pair reads containing the forward barcode in forward orientation
_R.bam  -> Pair reads containing the forward barcode in reverse complement orientation
_WT.bam -> Pair reads without barcode


If the fourth argument == True:
F and R outputs are merged, sorted by reference, indexed and removed.
2 outputs file
_F_plus_R.bam -> Pair reads containing the barcode
_WT.bam -> Pair reads without barcode


--python extract_barcodes_def.py bamFile barcode_A barcode_B 2--
Search for two barcodes, both in forward and reverse complement orientation
This is the case of overexpression library. One barcode identify the 
forward strand, the second barcode identify the reverse strand. 
The first argument is the barcode for the forward strand and the 
second argument is the barcode for the reverse strand.
Both barcodes are given in forward sequence.
The program will use the forward and reverse complement to search for both barcodes

4 outputs files expected 
_F.bam  -> Pair reads containing the forward barcode in forward orientation
_FR.bam -> Pair reads containing the forward barcode in reverse complement orientation
_R.bam  -> Pair reads containing the reverse barcode in forward orientation
_RR.bam -> Pair reads containing the reverse barcode in reverse complement orientation
"""

#import dill
#dill.detect.trace(True)
import pysam, argparse
from Bio.Seq import Seq
import sys
import tqdm
import os
import gc
#from multiprocessing import Pool
#import multiprocessing
#import bamnostic as bs
import functools
from subprocess import Popen, PIPE

def sort_and_index(f_name):
    """
    sort and index a bam file by reference
    new file replace old
    """
    pysam.sort("-o", f_name, f_name)
    pysam.index(f_name)

#https://www.biostars.org/p/1890/
def bam_read_count(bamfile):
    """ Return a tuple of the number of mapped and unmapped reads in a bam file
        Bam file reference sorted and indexed
    """
    p = Popen(['samtools', 'idxstats', bamfile], stdout=PIPE)
    mapped = 0
    unmapped = 0
    for line in p.stdout:
        rname, rlen, nm, nu = line.rstrip().split()
        mapped += int(nm)
        unmapped += int(nu)
    #p.close()
    return (mapped, unmapped)


#name of bam file
infile = sys.argv[1]
out_base = infile.replace('.bam','')
bamfile = pysam.Samfile(infile,"rb")

mapped, unmapped = bam_read_count(infile)
all_reads = mapped+unmapped

#store barcode Forward
f_barcode= Seq(sys.argv[2])
out_f = out_base+'_F.bam'
out_f = pysam.Samfile(out_f, "wb", template = bamfile)

#store barcode Reverse
r_barcode = Seq(sys.argv[3])
out_r = out_base+'_R.bam'
out_r = pysam.Samfile(out_r, "wb", template = bamfile)

#are we looking for one barcode (RIT, phenotyping with interference)
#or two barcodes (OEL overexpresiion libraries)?
barcode_type = int(sys.argv[4])

#dictonary to point to the output files
pointer_dict={
    'F':out_f,
    'R':out_r}

fcount = 0
rcount = 0

#setup to search for the reverse complement 
#sequences of the barcodes
if barcode_type == 2:
    print('two barcodes:' 'Remember: input both barcode seqences in forward orientation')
    #store barcode Forward (reverse complement)
    fr_barcode = f_barcode.reverse_complement()
    out_fr = out_base+'_FR.bam'
    out_fr = pysam.Samfile(out_fr, "wb", template = bamfile)
    #store barcode Reverse (reverse complement)
    rr_barcode = r_barcode.reverse_complement()
    out_rr = out_base+'_RR.bam'
    out_rr = pysam.Samfile(out_rr, "wb", template = bamfile)
    #update the pointer to the add two more output files
    pointer_dict['FR']=out_fr
    pointer_dict['RR']=out_fr
    frcount = 0
    rrcount = 0
elif barcode_type == 1:
    print('one barcode:' 'Remember: input the forward and reverse complement sequences')
else:
    print('something odd, i shuld not be here with barcode type')
    sys.exit(1)

tot=0
done = {}

#transform the Bio Seq object to string
f_barcode=str(f_barcode)
r_barcode=str(r_barcode)

if 'fr_barcode' in globals():
    fr_barcode=str(fr_barcode)
if 'rr_barcode' in globals():
    rr_barcode=str(rr_barcode)


for line in tqdm.tqdm(bamfile,total=all_reads,miniters=1000000):
    tot+=1
    #we catch a read as soon as we see a barcode sequence
    #and we store the fastq read name in done
    #we prioratize in the order F R FR RR
    if line.query_name not in done:
        if f_barcode in line.seq:
            fcount+=1
            #out_f.write(line)
            done[line.query_name]='F'
        
        elif r_barcode in line.seq:
            rcount+=1
            #out_r.write(line)
            done[line.query_name]='R'
        
        elif barcode_type == 2:
            if fr_barcode in line.seq:
                frcount+=1
                #out_fr.write(line)
                done[line.query_name]='FR'
            
            elif rr_barcode in line.seq:
                rrcount+=1
                #out_rr.write(line)
                done[line.query_name]='RR'
    else:
        #if here, we already saw the read/mate
        #by using read name we shuld catch:
        #any secondary aligment
        #any mate read where the other mate contains a barcode
        pass
    
    ##for testing, parsing a few alignments
    #if tot >=500:
    #    break



bamfile.close()


#done = set(done.keys())
print('I shell now write out files')
bamfile = pysam.Samfile(infile,"rb")
count_wo=0
out_wo = out_base+'_WO.bam'
out_wo = pysam.Samfile(out_wo, "wb", template = bamfile)
#tot=0
for line in tqdm.tqdm(bamfile,total=all_reads,miniters=1000000): 
    #tot+=1   
    if line.query_name not in done:
        count_wo+=1
        out_wo.write(line)
    else:
        out_file = done[line.query_name]
        pointer_dict[out_file].write(line)
    #if tot >=500:
    #    break


out_f.close()
out_r.close()
out_wo.close()        
bamfile.close()
print ('sorting - indexing reads without barcode')
sort_and_index(out_base+'_WO.bam')

if barcode_type == 2:
    out_fr.close()
    out_rr.close()
    print ('I found',tot,fcount,rcount,frcount,rrcount,count_wo)
    print ('sorting - indexing')
    sort_and_index(out_base+'_F.bam')
    sort_and_index(out_base+'_R.bam')
    sort_and_index(out_base+'_FR.bam')
    sort_and_index(out_base+'_RR.bam')

elif barcode_type == 1:
    if sys.argv[5] == 'True':
        print ('merging')
        pysam.merge("-f", out_base+'_F_plus_R.bam', out_base+'_F.bam', out_base+'_R.bam')
        print ('sorting - indexing')
        sort_and_index(out_base+'_F_plus_R.bam')
        print('delete startig files')
        os.remove(out_base+'_F.bam')
        os.remove(out_base+'_R.bam')
    
    elif sys.argv[5] == 'False':
        sort_and_index(out_base+'_F.bam')
        sort_and_index(out_base+'_R.bam')
    
    else:
        print('Merge is either True or False')
        sys.exit(1)
    print ('I found',tot,fcount,rcount,count_wo)
else:
    print('barcode_type is 1 or 2')
    sys.exit(1)

del done
gc.collect()
print('done')
sys.exit(0)
