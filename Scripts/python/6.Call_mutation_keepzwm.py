import sys
import os
import time
import re
import multiprocessing
import argparse
from Bio import SeqIO, pairwise2
sys.path.append("/mnt/data5/disk/yangwj/4.Analysis_PacBio/Scripts/")
from modules.SequenceParser import formFASTA

def getReads(seq):
    for record in seq:
        name = record.description
        seq_record = str(record.seq)
        yield name, seq_record

def call_mutations(refseq, readname, readseq):
    alignments = pairwise2.align.globalms(refseq, readseq, 2, -3, -2, -1)
    best_alignment = alignments[0]
    aligned_qual = best_alignment.score

    editMutation = []
    #indel_Ins = []
    #indel_Del = []
    aligned_seq = ""
    refpos = 0
    readpos = 0
    current_readpos = 0
    num_mut = 0
    num_indel_ins = 0
    num_indel_del = 0
    num_indel = 0
    num_same_readpos = 0
    maxNum_continuousDel = 0
    count_CT_GA = 0
    ratio_CT_GA = 0
    current_refbase = 0
    num_same_refpos = 0
    maxNum_continuousIns = 0
    for refbase, readbase in zip(best_alignment[0], best_alignment[1]):
        if refbase != "-":
            refpos += 1
        if readbase != "-":
            readpos += 1
            aligned_seq += readbase
            if refbase != readbase and refbase in "ATGC" and readbase in "ATGC":
                editMutation.append(f"{refbase}{refpos}{readbase}")
                num_mut += 1
                if (refbase, readbase) == ("C", "T") or (refbase, readbase) == ("G", "A") :
                    count_CT_GA +=1
                    
        if (refbase, readbase) == ("-", readbase):
            #indel_Ins.append(f"{readpos}I")
            num_indel_ins += 1
            num_indel += 1
            if current_refbase == refbase:
                num_same_refpos += 1
                if num_same_refpos >= maxNum_continuousIns:
                    maxNum_continuousIns = num_same_refpos
            else:
                num_same_refpos = 1
                current_refbase = refbase
                
        if (refbase, readbase) == (refbase, "-"):
            #indel_Del.append(f"{readpos}D")
            num_indel_del += 1
            num_indel += 1
            if current_readpos == readpos:
                num_same_readpos += 1
                if num_same_readpos >= maxNum_continuousDel:
                    maxNum_continuousDel = num_same_readpos
            else:
                num_same_readpos = 1
                current_readpos = readpos
    ratio_CT_GA = count_CT_GA/num_mut if num_mut !=0 else 0
                 
    return readname, aligned_qual, num_mut, num_indel_ins, num_indel_del, maxNum_continuousIns, maxNum_continuousDel, num_indel, ratio_CT_GA, editMutation 


def Parsers():
    parser = argparse.ArgumentParser(description="call_mutation base on pairwise2.align.globalms.")
    parser.add_argument('-r', '--reference', type=str, nargs='?', required=True, help="References file")
    parser.add_argument('-s', '--seq', type=str, nargs='?' ,required=True, help="Sequence file")
    #parser.add_argument('-ratio', '--ratio_CG_TA', type=str, nargs='?' ,required=True, help="filter_ratio_CG_TA")
    parser.add_argument('-o', '--outfile', type=str, nargs='?', required=True, help="Output editevents file")
    parser.add_argument('-c', '--cpu', type=int,  nargs='?',required=True, help="The Number of cpu needed")
    args = parser.parse_args()
    return args



if __name__ == "__main__":
    start = time.time()
    options = Parsers()
    ref_seq = open(options.reference,'r').read()
    ref = formFASTA(ref_seq)
    _, refseq = next(ref)
    query_seq_file = SeqIO.parse(options.seq, 'fasta')
    #Create a process pool with multiple processes
    pool = multiprocessing.Pool(options.cpu)
    results = []
    for readname, readseq in getReads(query_seq_file):
        pfq = pool.apply_async(call_mutations, args=(refseq, readname, readseq))
        results.append(pfq)
           
    print('wait')
    pool.close()
    #Wait for all processes to complete
    pool.join()
    print('done')
    
    #ratio = float(options.ratio_CG_TA)
    outFA = open(options.outfile, 'w')
    for align in results:
        result = align.get()
        outFA.write(f">{result[0]} qual={result[1]} num.Mut={result[2]} num.indel_ins={result[3]} num.indel_del={result[4]} num.Max_ConIns={result[5]} num.Max_ConDel={result[6]} num.indel={result[7]} ratio_CT_GA={result[8]}\neditMutation:{result[9]}\n")
        #if result[8] > ratio:
            #outFA.write(f">{result[0]} qual={result[1]} num.Mut={result[2]} num.indel_ins={result[3]} num.indel_del={result[4]} num.Max_ConIns={result[5]} num.Max_ConDel={result[6]} num.indel={result[7]} ratio_CT_GA={result[8]}\neditMutation:{result[9]}\n")

        #if result[5] <= 5 and result[6] <= 20:
            #outFA.write(f">{result[0]} qual={result[1]} num.Mut={result[2]} num.indel_ins={result[3]} num.indel_del={result[4]} num.Max_ConDel={result[5]} num.indel={result[6]}\neditMutation:{result[7]}\n")

    outFA.close()
     
    end = time.time()
    print(end - start)