import os
import sys
import re
import time
import argparse
from collections import OrderedDict, defaultdict
sys.path.append("/mnt/data5/disk/yangwj/4.Analysis_PacBio/Scripts/")
from modules.SequenceParser import formFASTA

def getFilterBClist(filterBC_file):
    if filterBC_file == "All":
        return None
    else:
        BClist = [i.strip() for i in open(filterBC_file, 'r')]
        return BClist
    
def identEdit(editString):
    # Use regular expressions to match the numbers and the bases before and after the numbers
    pattern = r'([ACGT])(\d+)([ACGT])'
    # Use the findall function to findall matching results
    matches = re.findall(pattern, editString)
    mutations = []
    for match in matches:
        refbase = match[0]
        refpos = int(match[1]) - 1
        mutbase = match[2]
        mutations.append((refbase, refpos, mutbase))
    return mutations

def convert(edits, ref):
    res = identEdit(edits)
    sequence_list = list(ref)
    for mut in res:
        oribase = mut[0]
        pos_cv = mut[1]
        replace_base = mut[2]
        if sequence_list[pos_cv] == oribase:
            sequence_list[pos_cv] = replace_base
        else:
            sequence_list[pos_cv] = "N"
    new_sequence = ''.join(sequence_list)
    return new_sequence

def count_edit(edits):
    freq = {}
    for edit in edits: 
        for e in edit.split("\t")[1].split(','):
            if e not in freq:
                freq[e] = 0
            freq[e] += 1
    return freq

def count_mutation(eventLineList):
    count_CG_TA = 0
    count_CG_GC = 0
    count_CG_AT = 0
    count_TA_CG = 0
    count_TA_AT = 0
    count_TA_GC = 0
    for line in eventLineList:
        spl = line.strip().split('\t')
        edits ='\t'.join(spl[1].split(','))
        pattern = r'([ACGT])(\d+)([ACGT])'
        matches = re.findall(pattern, edits)
        for match in matches:
            refbase = match[0]
            mutbase = match[2]
            if (refbase, mutbase) == ("C", "T") or (refbase, mutbase) == ("G", "A"):
                count_CG_TA +=1
            if (refbase, mutbase) == ("C", "G") or (refbase, mutbase) == ("G", "C"):
                count_CG_GC +=1
            if (refbase, mutbase) == ("C", "A") or (refbase, mutbase) == ("G", "T"):
                count_CG_AT +=1
            if (refbase, mutbase) == ("T", "C") or (refbase, mutbase) == ("A", "G"):
                count_TA_CG +=1
            if (refbase, mutbase) == ("T", "A") or (refbase, mutbase) == ("A", "T"):
                count_TA_AT +=1    
            if (refbase, mutbase) == ("T", "G") or (refbase, mutbase) == ("A", "C"):
                count_TA_GC +=1
    return [("count_CG_TA", count_CG_TA),("count_CG_GC", count_CG_GC),("count_CG_AT", count_CG_AT),
            ("count_TA_CG", count_TA_CG),("count_TA_AT", count_TA_AT),("count_TA_GC", count_TA_GC)]
            

def Parsers():
    parser = argparse.ArgumentParser(description="Convert reference to new sequence base on mutations.")
    parser.add_argument('-r', '--reference', type=str, nargs='?', required=True, help="References file")
    parser.add_argument('-e', '--event_file', type=str, nargs='?', required=True, help="Edit events file")
    parser.add_argument('-f', '--filterBC', type=str, nargs='?', required=True, default="All", help="filter BCs from scRNA-seq")
    parser.add_argument('--outDir', type=str, nargs='?', help="Output results PATH")
    parser.add_argument('-n', '--name', type=str, nargs='?', help="Output files prefix")
    #parser.add_argument('-o', '--outFA', type=str, nargs='?', required=True, help="Out put FASTA format file")
    #parser.add_argument('-c', '--cpu', type=int,  nargs='?',required=True, help="The Number of cpu needed")
    args = parser.parse_args()
    return args

    
if __name__ == "__main__":
    start = time.time()
    options = Parsers()
    ## input file
    ref_seq = open(options.reference,'r').read()
    #ref_seq = open("/mnt/data5/disk/yangwj/4.Analysis_PacBio/5.CallEvents/Reference_0830.fasta",'r').read()
    event_file = open(options.event_file, 'r')
    
    ## output files
    outFA = open(os.path.join(options.outDir, options.name + "_convertSeq.fasta"), 'w')
    outFreq = open(os.path.join(options.outDir, options.name + "_EditFreq.txt"), 'w')
    outMutType = open(os.path.join(options.outDir, options.name + "_MutType.txt"), 'w')
    # write headers
    outFreq.write("{}\t{}\t{}\n".format("Mutation", "Num", "Percent"))
    outMutType.write("{}\t{}\n".format("MutType","Num"))

    ## Run
    ref = formFASTA(ref_seq)
    _, refseq = next(ref)
    if options.filterBC == "All":
        BClist = False
    else:
        BClist = getFilterBClist(options.filterBC)
    #print(BClist)
    
    ## get eventLineList
    next(event_file)
    eventLineList = []
    if BClist:
        for line in event_file:
            BC = line.split('\t')[0].split("=")[1]
            #print(BC)
            if BC in BClist:
                eventLineList.append(line)
    else:
        for line in event_file:
            eventLineList.append(line)
    #print(eventLineList)

    ##convert sequence
    eventTOunique = defaultdict(list)
    for line in eventLineList:
        spl = line.strip().split('\t')
        BC = spl[0]
        edits ='\t'.join(spl[1].split(','))
        eventTOunique[edits].append(BC)

    for edits, BCs in eventTOunique.items():
        name = eventTOunique[edits][0]
        new_sequence = convert(edits, refseq)
        outFA.write(f">{name}\n{new_sequence}\n")
    outFA.close()

    ## count the edit frequency in every base
    num_BC = len(eventLineList)
    counts = count_edit(eventLineList)
    for pos_mut, num in counts.items():
        outFreq.write("{}\t{}\t{}\n".format(pos_mut, num, num/num_BC))
    outFreq.close()
    
    ## count number of every mutation type
    mutation_counts = count_mutation(eventLineList)
    for mutation_count in mutation_counts:
        outMutType.write("{}\t{}\n".format(mutation_count[0], mutation_count[1]))
    outMutType.close()
    
    end = time.time()
    print(end - start)
    
    