#!/usr/bin/env python
# coding: utf-8

import sys
from Bio import SeqIO
sys.path.append("/mnt/data5/disk/yangwj/4.Analysis_PacBio/Scripts/")
from modules.SequenceParser import getReference

def getReferences(ref):
    ref_id = next(ref)
    ref_seq = next(ref).upper()
    reference = ref_id + ref_seq
    return reference

ref = open("/mnt/data5/disk/yangwj/4.Analysis_PacBio/4.Align_to_Ref/Barcode.fasta",'r')
reference = getReferences(ref)

output_file = open("/mnt/data5/disk/yangwj/4.Analysis_PacBio/5.CallEvents/Reference_0830.fasta", "w")
output_file.write(reference)
output_file.close()

