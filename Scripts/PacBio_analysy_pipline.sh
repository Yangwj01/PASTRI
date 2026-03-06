##====== 1.blast and CCS ======
nohup python github_deal_big_AU_pacbio_file.py "/mnt/data5/disk/yangwj/2.Rawdata_PacBio/20230810/HSC_B7/HSC_B7.subreads.bcAd1001.bam" /mnt/data4/disk/YWJ/HSC_lineageTree/HSC_B7.logfile 10 "/mnt/data5/disk/yangwj/4.Analysis_PacBio/5.CallEvents/Reference_0830.fasta" 2>&1 | tee 5.log &
samtools view -bS "/mnt/data5/disk/yangwj/2.Rawdata_PacBio/20230810/HSC_B7/HSC_B7.subreads.bcAd1001.n.sam" -o /mnt/data5/disk/yangwj/2.Rawdata_PacBio/20230810/HSC_B7/HSC_B7.subreads.bcAd1001.n.bam
samtools view -bS "/mnt/data5/disk/yangwj/2.Rawdata_PacBio/20230810/HSC_B7/HSC_B7.subreads.bcAd1001.p.sam" -o /mnt/data5/disk/yangwj/2.Rawdata_PacBio/20230810/HSC_B7/HSC_B7.subreads.bcAd1001.p.bam
ccs -j=30 --min-passes=3 "/mnt/data5/disk/yangwj/2.Rawdata_PacBio/20230810/HSC_B7/HSC_B7.subreads.bcAd1001.n.bam" /mnt/data5/disk/yangwj/4.1.Pacbio_blastr/1.CCS_blasr/HSC_B7_blasr.n.ccs.bam
ccs -j=30 --min-passes=3 "/mnt/data5/disk/yangwj/2.Rawdata_PacBio/20230810/HSC_B7/HSC_B7.subreads.bcAd1001.p.bam" /mnt/data5/disk/yangwj/4.1.Pacbio_blastr/1.CCS_blasr/HSC_B7_blasr.p.ccs.bam

##====== 2.filter CCS and Extract BC/UMI ======
bedtools bamtofastq -i /mnt/data5/disk/yangwj/4.1.Pacbio_blastr/1.CCS_blasr/HSC_B7_blasr.n.ccs.bam -fq /mnt/data5/disk/yangwj/4.1.Pacbio_blastr/2.filterAndExtractBCUMI_CCS/bamtofastq/HSC_B7_blasr.n.ccs.fastq
bedtools bamtofastq -i /mnt/data5/disk/yangwj/4.1.Pacbio_blastr/1.CCS_blasr/HSC_B7_blasr.p.ccs.bam -fq /mnt/data5/disk/yangwj/4.1.Pacbio_blastr/2.filterAndExtractBCUMI_CCS/bamtofastq/HSC_B7_blasr.p.ccs.fastq
python /mnt/data5/disk/yangwj/4.Analysis_PacBio/Scripts/2.filterAndExtractBCUMI_CCS_PacBio.py -i "/mnt/data5/disk/yangwj/4.1.Pacbio_blastr/2.filterAndExtractBCUMI_CCS/bamtofastq/HSC_B7_blasr.n.ccs.fastq" -o /mnt/data5/disk/yangwj/4.1.Pacbio_blastr/2.filterAndExtractBCUMI_CCS/fasta_extractBC_UMI/HSC_B7.n.filter.fasta -s /mnt/data5/disk/yangwj/4.1.Pacbio_blastr/2.filterAndExtractBCUMI_CCS/fasta_extractBC_UMI/HSC_B7.n.staDF.txt --nopass /mnt/data5/disk/yangwj/4.1.Pacbio_blastr/2.filterAndExtractBCUMI_CCS/fasta_extractBC_UMI/noPass_HSC_B7.n.fa --format FASTQ --FP TGGACGAGCTGTACAAGTAA --RP AGATCGGAAGAGCGTCGTGTAG --inFP TCTCTATACGATCCGGACCT --inRP --pAw 12 --pMs 3 --mismatch 2 --numcpu 10
python /mnt/data5/disk/yangwj/4.Analysis_PacBio/Scripts/2.filterAndExtractBCUMI_CCS_PacBio.py -i "/mnt/data5/disk/yangwj/4.1.Pacbio_blastr/2.filterAndExtractBCUMI_CCS/bamtofastq/HSC_B7_blasr.p.ccs.fastq" -o /mnt/data5/disk/yangwj/4.1.Pacbio_blastr/2.filterAndExtractBCUMI_CCS/fasta_extractBC_UMI/HSC_B7.p.filter.fasta -s /mnt/data5/disk/yangwj/4.1.Pacbio_blastr/2.filterAndExtractBCUMI_CCS/fasta_extractBC_UMI/HSC_B7.p.staDF.txt --nopass /mnt/data5/disk/yangwj/4.1.Pacbio_blastr/2.filterAndExtractBCUMI_CCS/fasta_extractBC_UMI/noPass_HSC_B7.p.fa --format FASTQ --FP TGGACGAGCTGTACAAGTAA --RP AGATCGGAAGAGCGTCGTGTAG --inFP TCTCTATACGATCCGGACCT --inRP --pAw 12 --pMs 3 --mismatch 2 --numcpu 10
python /mnt/data5/disk/yangwj/4.Analysis_PacBio/Scripts/2.supp.MergeR1R2_filterCCS_PacBio.py -1 "/mnt/data5/disk/yangwj/4.1.Pacbio_blastr/2.filterAndExtractBCUMI_CCS/fasta_extractBC_UMI/HSC_B7.n.filter.fasta" -2 "/mnt/data5/disk/yangwj/4.1.Pacbio_blastr/2.filterAndExtractBCUMI_CCS/fasta_extractBC_UMI/HSC_B7.p.filter.fasta" -o /mnt/data5/disk/yangwj/4.1.Pacbio_blastr/2.filterAndExtractBCUMI_CCS/fasta_Merge.n.p/HSC_B7.union.filter.blasr.fasta

##====== 3.adjust BC base on 10X single cell BC ======
python /mnt/data5/disk/yangwj/4.Analysis_PacBio/Scripts/3.adjustBCandUMI_baseOn_scRNAseqBC_PacBio.py -i "/mnt/data5/disk/yangwj/4.1.Pacbio_blastr/2.filterAndExtractBCUMI_CCS/fasta_Merge.n.p/HSC_B7.union.filter.blasr.fasta" -o /mnt/data5/disk/yangwj/4.1.Pacbio_blastr/3.newFilter_BaseOn_sc/HSC_B7.adjust.base.sc.fa -f "/mnt/data5/disk/yangwj/Result_RData_Rds/seuratObj_Mt0.25/HSC_B7_filterBCs.txt"
python /mnt/data5/disk/yangwj/4.Analysis_PacBio/Scripts/4.1.GetConsensusSeqFromCCS_PacBio.py -i /mnt/data5/disk/yangwj/4.1.Pacbio_blastr/3.newFilter_BaseOn_sc/HSC_B7.adjust.base.sc.fa -o /mnt/data5/disk/yangwj/4.1.Pacbio_blastr/3.newFilter_BaseOn_sc/HSC_B7.union.UMIcon.sc.fa --numcpu 30

##====== 4.call Mutations ======
python /mnt/data5/disk/yangwj/4.Analysis_PacBio/Scripts/6.Call_mutation_keepzwm.py -r /mnt/data5/disk/yangwj/4.Analysis_PacBio/5.CallEvents/Reference_0830.fasta -s "/mnt/data5/disk/yangwj/4.1.Pacbio_blastr/3.newFilter_BaseOn_sc/HSC_B7.union.UMIcon.sc.fa" -o /mnt/data5/disk/yangwj/4.1.Pacbio_blastr/4.newCall_Mutations/HSC_B7.union.events.sc.txt -c 30

##====== 5.filter by numIndels <=20 and zws >=2 ======
python "/mnt/data5/disk/yangwj/4.Analysis_PacBio/Scripts/7.2.Events_filterIndels_baseRef.py" -i "/mnt/data5/disk/yangwj/4.1.Pacbio_blastr/4.newCall_Mutations/HSC_B7.union.events.sc.txt" -o /mnt/data5/disk/yangwj/4.1.Pacbio_blastr/4-6.filter_numIndel_zws/Indels20_zws2_7new/4.Mutations/HSC_B7.filter.union.events.sc.txt
cat HSC_B7.filter.union.events.sc.txt | grep '>' | awk '{print $2,$3,$4,$5}' | sed 's/zwCCS.num=//g' | awk '{a[$1]+=$3}END{for(i in a){print i,a[i]}}' | awk '$2>=2' > fi.zws.HSC_B7.bc.txt
awk 'NR==FNR{a[$1]; next} /^>/ {id=$2; if(id in a) {print; getline; print $0}}' fi.zws.HSC_B7.bc.txt HSC_B7.filter.union.events.sc.txt > fi.indel20.zws2.HSC_B7.union.Events.txt

##====== 6.get consensus Edit Events per cell======
python "/mnt/data5/disk/yangwj/4.Analysis_PacBio/Scripts/7.new_Get_consensus_Events_BC_comprehen_UMIandSCORE_with_filterBC.py" -i "/mnt/data5/disk/yangwj/4.1.Pacbio_blastr/4-6.filter_numIndel_zws/Indels20_zws2_7new/4.Mutations/fi.indel20.zws2.HSC_B7.union.Events.txt" --outDir "/mnt/data5/disk/yangwj/4.1.Pacbio_blastr/4-6.filter_numIndel_zws/Indels20_zws2_7new/5.ConEvents/" -n HSC_B7
cat HSC_B7_comSCOREandUMIconEvents.txt | awk -F"\t" -v OFS="\t" '$2!="" {print}' > HSC_B7_EditconEvents.txt

##====== 7. convert events to sequence for Iqtree input======
python "/mnt/data5/disk/yangwj/4.Analysis_PacBio/Scripts/8.2.ConvertSeq_unique_Mutation.py" -r /mnt/data5/disk/yangwj/4.Analysis_PacBio/5.CallEvents/Reference_0830.fasta -e "/mnt/data5/disk/yangwj/4.1.Pacbio_blastr/4-6.filter_numIndel_zws/Indels20_zws2_7new/5.ConEvents/HSC_B7_EditconEvents.txt" -f "/mnt/data5/disk/yangwj/Result_RData_Rds/seuratObj_Mt0.25/HSC_B7_filterBCs.txt" --outDir "/mnt/data5/disk/yangwj/4.1.Pacbio_blastr/4-6.filter_numIndel_zws/Indels20_zws2_7new/6.iqtree/" -n HSC_B7

##====== 8.get tree by Iqtree======
/home/yangwj/Download/iqtree-2.2.2.6-Linux/bin/iqtree2 -s "/mnt/data5/disk/yangwj/4.1.Pacbio_blastr/4-6.filter_numIndel_zws/Indels20_zws2_7new/6.iqtree/HSC_B7_convertSeq.fasta" -m GTR+FO+I+R10 -pre "/mnt/data5/disk/yangwj/4.1.Pacbio_blastr/4-6.filter_numIndel_zws/Indels20_zws2_7new/6.iqtree/HSC_B7/HSC_B7" -nt AUTO -czb -alrt 1000 -B 1000
python "/mnt/data5/disk/yangwj/4.Analysis_PacBio/Scripts/9.change_IQtree.py" -i "/mnt/data5/disk/yangwj/4.1.Pacbio_blastr/4-6.filter_numIndel_zws/Indels20_zws2_7new/6.iqtree/HSC_B7/HSC_B7.treefile" -d /mnt/data5/disk/yangwj/4.1.Pacbio_blastr/4-6.filter_numIndel_zws/Indels20_zws2_7new/6.iqtree/HSC_B7/ -n HSC_B7



