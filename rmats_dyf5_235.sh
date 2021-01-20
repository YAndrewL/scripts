#---star mapping for rMATS
STAR --readFilesIn /home/ldd/rnaseq/20181223/trimmed_data/N2-L1-1_TKR181200835_HTGWLCCXY_L8_1_val_1.fq.gz \
/home/ldd/rnaseq/20181223/trimmed_data/N2-L1-1_TKR181200835_HTGWLCCXY_L8_2_val_2.fq.gz  --chimSegmentMin 2 \
--outFilterMismatchNmax 3 --alignEndsType EndToEnd \
--runThreadN 12 --outSAMstrandField intronMotif --outSAMtype BAM SortedByCoordinate --alignSJDBoverhangMin 6 \
--alignIntronMax 299999 --genomeDir /home/data1/huangjingying/eff1/genome_235/star_index/ \
--sjdbGTFfile /home/data1/huangjingying/eff1/genome_235/annotation/Caenorhabditis_elegans.WBcel235.94.gtf \
--outFileNamePrefix starbam/N2-L1-1 --readFilesCommand zcat & 

STAR --readFilesIn /home/ldd/rnaseq/20181223/trimmed_data/N2-L1-2_TKR181200836_HTGWLCCXY_L8_1_val_1.fq.gz \
/home/ldd/rnaseq/20181223/trimmed_data/N2-L1-2_TKR181200836_HTGWLCCXY_L8_2_val_2.fq.gz  --chimSegmentMin 2 \
--outFilterMismatchNmax 3 --alignEndsType EndToEnd \
--runThreadN 12 --outSAMstrandField intronMotif --outSAMtype BAM SortedByCoordinate --alignSJDBoverhangMin 6 \
--alignIntronMax 299999 --genomeDir /home/data1/huangjingying/eff1/genome_235/star_index/ \
--sjdbGTFfile /home/data1/huangjingying/eff1/genome_235/annotation/Caenorhabditis_elegans.WBcel235.94.gtf \
--outFileNamePrefix starbam/N2-L1-2 --readFilesCommand zcat & 


STAR --readFilesIn /home/ldd/rnaseq/20180604/clean_data/cas501-1_HL2MNCCXY_L6_1_val_1.fq.gz \
/home/ldd/rnaseq/20180604/clean_data/cas501-1_HL2MNCCXY_L6_2_val_2.fq.gz  --chimSegmentMin 2 \
--outFilterMismatchNmax 3 --alignEndsType EndToEnd \
--runThreadN 12 --outSAMstrandField intronMotif --outSAMtype BAM SortedByCoordinate --alignSJDBoverhangMin 6 \
--alignIntronMax 299999 --genomeDir /home/data1/huangjingying/eff1/genome_235/star_index/ \
--sjdbGTFfile /home/data1/huangjingying/eff1/genome_235/annotation/Caenorhabditis_elegans.WBcel235.94.gtf \
--outFileNamePrefix starbam/cas501-1 --readFilesCommand zcat & 

STAR --readFilesIn /home/ldd/rnaseq/20180604/clean_data/cas501-2_HL2MNCCXY_L6_1_val_1.fq.gz \
/home/ldd/rnaseq/20180604/clean_data/cas501-2_HL2MNCCXY_L6_2_val_2.fq.gz  --chimSegmentMin 2 \
--outFilterMismatchNmax 3 --alignEndsType EndToEnd \
--runThreadN 12 --outSAMstrandField intronMotif --outSAMtype BAM SortedByCoordinate --alignSJDBoverhangMin 6 \
--alignIntronMax 299999 --genomeDir /home/data1/huangjingying/eff1/genome_235/star_index/ \
--sjdbGTFfile /home/data1/huangjingying/eff1/genome_235/annotation/Caenorhabditis_elegans.WBcel235.94.gtf \
--outFileNamePrefix starbam/cas501-2 --readFilesCommand zcat & 

#---rMATS
echo 'starbam/N2-L1-1Aligned.sortedByCoord.out.bam,starbam/N2-L1-2Aligned.sortedByCoord.out.bam' > N2-L1.txt
echo 'starbam/cas501-1Aligned.sortedByCoord.out.bam,starbam/cas501-2Aligned.sortedByCoord.out.bam' > cas501.txt

for files in $(find "./" -name '*.bam');do (samtools index $files);done &

../rMATS.4.0.2/rMATS-turbo-Linux-UCS4/rmats.py --b1 N2-L1.txt --b2 cas501.txt \
--gtf /home/data1/huangjingying/eff1/genome_235/annotation/Caenorhabditis_elegans.WBcel235.94.gtf \
--bi /home/data1/huangjingying/eff1/genome_235/star_index/ --od N2_vs_cas501_2 -t paired --nthread 6 --readLength 150 &

awk 'BEGIN {FS="\t" ; OFS="\t"}NR==1{print}$20<0.1{print}' RI.MATS.JCEC.txt > RI.MATS.JCEC.FDR.txt
sed 's/chr//g' RI.MATS.JCEC.FDR.txt > RI.MATS.JCEC.FDR.nochr.txt
rmats2sashimiplot --b1 ../starbam/N2-L1-1Aligned.sortedByCoord.out.bam,../starbam/N2-L1-2Aligned.sortedByCoord.out.bam \
--b2 ../starbam/cas501-1Aligned.sortedByCoord.out.bam,../starbam/cas501-2Aligned.sortedByCoord.out.bam \
-t RI -e ./RI.MATS.JCEC.FDR.nochr.txt --l1 N2 --l2 cas501 --exon_s 1 --intron_s 1 -o plot_RI_JCEC_cas501 >log.log &

awk 'BEGIN {FS="\t" ; OFS="\t"}NR==1{print}$20<0.1{print}' RI.MATS.JC.txt > RI.MATS.JC.FDR.txt
sed 's/chr//g' RI.MATS.JC.FDR.txt > RI.MATS.JC.FDR.nochr.txt
rmats2sashimiplot --b1 ../starbam/N2-L1-1Aligned.sortedByCoord.out.bam,../starbam/N2-L1-2Aligned.sortedByCoord.out.bam \
--b2 ../starbam/cas501-1Aligned.sortedByCoord.out.bam,../starbam/cas501-2Aligned.sortedByCoord.out.bam \
-t RI -e ./RI.MATS.JC.FDR.nochr.txt --l1 N2 --l2 cas501 --exon_s 1 --intron_s 1 -o plot_RI_JC_cas501 >log.log &

==========
Done processing each gene from dictionary to compile AS events
Found 3608 exon skipping events
Found 299 exon MX events
Found 1613 alt SS events
There are 953 alt 3 SS events and 660 alt 5 SS events.
Found 541 RI events
==========

#---
echo '/home/data3/lyf/edit_file/N2_1/all_combined.zz.sorted.bam,/home/data3/lyf/edit_file/N2_2/all_combined.zz.sorted.bam' > N2-L1.txt
echo '/home/data3/lyf/edit_file/cas501_1/all_combined.zz.sorted.bam,/home/data3/lyf/edit_file/cas501_2/all_combined.zz.sorted.bam' > cas501.txt

../rMATS.4.0.2/rMATS-turbo-Linux-UCS4/rmats.py --b1 N2-L1.txt --b2 cas501.txt \
--gtf /home/data1/huangjingying/eff1/genome_235/annotation/Caenorhabditis_elegans.WBcel235.94.gtf \
--bi /home/data1/huangjingying/eff1/genome_235/star_index/ --od N2_vs_cas501_editing -t paired --nthread 6 --readLength 150 &

awk 'BEGIN {FS="\t" ; OFS="\t"}NR==1{print}$20<0.05{print}' RI.MATS.JCEC.txt > RI.MATS.JCEC.FDR.txt
sed 's/chr//g' RI.MATS.JCEC.FDR.txt > RI.MATS.JCEC.FDR.nochr.txt
rmats2sashimiplot --b1 /home/data3/lyf/edit_file/N2_1/all_combined.zz.sorted.bam,/home/data3/lyf/edit_file/N2_2/all_combined.zz.sorted.bam \
--b2 /home/data3/lyf/edit_file/cas501_1/all_combined.zz.sorted.bam,/home/data3/lyf/edit_file/cas501_2/all_combined.zz.sorted.bam \
-t RI -e ./RI.MATS.JCEC.FDR.nochr.txt --l1 N2 --l2 cas501 --exon_s 1 --intron_s 1 -o plot_RI_JCEC_cas501 >log.log &

