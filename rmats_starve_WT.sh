STAR --readFilesIn /home/data3/ZZW/RNA-SEQ/WT-01/WT-01_1_val_1.fq.gz \
/home/data3/ZZW/RNA-SEQ/WT-01/WT-01_2_val_2.fq.gz  --chimSegmentMin 2 \
--outFilterMismatchNmax 3 --alignEndsType EndToEnd \
--runThreadN 12 --outSAMstrandField intronMotif --outSAMtype BAM SortedByCoordinate --alignSJDBoverhangMin 6 \
--alignIntronMax 299999 --genomeDir /home/data1/huangjingying/eff1/genome_235/star_index/ \
--sjdbGTFfile /home/data1/huangjingying/eff1/genome_235/annotation/Caenorhabditis_elegans.WBcel235.94.gtf \
--outFileNamePrefix starbam/WT-01 --readFilesCommand zcat &

#outSAMstrandField
#alignSJDBoverhangMin : minimum overhang for annotated junctions
#alignIntronMax : maximum intron length
#
STAR --readFilesIn /home/data3/ZZW/RNA-SEQ/WT-03/WT-03_1_val_1.fq.gz \
/home/data3/ZZW/RNA-SEQ/WT-03/WT-03_2_val_2.fq.gz  --chimSegmentMin 2 \
--outFilterMismatchNmax 3 --alignEndsType EndToEnd \
--runThreadN 12 --outSAMstrandField intronMotif --outSAMtype BAM SortedByCoordinate --alignSJDBoverhangMin 6 \
--alignIntronMax 299999 --genomeDir /home/data1/huangjingying/eff1/genome_235/star_index/ \
--sjdbGTFfile /home/data1/huangjingying/eff1/genome_235/annotation/Caenorhabditis_elegans.WBcel235.94.gtf \
--outFileNamePrefix starbam/WT-03 --readFilesCommand zcat &

STAR --readFilesIn /home/data3/ZZW/RNA-SEQ/WT-41/WT-41_1_val_1.fq.gz \
/home/data3/ZZW/RNA-SEQ/WT-41/WT-41_2_val_2.fq.gz  --chimSegmentMin 2 \
--outFilterMismatchNmax 3 --alignEndsType EndToEnd \
--runThreadN 12 --outSAMstrandField intronMotif --outSAMtype BAM SortedByCoordinate --alignSJDBoverhangMin 6 \
--alignIntronMax 299999 --genomeDir /home/data1/huangjingying/eff1/genome_235/star_index/ \
--sjdbGTFfile /home/data1/huangjingying/eff1/genome_235/annotation/Caenorhabditis_elegans.WBcel235.94.gtf \
--outFileNamePrefix starbam/WT-41 --readFilesCommand zcat &

STAR --readFilesIn /home/data3/ZZW/RNA-SEQ/WT-42/WT-42_1_val_1.fq.gz \
/home/data3/ZZW/RNA-SEQ/WT-42/WT-42_2_val_2.fq.gz  --chimSegmentMin 2 \
--outFilterMismatchNmax 3 --alignEndsType EndToEnd \
--runThreadN 12 --outSAMstrandField intronMotif --outSAMtype BAM SortedByCoordinate --alignSJDBoverhangMin 6 \
--alignIntronMax 299999 --genomeDir /home/data1/huangjingying/eff1/genome_235/star_index/ \
--sjdbGTFfile /home/data1/huangjingying/eff1/genome_235/annotation/Caenorhabditis_elegans.WBcel235.94.gtf \
--outFileNamePrefix starbam/WT-42 --readFilesCommand zcat &

for files in $(find "./" -name '*.bam');do (samtools index $files);done &

echo 'starbam/WT-01Aligned.sortedByCoord.out.bam,starbam/WT-03Aligned.sortedByCoord.out.bam' > WT_1.txt
echo 'starbam/WT-41Aligned.sortedByCoord.out.bam,starbam/WT-42Aligned.sortedByCoord.out.bam' > WT_4.txt

../rMATS.4.0.2/rMATS-turbo-Linux-UCS4/rmats.py --b1 WT_1.txt --b2 WT_4.txt \
--gtf /home/data1/huangjingying/eff1/genome_235/annotation/Caenorhabditis_elegans.WBcel235.94.gtf \
--bi /home/data1/huangjingying/eff1/genome_235/star_index/ --od WT_1_vs_4 -t paired --nthread 6 --readLength 150 &

#--intron retention
awk 'BEGIN {FS="\t" ; OFS="\t"}NR==1{print}$20<0.05{print}' RI.MATS.JCEC.txt | sed 's/chr//g' > RI.MATS.JCEC.FDR0.5.nochr.txt
rmats2sashimiplot --b1 ../starbam/WT-01Aligned.sortedByCoord.out.bam,../starbam/WT-03Aligned.sortedByCoord.out.bam \
--b2 ../starbam/WT-41Aligned.sortedByCoord.out.bam,../starbam/WT-42Aligned.sortedByCoord.out.bam \
-t RI -e ./RI.MATS.JCEC.FDR0.5.nochr.txt --l1 wt_1 --l2 wt_4 --exon_s 1 --intron_s 1 -o plot_RI_JCEC_wt_4 >log.log &

for i in ;do()

grep 'exon\|CDS\|five_prime_utr\|three_prime_utr' Caenorhabditis_elegans.WBcel235.94.gtf | awk 'BEGIN {FS="\t" ; OFS="\t"} {print $3,$4,$5,$7,$8} '
grep 'exon\|CDS\|five_prime_utr\|three_prime_utr' Caenorhabditis_elegans.WBcel235.94.gtf | awk 'BEGIN {FS="\t" ; OFS="\t"} {print $9,$10} ' |head

