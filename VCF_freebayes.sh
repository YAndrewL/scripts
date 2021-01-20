# project path (Server) /home/data3/lyf/chengh/seq1104/cleaned

# mapping with Chlamy. DO NOT use nohup
bwa-mem2 mem ../ref/chre.fa S1_FKDL202605949-1a_1.fq.gz S1_FKDL202605949-1a_2.fq.gz > s1.sam &
bwa-mem2 mem ../ref/chre.fa S2_FKDL202605950-1a_1.fq.gz S2_FKDL202605950-1a_2.fq.gz > s2.sam &
bwa-mem2 mem ../ref/chre.fa S3_FKDL202605951-1a_1.fq.gz S3_FKDL202605951-1a_2.fq.gz > s3.sam &
bwa-mem2 mem ../ref/chre.fa S4_FKDL202605952-1a_1.fq.gz S4_FKDL202605952-1a_2.fq.gz > s4.sam &
bwa-mem2 mem ../ref/chre.fa S5_FKDL202605953-1a_1.fq.gz S5_FKDL202605953-1a_2.fq.gz > s5.sam &
bwa-mem2 mem ../ref/chre.fa S6_FKDL202605954-1a_1.fq.gz S6_FKDL202605954-1a_2.fq.gz > s6.sam &
bwa-mem2 mem ../ref/chre.fa S7_FKDL202605955-1a_1.fq.gz S7_FKDL202605955-1a_2.fq.gz > s7.sam &

for i in {1..7};
cd s${i}
samtools view -bS s${i}.sam > s${i}.bam
samtools sort s${i}.bam > s${i}.sort.bam
~/tools/gatk-4.0.12.0/gatk MarkDuplicates -I s${i}.sort.bam -O s${i}.dup.bam -M dup.txt
cat s${i}.dup.bam | bamleftalign -f ../../ref/chre.fa -m 5 > s${i}.left.bam
freebayes -T 0.001 -p 1 -J -K -X -u -m 20 -q 30 -f ../../ref/chre.fa s${i}.dup.bam > s${i}.raw.vcf
cd ..
done