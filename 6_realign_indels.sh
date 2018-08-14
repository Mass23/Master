 java -jar GenomeAnalysisTK.jar \
   -T IndelRealigner \
   -R reference.fasta \
   -I input.bam \
   -known indels.vcf \
   -targetIntervals intervalListFromRTC.intervals \
   -o realignedBam.bam
