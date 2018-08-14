# SAM to BAM
java -jar picard.jar SortSam \ 
    INPUT=aligned_reads.sam \ 
    OUTPUT=sorted_reads.bam \ 
    SORT_ORDER=coordinate 

# Mark duplicates
java -Xmx1g -jar /apps/PICARD/1.95/MarkDuplicates.jar \
                            MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000\
                            METRICS_FILE=out.metrics \
                            REMOVE_DUPLICATES=true \
                            ASSUME_SORTED=true  \
                            VALIDATION_STRINGENCY=LENIENT \
                            INPUT=out.sorted.bam \
                            OUTPUT=out.dedupe.bam
