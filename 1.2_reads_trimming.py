# A script that processes alle the fastq files in a directory with trimmomatic.

# LEADING = removes the leading bases under the threshold
# TRAILING = removes the trailing bases under the threshold
# MINLEN = drops the reads shorter than the threshold
# SLIDINGWINDOW = removes the k-mers under a quality threshold

def trimmomatic_paired(leading, trailing, slidingwindow1, slidingwindow2, minlen):
    import os
    import subprocess
    import glob

    fastqlist = glob.glob("*/*/*.fastq.gz")
    print(fastqlist)
        
    for fastqfile in fastqlist:
        print(fastqfile, " in process...")

        if "_R1_" in fastqfile:
            r1dir = fastqfile
            r2dir = r1dir.replace("R1", "R2")
                    
        elif "_R2_" in fastqfile:
            r2dir = fastqfile
            r1dir = r2dir.replace("R2", "R1")

        no_paired = r1dir.replace("_R1", "")

        args = ["trimmomatic", "PE", "-phred33",
        # Input R1, R2
        r1dir, r2dir,
        # Output
        "/trimmed_reads/" + no_paired.replace(".fastq.gz", "") + "_forward_paired.fq.gz",
        "/trimmed_reads/" + no_paired.replace(".fastq.gz", "") + "_forward_unpaired.fq.gz",
        "/trimmed_reads/" + no_paired.replace(".fastq.gz", "") + "_reverse_paired.fq.gz",
        "/trimmed_reads/" + no_paired.replace(".fastq.gz", "") + "_reverse_unpaired.fq.gz",
        "LEADING:" + str(leading),
        "TRAILING:" + str(trailing),
        "SLIDINGWINDOW:" + str(slidingwindow1) + ":" + str(slidingwindow2),
        "MINLEN:" + str(minlen)]

        subprocess.call(" ".join(args))

        fastqlist.remove(r1dir)
        fastqlist.remove(r2dir)

trimmomatic_paired(3, 3, 4, 15, 36)
