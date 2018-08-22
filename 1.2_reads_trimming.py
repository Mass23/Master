# A script that processes alle the fastq files in a directory with trimmomatic.

# LEADING = removes the leading bases under the threshold
# TRAILING = removes the trailing bases under the threshold
# MINLEN = drops the reads shorter than the threshold
# SLIDINGWINDOW = removes the k-mers under a quality threshold

def trimmomatic_paired(leading, trailing, slidingwindow1, slidingwindow2, minlen):
    import subprocess
    import glob
    import os

    fastqlist = glob.glob("*.fastq.gz")
        
    for fastqfile in fastqlist:

        if "_R1" in fastqfile:
            r1dir = fastqfile
            r2dir = r1dir.replace("_R1", "_R2")
                    
        elif "_R2" in fastqfile:
            r2dir = fastqfile
            r1dir = r2dir.replace("_R2", "_R1")

        r1id = r1dir.split("/")[-1].replace(".fastq.gz","")
        r2id = r2dir.split("/")[-1].replace(".fastq.gz","")
        curr_dir = "/".join(r1dir.split("/")[:-1])

        no_paired = r1id.replace("_R1", "")
        no_paired = no_paired.split("_")[0] + "_" + no_paired.split("_")[2] + "_" + no_paired.split("_")[3]

        args = ["trimmomatic", "PE", "-threads 2", "-phred64",
        # Input R1, R2
        r1id + ".fastaq.gz", r2id + ".fastaq.gz",
        # Output forward/reverse, paired/unpaired
        no_paired + "_forward_paired.fq.gz",
        no_paired + "_forward_unpaired.fq.gz",
        no_paired + "_reverse_paired.fq.gz",
        no_paired + "_reverse_unpaired.fq.gz",
        "ILLUMINACLIP:TruSeq3-PE.fa:2:30:10",
        "LEADING:" + str(leading),
        "TRAILING:" + str(trailing),
        "SLIDINGWINDOW:" + str(slidingwindow1) + ":" + str(slidingwindow2),
        "MINLEN:" + str(minlen)]

        subprocess.call(" ".join(args), shell = True)
        print(no_paired, " done!")

trimmomatic_paired(3, 3, 4, 15, 36)
