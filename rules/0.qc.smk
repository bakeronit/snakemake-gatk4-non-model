#!/usr/bin/env python


rule fastqc:
    input:
        unpack(get_fastq)
    output:
        fwd = "analysis/qc/fastqc/{prefix}" + config["suffix"][0] + "_fastqc.zip",
        rev = "analysis/qc/fastqc/{prefix}" + config["suffix"][1] + "_fastqc.zip"
    threads: config["fastqc_threads"]
    shell:
        """
        fastqc --quiet --outdir analysis/qc/fastqc --noextract -f fastq {input} -t {threads}
        """

rule bamstats:
    input:
        "analysis/mapping/{prefix}_aligned_duplicates_marked_sorted.bam" 
    output:
        "analysis/qc/bamtools/{prefix}_bamtools.stats"
    shell:
        """
        bamtools stats -in {input} | grep -v "*" > {output}
        
        """

rule multiqc:
    input:
        expand("analysis/qc/fastqc/{prefix}{R}_fastqc.zip",prefix=samples.index,R=config["suffix"]),
        expand("analysis/qc/bamtools/{prefix}_bamtools.stats",prefix=samples.index),
        expand("analysis/mapping/{prefix}.duplicate_metrics",prefix=samples.index)
    output:
        "analysis/qc/multiqc.html"
    wrapper:
        "0.35.1/bio/multiqc"

rule genomeCov:
    input:
        "analysis/mapping/{prefix}_aligned_duplicates_marked_sorted.bam"
    output:
        "analysis/qc/coverage/{prefix}_coverage.txt"
    shell:
        """
        bamcov {input} -o {output}
        """

