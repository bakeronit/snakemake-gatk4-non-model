#!/usr/bin/env python
import pandas as pd

configfile: "config.yaml"

samples = pd.read_csv("metadata.tsv",sep='\t').set_index(["prefix"], drop=False)

REFERENCE=config["GENOME"]
scaffolds = []

with open(REFERENCE,'rt') as fh:
    for line in fh:
        line = line.strip()
        if line.startswith(">"):
            line = line.split(" ")[0]
            scaffolds.append(line[1:])

wildcard_constraints:
    prefix = "|".join(samples.index),
    sample = "|".join(samples["sample"]),
    scaffold = "|".join(scaffolds)


def get_fastq(wildcards):
    fastqs = samples.loc[(wildcards.prefix), ["fq1", "fq2"]].dropna()
    return {"r1": fastqs.fq1, "r2": fastqs.fq2}

trim = config["trim"]
def get_fastq2(wildcards):
    fastqs = samples.loc[(wildcards.prefix), ["fq1", "fq2"]].dropna()
    if trim:
        return {"r1": os.path.join("data/trimmed",wildcards.prefix + ".1P.fastq.gz"), "r2": os.path.join("data/trimmed",wildcards.prefix + ".2P.fastq.gz")}
    else:
        return {"r1": fastqs.fq1, "r2": fastqs.fq2}


include: "rules/0.qc.smk",
include: "rules/1.fastq2ubam.smk",
include: "rules/2.preprocessing.smk",
include: "rules/3.joinCalling.smk",
include: "rules/4.hardfilter.smk"

rule all:
    input:
        "analysis/qc/multiqc.html",
        expand("analysis/qc/coverage/{prefix}_coverage.txt",prefix=samples.index),
        "analysis/vcf/all_filtered_snp.vcf.gz",
        "analysis/vcf/all_filtered_indel.vcf.gz"

rule clean:
    shell:
        "rm -rf analysis data/ubam data/trimmed sample_map.txt logs"

rule dbclean:
    shell:
        "rm -rf analysis/genomicsDB/*.db"
