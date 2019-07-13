#!/usr/bin/env python
from datetime import datetime

"""
This workflow involves pre-processing the raw sequence data (uBAM format) to produce analysis-ready BAM files. Recalibrate Base Quality Scores is not included since our organism does not have know SNPs database.
"""

#today = datetime.today().strfmt('%Y-%m-%d')

rule WriteSampleMap:
    input:
        "metadata.tsv",
        expand("analysis/gvcf/{prefix}.g.vcf.gz",prefix=samples.index)
    output:
        "sample_map.txt"
    run:
        #import pandas as pd
        #samples = pd.read_csv("metadata.tsv",sep='\t').set_index(["prefix"], drop=False)
        files = ["analysis/gvcf/" + i + ".g.vcf.gz" for i in samples.index]
        sample = [os.path.basename(i).split(".")[0] for i in files]
        out = open(output[0],'w')
        for c1,c2 in zip(sample, files):
            print("%s\t%s"%(c1,c2),file = out)


rule GenomicsDBImport:
    input:
        "sample_map.txt"
    output:
        db = directory("analysis/genomicsDB/{scaffold}.db"),
        tar = "analysis/genomicsDB/{scaffold}.db.tar"
    shell:
        """
        gatk --java-options "-Xmx4g -Xms4g -Djava.io.tmpdir=/scratch/jc502059/tmp" \
        GenomicsDBImport \
        --genomicsdb-workspace-path {output.db} \
        --L {wildcards.scaffold} \
        --sample-name-map {input} \
        --reader-threads 5 \
        --batch-size 50

        tar -cf {output.tar} {output.db}
        """

rule GenotypeGVCFs:
    input:
        "analysis/genomicsDB/{scaffold}.db"
    output:
        "analysis/vcf/{scaffold}.vcf.gz"
    shell:
        """
        gatk --java-options "-Xmx8g -Xms8g -Djava.io.tmpdir=/scratch/jc502059/tmp" \
        GenotypeGVCFs \
        -R {REFERENCE} \
        -O {output} \
        --only-output-calls-starting-in-intervals \
        --use-new-qual-calculator \
        -V gendb://{input} \
        -L {wildcards.scaffold}
        """

rule GatherVcfs:
    input:
        expand("analysis/vcf/{scaffold}.vcf.gz",scaffold = scaffolds)
    output:
        file="analysis/vcf/all_raw.vcf.gz",
        index="analysis/vcf/all_raw.vcf.gz.tbi"
    params:
        " -I ".join("analysis/vcf/" + s + ".vcf.gz" for s in scaffolds)
    shell:
        """
        gatk --java-options "-Xmx6g -Xms6g" \
        GatherVcfsCloud \
        --ignore-safety-checks \
        --gather-type BLOCK \
        -I {params} \
        -O {output.file}

        gatk --java-options "-Xmx6g -Xms6g" \
        IndexFeatureFile \
        -F {output.file} \
        -O {output.index}
        """
    
