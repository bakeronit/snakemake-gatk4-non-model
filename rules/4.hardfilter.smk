#!/usr/bin/env python
from datetime import datetime

"""
This workflow involves pre-processing the raw sequence data (uBAM format) to produce analysis-ready BAM files. Recalibrate Base Quality Scores is not included since our organism does not have know SNPs database.
"""



rule selectSNPs:
    input:
        "analysis/vcf/all_raw.vcf.gz"
    output:
        "analysis/vcf/all_raw_snp.vcf.gz"
    shell:
        """
        gatk SelectVariants \
        -V {input} \
        -select-type SNP \
        -O {output}
        """


rule selectINDELs:
    input:
        "analysis/vcf/all_raw.vcf.gz"
    output:
        "analysis/vcf/all_raw_indel.vcf.gz"
    shell:
        """
        gatk SelectVariants \
        -V {input} \
        -select-type INDEL \
        -O {output}
        """

rule hardfilterSNPs:
    input:
        "analysis/vcf/all_raw_snp.vcf.gz"
    output:
        "analysis/vcf/all_filtered_snp.vcf.gz"
    shell:
        """
        gatk VariantFiltration \
        -V {input} \
        -filter "QD < 2.0" --filter-name "QD2" \
        -filter "QUAL < 30.0" --filter-name "QUAL30" \
        -filter "SOR > 3.0" --filter-name "SOR3" \
        -filter "FS > 60.0" --filter-name "FS60" \
        -filter "MQ < 40.0" --filter-name "MQ40" \
        -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
        -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
        -O {output}
        """

rule hardfilterINDELs:
    input:
        "analysis/vcf/all_raw_indel.vcf.gz"
    output:
        "analysis/vcf/all_filtered_indel.vcf.gz"
    shell:
        """
        gatk VariantFiltration \
        -V {input} \
        -filter "QD < 2.0" --filter-name "QD2" \
        -filter "QUAL < 30.0" --filter-name "QUAL30" \
        -filter "FS > 200.0" --filter-name "FS200" \
        -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
        -O {output}
        """


