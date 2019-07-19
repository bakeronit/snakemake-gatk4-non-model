#!/usr/bin/env python

def get_read_group(wildcards):
    flowcell_barcode = samples.loc[(wildcards.prefix),"flowcell_barcode"]
    lane = samples.loc[(wildcards.prefix),"lane"]
    sample_barcode = samples.loc[(wildcards.prefix),"sample_barcode"]
    sampleid = samples.loc[(wildcards.prefix),"sample"]
    return {"flowcell_barcode":flowcell_barcode,"lane":lane,"sample_barcode":sample_barcode,"sampleid":sampleid}

# this wrapper has a problem with multiple threads, it will use pigz if compressed file provided.Not available on this machine. Not sure about the reason. TODO:remove wrapper.
rule trim_pe:
    input:
        unpack(get_fastq)
    output:
        r1 = "data/trimmed/{prefix}.1P.fastq.gz",
        r2 = "data/trimmed/{prefix}.2P.fastq.gz",
        r1_unpaired="data/trimmed/{prefix}.1U.fastq.gz",
        r2_unpaired="data/trimmed/{prefix}.2U.fastq.gz"
    log:
        "logs/trimmomatic/{prefix}.log"
    params:
        trimmer=["LEADING:3","TRAILING:3","SLIDINGWINDOW:4:15","MINLEN:25"],
        compression_level="-9"
    wrapper:
        "0.35.2/bio/trimmomatic/pe"


rule fastq2ubam:
    input:
        unpack(get_fastq2)
    output:
        "data/ubam/{prefix}.unaligned_reads.bam"
    params:
        flowcell_barcode = lambda wildcards : get_read_group(wildcards)['flowcell_barcode'],
        lane=lambda wildcards : get_read_group(wildcards)['lane'],
        sample_barcode = lambda wildcards : get_read_group(wildcards)['sample_barcode'],
        sampleid = lambda wildcards : get_read_group(wildcards)['sampleid']
	platform = config['platform']
    shell:
       """ 
        picard FastqToSam \
         FASTQ={input.r1} \
         FASTQ2={input.r2} \
         OUTPUT={output} \
         READ_GROUP_NAME={params.flowcell_barcode}{params.lane} \
         PLATFORM=illumina \
         PLATFORM_UNIT={params.flowcell_barcode}.{params.lane}.{params.sample_barcode} \
         LIBRARY_NAME={params.sampleid}_{params.sample_barcode} \
         SAMPLE_NAME={params.sampleid} \
         PLATFORM_MODEL={params.platform}
        """
