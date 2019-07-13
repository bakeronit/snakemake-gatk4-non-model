# snakemake-gatk4-non-model
#### germline variant calling
Working with non-model organism means you don't have known SNPs and prepared interval list in [GATK bundle](https://software.broadinstitute.org/gatk/download/bundle). Alternatively I did:
- **Optional Base Quality Score Calibration(BQSR)** (working...)
- Apply hard filters to call sets.
- Split genome into scaffolds as intervals.

How to run:

0. `conda env create -f environment.yaml` and install [bamcov](https://github.com/fbreitwieser/bamcov) ...
1. prepare metadata.tsv.
2. modify config.yaml if needed.
3. `snakemake` 

To run on JCU HPC with PBSpro:

`snakemake -p --cluster-config jcu_hpc.json --cluster "qsub -j oe -l walltime={cluster.time} -l select=1:ncpus={cluster.ncpus}:mem={cluster.mem}" --jobs 100 --latency-wait 5`


Reference:

https://snakemake.readthedocs.io/en/stable/
https://github.com/gatk-workflows/gatk4-germline-snps-indels
https://github.com/snakemake-workflows/dna-seq-gatk-variant-calling
https://zhuanlan.zhihu.com/p/33891718
https://software.broadinstitute.org/gatk/documentation/article?id=11097
https://gatkforums.broadinstitute.org/gatk/discussion/12443/genomicsdbimport-run-slowly-with-multiple-samples
