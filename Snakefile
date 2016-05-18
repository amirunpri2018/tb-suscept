# Snakefile

# To run on RCC Midway:
# snakemake -j 100 --cluster-config config-rcc.json -c "sbatch --mem={cluster.mem} --nodes={cluster.n} --tasks-per-node={cluster.tasks}"

import glob
import os

# Configuration ----------------------------------------------------------------

# Paths must end with forward slash
scratch = "/scratch/midway/jdblischak/"
external = "/project/gilad/jdblischak/tb-suscept/"
fastq_dir = external + "fastq/"
kallisto_dir = external + "kallisto/"
code = "code/"
data = "data/"

# Input fastq files
fastq_files = glob.glob(external + "fastq/??-*fastq.gz")
samples = [os.path.basename(f).rstrip(".fastq.gz") for f in fastq_files]
#print(samples)

for d in [scratch, external, fastq_dir, kallisto_dir]:
    if not os.path.isdir(d):
        os.mkdir(d)

# Targets ----------------------------------------------------------------------

localrules: run_kallisto, prepare_kallisto

rule run_kallisto:
    input: data + "eff-counts.txt",
           data + "tpm.txt"

rule prepare_kallisto:
    input: scratch + "transcriptome-ensembl-GRCh38.idx"

# Rules ------------------------------------------------------------------------
 
rule download_transcriptome:
    output: scratch + "transcriptome-ensembl-GRCh38.fa.gz"
    shell: "wget -O {output} http://bio.math.berkeley.edu/kallisto/transcriptomes/Homo_sapiens.GRCh38.rel79.cdna.all.fa.gz"

rule kallisto_index:
    input:  scratch + "transcriptome-ensembl-GRCh38.fa.gz"
    output: scratch + "transcriptome-ensembl-GRCh38.idx"
    shell: "kallisto index -i {output} {input}"

rule kallisto_quant:
    input: read = fastq_dir + "{sample}.fastq.gz",
           index = scratch + "transcriptome-ensembl-GRCh38.idx"
    output: kallisto_dir + "{sample}/abundance.tsv"
    params: outdir = kallisto_dir + "{sample}", length = 300, sd_len = 50
    shell: "kallisto quant -i {input.index} -o {params.outdir} --single -l {params.length} -s {params.sd_len} {input.read}"

rule kallisto_collate:
    input: expand(kallisto_dir + "{sample}/abundance.tsv", sample = samples)
    output: eff_counts = data + "eff-counts.txt",
            tpm =  data + "tpm.txt",
    params: script = code + "kallisto-collate.R"
    shell: "Rscript {params.script} {output.eff_counts} {output.tpm} {input}"
