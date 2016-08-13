# Snakefile

# Experimental workflow to map to MTB genome.

# To run on RCC Midway:
# snakemake -s Snake-mtb.py -kp --ri -j 500 --cluster-config config-rcc.json -c "sbatch --mem={cluster.mem} --nodes={cluster.n} --tasks-per-node={cluster.tasks}"

import glob
import os

# Configuration ----------------------------------------------------------------

# Paths must end with forward slash
scratch = "/scratch/midway/jdblischak/"
external = "/project/gilad/jdblischak/tb-suscept/"
fastq_dir = external + "fastq/"
code = "code/"
data = "data/"
genome = scratch + "genome/"
bam_dir = external + "bam/"
counts_dir = external + "counts/"
mtb_dir = scratch + "mtb/"
genome_mtb = mtb_dir + "genome/"

# Input fastq files
fastq_files = glob.glob(external + "fastq/??-*fastq.gz")
samples = [os.path.basename(f).rstrip(".fastq.gz") for f in fastq_files]
#print(samples)
chromosomes = [str(x) for x in range(1, 23)] + ["X", "Y", "M"]

for d in [scratch, external, fastq_dir, genome, bam_dir, counts_dir,
          mtb_dir, genome_mtb]:
    if not os.path.isdir(d):
        os.mkdir(d)

# Targets ----------------------------------------------------------------------

localrules: run_mtb, prepare_mtb

rule run_mtb:
    input: expand(bam_dir + "{sample}-unmapped.bam", sample = samples)


rule prepare_mtb:
    input: genome_mtb + "h37rv.ASM19595v2.reads"

# Rules ------------------------------------------------------------------------

# MTB pipeline

# ftp://ftp.ensemblgenomes.org/pub/bacteria/release-32/fasta/bacteria_0_collection/mycobacterium_tuberculosis_h37rv/dna/
rule download_genome_mtb:
    output: genome_mtb + "h37rv.ASM19595v2.fa.gz"
    shell: "wget -O {output} ftp://ftp.ensemblgenomes.org/pub/bacteria/release-32/fasta/bacteria_0_collection/mycobacterium_tuberculosis_h37rv/dna/Mycobacterium_tuberculosis_h37rv.ASM19595v2.dna_rm.chromosome.Chromosome.fa.gz"

# Should test to see if this can replace unzip_genome.
rule unzip:
    input: "{g}.fa.gz"
    output: temp("{g}.fa")
    shell: "gunzip {input}"

rule subread_index_mtb:
    input: genome_mtb + "h37rv.ASM19595v2.fa"
    output: genome_mtb + "h37rv.ASM19595v2.reads"
    params: prefix = genome_mtb + "h37rv.ASM19595v2"
    shell: "subread-buildindex -o {params.prefix} {input}"

rule obtain_unmapped:
    input: bam_dir + "{sample}.bam"
    output: bam_dir + "{sample}-unmapped.sam"
    shell: "samtools view -f 4 {input} > {output}"

rule map_unmapped:
    input: read = bam_dir + "{sample}-unmapped.sam",
           index = genome_mtb + "h37rv.ASM19595v2.reads"
    output: bam_dir + "{sample}-unmapped.bam"
    params: prefix = genome_mtb + "h37rv.ASM19595v2"
    shell: "subread-align -i {params.prefix} -r {input.read} -t 1 -u --SAMinput > {output}"
