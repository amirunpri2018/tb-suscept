# Snakefile

# To run on RCC Midway:
# snakemake -kp --ri -j 500 --cluster-config config-rcc.json -c "sbatch --mem={cluster.mem} --nodes={cluster.n} --tasks-per-node={cluster.tasks}"

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
genome = scratch + "genome/"
bam_dir = external + "bam/"
counts_dir = external + "counts/"

# Input fastq files
fastq_files = glob.glob(external + "fastq/??-*fastq.gz")
samples = [os.path.basename(f).rstrip(".fastq.gz") for f in fastq_files]
#print(samples)
chromosomes = [str(x) for x in range(1, 23)] + ["X", "Y", "M"]

for d in [scratch, external, fastq_dir, kallisto_dir, genome, bam_dir, counts_dir]:
    if not os.path.isdir(d):
        os.mkdir(d)

# Targets ----------------------------------------------------------------------

localrules: run_kallisto, prepare_kallisto, run_subread, prepare_subread

rule run_kallisto:
    input: data + "eff-counts.txt",
           data + "tpm.txt"

rule prepare_kallisto:
    input: scratch + "transcriptome-ensembl-GRCh38.idx"

rule run_subread:
    input: data + "subread-counts-per-sample.txt"

rule prepare_subread:
    input: genome + "hg38.reads", genome + "exons.saf"

# Rules ------------------------------------------------------------------------

# Kallisto pipeline

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

# Subread pipeline

rule download_genome:
    output: genome + "chr{chr}.fa.gz"
    params: chr = "{chr}"
    shell: "wget -O {output} http://hgdownload.cse.ucsc.edu/goldenPath/hg38/chromosomes/chr{params.chr}.fa.gz"

rule unzip_chromosome_fasta:
    input: genome + "chr{chr}.fa.gz"
    output: temp(genome + "chr{chr}.fa")
    shell: "gunzip {input}"

rule subread_index:
    input: expand(genome + "chr{chr}.fa", chr = chromosomes)
    output: genome + "hg38.reads"
    params: prefix = genome + "hg38"
    shell: "subread-buildindex -o {params.prefix} {input}"

rule subread_align:
    input: read = fastq_dir + "{sample}.fastq.gz",
           index = genome + "hg38.reads"
    output: bam_dir + "{sample}.bam"
    params: prefix = genome + "hg38", threads = 8
    shell: "subread-align -i {params.prefix} -r {input.read} -t 0 -u -T {params.threads} > {output}"
    
rule create_exons_saf:
    output: genome + "exons.saf"
    params: script = code + "create-exons.R"
    shell: "Rscript {params.script} >  {output}"

rule subread_feauturecounts:
    input: bam = bam_dir + "{sample}.bam",
           exons = genome + "exons.saf"
    output: counts_dir + "{sample}.genecounts.txt"
    shell: "featureCounts -a {input.exons} -F SAF -R -o {output} {input.bam}"

rule subread_collate_per_lane:
    input: expand(counts_dir + "{sample}.genecounts.txt", sample = samples)
    output: data + "subread-counts-per-lane.txt"
    run:
        files = input
        files.sort()
        outfile = open(output[0], "w")
        
        # Get gene names from first file
        gene_list = []
        f = files[0]
        handle = open(f, "r")
        for line in handle:
            if line[0] == "#" or line[:6] == "Geneid":
                continue
            cols = line.strip("\n").split("\t")
            gene = cols[0]
            gene_list.append(gene)
        handle.close()

        # Create header row
        header = "individual\tstatus\ttreatment\tflow_cell\tlane\t" + \
                 "\t".join(gene_list) + "\n"
        outfile.write(header)

        # Process each input genecounts file
        for f in files:
            # Get counts from f
            g = 0 # iterator for indexing gene names
            counts = ["NA"] * len(gene_list)
            handle = open(f, "r")
            for line in handle:
                if line[0] == "#" or line[:6] == "Geneid":
                    continue
                cols = line.strip("\n").split("\t")
                gene = cols[0]
                assert gene == gene_list[g], \
                       "Gene names are not in correct order in file: %s"%(f)
                counts[g] = cols[6]
                g += 1
            handle.close()

            # Get meta data from filename
            path_parts = f.split("/")
            fname = path_parts[-1].rstrip(".genecounts.txt")
            fname_parts = fname.split("-")
            num, status, treatment = fname_parts[:3]
            if status == "tb":
                individual = "t" + num
            else:
                individual = "c" + num
            flow_cell, lane = fname_parts[3:5]

            # Write line
            outfile.write(individual + "\t" + status + "\t" + \
                          treatment + "\t" + flow_cell + "\t" + \
                          lane + "\t" + "\t".join(counts) + "\n")

        outfile.close()

rule subread_collate_per_sample:
    input: data + "subread-counts-per-lane.txt"
    output: data + "subread-counts-per-sample.txt"
    run:
        import pandas as pd

        lane = pd.read_table(input[0])
        sample = lane.groupby(["individual", "status", "treatment"],
                              as_index = False).sum()
        sample.to_csv(output[0], sep = "\t", na_rep = "NA", index = False)
