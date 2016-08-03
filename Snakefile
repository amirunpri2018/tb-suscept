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
mtb_dir = scratch + "mtb/"
genome_mtb = mtb_dir + "genome/"

# Input fastq files
fastq_files = glob.glob(external + "fastq/??-*fastq.gz")
samples = [os.path.basename(f).rstrip(".fastq.gz") for f in fastq_files]
#print(samples)
chromosomes = [str(x) for x in range(1, 23)] + ["X", "Y", "M"]

for d in [scratch, external, fastq_dir, kallisto_dir, genome, bam_dir, counts_dir,
          mtb_dir, genome_mtb]:
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
    input: data + "subread-counts-per-sample.txt", data + "total-counts.txt"

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
    shell: "featureCounts -a {input.exons} -F SAF -o {output} {input.bam}"

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

rule count_fastq:
    input: fastq_dir + "{sample}.fastq.gz"
    output: fastq_dir + "{sample}.count.txt"
    shell: "bioawk -c fastx 'END{{print NR}}' {input} > {output}"

rule count_bam:
    input: bam_dir + "{sample}.bam"
    output: bam_dir + "{sample}.count.txt"
    shell: "samtools view -c -q 1 {input} > {output}"

rule count_gather:
    input: fastq = expand(fastq_dir + "{sample}.count.txt", sample = samples),
           bam = expand(bam_dir + "{sample}.count.txt", sample = samples),
           featureCounts = expand(counts_dir + "{sample}.genecounts.txt.summary",
                                  sample = samples)
    output: data + "total-counts.txt"
    run:
        import glob
        import sys

        files = input.fastq + input.bam + input.featureCounts

        # Output file
        outfile = open(output[0], "w")

        # Output header:
        #   stage - the stage of the processing pipeline
        outfile.write("stage\tid\tstatus\ttreatment\tflow_cell\tlane\tcounts\n")

        # Function definitions ---------------------------------------------------------

        def read_lines(f):
            """
            Input: Path to file
            Output: Lines of the file (list)
            """
            handle = open(f, "r")
            lines = handle.readlines()
            handle.close()
            return lines

        def read_count_file(f):
            """
            Input: Path to file
            Output: The number contained in the file (str)
            Explanation: The count file only contains one number, the number of
              sequences at that stage in the processing pipeline.
            """
            lines = read_lines(f)
            assert len(lines) == 1, "Count file has only one line"
            counts = lines[0].strip("\n")
            return counts

        def read_featureCounts_summary(f):
            """
            Input: Path to file
            Output: The number of sequences labeled Assigned (str)
            Explanation: featureCounts outputs a file with the extension .summmary that
                 details the number of sequences per result category. The
                 category Assigned is for sequences that map unambiguously to a
                 gene.
            """
            assert  f[-8:] == ".summary", \
                "featureCounts summary file has correct extension"
            lines = read_lines(f)
            assert len(lines) == 12, \
                "featureCounts summary file has 12 lines"
            assert lines[1].split("\t")[0] == "Assigned", \
                "The Assigned category is the first entry after the header"
            counts = lines[1].strip("\n").split("\t")[1]
            return counts

        # Process each file ------------------------------------------------------------

        for f in files:
            # Set default values
            path = f.split("/")
            fname = path[-1]
            dir = path[-2]

            if dir == "fastq":
                stage = "raw"
                counts = read_count_file(f)
            elif dir == "bam":
                stage = "mapped to genome"
                counts = read_count_file(f)
            elif dir == "counts":
                stage = "mapped to exons"
                counts = read_featureCounts_summary(f)

            # Get meta data from filename
            fname_parts = fname.rstrip("genecounts.txt.summary").split("-")
            num, status, treatment, flow_cell, lane = fname_parts[:5]
            id = "-".join([num, status, treatment])

            outfile.write("\t".join([stage, id, status, treatment,
                                     flow_cell, lane, counts]) + "\n")

        outfile.close()

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

rule prepare_mtb:
    input: genome_mtb + "h37rv.ASM19595v2.reads"
