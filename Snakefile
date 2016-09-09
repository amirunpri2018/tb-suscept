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
fastqc_dir = external + "fastqc/"
code = "code/"
data = "data/"
genome = scratch + "genome/"
bam_dir = external + "bam/"
counts_dir = external + "counts/"
figure = "figure/"

# Input fastq files
fastq_files = glob.glob(external + "fastq/??-*fastq.gz")
samples = [os.path.basename(f).rstrip(".fastq.gz") for f in fastq_files]
#print(samples)
chromosomes = [str(x) for x in range(1, 23)] + ["X", "Y", "M"]

for d in [scratch, external, fastq_dir, fastqc_dir, genome, bam_dir, counts_dir,
          code, data, figure]:
    if not os.path.isdir(d):
        os.mkdir(d)

# Targets ----------------------------------------------------------------------

localrules: run_analysis, run_subread, prepare_subread, qc

rule run_analysis:
    input: figure + "limma.eps",
           figure + "gwas.eps",
           figure + "classifier-svm.eps",
           data + "Supplementary_Data_S1.tds",
           data + "Supplementary_Data_S2.xlsx",
           data + "Supplementary_Data_S3.xlsx",
           data + "Supplementary_Data_S4.xlsx"

rule run_subread:
    input: data + "subread-counts-per-sample.txt", data + "total-counts.txt"

rule prepare_subread:
    input: genome + "hg38.reads", genome + "exons.saf"

rule qc:
    input: external + "multiqc_report.html"

# Rules ------------------------------------------------------------------------

# Sequence quality control

rule fastqc:
    input: fastq_dir + "{sample}.fastq.gz"
    output: fastqc_dir + "{sample}_fastqc.html"
    shell: "fastqc --outdir {fastqc_dir} {input}"

rule mutliqc:
    input: expand(fastqc_dir + "{sample}_fastqc.html", sample = samples)
    output: external + "multiqc_report.html"
    shell: "multiqc --outdir {external} {external}"

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

# Analysis

# Filter lowly expressed genes
rule filter_genes:
    input: counts = data + "subread-counts-per-sample.txt",
           info = data + "experiment-info.txt",
           script = code + "qc-genes.R"
    output: counts = data + "counts.txt",
            cpm = data + "cpm-all.rds"
    shell: "Rscript {input.script} {data}"

# Download gene annotation
rule download_gene_annotation:
    input: counts = data + "counts.txt",
           script = code + "download-gene-annotation.R"
    output: gene_anno = data + "gene-annotation.txt"
    shell: "Rscript {input.script} {data}"

# Remove outliers based on PCA
rule remove_outliers:
    input: counts = data + "counts.txt",
           info = data + "experiment-info.txt",
           script = code + "qc-outliers.R"
    output: counts = data + "counts-filtered.txt",
            info = data + "experiment-info-filtered.txt",
            pca_outliers = data + "pca-outliers.txt",
            explained_outliers = data + "explained-outliers.rds",
            outliers = data + "outliers.txt"
    shell: "Rscript {input.script} {data}"

# Check for technical batch effects
rule batch:
    input: counts = data + "counts-filtered.txt",
           info = data + "experiment-info-filtered.txt",
           script = code + "qc-batch.R"
    output: pca = data + "results-pca.txt",
            explained = data + "results-pca-explained.rds",
            covariates = data + "results-pca-covariates.txt"
    shell: "Rscript {input.script} {data}"

# Differential expression analysis with limma
rule limma:
    input: counts = data + "counts-filtered.txt",
           info = data + "experiment-info-filtered.txt",
           script = code + "main-limma.R"
    output: voom = data + "results-limma-voom.rds",
            fit = data + "results-limma-fit.rds",
            results = data + "results-limma-stats.rds"
    shell: "Rscript {input.script} {data}"

# Compare to GWAS results
rule gwas:
    input: fit = data + "results-limma-fit.rds",
           script = code + "main-gwas.R"
    output: results = data + "results-gwas.txt",
            lm_stats = data + "results-gwas-lm.txt"
    shell: "Rscript {input.script} {data}"

# Combine with lbb2012
rule combine_studies:
    input: counts = data + "counts-filtered.txt",
           info = data + "experiment-info-filtered.txt",
           lbb = data + "Exp_final_Batch_corrected.Rdata",
           script = code + "combine-studies.R"
    output: combined_anno = data + "combined-annotation.txt",
            combined_raw = data + "combined-raw.txt",
            combined_norm = data + "combined-normalized.txt",
            combined_regr = data + "combined-regressed.txt",
            combined_pca = data + "combined-pca.txt",
            combined_pca_exp = data + "combined-pca-explained.txt",    
            combined_pca_regr = data + "combined-pca-regressed.txt",
            combined_pca_regr_exp = data + "combined-pca-explained-regressed.txt",
            combined_limma = data + "combined-limma.rds"
    shell: "Rscript {input.script} {data}"

# Build classifier with caret
rule classifier:
    input: combined_anno = data + "combined-annotation.txt",
           combined_regr = data + "combined-regressed.txt",
           combined_limma = data + "combined-limma.rds",
           script = code + "main-classifier.R"
    output: predictions = data + "classifier-predictions.rds",
            predictions_lbb = data + "classifier-predictions-lbb.rds",
            training = data + "training-input.txt",
            testing = data + "testing-input.txt"
    shell: "Rscript {input.script} {data}"

rule create_figures:
    input: data +  "cpm-all.rds",
           data + "counts.txt",
           data + "counts-filtered.txt",
           pca_outliers = data + "pca-outliers.txt",
           explained_outliers = data + "explained-outliers.rds",
           outliers = data + "outliers.txt",
           pca = data + "results-pca.txt",
           explained = data + "results-pca-explained.rds",
           covariates = data + "results-pca-covariates.txt",
           gwas = data + "results-gwas.txt",
           gwas_lm = data + "results-gwas-lm.txt",
           combined_anno = data + "combined-annotation.txt",
           combined_raw = data + "combined-raw.txt",
           combined_norm = data + "combined-normalized.txt",
           combined_regr = data + "combined-regressed.txt",
           combined_pca = data + "combined-pca.txt",
           combined_pca_exp = data + "combined-pca-explained.txt",    
           combined_pca_regr = data + "combined-pca-regressed.txt",
           combined_pca_regr_exp = data + "combined-pca-explained-regressed.txt",
           predictions = data + "classifier-predictions.rds",
           predictions_lbb = data + "classifier-predictions-lbb.rds",
           script = code + "create-figures.R"
    output: figure + "gene-exp-distribution.pdf",
            figure + "gene-exp-distribution.png",
            figure + "outliers.pdf",
            figure + "outliers.png",
            figure + "heatmap-all-samples.pdf",
            figure + "heatmap-all-samples.png",
            figure + "heatmap-no-outliers.pdf",
            figure + "heatmap-no-outliers.png",
            figure + "batch-pca.pdf",
            figure + "batch-pca.png",
            figure + "batch-infection.pdf",
            figure + "batch-infection.png",
            figure + "limma.eps",
            figure + "limma.pdf",
            figure + "limma.png",
            figure + "limma-supp.pdf",
            figure + "limma-supp.png",
            figure + "gwas.eps",
            figure + "gwas.pdf",
            figure + "gwas.png",
            figure + "gwas-n-snps.pdf",
            figure + "gwas-n-snps.png",
            figure + "combined-distributions.pdf",
            figure + "combined-distributions.png",
            figure + "combined-pca.pdf",
            figure + "combined-pca.png",
            figure + "classifier-compare.pdf",
            figure + "classifier-compare.png",
            figure + "classifier-en.pdf",
            figure + "classifier-en.png",
            figure + "classifier-svm.eps",
            figure + "classifier-svm.pdf",
            figure + "classifier-svm.png",
            figure + "classifier-rf.pdf",
            figure + "classifier-rf.png",
            figure + "classifier-exp.pdf",
            figure + "classifier-exp.png"
    shell: "Rscript {input.script} {data} {figure}"

rule create_supp_data:
    input: data + "counts.txt",
           data + "results-limma-stats.rds",
           data + "gene-annotation.txt",
           data + "results-gwas.txt",
           data + "results-gwas-lm.txt",
           data + "training-input.txt",
           data + "testing-input.txt",
           data + "combined-limma.rds",
           data + "classifier-predictions.rds",
           data + "classifier-predictions-lbb.rds",
           script = code + "create-supp-data.R"
    output: data + "Supplementary_Data_S1.tds",
            data + "Supplementary_Data_S2.xlsx",
            data + "Supplementary_Data_S3.xlsx",
            data + "Supplementary_Data_S4.xlsx"
    shell: "Rscript {input.script} {data}"
