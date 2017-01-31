# Supplementary data

The text below is copied from the manuscript (see `paper/`):

## Supplementary Data S1

Supplementary Data S1 contains information on the 50 samples. Most variables
describe the batch processing steps outlined in Supplementary Fig. S1. “id” is a
unique identifier for each sample, “individual” is the individual identifier
(“s” = susceptible, “r” = resistant), “status” is the susceptibility status,
“treatment” is if the sample was infected or non-infected, “infection” is the 
date of the infection experiment (12 total), “arrival” is the identifier for the
arrival batch (4 total), “extraction” is the batch for RNA extraction (5 total),
“master mix” is the batch for library preparation (3 total), “rin” is the RNA
Integrity Number from the Agilent Bioanalyzer, and “outlier” is a Boolean
variable indicating if the sample was identified as an outlier (Supplementary 
Fig. S5) and removed from the analysis. (tds)

## Supplementary Data S2

Supplementary Data S2 contains the gene expression counts for the 11,336 genes
after filtering lowly expressed genes for all 50 samples (Supplementary Fig.
S2). Each row is a gene labeled with its Ensembl gene ID. Each column is a
sample. Each sample is labeled according to the pattern “x##-status-treatment”,
where x is “r” for resistant or “s” for susceptible, ## is the ID number, status
is “resist” for resistant or “suscep” for susceptible, and treatment is “noninf”
for non-infected or “infect” for infected. (tds)

## Supplementary Data S3

Supplementary Data S3 contains the results of the differential expression
analysis with limma (Fig. 1). The workbook contains 4 sheets corresponding to
the 4 tests performed. “status ni” is the test between resistant and susceptible
individuals in the non-infected state, “status ii” is the test between resistant
and susceptible individuals in the infected state, “treat resist” is the test
between the non-infected and infected states for resistant individuals, and
“treat suscep” is the test between the non-infected and infected states for
susceptible individuals. Each sheet has the same columns. “id” is the Ensembl
gene ID, “gene” is the gene name, “logFC” is the log fold change from limma,
“AveExpr” is the average log expression from limma, “t” is the t-statistic from
limma, “P.Value” is the p-value from limma, “adj.P.Val” is the adjusted p-value
from limma, “qvalue” is the q-value calculated with adaptive shrinkage, “chr” is
the chromosome where the gene is located, “description” is the description of
the gene from Ensembl, “phenotype” is the associated phenotype(s) assigned my
Ensembl, “go id” is the associated GO term(s) assigned by Ensembl, and “go
description” is the corresponding name(s) of the GO term(s). (xlsx)

## Supplementary Data S4

Supplementary Data S4 contains the results of the GWAS comparison analysis (Fig.
2). The first sheet “input-data” contains the p-values for the GWAS SNP assigned
to each gene from each study. The columns “gwas p russia”, “gwas p gambia”, 
“gwas p ghana”, “gwas p uganda”, “gwas p height” contain the p-values from the
TB susceptibility GWAS in Russia, The Gambia, Ghana, Uganda and Tanzania, and
the height GWAS in Europeans, respectively. The columns “status ni”, “status
ii”, “treat resist”, and “treat suscep” refer to the tests described for
Supplementary Data S3 and contain the log fold changes for each comparison. All
the other gene annotation columns are the same as described for Supplementary
Data S3. The second sheet “top-genes” contains a subset of the full results to
highlight those genes which had an absolute log fold change greater than 2
between resistant and susceptible individuals in the non-infected state (“status
ni”). (xlsx)

## Supplementary Data S5

Supplementary Data S5 contains the results of the classifier analysis.
Specifically it contains the results from the support vector machine using the
genes with a q-value less than 0.05 (Fig. 3). The sheet “gene-list” contains
information about the genes used for the classifier (the columns are described
in the section for Supplementary Data S3). The sheet “training-input” contains
the input gene expression data for training the model. The sheet
“training-results” contains the results of the leave-one-out-cross- validation
when training the model on the samples from the current study. The sheet
“testing-input” contains the input gene expression data for testing the model.
The sheet “testing-results” contains the results from testing the model on the
samples from Barreiro et al., 2012 24 . The column “prob tb suscep” is the
probability of being susceptible to TB assigned by the model. (xlsx)
