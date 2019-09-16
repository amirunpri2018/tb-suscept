# Predicting susceptibility to tuberculosis based on gene expression profiling in dendritic cells
[![Anaconda-Server Badge](https://anaconda.org/jdblischak/tb-suscept/badges/version.svg)](https://anaconda.org/jdblischak/tb-suscept)

John D. Blischak, Ludovic Tailleux, Marsha Myrthil, Cécile Charlois, Emmanuel
Bergot, Aurélien Dinh, Gloria Morizot, Olivia Chény, Cassandre Von Platen,
Jean-Louis Herrmann, Roland Brosch, Luis B. Barreiro, Yoav Gilad

## Abstract

Tuberculosis (TB) is a deadly infectious disease, which kills millions of people
every year. The causative pathogen, Mycobac- terium tuberculosis (MTB), is
estimated to have infected up to a third of the world’s population; however,
only approximately 10% of infected healthy individuals progress to active TB.
Despite evidence for heritability, it is not currently possible to predict who
may develop TB. To explore approaches to classify susceptibility to TB, we
infected with MTB dendritic cells (DCs) from putatively resistant individuals
diagnosed with latent TB, and from susceptible individuals that had recovered
from active TB. We measured gene expression levels in infected and non-infected
cells and found hundreds of differentially expressed genes between susceptible
and resistant individuals in the non-infected cells. We further found that
genetic polymorphisms nearby the differentially expressed genes between
susceptible and resistant individuals are more likely to be associated with TB 
susceptibility in published GWAS data. Lastly, we trained a classifier based on
the gene expression levels in the non-infected cells, and demonstrated decent
performance on our data and an independent data set. Overall, our promising
results from this small study suggest that training a classifier on a larger
cohort may enable us to accurately predict TB susceptibility.

## Important links

* [Final publication][scirep] ([PubMed][pmd])
* [Preprint][biorxiv] on bioRxiv
* [GSE94116][geo] - The raw FASTQ files as well as some supplementary data files are available via at the Gene Expression Omnibus.
  To download all the raw FASTQ files, run the script [code/download-fastq-sra.sh](./code/download-fastq-sra.sh).

[scirep]: https://www.nature.com/articles/s41598-017-05878-w
[pmd]: https://www.ncbi.nlm.nih.gov/pubmed/28720766
[biorxiv]: http://biorxiv.org/content/early/2017/02/03/104729
[geo]: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE94116

## Organization

For accessing summary data, please see [data/](./data).

For reading and/or re-building the manuscript, please see [paper/](./paper).

For the pipeline to convert raw sequencing data (fastq) to gene counts, please
see [Snakefile](./Snakefile).

For the analsysis that produced Fig. 1, please see
[code/main-limma.R](./code/main-limma.R).

For the analysis that produced Fig. 2, please see
[pkg/tbsuscept/R/gwas.R](./pkg/tbsuscept/R/gwas.R) and
[code/test-package.R](./code/test-package.R).

For the analsysis that produced Fig. 3, please see
[code/main-classifier.R](./code/main-classifier.R).

For accessing the figures, please see [figure/](./figure). For creating the
figures, please see [code/create-figures.R](./code/create-figures.R).

## Computational environment

To recreate a similar\* computational environment, follow these
instructions to use the [conda][] package manager. Note that this will
only work on a Linux-64 based architecture.

1. Download and install Miniconda ([instructions](https://conda.io/miniconda.html))
2. Download environment file. If you cloned the Git repo, you already have the environment file.
    - `conda install anaconda-client; anaconda download jdblischak/tb-suscept`
3. Create conda environment:
    - `conda env create -n tb-suscept --file environment.yaml`
4. Activate conda environment:
    - To activate: `source activate tb-suscept`
    - To deactivate: `source deactivate tb-suscept`

\* I started this project before I started heavily relying on
[conda][] to manage my computational environment to make it more
reproducible. The version of R on the cluster I used was R
3.2.5. Unfortunately the available R packages from the official
[conda][] repositories, as well as community provided packages from
[bioconda][] and [conda-forge][], are not available for R 3.2.5. Thus
the file `environment.yaml` creates an environment with R 3.3.2.

[conda]: https://conda.io/docs/
[bioconda]: https://bioconda.github.io
[conda-forge]: https://conda-forge.org/

## License

The code is available under the [MIT license][mit]. The text and data are
availabe under the [CC-BY license][cc-by]. See the file [LICENSE](./LICENSE) for
full details.

[mit]: https://choosealicense.com/licenses/mit/
[cc-by]: https://choosealicense.com/licenses/cc-by-4.0/

## Questions

Please feel free to open an Issue or send me an email directly if you need any
assistantance using the code or data in this repository.
