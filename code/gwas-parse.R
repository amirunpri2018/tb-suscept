
library("dplyr")

if(interactive()) {
  data_dir <- "../data"
} else {
  args <- commandArgs(trailingOnly = TRUE)
  data_dir <- args[1]
}
stopifnot(dir.exists(data_dir))

# Load GWAS enrichment data -----------------------------------------------

result_ghana <- readRDS(file.path(data_dir, "gwas-ghana.rds"))
result_gambia <- readRDS(file.path(data_dir, "gwas-gambia.rds"))
result_height <- readRDS(file.path(data_dir, "gwas-height.rds"))
result_russia <- readRDS(file.path(data_dir, "gwas-russia.rds"))
result_sobota <- readRDS(file.path(data_dir, "gwas-sobota.rds"))

anno_gene <- read.delim(file.path(data_dir, "gene-annotation.txt"),
                        stringsAsFactors = FALSE)

# Collate -----------------------------------------------------------------

d_ghana <- data.frame(gene = names(result_ghana$status_ni$input$y),
                       de = result_ghana$status_ni$input$x,
                       p = result_ghana$status_ni$input$y)
rownames(d_ghana) <- d_ghana$gene
d_gambia <- data.frame(gene = names(result_gambia$status_ni$input$y),
                       de = result_gambia$status_ni$input$x,
                       p = result_gambia$status_ni$input$y)
rownames(d_gambia) <- d_gambia$gene

d_russia <- data.frame(gene = names(result_russia$status_ni$input$y),
                       de = result_russia$status_ni$input$x,
                       p = result_russia$status_ni$input$y)
rownames(d_russia) <- d_russia$gene
d_height <- data.frame(gene = names(result_height$status_ni$input$y),
                       de = result_height$status_ni$input$x,
                       p = result_height$status_ni$input$y)
rownames(d_height) <- d_height$gene

d_sobota <- data.frame(gene = names(result_sobota$status_ni$input$y),
                       de = result_sobota$status_ni$input$x,
                       p = result_sobota$status_ni$input$y)
rownames(d_sobota) <- d_sobota$gene


d1 <- merge(d_ghana, d_gambia, by = "gene", suffixes = c("_ghana", "_gambia"),
            all = TRUE)
d2 <- merge(d_russia, d_height, by = "gene", suffixes = c("_russia", "_sobota"),
            all = TRUE)
d <- merge(d1, d2, by = "gene", all = TRUE)
rownames(d) <- d$gene

# Rank genes --------------------------------------------------------------

d <- d %>% mutate(rank_ghana = rank(p_ghana) / n(),
                  rank_gambia = rank(p_gambia) / n(),
                  rank_russia = rank(p_russia) / n(),
                  rank_sobota = rank(p_sobota) / n(),
                  rank_de = rank(de_ghana) / n())
rownames(d) <- d$gene
d["ENSG00000108702", ]
d["ENSG00000130477", ]

d %>% filter(de_ghana > 2, p_ghana < .1, p_gambia < .1, p_russia < .1)

d %>% filter(de_ghana > 2, rank_ghana < .1, rank_gambia < .1, rank_russia < .1,
             rank_height > .1)

d %>% filter(de_ghana > 2, rank_ghana < .25, rank_gambia < .25, rank_height > .1)

# ENSG00000130477
d %>% filter(de_ghana > 2, p_ghana < .01, p_gambia < .01, p_russia < .01)

# ENSG00000108702
d %>% filter(de_ghana > 2, p_ghana < .01, p_gambia < .01, p_russia < .1)

d %>% filter(de_ghana > 2) %>% arrange(p_ghana) %>% select(gene) %>% head(1)
d %>% filter(de_ghana > 2) %>% arrange(p_gambia) %>% select(gene) %>% head(1)
d %>% filter(de_ghana > 2) %>% arrange(p_russia) %>% select(gene) %>% head(1)
d %>% filter(de_ghana > 2) %>% arrange(p_russia, p_gambia, p_ghana) %>% select(gene, starts_with("p_"))

# Number of GWAS with p < .01
d_count <- d %>%
  filter(de_ghana > 2) %>%
  select(gene, p_russia, p_gambia, p_ghana, p_sobota)

rownames(d_count) <- d_count$gene
d_count <- d_count[, -1]
d_count <- d_count < .05
sort(rowSums(d_count))


