# install the packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!requireNamespace("DESeq2", quietly = TRUE))
  BiocManager::install("DESeq2")
if (!requireNamespace("getopt", quietly = TRUE))
  install.packages("getopt")
if (!requireNamespace("ggplot2", quietly = TRUE))
  install.packages("ggplot2")
if (!requireNamespace("tools", quietly = TRUE))
  install.packages("tools")
############################################################
# get the options
get_opts <- function() {
  suppressMessages(require(getopt, quietly = T))
  spec <- matrix(
    c("input", "i", 1, "character", "Input dataset of gene expression",
      "design", "d", 1, "character", "Design file which contian samples and groups",
      "suffix", "s", 2, "character", "suffix of samples ['utr5','cds','utr3', 'NA'] (default NA)",
      "threshold", "t", 2, "double", "Average threshold of each group (default 1)",
      "log2FC", "l", 2, "double", "filter the log2FoldChange (default 1)",
      "padj", "p", 2, "double", "filter the adjust p_value (default 0.05)",
      "help", "h", 0, "logical", "DESeq2 analysis."),
    byrow = TRUE, ncol = 5)

  args <- getopt(spec = spec)

  # check the options
  if (!is.null(args$help) ||
    is.null(args$input) ||
    is.null(args$design)) {
    # print useage
    cat(paste("\n", getopt(spec = spec, usage = T), "\n"))
    quit(status = 1)
  }

  # set the default options
  if (is.null(args$threshold)) { args$threshold <- 1 }
  if (is.null(args$log2FC)) { args$log2FC <- 2 }
  if (is.null(args$padj)) { args$padj <- 0.05 }
  if (is.null(args$suffix) || args$suffix == 'NA') { args$suffix <- NA }

  # print all options
  cat("===================================\n")
  cat(paste0(as.character(Sys.time()), "\n"))
  cat(paste0("Step1: Get options", "\n\n"))
  cat(paste("input:     ", args$input, "\n"))
  cat(paste("design:    ", args$design, "\n"))
  cat(paste("threshold: ", args$threshold, "\n"))
  cat(paste("log2FC:    ", args$log2FC, "\n"))
  cat(paste("padj:      ", args$padj, "\n\n"))

  return(args)
}

############################################################
# import gene expression table
import_gene <- function(args) {
  suppressMessages(require(tools, quietly = T))
  # read results
  cat(paste0(as.character(Sys.time()), ":\t", "Read expression table", "\n"))

  if (file_ext(args$input) == "csv") {
    rpf_data <- read.csv(file = args$input, header = T, row.names = 1, comment.char = '#')

  }else if (file_ext(args$input) == "txt" || file_ext(args$input) == "TXT") {
    rpf_data <- read.table(file = args$input, sep = "\t", header = T, row.names = 1, comment.char = '#')
  }

  return(rpf_data)
}

############################################################
# make design table
make_design <- function() {

  # read sample list and class, like this
  # name  group
  # sample1 ck
  # sample2 ck
  # sample3 treat
  # sample4 treat

  sample <- read.table(file = args$design, sep = "\t", header = T, row.names = NULL)
  keys <- as.matrix(summary.factor(sample[2]))
  sampName1 <- rownames(keys)[1]
  sampName2 <- rownames(keys)[2]
  sampNu1 <- keys[1]
  sampNu2 <- keys[2]


}

############################################################
# make design table
make_design <- function() {
  cat("===================================\n")
  cat(paste0(as.character(Sys.time()), "\n"))
  cat(paste0("Step2: Run DESeq2", "\n\n"))
  suppressMessages(require(DESeq2, quietly = T))
}

# make the phenotype list
cname <- as.data.frame(colnames(data))
pheno <- merge(cname, sample, by.x = "colnames(data)", by.y = "sample", all.x = T, sort = F)
colnames(pheno) <- c("sample", "class")

# split the data by phenotype
id_1 <- grep(sampName1, pheno[, 2])
samp1_data <- data[, as.numeric(id_1)]
id_2 <- grep(sampName2, pheno[, 2])
samp2_data <- data[, as.numeric(id_2)]
df <- cbind(samp1_data, samp2_data)

# filter the data by the threshold
cat(paste0(as.character(Sys.time()), ":\t", "Filter the data by threshold", "\n"))
FilteredMat <- subset(df, rowMeans(df[1:sampNu1]) >= args$threshold |
  rowMeans(df[(sampNu1 + 1):(sampNu1 + sampNu2)]) >= args$threshold)
mat <- round(FilteredMat)

# make the factors by phenotype
cat(paste0(as.character(Sys.time()), ":\t", "Make the phenotype table", "\n"))

type <- factor(c(rep(sampName1, sampNu1), rep(sampName2, sampNu2)))
condition <- factor(c(rep(sampName1, sampNu1), rep(sampName2, sampNu2)), levels = c(sampName1, sampName2))
colData <- data.frame(row.names = colnames(mat), condition)

# normalized the dataset
cat(paste0(as.character(Sys.time()), ":\t", "Run the diff analysis", "\n"))
dds <- DESeqDataSetFromMatrix(mat, colData, design = ~condition)
dds <- DESeq(dds)

# differential analysis and retrieve DESeq2 results from
res <- results(dds, contrast = c("condition", sampName2, sampName1))
out <- data.frame(res)

# make input data frame
names <- rownames(mat)
mat <- data.frame(names, mat)

# make normalized input data frame
newData <- counts(dds, normalized = TRUE)
names <- rownames(newData)
newdata_cnames <- colnames(newData)
newdata_cnames <- paste0(newdata_cnames, "_norm")
colnames(newData) <- newdata_cnames
newData <- data.frame(names, newData)

# make output data frame
names <- rownames(res)
out <- data.frame(names, out)

# merge DESeq2 results and raw reads and normalized reads
mer_list <- list(out, mat, newData)
mer <- Reduce(function(x, y) merge(x, y, all = T), mer_list)

# output all results
cat(paste0(as.character(Sys.time()), ":\t", "Ouput results table", "\n"))

sample_basename <- unlist(strsplit(args$input, "[.]"))[1]
FileName <- paste0(sample_basename, "_DESeq2_result_th", args$threshold, ".txt")
write.table(mer, file = FileName, sep = "\t", quote = F, col.names = T, row.names = F)

cat(paste0(as.character(Sys.time()), ":\t", "DESeq2 done!", "\n\n"))
############################################################
# import ggplot2 package
cat("===================================\n")
cat(paste0(as.character(Sys.time()), "\n"))
cat(paste0("Step3: Draw volcano", "\n\n"))
suppressMessages(require(ggplot2, quietly = T))

pdf(paste0(sample_basename, "_DESeq2_Plots.pdf"), width = 6, height = 5)

#filter the different factors
cat(paste0(as.character(Sys.time()), ":\t", "Filter data by log2FC and pvalue", "\n"))

for (i in seq_len(nrow(out))) {
  if (abs(out[i, 'log2FoldChange']) >= args$log2FC)
    out[i, 'select_change'] <- 'f'
  else
    out[i, 'select_change'] <- 'n'
  if (out[i, 'padj'] %in% NA | abs(out[i, 'padj']) >= args$padj)
    out[i, 'select_padj'] <- 'n'
  else
    out[i, 'select_padj'] <- 'p'
  out[i, 'select'] <- paste0(out[i, 'select_change'], out[i, 'select_padj'])
}

# draw the volcano plot
nn <- paste0("FC < ", args$log2FC, "p >= ", args$padj, ",")
np <- paste0("FC < ", args$log2FC, "p < ", args$padj, ",")
fn <- paste0("FC >= ", args$log2FC, "p >= ", args$padj, ",")
fp <- paste0("FC >= ", args$log2FC, "p < ", args$padj)
out$select <- factor(out$select, levels = c('nn', 'np', 'fn', 'fp'),
                     labels = c(nn, np, fn, fp))

# draw volcano plot by the pvalue
cat(paste0(as.character(Sys.time()), ":\t", "Draw volcano by pvalue", "\n"))

volcano_plot_pvalue <- ggplot(out, aes(log2FoldChange, -log(padj, 10))) +
  # set the color
  geom_point(aes(color = select), alpha = 0.6) +
  scale_color_manual(values = c('gray30', 'green4', 'red2', 'blue2')) +
  # set the theme of background
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent')) +
  theme(legend.title = element_blank(), legend.key = element_rect(fill = 'transparent'),
        legend.background = element_rect(fill = 'transparent')) +
  # set the theme of line
  geom_vline(xintercept = c(-args$log2FC, args$log2FC), color = 'gray', size = 0.5) +
  geom_hline(yintercept = -log(args$padj, 10), color = 'gray', size = 0.5) +
  # set the labels
  labs(x = 'log2 Fold Change', y = '-log10 padj', size = 10) +
  theme(title = element_text(size = 12, color = "black", hjust = 0.5, lineheight = 0.2))

ggsave(paste0(sample_basename, "_DESeq2_volcano1_th", args$threshold, "_FC", args$log2FC, "_padj", args$padj, "-pvalue.png"),
       volcano_plot_pvalue, width = 6, height = 5)

plot(volcano_plot_pvalue)

#########################################################
# main function is here
cat("===================================\n")
cat(paste0(as.character(Sys.time()), "\n"))
cat(paste0("Step1: Check input files", "\n"))
args <- get_opts()

cat(paste0("Step2: import the gene expression table", "\n"))
rpf_data <- import_gene(args)

cat(paste0("Step3: make design table", "\n"))
rpf_data <- import_gene(args)

cat(paste0("Step4: run DESeq2", "\n"))
rpf_data <- import_gene(args)

cat(paste0("Step5: output results", "\n"))
rpf_data <- import_gene(args)

cat(paste0(as.character(Sys.time()), ":\n", "All done!", "\n"))
cat("===================================\n")




