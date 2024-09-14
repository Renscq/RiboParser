################################################################
# Draw the volcano of different expression genes in DESeq2 results
# Rensc
################################################################
# install the packages
if (!requireNamespace("getopt", quietly = TRUE))
  install.packages("getopt")
if (!requireNamespace("ggplot2", quietly = TRUE))
  install.packages("ggplot2")
if (!requireNamespace("ggplotify", quietly = TRUE))
  install.packages("ggplotify")
############################################################
# get the options
suppressMessages(require(getopt, quietly = T))
spec <- matrix(
  c("input", "i", 1, "character", "Input the DESeq2 results file.",
    "design", "d", 1, "character", "Input the DESeq2 design file.",
    "volcano", "v", 1, "character", "output the volcano file name)",
    "log2FC", "l", 1, "double", "Specify the threshold of log2FC (default 1).",
    "padj", "p", 2, "double", "Specify the threshold of padj (default 0.05).",
    "help", "h", 0, "logical", "Volcano plot of DESeq2 results."),
  byrow = TRUE, ncol = 5)

args <- getopt(spec = spec)

# check the options
if (!is.null(args$help) || is.null(args$input)) {
  # print usage
  cat(paste("\n", getopt(spec = spec, usage = T), "\n"))
  quit(status = 1)
}

# set the default options
if (is.null(args$log2FC)) { args$log2FC <- 1 }
if (is.null(args$padj)) { args$padj <- 0.05 }

# print all options
cat("=================================== \n")
cat(paste0(as.character(Sys.time()), "\n"))
cat(paste0("Step1: Get options", "\n"))
cat(paste("input:", args$input, "\n"))
cat(paste("design:", args$design, "\n"))
cat(paste("log2FC: ", args$log2FC, "\n"))
cat(paste("padj: ", args$padj, "\n"))
cat(paste("output:", args$output, "\n"))

############################################################
# import ggplot2 package
cat("=================================== \n")
cat(paste0(as.character(Sys.time()), "\n"))
cat(paste0("Step2: import ggplot2 packages", "\n"))
options(warn = -1)
suppressMessages(require(ggplot2, quietly = T))
suppressMessages(require(ggplotify, quietly = T))
suppressMessages(require(pheatmap, quietly = T))

############################################################
# check the files
if (is.null(args$input)) {
  cat(paste("\n", "Miss the DESeq2 results", "\n"))
  quit(status = 2)
}

############################################################
# impor the DESeq2 results
import_deseq <- function(file_name, log2fc, padj) {
  # now group
  prefix_name <- unlist(strsplit(basename(file_name), '[.csv]'))[1]
  cat(paste("Now group is: ", file_name, "\n", sep = ""))

  # read data
  data <- read.csv(file_name)
  deseq_res <- data.frame(row.names = data[, 1], data[, -1])

  # filter significant genes
  deseq_res$anno <- ifelse(deseq_res$padj >= padj |
                             deseq_res$padj %in% NA |
                             abs(deseq_res$log2FoldChange) < log2fc, 'NS',
                           ifelse(deseq_res$padj < padj & deseq_res$log2FoldChange >= log2fc, 'Up', 'Down'))
  deseq_res$anno <- factor(deseq_res$anno, levels = c("Up", "Down", "NS"))

  min_padj <- min(deseq_res$padj[deseq_res$padj > 0], na.rm = T)
  deseq_res$padj[deseq_res$padj == 0] <- min_padj / 100

  return(list(prefix_name = prefix_name,
              deseq_res = deseq_res))
}

############################################################
# draw the volcano
draw_volcano <- function(out, prefix_name, volcano_plot) {
  up_num <- length(which(deseq_res$anno == "Up"))
  down_num <- length(which(deseq_res$anno == "Down"))
  xlim_max <- max(deseq_res$log2FoldChange)

  volcano_plot_pvalue <- ggplot(out, aes(log2FoldChange, -log(padj, 10))) +
    geom_point(size = 1,
               aes(color = anno),
               alpha = 1) +
    scale_color_manual(values = c('Up' = '#bc3c28',
                                  'Down' = '#0072b5',
                                  'NS' = '#bebebe'),
                       labels = c(paste0("Up ", '(', up_num, ')'),
                                  paste0("Down ", '(', down_num, ')'),
                                  paste("NS"))) +
    theme(axis.title.x = element_text(size = 14, color = "black"),
          axis.title.y = element_text(size = 14, colour = 'black'),
          axis.text.x = element_text(size = 14, colour = 'black'),
          axis.text.y = element_text(size = 14, colour = 'black'),
          panel.grid = element_blank(),
          panel.background = element_rect(color = 'black', fill = 'transparent'),
          legend.title = element_blank(),
          legend.key = element_rect(fill = 'transparent'),
          legend.text = element_text(size = 12),
          legend.background = element_rect(fill = 'transparent'),
          plot.title = element_text(hjust = 0.5)) +
    # set the x/y range
    coord_cartesian(xlim = c(-xlim_max * 1.1,
                             xlim_max * 1.1)) +
    # plot the line to mark log2FC
    geom_vline(xintercept = c(-args$log2FC, args$log2FC),
               color = 'black',
               size = 0.4,
               linetype = 4) +
    # plot the line to mark padj
    geom_hline(yintercept = -log(args$padj, 10),
               color = 'black',
               size = 0.4,
               linetype = 4) +
    labs(title = prefix_name,
         x = 'log2 Fold Change',
         y = '-log10 padj',
         size = 16)

  # save the figure
  ggsave(volcano_plot, volcano_plot_pvalue, width = 6, height = 5)
  # ggsave(volcano_plot, volcano_plot_pvalue, width = 6, height = 5, dpi = 600)
}

############################################################
# plot the volcano
# main function is here
cat("=================================== \n")
cat(paste0(as.character(Sys.time()), "\n"))
cat(paste0("Step3: plot the volcano", "\n\n"))

for (file_name in files) {
  # impor the DESeq2 results
  res <- import_deseq(file_name, args$log2FC, args$padj)
  prefix_name <- res['prefix_name']
  deseq_res <- res['deseq_res']

  ############################################################
  # draw volcano plot by the padj
  draw_volcano(deseq_res, prefix_name, args$output)
}
############################################################
# End here
cat("=================================== \n")
cat(paste0(as.character(Sys.time()), "\n"))
cat(paste0("All done.", "\n\n"))
