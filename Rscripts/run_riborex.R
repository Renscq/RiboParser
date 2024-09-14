# install the packages
# options(repos=structure( c(CRAN = "https://cloud.r-project.org/")))
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}
if (!requireNamespace("riborex", quietly = TRUE)) {
    library(devtools)
    options(unzip = "internal")
    devtools::install_github("smithlabcode/riborex")
}
if (!requireNamespace("EnhancedVolcano", quietly = TRUE)) {
    BiocManager::install("EnhancedVolcano")
}
if (!requireNamespace("getopt", quietly = TRUE)) {
    install.packages("getopt")
}
if (!requireNamespace("dplyr", quietly = TRUE)) {
    install.packages("dplyr")
}
################################################################
# get the options
get_opts <- function() {
    suppressMessages(require(getopt, quietly = T))
    spec <- matrix(
        c(
            "ribo", "r", 1, "character", "Input gene expression table of Ribo-seq",
            "mrna", "m", 1, "character", "Input gene expression table of RNA-seq",
            "engine", "e", 2, "character", "Specify the engine for riborex, default edgeR. [DESeq2, edgeR, Voom]",
            "output", "o", 1, "character", "Output riborex results",
            "design", "d", 1, "character", "Design file which contian samples and groups [name groups seqtype file]",
            "compare", "c", 1, "character", "Group name of control_vs_treat, multiple groups combine with '+', eg. ck1+wt2_vs_treat1+treat3",
            "threshold", "t", 2, "double", "Average threshold of each group (default 1)",
            "log2FC", "l", 2, "double", "filter the log2FoldChange (default 1)",
            "padj", "p", 2, "double", "filter the adjust p_value (default 0.05)",
            "help", "h", 0, "logical", "Pipeline of edgeR analysis."
        ),
        byrow = TRUE, ncol = 5
    )

    args <- getopt(spec = spec)

    # check the options
    if (!is.null(args$help) ||
        is.null(args$ribo) ||
        is.null(args$ribo) ||
        is.null(args$design)) {
        # print useage
        cat(paste("\n", getopt(spec = spec, usage = T), "\n"))
        quit(status = 1)
    }

    # set the default options
    if (is.null(args$engine)) {
        args$engine <- "edgeR"
    }
    if (is.null(args$threshold)) {
        args$threshold <- 1
    }
    if (is.null(args$log2FC)) {
        args$log2FC <- 1
    }
    if (is.null(args$padj)) {
        args$padj <- 0.05
    }

    # print all options
    cat("===================================\n")
    cat(paste0(as.character(Sys.time()), "\n"))
    cat(paste0("Step1: Get options", "\n\n"))
    cat(paste("ribo:     ", args$ribo, "\n"))
    cat(paste("mrna:     ", args$mrna, "\n"))
    cat(paste("design:    ", args$design, "\n"))
    cat(paste("compare:    ", args$compare, "\n"))
    cat(paste("threshold: ", args$threshold, "\n"))
    cat(paste("log2FC:    ", args$log2FC, "\n"))
    cat(paste("padj:      ", args$padj, "\n\n"))

    return(args)
}

################################################################
# make design table
make_design <- function(args, seqtype) {
    # import design
    sample <- read.table(file = args$design, sep = "\t", header = T, row.names = NULL)
    sample[, 3] <- toupper(sample[, 3])
    sample <- sample[sample[, 3] %in% seqtype,]

    # make the comparasion
    comp <- unlist(strsplit(args$compare, split = "_vs_"))

    ck_comp <- unlist(strsplit(comp[1], split = "\\+"))
    ck <- sample[sample[, 2] %in% ck_comp,]

    treat_comp <- unlist(strsplit(comp[2], split = "\\+"))
    treat <- sample[sample[, 2] %in% treat_comp,]

    ck_list <- ck[, 1]
    treat_list <- treat[, 1]

    return(list(ck_list, treat_list))
}

################################################################
# clean the colnames
clean_colnames <- function(data) {
    cname <- colnames(data)

    for (num in c(1:length(cname))) {
        cname[num] <- basename(cname[num])
    }
    cname <- gsub("Aligned.sortedByCoord.out.bam", "", cname)
    cname <- gsub("_trim_cds_rpf", "", cname)
    cname <- gsub("_cds_rpf", "", cname)
    cname <- gsub("_utr5_rpf", "", cname)
    cname <- gsub("_utr3_rpf", "", cname)
    cname <- gsub(".genes.results", "", cname)
    cname <- gsub(".isoforms.results", "", cname)
    colnames(data) <- cname

    return(data)
}

################################################################
# retrieve the gene expression from original dataframe
import_gene <- function(args) {
    # import gene expression
    rna_df <- read.table(args$mrna, header = T, row.names = 1, comment.char = '#', check.names = FALSE)
    ribo_df <- read.table(args$ribo, header = T, row.names = 1, comment.char = '#', check.names = FALSE)

    # clean the suffix of colnames
    rna_df <- clean_colnames(rna_df)
    ribo_df <- clean_colnames(ribo_df)

    # make comparasion of rna and ribo data
    rna_sample <- make_design(args, "RNA")
    ck_rna <- unlist(rna_sample[1])
    treat_rna <- unlist(rna_sample[2])

    ribo_sample <- make_design(args, "RIBO")
    ck_ribo <- unlist(ribo_sample[1])
    treat_ribo <- unlist(ribo_sample[2])

    rna_class <- rep(c("control", "treated"), c(length(ck_rna), length(treat_rna)))
    ribo_class <- rep(c("control", "treated"), c(length(ck_ribo), length(treat_ribo)))
    
    # retrieve the gene expression
    rna <- rna_df[c(ck_rna, treat_rna)]
    ribo <- ribo_df[c(ck_ribo, treat_ribo)]
    
    # filter out low expression genes
    rna <- subset(rna, rowMeans(rna) >= args$threshold)
    ribo <- subset(ribo, rowMeans(ribo) >= args$threshold)

    # retrieve the gene intersection
    rna_gene <- rownames(rna)
    ribo_gene <- rownames(ribo)
    gene <- intersect(rna_gene, ribo_gene)
    
    rna <- round(rna[gene, ])
    ribo <- round(ribo[gene, ])
    
    return(list(rna, ribo, rna_class, ribo_class))
}

################################################################
# filter the significant differential genes
flt_sig <- function(args, res){
    if (args$engine == "DESeq2") {
        up_gene <- subset(res, padj < args$padj & log2FoldChange >= args$log2FC)
        down_gene <- subset(res, padj < args$padj & log2FoldChange <= -args$log2FC)
    } else if (args$engine == "edgeR") {
        up_gene <- subset(res, FDR < args$padj & logFC >= args$log2FC)
        down_gene <- subset(res, FDR < args$padj & logFC <= -args$log2FC)
    }
    
    up_file <- paste0(args$output, "-riborex-", args$engine, "-up.txt")
    down_file <- paste0(args$output, "-riborex-", args$engine, "-down.txt")
    write.table(up_gene, file = up_file, sep = "\t", row.names = F, col.names = TRUE, quote = F)
    write.table(down_gene, file = down_file, sep = "\t", row.names = F, col.names = TRUE, quote = F)
}

################################################################
# run the riborex pipeline

run_riborex <- function(args, rna, ribo, rna_class, ribo_class) {
    # run the riborex
    # DESeq2, edgeR, edgeRD, Voom
    res.riborex <- riborex(rna, ribo, rna_class, ribo_class, args$engine)
    res.riborex.tab <- data.frame(res.riborex)

    # merge output data
    colnames(rna) <- paste0(colnames(rna), "_rna")
    colnames(ribo) <- paste0(colnames(ribo), "_ribo")

    merge_res <- merge(res.riborex.tab, rna, by = "row.names", all.x = TRUE)
    row.names(merge_res) <- merge_res[, 1]
    merge_res <- merge_res[, -1]

    merge_res <- merge(merge_res, ribo, by = "row.names", all.x = TRUE)
    
    # filter sig diff genes
    flt_sig(args, merge_res)
    
    # output the results
    colnames(merge_res)[1] <- "name"

    FileName <- paste(args$output, "-riborex-", args$engine, "-result.txt", sep = "")
    write.table(merge_res, file = FileName, sep = "\t", row.names = F, col.names = TRUE, quote = F)

    return(merge_res)
}

############################################
draw_volcano <- function(args, merge_res) {
    if (args$engine == "DESeq2") {
        logfc <- "log2FoldChange"
        pvalue <- "padj"
    } else if (args$engine == "edgeR") {
        logfc <- "logFC"
        pvalue <- "FDR"
    }

    merge_res <- merge_res[order(merge_res[, pvalue]), ]
    labels <- merge_res[1:10, "name"]

    vol <- EnhancedVolcano(merge_res,
                           lab = merge_res$name,
                           selectLab = labels,
                           x = logfc,
                           y = pvalue,
                           title = args$compare,
                           subtitle = "",
                           pCutoff = 0.05,
                           FCcutoff = 1,
                           pointSize = 3.0,
                           labSize = 2.0,
                           labCol = "black",
                           legendLabSize = 8,
                           legendPosition = "right",
                           # col = c("grey30", "royalblue", "#0072b5", "#bc3c28"),
                           legendLabels = c("NS", "logFC", "FDR", "FDR & logFC")) +
            theme(plot.title = element_text(vjust = 0.5, hjust = 0.5))
    ggsave(paste0(args$output, "-", args$engine, "-te-volcano.pdf"), vol, width = 8, height = 7)
    ggsave(paste0(args$output, "-", args$engine, "-te-volcano.png"), vol, width = 8, height = 7, dpi = 600)
}

#############################################
# main function is here
cat("===================================\n")
cat(paste0(as.character(Sys.time()), "\n"))
cat(paste0("Step1: Check input files", "\n"))
args <- get_opts()

# import the packages
suppressMessages(require(riborex, quietly = T))
suppressMessages(require(dplyr, quietly = T))
suppressMessages(require(EnhancedVolcano, quietly = T))

cat(paste0("Step2: import the gene expression table", "\n"))
expr_data <- import_gene(args)
rna <- data.frame(expr_data[1])
ribo <- data.frame(expr_data[2])
rna_class <- unlist(expr_data[3])
ribo_class <- unlist(expr_data[4])

cat(paste0("Step3: run riborex analysis", "\n"))
merge_res <- run_riborex(args, rna, ribo, rna_class, ribo_class)

cat(paste0("Step4: draw volcano", "\n"))
draw_volcano(args, merge_res)

cat(paste0(as.character(Sys.time()), ":\n", "All done!", "\n"))
cat("===================================\n")