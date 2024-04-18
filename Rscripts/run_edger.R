#!/usr/bin/Rscript

# install the packages

options(repos = structure(c(CRAN = "https://cloud.r-project.org/")))
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}
if (!requireNamespace("edgeR", quietly = TRUE)) {
    BiocManager::install("edgeR")
}
if (!requireNamespace("EnhancedVolcano", quietly = TRUE)) {
    BiocManager::install("EnhancedVolcano")
}
if (!requireNamespace("getopt", quietly = TRUE)) {
    install.packages("getopt")
}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
    install.packages("ggplot2")
}
if (!requireNamespace("tools", quietly = TRUE)) {
    install.packages("tools")
}
if (!requireNamespace("dplyr", quietly = TRUE)) {
    install.packages("dplyr")
}

############################################################
# get the options
get_opts <- function() {
    suppressMessages(require(getopt, quietly = T))
    spec <- matrix(
        c(
            "input", "i", 1, "character", "Input dataset of gene expression",
            "output", "o", 1, "character", "Output edgeR results",
            "design", "d", 1, "character", "Design file which contian samples and groups",
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
        is.null(args$input) ||
        is.null(args$design)) {
        # print useage
        cat(paste("\n", getopt(spec = spec, usage = T), "\n"))
        quit(status = 1)
    }

    # set the default options
    if (is.null(args$threshold)) {
        args$threshold <- 1
    }
    if (is.null(args$log2FC)) {
        args$log2FC <- 1
    }
    if (is.null(args$padj)) {
        args$padj <- 0.05
    }
    if (is.null(args$suffix) || args$suffix == "NA") {
        args$suffix <- NA
    }

    # print all options
    cat("===================================\n")
    cat(paste0(as.character(Sys.time()), "\n"))
    cat(paste0("Step1: Get options", "\n\n"))
    cat(paste("input:     ", args$input, "\n"))
    cat(paste("design:    ", args$design, "\n"))
    cat(paste("compare:    ", args$compare, "\n"))
    cat(paste("threshold: ", args$threshold, "\n"))
    cat(paste("log2FC:    ", args$log2FC, "\n"))
    cat(paste("padj:      ", args$padj, "\n\n"))

    return(args)
}

############################################################
# make design table
make_design <- function(args) {
    # import design
    sample <- read.table(file = args$design, sep = "\t", header = T, row.names = NULL)
    # make the comparasion
    comp <- unlist(strsplit(args$compare, split = "_vs_"))

    ck_comp <- unlist(strsplit(comp[1], split = "\\+"))
    ck <- sample[sample[, 2] %in% ck_comp,]

    treat_comp <- unlist(strsplit(comp[2], split = "\\+"))
    treat <- sample[sample[, 2] %in% treat_comp,]

    ck_list <- ck[, 1]
    treat_list <- treat[, 1]

    # output the design file
    FileName <- paste(args$output, "-design.txt", sep = "")
    ck[, 2] <- "ck"
    treat[, 2] <- "treat"

    now_design <- rbind(ck, treat)
    now_design[, 1] <- paste0(now_design[, 1], "_norm")

    write.table(now_design, file = FileName, sep = "\t", row.names = F, col.names = TRUE, quote = F)

    return(list(ck_list, treat_list))
}

############################################################
# import gene expression table
import_gene <- function(args, ck_list, treat_list) {

    # read results
    if (file_ext(args$input) == "csv") {
        data <- read.csv(file = args$input, header = T, row.names = 1, comment.char = "#", check.names=FALSE)
    } else if (file_ext(args$input) == "txt" || file_ext(args$input) == "TXT") {
        data <- read.table(file = args$input, sep = "\t", header = T, row.names = 1, comment.char = "#", check.names=FALSE)
    } else {
        data <- read.table(file = args$input, sep = "\t", header = T, row.names = 1, comment.char = "#", check.names=FALSE)
    }

    head(data)

    # trim the suffix message
    cname <- colnames(data)
    cname <- gsub("_trim_cds_rpf", "", cname)
    cname <- gsub("_cds_rpf", "", cname)
    cname <- gsub("_utr5_rpf", "", cname)
    cname <- gsub("_utr3_rpf", "", cname)
    cname <- gsub(".genes.results", "", cname)
    cname <- gsub(".isoforms.results", "", cname)
    cname <- gsub("Aligned.sortedByCoord.out.bam", "", cname)
    for (ids in c(1:length(cname))) {
        cname[ids] <- basename(cname[ids])
    }

    colnames(data) <- cname

    # retrieve the groups for compararsion

    sample_list <- c(ck_list, treat_list)

    expr_data <- subset(data, select = sample_list)

    return(expr_data)
}

#########################################################
# run edgeR analysis

run_edgeR <- function(args, expr_data, ck_list, treat_list) {
    # filter the gene expression
    ck_num <- length(ck_list)
    treat_num <- length(treat_list)

    FilteredMat <- subset(expr_data, rowMeans(expr_data[1:ck_num]) >= args$threshold |
        rowMeans(expr_data[(ck_num + 1):(ck_num + treat_num)]) >= args$threshold)

    # run the edgeR
    group <- factor(rep(c(1, 2), times = c(ck_num, treat_num)))
    design <- model.matrix(~group)
    y <- DGEList(counts = FilteredMat, group = group)
    y <- estimateGLMCommonDisp(y, design)
    y <- estimateGLMTrendedDisp(y, design)
    y <- estimateGLMTagwiseDisp(y, design)
    fit <- glmFit(y, design)
    lrt <- glmLRT(fit, coef = 2)

    # output sig gene list
    gene <- data.frame(topTags(lrt, nrow(FilteredMat)))
    sig <- subset(gene, FDR < args$padj)
    up_gene <- subset(sig, logFC >= args$log2FC )
    down_gene <- subset(sig, logFC <= -args$log2FC)

    FileName <- paste(args$output, "-edger-up.txt", sep = "")
    write.table(up_gene, file = FileName, sep = "\t", row.names = F, col.names = F, quote = F)
    FileName <- paste(args$output, "-edger-down.txt", sep = "")
    write.table(down_gene, file = FileName, sep = "\t", row.names = F, col.names = F, quote = F)

    # merge the results
    fitted_value <- data.frame(lrt$fitted.values)
    colnames(fitted_value) <- paste0(colnames(FilteredMat), "_norm")

    edger_res <- merge(gene, fitted_value, by = "row.names", all.x = TRUE)
    row.names(edger_res) <- edger_res[, 1]
    edger_res <- edger_res[, -1]

    merge_res <- merge(edger_res, FilteredMat, by = "row.names", all.x = TRUE)
    colnames(merge_res)[1] <- "name"

    # output the results
    FileName <- paste(args$output, "-edger-result.txt", sep = "")
    write.table(merge_res, file = FileName, sep = "\t", row.names = F, col.names = TRUE, quote = F)

    return(merge_res)
}

#########################################################
# run edgeR analysis

run_bcv_edgeR <- function(args, expr_data, ck_list, treat_list) {
    ck_num <- length(ck_list)
    treat_num <- length(treat_list)

    FilteredMat <- subset(expr_data, expr_data[, ck_num] >= args$threshold |
        expr_data[, (ck_num + treat_num)] >= args$threshold)
    group <- factor(rep(c(1, 2), times = c(ck_num, treat_num)))

    # run the edger with bcv=0.2
    exprSet <- DGEList(counts = FilteredMat, group = group)
    keep <- rowSums(cpm(exprSet) > 1) >= 1
    exprSet <- exprSet[keep, , keep.lib.sizes = FALSE]
    exprSet <- calcNormFactors(exprSet)

    bcv <- 0.2
    et <- exactTest(exprSet, dispersion = bcv^2)

    # output sig gene list
    gene <- data.frame(topTags(et, nrow(FilteredMat)))
    sig <- subset(gene, logFC >= args$log2FC | logFC <= -args$log2FC)
    sig <- subset(sig, FDR < args$padj)

    up_gene <- rownames(sig[sig$logFC > 0,])
    down_gene <- rownames(sig[sig$logFC < 0,])

    FileName <- paste(args$output, "-edger-up.txt", sep = "")
    write.table(up_gene, file = FileName, sep = "\t", row.names = F, col.names = F, quote = F)
    FileName <- paste(args$output, "-edger-down.txt", sep = "")
    write.table(down_gene, file = FileName, sep = "\t", row.names = F, col.names = F, quote = F)

    # merge the results
    edger_res <- data.frame(topTags(et, nrow(FilteredMat)))

    merge_res <- merge(edger_res, FilteredMat, by = "row.names", all.x = TRUE)
    colnames(merge_res)[1] <- "name"

    # output the results
    FileName <- paste0(args$output, "-edger-result.txt")
    write.table(merge_res, file = FileName, sep = "\t", row.names = F, col.names = TRUE, quote = F)

    return(merge_res)
}

#########################################################
# draw volcano

draw_volcano <- function(args, edger_out) {
    edger_sort <- edger_out[order(edger_out$FDR),]
    # labels <- edger_sort$name[1:10]
    labels <- rownames(edger_sort[1:10,])

    vol <- EnhancedVolcano(edger_out,
                           lab = edger_out$name,
                           selectLab = labels,
                           x = "logFC",
                           y = "FDR",
                           title = args$compare,
                           subtitle = "",
                           pCutoff = 0.05,
                           FCcutoff = 1,
                           pointSize = 3.0,
                           labSize = 2.0,
                           labCol = "black",
                           legendLabSize = 8,
                           legendPosition = "right", ,
                           # col = c("grey30", "royalblue", "#0072b5", "#bc3c28"),
                           legendLabels = c("NS", "logFC", "FDR", "FDR & logFC")) +
      theme(plot.title = element_text(hjust = 0.5, vjust = 0.5))

    ggsave(paste0(args$output, "-edger-volcano.pdf"), vol, width = 8, height = 7)
    ggsave(paste0(args$output, "-edger-volcano.png"), vol, width = 8, height = 7, dpi = 600)
}

#########################################################
# main function is here
cat("===================================\n")
cat(paste0(as.character(Sys.time()), "\n"))
cat(paste0("Step1: Check input files", "\n"))
args <- get_opts()

# import the packages
suppressMessages(require(tools, quietly = T))
suppressMessages(require(edgeR, quietly = T))
suppressMessages(require(ggplot2, quietly = T))
suppressMessages(require(dplyr, quietly = T))
suppressMessages(require(EnhancedVolcano, quietly = T))

cat(paste0("Step2: import design table", "\n"))
now_design <- make_design(args)
ck_list <- unlist(now_design[1])
treat_list <- unlist(now_design[2])

cat(paste0("Step3: import the gene expression table", "\n"))
expr_data <- import_gene(args, ck_list, treat_list)

cat(paste0("Step4: run edgeR analysis", "\n"))
if (length(c(ck_list, treat_list)) > 2) {
    edger_out <- run_edgeR(args, expr_data, ck_list, treat_list)
} else {
    edger_out <- run_bcv_edgeR(args, expr_data, ck_list, treat_list)
}

cat(paste0("Step5: draw volcano", "\n"))
draw_volcano(args, edger_out)

cat(paste0(as.character(Sys.time()), ":\n", "All done!", "\n"))
cat("===================================\n")
