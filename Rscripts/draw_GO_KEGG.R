################################################################
# Draw the go and kegg enrichment
# Rensc
# 2021-8-9
################################################################
# install the packages
# if (!requireNamespace("getopt", quietly = TRUE))
#   install.packages("getopt")
# if (!requireNamespace("ggplot2", quietly = TRUE))
#   install.packages("ggplot2")
# if (!requireNamespace("clusterProfiler", quietly = TRUE))
#   BiocManager::install("clusterProfiler")
# if (!requireNamespace("enrichplot", quietly = TRUE))
#   BiocManager::install("enrichplot")
# if (!requireNamespace("AnnotationHub", quietly = TRUE))
#   BiocManager::install("AnnotationHub")
# if (!requireNamespace("pathview", quietly = TRUE))
#   BiocManager::install("pathview")
# if (!requireNamespace("topGO", quietly = TRUE))
#   BiocManager::install("topGO")
# if (!requireNamespace("Rgraphviz", quietly = TRUE))
#   BiocManager::install("Rgraphviz")
# if (!requireNamespace("mygene", quietly = TRUE))
#   BiocManager::install("mygene")
# if (!requireNamespace("stringr", quietly = TRUE))
#   BiocManager::install("stringr")

################################################################
# get the options
suppressMessages(require(getopt, quietly = T))
spec <- matrix(
  c("input", "i", 1, "character", "Input the gene list.",
    "output", "o", 1, "character", "output file (path) prefix name",
    "log2FC", "l", 1, "double", "Specify the log2FC threshold of sig. diff. genes (default 1).",
    "padj", "p", 2, "double", "Specify the padj threshold of sig. diff. genes (default 0.05).",
    "enrich", "e", 2, "double", "Specify the padj threshold of enrichment items (default 0.1).",
    "number", "n", 2, "double", "Specify the number of terms to show (default 30).",
    "species", "s", 1, "character", "Specify the kegg abbreviations of species.",
    "taxon", "t", 1, "double", "Specify the taxon id.",
    "fontsize", "f", 2, "double", "Specify the fontsize (default 12).",
    "color", "c", 2, "character", "Specify the color, (default r_b); eg: r_g / r_b / b_p.",
    "help", "h", 0, "logical", "Draw the figure of GO, KEGG and GSEA enrichment."),
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
if (is.null(args$enrich)) { args$enrich <- 0.1 }
if (is.null(args$number)) { args$number <- 30 }
if (is.null(args$fontsize)) { args$fontsize <- 12 }
if (is.null(args$color)) { args$color <- 'r_b' }

default_color <- data.frame(r = "red", o = "orange",
                            y = "yellow", g = "green",
                            c = "cyan", b = "blue", p = "purple")

now_color <- default_color[strsplit(args$color, '_')[[1]]]

# print all options
cat("=================================== \n")
cat(paste0(as.character(Sys.time()), "\n"))
cat(paste0("Step1: Get options", "\n\n"))
cat(paste("input:", args$input, "\n"))
cat(paste("output:", args$output, "\n"))
cat(paste("log2FC: ", args$log2FC, "\n"))
cat(paste("padj: ", args$padj, "\n"))
cat(paste("enrich: ", args$enrich, "\n"))
cat(paste("number: ", args$number, "\n"))
cat(paste("fontsize: ", args$fontsize, "\n"))
cat(paste("color: ", args$color, "\n"))
cat(paste("species: ", args$species, "\n"))
cat(paste("taxon: ", args$taxon, "\n\n"))

############################################################
# import the packages
cat("=================================== \n")
cat(paste0(as.character(Sys.time()), "\n"))
cat(paste0("Step2: import the packages", "\n\n"))
options(warn = -1)
suppressMessages(require(ggplot2, quietly = T))
suppressMessages(require(AnnotationHub, quietly = T))
suppressMessages(require(pathview, quietly = T))
suppressMessages(require(clusterProfiler, quietly = T))
suppressMessages(require(enrichplot, quietly = T))
suppressMessages(require(topGO, quietly = T))
suppressMessages(require(Rgraphviz, quietly = T))
suppressMessages(require(mygene, quietly = T))
suppressMessages(require(stringr, quietly = T))
suppressMessages(require(tools, quietly = T))
suppressMessages(require(dylyr, quietly = T))
############################################################
# check the files
check_files <- function(input_file, output_path) {
  # check the input file
  if (!file.exists(input_file)) {
    cat(paste("\n", "Miss DESeq2/edgeR results file", "\n"))
    quit(status = 2)
  }

  # set the output path
  if (file.exists(output_path)) {
    setwd(output_path)
  }else {
    dir.create(output_path)
    setwd(output_path)
  }

}

############################################################
# output the orgdb file
output_orgdb <- function(my_OrgDb, out_path) {
  org_colname <- columns(my_OrgDb)
  total_gene_num <- length(keys(my_OrgDb))
  orgdb_file <- paste0(out_path, '/', 'my_OrgDb.Rdata')
  if (!file.exists(orgdb_file)) {
    cat(paste0("Save the OrgDb.", "\n\n"))
    saveDb(my_OrgDb, file = orgdb_file)
  }
}

# output the gene list
output_gene <- function(my_OrgDb, out_path) {
  gene_list_file <- paste0(out_path, '/', 'gene_list.txt')
  if (!file.exists(gene_list_file)) {
    cat(paste0("Reterieve gene message.", "\n"))
    gene_list <- clusterProfiler::bitr(OrgDb = my_OrgDb,
                                       geneID = keys(my_OrgDb, 'ENTREZID'),
                                       fromType = 'ENTREZID',
                                       toType = c('SYMBOL', 'GENENAME'))
    cat(paste0("Save the gene list.", "\n"))
    write.table(x = gene_list,
                file = gene_list_file,
                sep = '\t',
                quote = F,
                row.names = F)
  }
}

# output the go list
output_go <- function(my_OrgDb, out_path) {
  go_list <- paste0(out_path, '/', 'GO_list.txt')
  if (!file.exists(go_list)) {
    go_list <- bitr(OrgDb = my_OrgDb,
                    geneID = keys(my_OrgDb, 'ENTREZID'),
                    fromType = 'ENTREZID',
                    toType = c("REFSEQ", "GO", "ONTOLOGY", 'SYMBOL', 'GENENAME'))
    write.table(x = go_list,
                file = paste0(out_path, '/', 'GO_list.txt'),
                sep = '\t',
                quote = F,
                row.names = F)
  }
}

# make the gene database
make_orgdb <- function(taxon, out_path) {

  # hub <- AnnotationHub(localHub=TRUE)
  hub <- AnnotationHub()

  cat(paste0("query taxon id: ", taxon, "\n"))
  my_hub <- AnnotationHub::query(hub, taxon)

  org_num <- grep('org.', my_hub$title)
  org_ID <- my_hub$ah_id[org_num]
  org_ID <- tail(org_ID, 1)
  cat(paste0('AnnotationHub ID: ', org_ID, '\n'))
  my_OrgDb <- my_hub[[org_ID]]

  output_gene(my_OrgDb, out_path)
  output_orgdb(my_OrgDb, out_path)
  # output_go(my_OrgDb, out_path)

  return(my_OrgDb)
}


############################################################
# import the gene list
import_diff <- function(now_filename, now_padj, now_log2fc, now_taxonomy) {

  # read data
  if (file_ext(now_filename) == "csv") {
    data <- read.csv(now_filename)

  }else if (file_ext(now_filename) == "txt" || file_ext(now_filename) == "TXT") {
    data <- read.table(file = now_filename, sep = "\t", header = T, row.names = NULL, comment.char = '#')
  }

  diff_res <- data.frame(row.names = data[, 1], data[, -1])

  if ("FDR" %in% colnames(diff_res)) {
    cat("Results generate from edgeR")
    padj <- "FDR"
    logFC <- "logFC"
  } else {
    cat("Results generate from DESeq2")
    padj <- "padj"
    logFC <- "log2FoldChange"
  }

  # filter significant genes
  diff_res$anno <- ifelse(diff_res[, padj] >= now_padj |
                            diff_res[, padj] %in% NA |
                            abs(diff_res[, logFC]) < now_log2fc, 'NS',
                          ifelse(diff_res[, padj] < now_padj & diff_res[, logFC] >= now_log2fc, 'Up', 'Down'))
  diff_res$anno <- factor(diff_res$anno,
                          levels = c("Up", "Down", "NS"))

  up_num <- length(which(diff_res$anno == "Up"))
  down_num <- length(which(diff_res$anno == "Down"))

  # repalce the padj=0 to min padj
  min_padj <- min(diff_res[diff_res[, padj] > 0, padj], na.rm = T)
  diff_res[diff_res[, padj] > 0, padj] <- min_padj / 100

  # make the gene dataframe
  rnames <- rownames(diff_res)
  diff_anno <- data.frame(gene = rnames,
                          logfc = diff_res[, logFC],
                          anno = diff_res$anno)

  # get the entrez gene id
  cat(paste0("Reterieve entrez gene id.", "\n"))
  # scope = c('refseq', 'ensembl.transcript'),
  converter <- queryMany(rnames, fields = c('entrezgene'), species = now_taxonomy)
  converter_df <- data.frame(converter)[c('query', 'entrezgene')]
  diff_anno <- merge(diff_anno, converter_df, by.x = 'gene', by.y = 'query')

  is.na(diff_anno) <- sapply(diff_anno, is.infinite)
  diff_anno <- na.omit(diff_anno)

  # diff_anno <- data.frame(gene = rnames,
  #                          logfc = diff_res[, logFC],
  #                          anno = diff_res$anno,
  #                          entrezgene = diff_res$entrezgene,
  #                          genename = diff_res$symbol,
  #                          description = diff_res$description)

  # make the list for GSEA
  gene_logFC <- diff_anno$logfc
  names(gene_logFC) <- diff_anno$entrezgene
  gsea_gene <- sort(gene_logFC, decreasing = TRUE)

  # make the list for GO and KEGG enrichment
  # go_up_gene <- diff_anno[diff_anno$anno == "Up", 'geneid']
  # go_down_gene <- diff_anno[diff_anno$anno == "Down", 'geneid']
  sig_gene <- diff_anno[diff_anno$anno == "Up" | diff_anno$anno == "Down", 'entrezgene']
  up_gene <- diff_anno[diff_anno$anno == "Up", 'entrezgene']
  down_gene <- diff_anno[diff_anno$anno == "Down", 'entrezgene']
  write.csv(file = 'sig_gene.txt', sig_gene, sep = '\t', quote = FALSE)
  write.csv(file = 'gsea_gene.txt', gsea_gene, sep = '\t', quote = FALSE)
  write.csv(file = 'diff_anno.txt', diff_anno, sep = '\t', quote = FALSE)

  rm(data)
  rm(diff_res)
  rm(converter_df)

  return(list(up_gene = up_gene,
              down_gene = down_gene,
              sig_gene = sig_gene,
              gsea_gene = gsea_gene,
              diff_anno = diff_anno))
}

############################################################
# GO
get_go_enrich <- function(my_OrgDb, gene_id, number, now_color, font_num, prefix_path, prefix_name) {
  ego <- enrichGO(gene_id, OrgDb = my_OrgDb,
                  keyType = "ENTREZID",
                  ont = "ALL",
                  pvalueCutoff = 1,
                  qvalueCutoff = 1)

  if (is.null(ego)) {
    cat(paste('group have none enrichment terms: ', prefix_name, '\n', sep = ''))
    file.create(paste(prefix_path, '/1_GO/', prefix_name, '_GO.csv', sep = ''))
  }else if (min(ego$p.adjust) >= 0.05) {
    cat(paste0('group have none enrichment terms with pvalue < 0.05: ', prefix_name, '\n'))
    file.create(paste0(prefix_path, '/1_GO/', prefix_name, '_GO_dot.pdf'))
  } else {
    ego_res <- data.frame(ego@result)
    write.csv(x = ego_res,
              file = paste(prefix_path, '/1_GO/', prefix_name, '_GO.csv', sep = ''),
              row.names = F)

    ego_dotplot <- dotplot(ego,
                           orderBy = "x",
                           x = "GeneRatio",
                           color = "p.adjust",
                           font.size = font_num,
                           showCategory = number) +
      scale_y_discrete(labels = function(x) str_wrap(x, width = 80)) +
      scale_colour_gradient(low = now_color[1], high = now_color[2])

    # scale_color_distiller(palette = "Spectral")

    ggsave(plot = ego_dotplot,
           filename = paste(prefix_path, '/1_GO/', prefix_name, '_GO_dot.pdf', sep = ''),
           width = 9,
           height = 7 + 0.1 * number)
    ggsave(plot = ego_dotplot,
           filename = paste(prefix_path, '/1_GO/', prefix_name, '_GO_dot.png', sep = ''),
           width = 9,
           height = 7 + 0.1 * number,
           dpi = 300)

    ego_barplot <- dotplot(ego,
                           split = "ONTOLOGY",
                           font.size = font_num,
                           showCategory = round(number / 3)) +
      scale_colour_gradient(low = now_color[1], high = now_color[2]) +
      facet_grid(ONTOLOGY ~ ., scale = "free") +
      scale_y_discrete(labels = function(x) str_wrap(x, width = 80))

    # scale_color_distiller(palette = "Spectral")

    ggsave(plot = ego_barplot,
           filename = paste0(prefix_path, '/1_GO/', prefix_name, '_GO_dot_split.pdf'),
           width = 9,
           height = 5)
    ggsave(plot = ego_barplot,
           filename = paste0(prefix_path, '/1_GO/', prefix_name, '_GO_dot_split.png'),
           width = 9,
           height = 7 + 0.1 * number,
           dpi = 300)

    rm(ego_dotplot)
    rm(ego_barplot)

    return(ego)
  }
}

############################################################
# group GO enrich
group_go <- function(my_OrgDb,
                     gene_list,
                     keyType = "ENTREZID",
                     ont = "BP",
                     number = 30,
                     now_color,
                     font_num = 12,
                     prefix_path,
                     prefix_name) {
  tryCatch(
    res <- compareCluster(geneCluster = gene_list,
                          fun = enrichGO,
                          OrgDb = my_OrgDb,
                          keyType = keyType,
                          ont = ont,
                          pvalueCutoff = 1,
                          qvalueCutoff = 1),
    error = function(e) { cat("error!") }
  )
  if (is.null(res)) {
    cat(paste('group have none enrichment terms: ', prefix_name, '\n', sep = ''))
    file.create(paste(prefix_path, '/1_GO/', prefix_name, '_', ont, '_group.csv', sep = ''))
  }else if (min(res@compareClusterResult$p.adjust) >= 0.05) {
    cat(paste0('group have none enrichment terms with pvalue < 0.05: ', prefix_name, '\n'))
    file.create(paste0(prefix_path, '/1_GO/', prefix_name, '_', ont, '_group_dot.pdf'))
  } else {
    # save the enrichment results
    write.csv(x = res@compareClusterResult,
              file = paste0(prefix_path, '/1_GO/', prefix_name, '_', ont, '_group_enrich.csv'),
              row.names = F)
    # draw the figure
    go_plot <- dotplot(res,
                       color = "p.adjust",
                       font.size = font_num,
                       showCategory = round(number / 3)) +
      scale_y_discrete(labels = function(x) str_wrap(x, width = 80)) +
      scale_colour_gradient(low = now_color[1], high = now_color[2])
    # save the figure
    ggsave(
      plot = go_plot,
      filename = paste0(prefix_path, '/1_GO/', prefix_name, '_', ont, '_group_dot.pdf'),
      width = 10,
      height = 7 + 0.1 * number
    )
  }
}

############################################################
# KEGG
get_kegg_enrich <- function(gsea_gene, gene_id, number, now_color, font_num, prefix_path, prefix_name, species) {
  my_kegg_db <- clusterProfiler::download_KEGG(species)
  # head(my_kegg_db$KEGGPATHID2NAME)
  # head(my_kegg_db$KEGGPATHID2EXTID)

  term2name <- my_kegg_db$KEGGPATHID2NAME
  term2gene <- my_kegg_db$KEGGPATHID2EXTID

  EntrezID <- bitr_kegg(my_kegg_db$KEGGPATHID2EXTID$to,
                        fromType = 'kegg',
                        toType = 'ncbi-geneid',
                        organism = species)
  mergedID <- merge(term2gene, EntrezID,
                    by.x = 'to', by.y = 'kegg')
  term2gene <- data.frame(from = mergedID$from,
                          to = mergedID[, 3])

  egg <- enricher(gene = gene_id,
                  TERM2GENE = term2gene,
                  TERM2NAME = term2name,
                  pvalueCutoff = 1,
                  qvalueCutoff = 1)

  if (is.null(egg)) {
    cat(paste0('group have none enrichment terms: ', prefix_name, '\n'))
    file.create(paste0(prefix_path, '/2_KEGG/', prefix_name, '_KEGG.csv'))
  } else if (min(egg$p.adjust) >= 0.05) {
    cat(paste0('group have none enrichment terms with pvalue < 0.05: ', prefix_name, '\n'))
    file.create(paste0(prefix_path, '/2_KEGG/', prefix_name, '_KEGG_dot.pdf'))
  } else {
    egg_res <- data.frame(egg@result)
    write.csv(x = egg_res,
              file = paste0(prefix_path, '/2_KEGG/', prefix_name, '_KEGG.csv'),
              row.names = F)

    egg_dotplot <- dotplot(egg,
                           orderBy = "x",
                           x = "GeneRatio",
                           color = "p.adjust",
                           font.size = font_num,
                           showCategory = number) +
      scale_y_discrete(labels = function(x) str_wrap(x, width = 80)) +
      scale_colour_gradient(low = now_color[1], high = now_color[2])

    ggsave(plot = egg_dotplot,
           filename = paste(prefix_path, '/2_KEGG/', prefix_name, '_KEGG_dot.pdf', sep = ''),
           width = 8,
           height = 5 + 0.1 * number)
    ggsave(plot = egg_dotplot,
           filename = paste(prefix_path, '/2_KEGG/', prefix_name, '_KEGG_dot.png', sep = ''),
           width = 8,
           height = 5 + 0.1 * number,
           dpi = 300)

    rm(egg_dotplot)

    return(egg)
  }
}

############################################################
# group kegg enrich
group_kegg <- function(gene_list,
                       keyType = "ncbi-geneid",
                       organism = organism,
                       number = 30,
                       now_color,
                       font_num = 12,
                       prefix_path,
                       prefix_name) {
  tryCatch(
    res <- compareCluster(geneCluster = gene_list,
                          fun = enrichKEGG,
                          organism = organism,
                          keyType = keyType,
                          pvalueCutoff = 1,
                          qvalueCutoff = 1),
    error = function(e) { cat("error!") }
  )
  if (is.null(res)) {
    cat(paste('group have none enrichment terms: ', prefix_name, '\n', sep = ''))
    file.create(paste(prefix_path, '/2_KEGG/', prefix_name, '_kegg_group.csv', sep = ''))
  }else if (min(res@compareClusterResult$p.adjust) >= 0.05) {
    cat(paste0('group have none enrichment terms with pvalue < 0.05: ', prefix_name, '\n'))
    file.create(paste0(prefix_path, '/2_KEGG/', prefix_name, '_kegg_group_dot.pdf'))
  } else {
    # save the enrichment results
    write.csv(x = res@compareClusterResult,
              file = paste0(prefix_path, '/2_KEGG/', prefix_name, '_kegg_group_enrich.csv'),
              row.names = F)
    # draw the figure
    kegg_plot <- dotplot(res,
                         color = "p.adjust",
                         font.size = font_num,
                         showCategory = round(number / 3)) +
      scale_y_discrete(labels = function(x) str_wrap(x, width = 80)) +
      scale_colour_gradient(low = now_color[1], high = now_color[2])
    # save the figure
    ggsave(
      plot = kegg_plot,
      filename = paste0(prefix_path, '/2_KEGG/', prefix_name, '_kegg_group_dot.pdf'),
      width = 10,
      height = 7 + 0.1 * number
    )
  }
}

############################################################
# GSEA
get_gsea_enrich <- function(my_OrgDb, gsea_gene, number, font_num, prefix_path, prefix_name, species) {
  # go gsea enrichment
  go_gsea <- gseGO(gsea_gene,
                   OrgDb = my_OrgDb,
                   keyType = "ENTREZID",
                   minGSSize = 10,
                   maxGSSize = 1000,
                   pvalueCutoff = 0.05)

  if (is.null(go_gsea)) {
    cat(paste0('group have none enrichment terms: ', prefix_name, '\n'))
    file.create(paste0(prefix_path, '/3_GO_GSEA/', prefix_name, '_go_gsea.png'))
  } else if (min(go_gsea$p.adjust) >= 0.05) {
    cat(paste0('group have none enrichment terms with pvalue < 0.05: ', prefix_name, '\n'))
  } else {
    write.csv(x = go_gsea@result,
              file = paste0(prefix_path, '/3_GO_GSEA/', prefix_name, '_GO_GSEA.csv'),
              row.names = F)

    go_gsea_ID <- head(go_gsea@result$ID, number)
    for (id_name in go_gsea_ID) {
      go_gsea_plot <- gseaplot2(x = go_gsea,
                                geneSetID = id_name,
                                base_size = font_num - 1,
                                pvalue_table = TRUE)

      ggsave(plot = go_gsea_plot,
             filename = paste0(prefix_path, '/3_GO_GSEA/', str_replace(id_name, ':', ''), '_go_gsea.pdf'),
             width = 8, height = 5)
      ggsave(plot = go_gsea_plot,
             filename = paste0(prefix_path, '/3_GO_GSEA/', str_replace(id_name, ':', ''), '_go_gsea.png'),
             width = 8, height = 5, dpi = 300)
    }
  }

  # kegg gsea enrichment
  kegg_gsea <- gseKEGG(gsea_gene,
                       minGSSize = 10,
                       maxGSSize = 1000,
                       pvalueCutoff = 0.05,
                       organism = species,
                       keyType = "ncbi-geneid")

  if (is.null(kegg_gsea)) {
    cat(paste('group have none enrichment terms: ', prefix_name, '\n', sep = ''))
    file.create(paste(prefix_path, '/4_KEGG_GSEA/', prefix_name, '_kegg_gsea.png', sep = ''))
  } else if (min(kegg_gsea$p.adjust) >= 0.05) {
    cat(paste0('group have none enrichment terms with pvalue < 0.05: ', prefix_name, '\n'))
  } else {
    write.csv(x = kegg_gsea@result,
              file = paste(prefix_path, '/4_KEGG_GSEA/', prefix_name, '_KEGG_GSEA.csv', sep = ''),
              row.names = F)

    kegg_gsea_ID <- head(kegg_gsea@result$ID, number)
    for (id_name in kegg_gsea_ID) {
      kegg_gsea_plot <- gseaplot2(x = kegg_gsea,
                                  geneSetID = id_name,
                                  base_size = font_num - 1,
                                  pvalue_table = TRUE)
      ggsave(plot = kegg_gsea_plot,
             filename = paste(prefix_path, '/4_KEGG_GSEA/', id_name, '_kegg_gsea.pdf', sep = ''),
             width = 8, height = 5)
      ggsave(plot = kegg_gsea_plot,
             filename = paste(prefix_path, '/4_KEGG_GSEA/', id_name, '_kegg_gsea.png', sep = ''),
             width = 8, height = 5, dpi = 300)
    }
  }
}

############################################################
# download pathway
down_pathview <- function(egg_sig, gsea_gene, species) {
  # draw the pathway
  path_ID <- head(egg_sig@result$ID, 5)
  for (id_name in path_ID) {
    tryCatch(
    { pathview(gene.data = gsea_gene,
               pathway.id = id_name,
               species = species,
               limit = list(gene = max(abs(gsea_gene)), cpd = 1))
    },
      error = function(e) {
        message('\nError @ ', id_name)
        return(NA) }
    )
  }
}

############################################################
# main program is here
cat("=================================== \n")
cat(paste(as.character(Sys.time()), "\n", sep = ""))
cat(paste("Step3: Check the input files: ", "\n", sep = ""))
check_files(args$input, args$output)
prefix_path <- getwd()
prefix_name <- unlist(strsplit(basename(args$input), "[.]"))[1]
cat(paste("prefix_path: ", prefix_path, "\n", sep = ""))
cat(paste("prefix_name: ", prefix_name, "\n", sep = ""))

cat("=================================== \n")
cat(paste(as.character(Sys.time()), "\n", sep = ""))
cat(paste("Step4: Build the species OrgDb: ", args$species, "\n", sep = ""))
my_OrgDb <- make_orgdb(args$taxon, prefix_path)

cat("=================================== \n")
cat(paste(as.character(Sys.time()), "\n", sep = ""))
cat(paste("Step5: Import the DESeq2/edgeR results: ", args$input, "\n", sep = ""))
diff_res <- import_diff(args$input, args$padj, args$log2FC, args$taxon)
up_gene <- diff_res$up_gene
down_gene <- diff_res$down_gene
sig_gene <- diff_res$sig_gene
gsea_gene <- diff_res$gsea_gene
diff_anno <- diff_res$diff_anno

cat("=================================== \n")
cat(paste(as.character(Sys.time()), "\n", sep = ""))
cat(paste("Step6: Run the go enrichment: ", "\n", sep = ""))
dir.create(paste(prefix_path, '/1_GO', sep = ''))
# run the whole gene go enrichment
if (!is.null(sig_gene)) {
  ego_sig <- get_go_enrich(my_OrgDb, sig_gene, args$number, now_color, args$fontsize, prefix_path, prefix_name)

}else {
  cat(paste('group have none significant genes: ', prefix_name, '\n', sep = ''))
  file.create(paste(prefix_path, '/', prefix_name, '_no_sig_genes.txt', sep = ''))
}
# run the down/up gene go enrichment
if (!is.null(up_gene) | !is.null(down_gene)) {
  gene_list <- list(up = up_gene, down = down_gene)
  group_go(my_OrgDb, gene_list, "ENTREZID", "BP", args$number, now_color, args$fontsize, prefix_path, prefix_name)
  group_go(my_OrgDb, gene_list, "ENTREZID", "CC", args$number, now_color, args$fontsize, prefix_path, prefix_name)
  group_go(my_OrgDb, gene_list, "ENTREZID", "MF", args$number, now_color, args$fontsize, prefix_path, prefix_name)
}

cat("=================================== \n")
cat(paste(as.character(Sys.time()), "\n", sep = ""))
cat(paste("Step7: Run the kegg enrichment: ", "\n", sep = ""))
dir.create(paste0(prefix_path, '/2_KEGG'))
# run the whole gene kegg enrichment
if (!is.null(sig_gene)) {
  egg_sig <- get_kegg_enrich(gsea_gene, sig_gene, args$number, now_color, args$fontsize, prefix_path, prefix_name, args$species)
}else {
  cat(paste0('group have none significant genes: ', prefix_name, '\n'))
  file.create(paste0(prefix_path, '/', prefix_name, '_no_sig_genes.txt'))
}
# run the down/up gene kegg enrichment
if (!is.null(up_gene) | !is.null(down_gene)) {
  gene_list <- list(up = up_gene, down = down_gene)
  group_kegg(gene_list, "ncbi-geneid", args$species, args$number, now_color, args$fontsize, prefix_path, prefix_name)
  group_kegg(gene_list, "ncbi-geneid", args$species, args$number, now_color, args$fontsize, prefix_path, prefix_name)
  group_kegg(gene_list, "ncbi-geneid", args$species, args$number, now_color, args$fontsize, prefix_path, prefix_name)
}

cat("=================================== \n")
cat(paste0(as.character(Sys.time()), "\n"))
cat(paste0("Step8: Run the gsea enrichment: ", "\n"))
dir.create(paste0(prefix_path, '/3_GO_GSEA'))
dir.create(paste0(prefix_path, '/4_KEGG_GSEA'))
get_gsea_enrich(my_OrgDb, gsea_gene, number, args$fontsize, prefix_path, prefix_name, args$species)
# run_GSEA(gsea_gene, number, prefix_name, args$species)

cat("=================================== \n")

cat(paste0(as.character(Sys.time()), "\n"))
cat(paste0("Step9: download the pathway figure: ", "\n"))
if (!is.logical(egg_sig)) {
  setwd(paste0(prefix_path, '/2_KEGG'))
  down_pathview(egg_sig, gsea_gene, args$species)
}
############################################################
# DONE
cat(paste0("\n", "All enrichment pipeline is done!", "\n"))
