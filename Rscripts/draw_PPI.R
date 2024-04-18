################################################################
# Draw the gene interaction network (PPI)
# Rensc
# 2021-8-12
################################################################
# install the packages
# if (!requireNamespace("getopt", quietly = TRUE))
#   install.packages("getopt")
# if (!requireNamespace("tidyverse", quietly = TRUE))
#   install.packages("tidyverse")
# if (!requireNamespace("dplyr", quietly = TRUE))
#   install.packages("dplyr")
# if (!requireNamespace("igraph", quietly = TRUE))
#   install.packages("igraph")
# if (!requireNamespace("ggraph", quietly = TRUE))
#   install.packages("ggraph")
# if (!requireNamespace("STRINGdb", quietly = TRUE))
#   BiocManager::install("STRINGdb")
################################################################
# get the options
suppressMessages(require(getopt, quietly = T))
spec <- matrix(
  c("input", "i", 2, "character", "Input the gene list.",
    "output", "f", 2, "character", "output figure file name",
    "taxon", "t", 2, "double", "Specify the taxon id  (default same as string_db).",
    "string_db", "d", 2, "double", "Specify the string_db taxon id.",
    "score", "s", 2, "double", "Specify the interaction score (default 400).",
    "edge", "e", 2, "double", "Specify the minimum interaction edge (default 1).",
    "help", "h", 0, "logical", "Draw the figure of string PPI."),
  byrow = TRUE, ncol = 5)

args <- getopt(spec = spec)

# check the options
if (!is.null(args$help) || is.null(args$input)) {
  # print usage
  cat(paste("\n", getopt(spec = spec, usage = T), "\n"))
  quit(status = 1)
}

# set the default options
if (is.null(args$score)) { args$score <- 400 }
if (is.null(args$edge)) { args$edge <- 1 }
if (is.null(args$taxon)) { args$taxon <- args$string_db }
if (is.null(args$string_db)) { args$string_db <- args$taxon }

# print all options
cat("=================================== \n")
cat(paste(as.character(Sys.time()), "\n", sep = ""))
cat(paste("Step1: Get options", "\n\n", sep = ""))
cat(paste("input:", args$input, "\n"))
cat(paste("output:", args$output, "\n"))
cat(paste("taxon: ", args$taxon, "\n"))
cat(paste("string_db: ", args$string_db, "\n"))
cat(paste("score: ", args$score, "\n"))
cat(paste("edge: ", args$edge, "\n\n"))

############################################################
# import the packages
cat("=================================== \n")
cat(paste(as.character(Sys.time()), "\n", sep = ""))
cat(paste("Step2: import the packages", "\n\n", sep = ""))
options(warn = -1)
suppressMessages(require(tidyverse, quietly = T))
suppressMessages(require(dplyr, quietly = T))
suppressMessages(require(STRINGdb, quietly = T))
suppressMessages(require(igraph, quietly = T))
suppressMessages(require(ggraph, quietly = T))
suppressMessages(require(tools, quietly = T))
############################################################
# check the files
check_files <- function(input_file, figure) {
  # check the input file
  if (!file.exists(input_file)) {
    cat(paste("\n", "Miss DESeq2 results file", "\n"))
    quit(status = 2)
  }

  # set the output path
  output_path <- dirname(figure)
  if (file.exists(output_path)) {
    setwd(output_path)
  }else {
    dir.create(output_path)
    setwd(output_path)
  }

}

############################################################
# import the deseq
import_deseq <- function(now_filename, now_padj, now_log2fc) {
  
  # read data
  if (file_ext(now_filename) == "csv"){
    data <- read.csv(now_filename)
    
  }else if (file_ext(now_filename) == "txt" || file_ext(now_filename) == "TXT"){
    data <- read.table(file = now_filename, sep = "\t", header = T, row.names = NULL, comment.char = '#')
  }
  
  deseq_res <- data.frame(row.names = data[, 1], data[, -1])
  
  # filter significant genes
  deseq_res$anno <- ifelse(deseq_res$padj >= now_padj |
                             deseq_res$padj %in% NA |
                             abs(deseq_res$log2FoldChange) < now_log2fc, 'NS',
                           ifelse(deseq_res$padj < now_padj & deseq_res$log2FoldChange >= now_log2fc, 'Up', 'Down'))
  deseq_res$anno <- factor(deseq_res$anno,
                           levels = c("Up", "Down", "NS"))
  
  up_num <- length(which(deseq_res$anno == "Up"))
  down_num <- length(which(deseq_res$anno == "Down"))
  
  min_padj <- min(deseq_res$padj[deseq_res$padj > 0], na.rm = T)
  deseq_res$padj[deseq_res$padj == 0] <- min_padj / 100
  
  rnames <- rownames(deseq_res)
  
  deseq_anno <- data.frame(gene = rnames,
                           logfc = deseq_res$log2FoldChange,
                           padj = deseq_res$padj,
                           anno = deseq_res$anno,
                           entrezgene = deseq_res$entrezgene,
                           genename = deseq_res$symbol,
                           description = deseq_res$description)
  
  # make the list for network
  up_gene <- deseq_anno[deseq_anno$anno == "Up", ]
  if (nrow(up_gene) > 200){
    up_gene <- up_gene[order(up_gene$padj), ]
    up_gene <- head(up_gene, 200)
  }
  down_gene <- deseq_anno[deseq_anno$anno == "Down", ]
  if (nrow(down_gene) > 200){
    down_gene <- down_gene[order(down_gene$padj), ]
    down_gene <- head(down_gene, 200)
  }
  sig_gene <- rbind(up_gene, down_gene)
  
  return(sig_gene)
}

############################################################
# import the gene list
import_gene <- function(now_filename) {
  gene_list <- read.table(file = now_filename, sep = "\t", header = T, row.names = NULL, comment.char = '#')
  
  
}




############################################################
make_string_db <- function(taxon, sig_gene, ppi_score, nodes_name, edges_name, figure) {
  string_db <- STRINGdb$new(version="11",
                            species=taxon,
                            score_threshold=ppi_score,
                            input_directory="")

  data_mapped <- string_db$map(my_data_frame = sig_gene,
                               my_data_frame_id_col_names = "entrezgene",
                               removeUnmappedRows = TRUE)
  hit <- data_mapped$STRING_id
  net_info <- string_db$get_interactions(hit)

  write.table(data_mapped, file = nodes_name, sep = '\t', quote = FALSE, col.names = T, row.names = F)
  write.table(net_info, file = edges_name, sep = '\t', quote = FALSE, col.names = T, row.names = F)

  pdf(file=figure, width = 6, height = 6)
  network_plot <- string_db$plot_network(hit)
  dev.off()

  return(list(string_db=string_db,
              data_mapped=data_mapped,
              net_info=net_info))
}

############################################################
draw_ggraph <- function(gene_res, edge_num, prefix_name) {
  data_mapped <- gene_res['data_mapped']
  net_info <- gene_res['net_info']

  links <- net_info %>%
    mutate(from = data_mapped[match(from, data_mapped$STRING_id), "symbol"]) %>%
    mutate(to = data_mapped[match(to, data_mapped$STRING_id), "symbol"]) %>%
    dplyr::select(from, to , last_col()) %>%
    dplyr::rename(weight = combined_score)

  links_2 <- links %>% mutate(from_c = count(., from)$n[match(from, count(., from)$from)]) %>%
    mutate(to_c = count(., to)$n[match(to, count(., to)$to)]) %>%
    filter(!(from_c >= edge_num & to_c >= edge_num)) %>%
    dplyr::select(1,2,3)

  if (nrow(links_2) > 0){
    # get the nodes
    nodes_2 <- links_2 %>% { data.frame(gene = c(.$from, .$to)) } %>% distinct()

    # create the network
    net_2 <- igraph::graph_from_data_frame(d=links_2,vertices=nodes_2,directed = F)

    # get the details of network
    igraph::V(net_2)$deg <- igraph::degree(net_2)
    igraph::V(net_2)$size <- igraph::degree(net_2)/5
    igraph::E(net_2)$weight <- igraph::E(net_2)$weight/10

    linear_fig <- ggraph(net_2, layout = "linear", circular = TRUE)+
                    geom_edge_arc(aes(edge_width=width), color = "lightblue", show.legend = F)+
                    geom_node_point(aes(size=size), color="orange", alpha=0.7)+
                    geom_node_text(aes(filter=deg>5, label=name), size = 5, repel = F)+
                    scale_edge_width(range = c(0.2,1))+
                    scale_size_continuous(range = c(1,10) )+
                    guides(size=F)+
                    theme_graph()

    kk_fig <- ggraph(net_2, layout = "kk")+
                geom_edge_fan(aes(edge_width=width), color = "lightblue", show.legend = F)+
                geom_node_point(aes(size=size), color="orange", alpha=0.7)+
                geom_node_text(aes(filter=deg>5, label=name), size = 5, repel = T)+
                scale_edge_width(range = c(0.2,1))+
                scale_size_continuous(range = c(1,10) )+
                guides(size=F)+
                theme_graph()

    ggsave(paste(prefix_name, "-linear.pdf", sep = ""),
           linear_fig, width = 6, height = 6)
    ggsave(paste(prefix_name, "-kk.pdf", sep = ""),
           kk_fig, width = 6, height = 6)
  }

}
############################################################
# main program is here
cat("=================================== \n")
cat(paste(as.character(Sys.time()), "\n", sep = ""))
cat(paste("Step3: Check the input files: ", "\n", sep = ""))
check_files(args$input, args$figure)
prefix_name <- paste(dirname(args$figure), unlist(strsplit(basename(args$figure), '[.pdf]')[1]),sep='')

cat("=================================== \n")
cat(paste(as.character(Sys.time()), "\n", sep = ""))
cat(paste("Step4: Import the deseq2 results: ", args$input, "\n", sep = ""))
sig_gene <- import_deseq(args$input, args$padj, args$log2FC)
# 
# cat("=================================== \n")
# cat(paste(as.character(Sys.time()), "\n", sep = ""))
# cat(paste("Step5: Convert the gene name to entrezID: ", "\n", sep = ""))
# sig_gene <- convert_geneID(deseq_res, args$taxon)

cat("=================================== \n")
cat(paste(as.character(Sys.time()), "\n", sep = ""))
cat(paste("Step5: Build the PPI network: ", args$string_db, "\n", sep = ""))
if (!is.null(sig_gene)){
  sig_gene_res <- make_string_db(args$string_db, sig_gene, args$score, args$nodes, args$edges, args$figure )
}else{
  cat(paste('group have none significant genes: ', '\n', sep = ''))
  file.create(paste(prefix_name, '_no_sig_up_genes.txt', sep = ''))
}

# cat("=================================== \n")
# cat(paste(as.character(Sys.time()), "\n", sep = ""))
# cat(paste("Step6: Draw the PPI network with ggraph: ", args$taxon, "\n", sep = ""))
# if (!is.null(sig_gene)){
#   draw_ggraph(sig_gene_res, args$edgenum, prefix_name)
# }else{
#   cat(paste('group have none significant genes: ', '\n', sep = ''))
#   file.create(paste(prefix_name, '_no_sig_up_genes.txt', sep = ''))
# }


############################################################
# DONE
cat(paste("\n", "ALL DONE!", "\n"))

