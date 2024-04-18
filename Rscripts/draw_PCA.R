################################################################
# Draw the PCA plot
# Rensc
# 2021-8-3
################################################################
# install the packages
# if (!requireNamespace("getopt", quietly = TRUE))
#   install.packages("getopt")
# if (!requireNamespace("FactoMineR", quietly = TRUE))
#   install.packages("FactoMineR")
# if (!requireNamespace("factoextra", quietly = TRUE))
#   install.packages("factoextra")
# if (!requireNamespace("cowplot", quietly = TRUE))
#   install.packages("cowplot")
# if (!requireNamespace("dplyr", quietly = TRUE))
#   install.packages("dplyr")

############################################################
# get the options
suppressMessages(require(getopt, quietly = T))
spec <- matrix(
  c("input", "i", 1,"character", "Input dataset of gene expression",
    "output", "o", 1, "character", "output path name",
    "design", "d", 1, "character", "Phenotype file which contian [samples, groups]",
    "threshold", "t", 2, "double", "Average threshold of all groups (default 1)",
    "log2", 'l', 2, "character","log the data or not (y/n default n)",
    "help", "h", 0, "logical", "PCA analysis and plot!"),
  byrow = TRUE, ncol = 5)

args <- getopt(spec = spec)

# check the options
if ( !is.null(args$help) || is.null(args$input) || is.null(args$design)) {
  # print useage 
  cat(paste("\n",getopt(spec = spec, usage = T), "\n"))
  quit(status = 1)
}

if (!file.exists(args$output)) {
  dir.create(args$output)
}

# set the default options
if ( is.null(args$output ) ) { args$output <- args$input }
if ( is.null(args$threshold) ) { args$threshold <- 1 }
if ( is.null(args$log2) ) { args$log2 <- "n" }

# print all options
cat("=================================== \n")
cat(paste(as.character(Sys.time()), "\n", sep = ""))
cat(paste("Step1: Get options", "\n\n", sep = ""))
cat(paste("input: ", args$input, "\n"))
cat(paste("output: ", args$output, "\n"))
cat(paste("design: ", args$design, "\n"))
cat(paste("threshold: ", args$threshold, "\n"))
cat(paste("log2: ", args$log2, "\n\n"))

############################################################
# import PCA package
cat("=================================== \n")
cat(paste(as.character(Sys.time()), "\n", sep = ""))
cat(paste("Step2: import PCA packages", "\n\n", sep = ""))
options(warn  = -1)
require(ggplot2, quietly = T)
suppressMessages(require(dplyr, quietly = T))
suppressMessages(require(FactoMineR, quietly = T))
suppressMessages(require(factoextra, quietly = T))
suppressMessages(require(cowplot, quietly = T))
suppressMessages(require(ggforce, quietly = T))
suppressMessages(require(tools, quietly = T))
############################################################
# read results
cat(paste(as.character(Sys.time()), ":\t", "Read expresstion table", "\n", sep = ""))

if (file_ext(args$input) == "csv"){
  data <- read.csv(file = args$input, header = T, row.names = NULL, comment.char = '#')
  
}else if (file_ext(args$input) == "txt" || file_ext(args$input) == "TXT"){
  data <- read.table(file = args$input, sep = "\t", header = T, row.names = NULL, comment.char = '#')
}

rowns <- as.character(data[, 1])
rownames(data) <- rowns
data <- data[,-1:-6]
colnames(data) <- gsub('[_-]', '.', colnames(data))
# col_name <- gsub('.*[.]bam.', '', colnames(data))
# col_name <- gsub('.bam', '', col_name, fixed = TRUE)
# colnames(data) <- col_name

############################################################
# make the phenotype list
cat(paste(as.character(Sys.time()), ":\t", "Make the phenotype table", "\n", sep = ""))

prefix_name <- unlist(strsplit(basename(args$input), "[.]"))[1]
sample <- read.csv(file = args$design, sep = "\t", header = T, row.names = NULL)
sample_name <- gsub('[_-]', '.', sample[, 1])
sample[, 1] <- sample_name

class_list <- unique(sample[, 2])
sample_mat <- data.frame(table(sample[, 2])[class_list])

cname <- data.frame(name = colnames(data))
colnames(sample)[1] <- 'name'
pheno <- left_join(cname, sample, copy = T, by = "name")
# write.table(sample_mat,file='sample_mat.txt')
############################################################
# split the data by phenotype
dat <- rowns
for (s in sample_mat$Var1) {
  
  ids <- which(pheno[, 2] ==  s)
  # ids <- grep(s,pheno[, 2])
  samp_data <- data[, as.numeric(ids)]
  dat <- cbind(dat, samp_data)
  
}
dat <- dat[, -1]

############################################################
# set the colors of samples

color_set <- c(RColorBrewer::brewer.pal(name = 'Set2',n = 7),
               RColorBrewer::brewer.pal(name = 'Set1',n = 8),
               RColorBrewer::brewer.pal(name = 'Set3',n = 8))

class_nums <- length(class_list)

anno_color <- data.frame(row.names = unique(sample[, 2]), gorups = color_set[1:class_nums])

anno_color <- t(anno_color)[, ]

############################################################
# filter the data by the threshold
cat(paste(as.character(Sys.time()), ":\t", "Filter data by threshold", " ",args$threshold,"\n", sep = ""))

FilteredMat <- subset(dat, rowMeans(dat) >= args$threshold)

if ( args$log2 ==  "n") {
  Mat <- FilteredMat
  write.table(Mat, file = paste(args$output, '/', prefix_name, "_PCA_th_", args$threshold, ".txt", sep = ''),
              col.names = NA, sep = "\t", quote = F)
  
} else if (args$log2 ==  "y") {
  Mat <- log2(FilteredMat + 1)
  rawname <- colnames(FilteredMat)
  logname <- paste(rawname, "_log2",  sep = "")
  newdata_cnames <- c(rawname,logname)
  
  mer <- cbind(FilteredMat, Mat)
  colnames(mer) <- newdata_cnames
  write.table(mer, file = paste(args$output, '/', prefix_name, "_PCA_th_", args$threshold, "_log2.txt", sep = ''),
              col.names = NA, sep = "\t", quote = F)
}

############################################################
# PCA main analysis
cat(paste(as.character(Sys.time()), ":\t", "Run the PCA analysis", "\n", sep = ""))

res.pca <- PCA(Mat, graph = FALSE, scale = T)
res.f.pca <- PCA(Mat, graph = FALSE, scale = F)

# write the PCA eigenvalues
pc_eig <- get_eig(res.pca)

write.table(pc_eig, file = paste(args$output, '/',prefix_name, "_PCA_eigenvalues_th", args$threshold, ".txt", sep = ''),
            col.names = NA, sep = "\t", quote = F)

############################################################
# Draw plot1: get pca eigenvalues
cat("=================================== \n")
cat(paste(as.character(Sys.time()), "\n", sep = ""))
cat(paste("Step3: Draw the plots of PCA", "\n\n", sep = ""))
cat(paste(as.character(Sys.time()), ":\t", "Plot1: PCA eigenvalues", "\n", sep = ""))
# save the default Rplots.pdf
pdf(file=paste(args$output,'Rplots.pdf',sep=''))

plot1 <- fviz_eig(res.pca, addlabels = TRUE) +
  theme(text = element_text(size = 13),
        panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white"),
        plot.title = element_text(size = 16, hjust = 0.5))

ggsave(paste(args$output, '/', prefix_name, "_PCA_plot1_eigenvalues_th", args$threshold, ".png", sep = ""),
       plot1, width = 5, height = 4, dpi = 300)
ggsave(paste(args$output, '/', prefix_name, "_PCA_plot1_eigenvalues_th", args$threshold, ".pdf", sep = ""),
       plot1, width = 5, height = 4)

# write the pca contribution
var <- get_pca_var(res.pca)

write.table(var$contrib, file = paste(args$output, '/', prefix_name, "_PCA-contrib-samples_th", args$threshold,".txt",sep = ''),
            col.names = NA, sep = "\t", quote = F)

ind <- get_pca_ind(res.pca)

write.table(ind$contrib, file = paste(args$output, '/', prefix_name, "_PCA-contrib-factors_th",args$threshold,".txt",sep = ''),
            col.names = NA, sep = "\t", quote = F)

############################################################
# Draw plot2: top contrib factors
cat(paste(as.character(Sys.time()), ":\t", "Plot2: top PC contribution factors by scale data", "\n", sep = ""))

p1 <- fviz_pca_var(res.pca, title = "top10-contribution-samples(scale)", col.var = "contrib",
                   geom = c("point", "text"), gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                   select.var = list(contrib = 10), labelsize = 3, pointsize = 3, repel = T) +
  theme(text = element_text(size = 13),
        panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white"),
        plot.title = element_text(size = 16, hjust = 0.5))

p2 <- fviz_pca_var(res.pca, title = "all-contribution-samples(scale)", col.var = "contrib",
                   geom = "point", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                   pointsize = 2) +
  theme(text = element_text(size = 13),
        panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white"),
        plot.title = element_text(size = 16, hjust = 0.5))

p3 <- fviz_pca_var(res.f.pca, title = "top10-contribution-samples(no scale)", col.var = "contrib",
                   geom = c("point", "text"), gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                   select.var = list(contrib = 10), labelsize = 3, pointsize = 3, repel = T) +
  theme(text = element_text(size = 13),
        panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white"),
        plot.title = element_text(size = 16, hjust = 0.5))

p4 <- fviz_pca_var(res.f.pca, title = "all-contribution-samples(no scale)", col.var = "contrib",
                   geom = "point", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                   pointsize = 2) +
  theme(text = element_text(size = 13),
        panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white"),
        plot.title = element_text(size = 16, hjust = 0.5))

plot2 <- plot_grid(p1, p2, p3, p4, ncol = 2)

ggsave(paste(args$output, '/', prefix_name, "_PCA_plot2_top_contribution_samples_th",args$threshold,".png", sep = ""),
       plot2, width = 10, height = 8, dpi = 300)
ggsave(paste(args$output, '/', prefix_name, "_PCA_plot2_top_contribution_samples_th",args$threshold,".pdf", sep = ""),
       plot2, width = 10, height = 8)
print(plot2)

############################################################
# Draw plot3: top contrib factors
cat(paste(as.character(Sys.time()), ":\t", "Plot3: top PC contribution factors", "\n", sep = ""))

p5 <- fviz_pca_ind(res.pca, title = "top10-contribution-factors", col.ind = "contrib",
                   geom = c("point", "text"), gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                   select.ind = list(contrib = 10), labelsize = 3, pointsize = 3, repel = TRUE) +
  theme(text = element_text(size = 13),
        panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white"),
        plot.title = element_text(size = 16, hjust = 0.5))

p6 <- fviz_pca_ind(res.pca, title = "all-contribution-factors", col.ind = "contrib",
                   geom = "point", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                   pointsize = 2) +
  theme(text = element_text(size = 13),
        panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white"),
        plot.title = element_text(size = 16, hjust = 0.5))

plot3<- plot_grid(p5, p6, ncol = 2)

ggsave(paste(args$output, '/', prefix_name, "_PCA_plot3_top_contribution_factors_th",args$threshold,".png", sep = ""),
       plot3, width = 10, height = 4, dpi = 300)
ggsave(paste(args$output, '/', prefix_name, "_PCA_plot3_top_contribution_factors_th",args$threshold,".pdf", sep = ""),
       plot3, width = 10, height = 4)
print(plot3)

############################################################
# Draw plot4: cos2/contrib on all the dimensions
cat(paste(as.character(Sys.time()), ":\t", "Plot4: PCA factors clusters", "\n", sep = ""))

r = nrow(sample_mat)
res.km <- kmeans(res.pca$var$coord, centers = r, nstart = 25)
grp <- as.factor(res.km$cluster)

p7 <- fviz_pca_var(res.pca, col.var = grp, title = "PCA samples clusters(scale)",
                   legend.title = "Cluster",repel = TRUE) +
  theme(text = element_text(size = 12), 
        panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white"),
        plot.title = element_text(size = 14, hjust = 0.5))

p8 <- fviz_pca_var(res.f.pca, col.var = grp, title = "PCA samples clusters(no scale)",
                   legend.title = "Cluster",repel = TRUE) +
  theme(text = element_text(size = 13),
        panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white"),
        plot.title = element_text(size = 16, hjust = 0.5))

plot4 <- plot_grid(p7, p8, ncol = 2)

ggsave(paste(args$output, '/', prefix_name, "_PCA_plot4_samples_clusters_th", args$threshold, ".png", sep = ""),
       plot4, width = 10, height = 4, dpi = 300)
ggsave(paste(args$output, '/', prefix_name, "_PCA_plot4_samples_clusters_th", args$threshold, ".pdf", sep = ""), 
       plot4, width = 10, height = 4)
print(plot4)

############################################################
# Draw plot5: draw the PCA by design
cat(paste(as.character(Sys.time()), ":\t", "Plot5: PCA factors design", "\n", sep = ""))

k = rep(sample_mat$Var1,sample_mat$Freq)

x_left <- min(res.f.pca$var$coord[,1]) - abs(min(res.f.pca$var$coord[,1]) * 0.15)
y_down <- min(res.f.pca$var$coord[,2]) - abs(min(res.f.pca$var$coord[,2]) * 0.15)
x_right <- max(res.f.pca$var$coord[,1]) + abs(max(res.f.pca$var$coord[,1]) * 0.15)
y_up <- max(res.f.pca$var$coord[,2]) + abs(max(res.f.pca$var$coord[,2]) * 0.15)

plot5 <- fviz_pca_var(res.f.pca, title = "PCA",
                      geom = c("text", "point"), 
                      pointshape = 19, 
                      pointsize = 3, 
                      col.var = k,
                      palette = anno_color,
                      label = "all",
                      # addEllipses = TRUE, ellipse.type = "norm", mean.point = F,
                      legend.title = "design",
                      repel = TRUE) +
  coord_cartesian(clip = "off")+
  theme(text = element_text(size = 13), 
        panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white"),
        # panel.border = element_rect(colour = 'black', fill = NA),
        plot.title = element_text(size = 16, hjust = 0.5)) +
  # geom_text(aes(label = rownames(data))) +
  xlim(x_left, x_right) + ylim(y_down, y_up)

ggsave(paste(args$output, '/', prefix_name, "_PCA_plot5_samples_eig_th", args$threshold, ".png", sep = ""), 
       plot5, width = 5, height = 4, dpi = 600)

ggsave(paste(args$output, '/', prefix_name, "_PCA_plot5_samples_eig_th", args$threshold, ".pdf", sep = ""), 
       plot5, width = 5, height = 4)
# close and delete the default Rplots.pdf
dev.off()
file.remove(paste(args$output,'Rplots.pdf',sep=''))
############################################################
# End here
cat(paste(as.character(Sys.time()), ":\t", "All done!", "\n", sep = ""))
cat("=================================== \n")

