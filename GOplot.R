### GO Plot for DEGS 
#Installation of the latest released version
install.packages('GOplot')
# Installation of the latest development version
#install_github('wencke/wencke.github.io')

# load GO
library(GOplot)

# read degs file and handling of it to be appopriate to be used for DAVID database
exp1.degs = read.table(file, header = FALSE, sep = "", dec = ".")
AveExpr <- rowMeans(exp1.degs)
res=cbind(res.degs, AveExpr)
res=cbind(ID=rownames(res),res)

cols = c("ID", "AveExpr", "logFC", "P.Value", "adj_pval")
colnames(res)=cols
res = as.data.frame(res[-6])
head(res)

# Load the dataset
david <- read.delim(
  "chart_5332916E65FB1685717161410.txt",
  sep="\t", header=TRUE)

genelist <- res
# Get a glimpse of the data format of the results of the functional analysis... 
tail(david)

# handling of david table
david$ID = unlist(sapply(strsplit(david$Term, "~"),c)[1,])
david$Term = unlist(sapply(strsplit(david$Term, "~"),c)[2,])
david["adj_pval"] = david["Benjamini"]
david$Category = unlist(sapply(strsplit(david$Category, "_"),c)[2,])
head(david)

david = david[, c("Category", "ID", "Term", "Genes", "adj_pval")]
head(david)

# Generate the plotting object
circ <- circle_dat(david, genelist)

# Generate a simple barplot for BP category
GOBar(subset(circ, category == 'BP'))

# Facet the barplot according to the categories of the terms 
GOBar(circ, display = 'multiple')

# Facet the barplot, add a title and change the colour scale for the z-score
GOBar(circ, display = 'multiple', title = 'Z-score coloured barplot', zsc.col = c('yellow', 'black', 'cyan'))

# Generate the bubble plot with a label threshold of 3
GOBubble(circ, labels = 3)

# Add a title, change the colour of the circles, facet the plot according to the categories and change the label threshold
GOBubble(circ, title = 'Bubble plot', colour = c('orange', 'darkred', 'gold'), display = 'multiple', labels = 3)

# Colour the background according to the category
GOBubble(circ, title = 'Bubble plot with background colour', display = 'multiple', bg.col = T, labels = 3)

# Reduce redundant terms with a gene overlap >= 0.75...
reduced_circ <- reduce_overlap(circ, overlap = 0.75)
# ...and plot it
GOBubble(reduced_circ, labels = 2.8)

#Generate a circular visualization of the results of gene- annotation enrichment analysis
GOCircle(circ)


# Generate a circular visualization of selected terms
IDs <- c('GO:0007507', 'GO:0001568', 'GO:0001944', 'GO:0048729', 'GO:0048514', 'GO:0005886', 'GO:0008092', 'GO:0008047')
GOCircle(circ, nsub = IDs)

# Generate a circular visualization for 10 terms
GOCircle(circ, nsub = 10)

c(1:100)
dim(genes)[1]

# Define a list of genes and sample contains the data frame of selected genes and their logFC.
genes <- genelist[, c("ID", "logFC")]
rownames(genes)= c(1:dim(genes)[1])
process <- david$Term[1:10]
head(genes)
head(process)


# generate the binary matrix
chord <- chord_dat(circ, genes[1:100,], david$Term[1:10])

# Generate the matrix with a list of selected genes
chord <- chord_dat(data = circ, genes = genes[1:100,])

# Generate the matrix with selected processes
chord <- chord_dat(data = circ, process = david$Term[1:10])

# Create chordplot
chord <- chord_dat(data = circ, genes = genes[1:100,], process = david$Term[1:10])
head(chord)
GOChord(chord, space = 0.02, gene.order = 'logFC', gene.space = 0.25, gene.size = 5)

# Display only genes which are assigned to at least three processes
GOChord(chord, limit = c(0, 0), gene.order = 'logFC')

# First, we use the chord object without logFC column to create the heatmap
GOHeat(chord[,-11], nlfc = 0)

# Now we create the heatmap with logFC values and user-defined colour scale
GOHeat(chord, nlfc = 1, fill.col = c('red', 'yellow', 'green'))

GOCluster(circ,david$Term[1:10] , clust.by = 'logFC', term.width = 2)


GOCluster(circ, david$Term[1:10], clust.by = 'term', lfc.col = c('darkgoldenrod1', 'black', 'cyan1'))

# Venndiagram
l1 <- subset(circ, term == 'cell-cell adhesion', c(genes,logFC))
l2 <- subset(circ, term == 'cell division', c(genes,logFC))
l3 <- subset(circ, term == 'angiogenesis', c(genes,logFC))
GOVenn(l1,l2,l3, label = c('cell-cell adhesion', 'cell division', 'angiogenesis'))






