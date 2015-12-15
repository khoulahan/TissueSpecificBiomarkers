### create_PCA_scatterplot.R ######################################################################
# generates scatterplot of PC1 versus PC2 for inputted expression data 
###################################################################################################
library(BoutrosLab.plotting.general);
library(plyr);
library(getopt);

date <- Sys.Date();

### OBTAIN COMMAND LINE ARGUMENTS #################################################################
spec = matrix(c(
	# required parameters
	'expression',	'e',	1,	"character",
	'filename',		'f',	1,	"character",
	'groups',		'g',	1,	"character"
	), byrow = TRUE, ncol = 4);
opt = getopt(spec);,

###################################################################################################
### READ IN DATA ##################################################################################
# read in gene expression 
gene.expression <- read.delim(
	opt$expression,
	sep = '\t',
	header = TRUE
	);

### FORMAT AND NORMALIZE ##########################################################################
# remove any genes with NA expression 
gene.expression <- t(na.omit(t(gene.expression)))
# exponential normalization 
source("../R/exponential_normalization.R");
normalized.expression <- exponential.normalization(gene.expression);
# keep only the most variable genes
magnitude <- 10^(nchar(as.character(nrow(normalized.expression)))-1);
number.features <- round_any(nrow(normalized.expression), magnitude, f = floor);
normalized.expression <- normalized.expression[,order(apply(normalized.expression, 2, var), decreasing = TRUE)][,1:number.features];

### PCA ###########################################################################################
# run pca
pca <- princomp(~., data = as.data.frame(normalized.expression));

### PLOT ##########################################################################################
# format plot data
plot.data <- data.frame(
	Comp.1 = pca$scores[,'Comp.1'],
	Comp.2 = pca$scores[,'Comp.2'],
	Sample = rownames(pca$scores)
	);
plot.data$Group <- rep(NA, nrow(plot.data));

# add group
groups <- unlist(strsplit(opt$groups, split = ','));
if (!all(unlist(sapply(groups, function(x) { grep(x, plot.data$Sample)})))) {
	stop("Please append dataset name to sample IDs ...")
} else {
	for (group in groups) {
		plot.data[grep(group,plot.data$Sample),'Group'] <- group;
		}
	}

# create scatterplot
create.scatterplot(
	Comp.1 ~ Comp.2,
	data = plot.data,
	filename = opt$filename,
	groups = plot.data$Group,
	col = default.colours(length(groups)),
	ylab.label = 'Principal Component 1',
	xlab.label = 'Principal Component 2',
	xaxis.cex = 1.5,
	yaxis.cex = 1.5,
	xlab.cex = 1.5,
	ylab.cex = 1.5,
	key = list(
		text = list(
			lab = groups[order(groups)],
			cex = 1,
			col = 'black'
			),
		points = list(
			pch = 19,
			col = default.colours(length(groups)),
			cex = 1
			),
		x = 0.04,
		y = 0.5,
		padding.text = 2
		)
	);