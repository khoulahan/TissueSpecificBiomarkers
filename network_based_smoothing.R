#### network-based smoothing ######################################################################
# implementation of network based smoothing in R
###################################################################################################
suppressWarnings(suppressMessages(library(org.Hs.eg.db)));
suppressWarnings(suppressMessages(library(annotate)));

# read in human net network 
humannet <- read.delim(
	'~/Documents/Pathway/HumanNet.v1.top10.txt', 
	sep = '\t', 
	header = FALSE
	);

# convert Entrez gene ID to Hugo symbol
humannet[,1] <- getSYMBOL(as.character(humannet[,1]), data='org.Hs.eg');
humannet[,2] <- getSYMBOL(as.character(humannet[,2]), data='org.Hs.eg');

# convert to an adjacency matrix 
humannet.adjacency <- as.matrix(table(humannet[,1], humannet[,2]));

# read in mutations 
mutations <- read.delim(
	'~/Documents/Mutations/CCLE_mutation_featurematrix.tab',
	sep = '\t',
	header = TRUE
	);
# convert to matrix
mutations <- as.matrix(mutations);

# ensure rows (genes) of network are the same as genes in mutations matrix
humannet.adjacency <- humannet.adjacency[rownames(humannet.adjacency) %in% colnames(mutations),];
# ensure genes in rows equal genes in columns of network
humannet.adjacency <- humannet.adjacency[rownames(humannet.adjacency) %in% colnames(humannet.adjacency), 
	colnames(humannet.adjacency) %in% rownames(humannet.adjacency)];
# keep only genes in network in mutation matrix
mutations <- mutations[,colnames(mutations) %in% rownames(humannet.adjacency)];

# convert adjacency matrix to transition matrix of probabilites
humannet.transition <- humannet.adjacency/rowSums(humannet.adjacency);
humannet.transition[is.na(humannet.transition)] <- 0;

#### NETWORK PROPOGATION ##########################################################################
# initialize mat.norm at 1
mat.norm <- 1;
# initalize F matrix
F.mat <- mutations;
F.mat.previous <- mutations;
# set alpha 
alpha <- 0.7;
while (mat.norm > 1e10^-6) {
	F.mat <- alpha*F.mat%*%humannet.transition+(1-alpha)*mutations;
	# calculate normal of difference between updated matrix and previous matrix
	mat.norm <- norm(F.mat-F.mat.previous);
	# update previous matrix
	F.mat.previous <- F.mat;
}

# write to file
write.table(
	x = F.mat,
	file = '~/Documents/Mutations/2015-10-28_DiffusedMutationsCCLETop10Humannet.txt',
	sep = '\t',
	quote = FALSE
	);

