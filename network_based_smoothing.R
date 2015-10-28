#### network-based smoothing ######################################################################
# implementation of network based smoothing in R
###################################################################################################
library(org.Hs.eg.db)
library(annotate)
# read in human net network 
humannet <- read.delim(
	'~/Documents/Pathway/HumanNet.v1.benchmark.txt', 
	sep = '\t', 
	header = FALSE
	);

# convert Entrez gene ID to Hugo symbol
humannet[,1] <- getSYMBOL(as.character(humannet[,1]), data='org.Hs.eg');
humannet[,2] <- getSYMBOL(as.character(humannet[,2]), data='org.Hs.eg');

# convert to an adjacency matrix 
humannet.adjacency <- as.matrix(table(humannet[,1], humannet[,2]));

# convert adjacency matrix to transition matrix of probabilites
humannet.transition <- humannet.adjacency/rowSums(humannet.adjacency);

# read in mutations 
mutations <- read.delim(
	'~/Documents/Mutations/CCLE_mutation_featurematrix.tab',
	sep = '\t',
	header = TRUE
	);
# convert to matrix
mutations <- as.matrix(mutations);

# ensure genes in rows equal genes in columns of network
humannet.transition <- humannet.transition[rownames(humannet.transition) %in% colnames(humannet.transition), 
	colnames(humannet.transition) %in% rownames(humannet.transition)];
# ensure rows (genes) of network are the same as genes in mutations matrix
humannet.transition <- humannet.transition[rownames(humannet.transition) %in% colnames(mutations),];
# keep only genes in network in mutation matrix
mutations 			<- mutations[,colnames(mutations) %in% rownames(humannet.transition)];

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
	file = '~/Documents/Mutations/2015-10-28_DiffusedMutationsCCLE.txt',
	sep = '\t',
	quote = FALSE
	);

