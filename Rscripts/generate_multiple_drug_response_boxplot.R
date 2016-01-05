### generate_multiple_drug_response_boxplot.R #####################################################
# generate boxplots showing tissue vs non tissue for all three datasets 
###################################################################################################
library(BoutrosLab.plotting.general)
library(getopt)

date <- Sys.Date();

### OBTAIN COMMAND LINE ARGUMENTS #################################################################
spec = matrix(c(
	# required parameters
	'literature',			'l',	1,	"character",
	'variant',				'v',	1,	"character",
	'filename',				'f',	1,	"character",
	'mutations.ccle',		'm',	2,	"character",
	'drug.response.ccle',	'd',	2,	"character",
	'drug.response.ctdd',	'e',	2,	"character",
	'mutations.cgp',		'p',	2,	"character",
	'drug.response.cgp',	'r',	2,	"character",
	'additional.variants',	'a',	2,	"character"
	), byrow = TRUE, ncol = 4);
opt = getopt(spec);,

if (is.null(opt$width.subset))	{opt$width.subset = 8}
if (is.null(opt$width.all)) {opt$width.all = 8}
### READ IN MAPPING FILES #########################################################################
# read in mapping files: tissue and therapeutic terms and drug names
map.tissue <- read.delim(
	'../MappingFiles/MapDiseaseToTissue.txt',
	sep = '\t',
	header = TRUE
	);
therapeutic.terms <- read.delim(
	'../MappingFiles/TherapeuticTermsMappedDrugs.txt',
	sep = '\t',
	header = TRUE
	);
map.drug.names <- read.delim(
	'../MappingFiles/MapDrugNames.tab',
	sep = '\t',
	header = TRUE,
	stringsAsFactors = FALSE
	);
### READ IN LITERATURE BIOMARKERS FILE ############################################################
# read in biomarkers files
biomarkers <- read.delim(
	opt$literature,
	sep = '\t',
	header = TRUE
	);
### NORMALIZE DRUG RESPONSE #######################################################################
# normalize drug response to be between 0 and 1
normalize.range <- function(x) {
	max.x <- max(x, na.rm = TRUE);
	min.x <- min(x, na.rm = TRUE);
	(x - min.x)/(max.x - min.x);
	}
## CCLE DATASET ###################################################################################
# generate vector of datasets to test 
datasets <- c();
if (!is.null(opt$mutations.ccle)) {
	# read in mutation profile
	mutations.ccle <- read.delim(
		opt$mutations.ccle,
		sep = '\t',
		header = TRUE
		);
}
if (!is.null(opt$drug.response.ccle)) {
	# read in drug response
	drug.response.ccle <- read.delim(
		opt$drug.response.ccle,
		sep = '\t',
		header = TRUE
		);
	# normalize drug response
	drug.response.ccle 	<- t(apply(drug.response.ccle, 1, normalize.range));
	# because ccle is the area above the curve calculate 1 - the value
	drug.response.ccle <- 1 - drug.response.ccle;
	# add CCLE to datasets
	datasets <- c(datasets, 'CCLE');
}

## CGP DATASET ####################################################################################
if (!is.null(opt$mutations.cgp) & !is.null(opt$drug.response.cgp)) {
	# read in mutation profile
	mutations.cgp <- read.delim(
		opt$mutations.cgp,
		sep = '\t',
		header = TRUE
		);
	# read in drug response
	drug.response.cgp <- read.delim(
		opt$drug.response.cgp,
		sep = '\t',
		header = TRUE
		);
	drug.response.cgp <- drug.response.cgp[!apply(drug.response.cgp, 1, function(x) {all(is.na(x))}),];
	# normalize drug response
	drug.response.cgp 	<- t(apply(drug.response.cgp, 1, normalize.range));
	# add CGP to datasets
	datasets <- c(datasets, 'CGP');
} 

## CTDD DATASET ###################################################################################
# same mutations as ccle
if (!is.null(opt$drug.response.ctdd)) {
	# read in drug response
	drug.response.ctdd <- read.delim(
		opt$drug.response.ctdd,
		sep = '\t',
		header = TRUE
		);
	# remove any cell lines that all have NA
	drug.response.ctdd <- drug.response.ctdd[!apply(drug.response.ctdd, 1, function(x) {all(is.na(x))}),];
	# normalize drug response
	drug.response.ctdd 	<- t(apply(drug.response.ctdd, 1, normalize.range));
	# add CCLE to datasets 
	datasets <- c(datasets,'CCLE');
}

### EXTRACT VARIANTS ##############################################################################
grid <- expand.grid(
	c('AZD6244','Erlotinib','PD.0325901','PD.0332991','PLX4720','RAF265','L.685458',
		'Lapatinib','Nutlin.3','Sorafenib','X17.AAG','Paclitaxel','Panobinostat','TKI258','Irinotecan'),
	datasets)
# extract variant information
source("../R/extract_variant_multiple_datasets.R");
plot.data <- extract.variant.multiple.datasets(grid, additional.variant = opt$additional.variant);

## BOXPLOT ########################################################################################
# generate boxplot
source("../R/create_multiple_drug_response_boxplot.R");
create.multiple.drug.response.boxplot(
	plot.data[plot.data$gene == opt$variant,],
	file = opt$filename,
	xaxis.labels = unique(plot.data[plot.data$gene == opt$variant,'compound'])[order(
		as.character(unique(plot.data[plot.data$gene == opt$variant,'compound'])))],
	variant = opt$variant,
	width = 10
	)