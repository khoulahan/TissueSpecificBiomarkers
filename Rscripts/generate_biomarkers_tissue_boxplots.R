### generate_biomarkers_tissue_boxplots.R #########################################################
# generates boxplots comparing drug response differentiation between samples with specified variant
# and without in both associated (as definied by literature) and tissues
###################################################################################################
library(BoutrosLab.plotting.general)
library(BEST)
library(getopt)

date <- Sys.Date();

### OBTAIN COMMAND LINE ARGUMENTS #################################################################
spec = matrix(c(
	# required parameters
	'compound',			'c',	1,	"character",
	'mutations',		'm',	1,	"character",
	'drug.response',	'd',	1,	"character",
	'therapeutic.terms','t',	1,	"character",
	'map.tissue',		'p',	1,	"character",
	'literature',		'l',	1,	"character",
	'output.dir',		'o',	2,	"character",
	'width.subset',		's',	2,	"integer",
	'width.all',		'a',	2,	"integer"
	), byrow = TRUE, ncol = 4);
opt = getopt(spec);,

if (is.null(opt$width.subset))	{opt$width.subset = 8}
if (is.null(opt$width.all)) {opt$width.all = 8}

### READ IN DATA ##################################################################################
# read in literature data 
literature <- read.delim(
	opt$literature,
	sep = '\t',
	header = TRUE,
	strip.white = TRUE
	);
# read in file that maps therapeutic terms in biomarkers file to drugs in dataset
therapeutic.terms <- read.delim(
	opt$therapeutic.terms,
	sep = '\t',
	header = TRUE
	);
# read in map tissue 
map.tissue <- read.delim(
	opt$map.tissue,
	sep = '\t',
	header = TRUE
	);

### DATASET #######################################################################################
drug.response <- read.delim(
	opt$drug.response,
	sep = '\t',
	header = TRUE
	);
mutations <- read.delim(
	opt$mutations,
	sep = '\t',
	header = TRUE
	);

### GENERATE BOXPLOTS #############################################################################
source("../R/extract_biomarkers_drug_reponse.R")
plot.data <- extract.biomarkers.drug.response(
	therapeutic.terms = therapeutic.terms, 
	biomarkers = biomarkers, 
	mutations = mutations, 
	drug.response = drug.response,
	compound = opt$compound, 
	map.tissue = map.tissue
	);
source("../R/generate_associated_tissue_boxplot.R")
tissue.subset.hypothesis.results <- generate.associated.tissue.boxplot(
	plot.data = plot.data,
	output.dir = opt$output.dir,
	width = opt$width.subset
	);
source("../R/generate_associated_vs_nonassocaited_boxplot.R")
aggregate.tissue.hypothesis.results <- generate.associated.vs.nonassociated.boxplot(
	plot.data = plot.data,
	output.dir = opt$output.dir,
	width = opt$width.all
	);

### WRITE TO FILE #################################################################################
write.table(
	x = tissue.subset.hypothesis.results,
	file = paste0(
		opt$output.dir,
		date,
		"_",
		opt$compound,
		"_TissueSubsetHypothesisTesting.txt"
		),
	sep = '\t',
	quote = FALSE
	);
write.table(
	x = aggregate.tissue.hypothesis.results,
	file = paste0(
		opt$output.dir,
		date,
		"_",
		opt$compound,
		"_AggregateTissueHypothesisTesting.txt"
		),
	sep = '\t',
	quote = FALSE
	);