### test_variant_drug_combinations.R ##############################################################
library(rjags);
library(BEST);
library(doParallel);
library(getopt);

date <- Sys.Date();

### OBTAIN COMMAND LINE ARGUMENTS #################################################################
spec = matrix(c(
	# required parameters
	'gene',				'g',	1,	"character",
	'mutations',		'm',	1,	"character",
	'drug.response',	'd',	1,	"character",
	'output.dir',		'o',	2,	"character",
	'parallel',			'p',	2,	"logicial",
	'power',			'w',	2,	"logical",
	'rep',				'r',	2,	"integer",
	'cores',			'c',	2,	"integer"	
	), byrow = TRUE, ncol = 4);
opt = getopt(spec);

if (is.null(opt$parallel)) {opt$parallel = TRUE}
if (is.null(opt$power)) {opt$power = FALSE}
if (is.null(opt$rep)) {opt$rep = 200}
if (is.null(opt$cores)) {opt$cores = 3}

### READ IN DATA ##################################################################################
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

#### TEST ASSOCIATIONS ############################################################################
# set parallel background if parallel is set to true
if (opt$parallel) {
	cl <- makeCluster(opt$cores);
	registerDoParallel(cl, cores = opt$cores);
	}
# create function to be run in parallel
calculate.drug.response.differences <- function(compound, variant, mutations, drug.response, repetitions = 200, power = FALSE) {

	print(paste("Comparing response to", compound, "..."));

	# find drug response of cell lines with variant and without
	drug.response.variant <- drug.response[rownames(drug.response) %in% rownames(mutations)[which(
		mutations[,variant] == 1)],compound];
	drug.response.no.variant <- drug.response[rownames(drug.response) %in% rownames(mutations)[which(
		mutations[,variant] == 0)],compound];

	# ensure no NAs in drug response
	drug.response.variant <- drug.response.variant[!is.na(drug.response.variant)];
	drug.response.no.variant <- drug.response.no.variant[!is.na(drug.response.no.variant)];

	# ensure there are more than 5 cell lines with mutation 
	if (length(drug.response.variant) < 5 | length(drug.response.no.variant) < 5) {
		print("Not enough cell lines with/without variant to test ...")
		return(NULL);
		}

	# bayesian estimation of posterior distributions
	BESTout 		<- BESTmcmc(drug.response.variant, drug.response.no.variant);
	summary.BESTout <- summary(BESTout);

	# Wilcox rank sum test
	# not sure the direction therefore testing two sided
	p.values <- wilcox.test(
		drug.response.variant,
		drug.response.no.variant,
		alternative = 'two.sided'
		)$p.value;

	if (power) {
		# calculate power of difference in median being outside ROPE of (-0.01, 0.01)
		BESTPower <- BESTpower(
			BESTout,
			N1 = length(drug.response.variant),
			N2 = length(drug.response.no.variant),
			ROPEm = c(-0.05,0.05),
			nRep = repetitions
			);
		# return data frame with power analysis results
		data.frame(
			compound = compound,
			variant = variant,
			p.value = round(p.values, digits=4),
			median.diff = round(summary.BESTout['muDiff','median'], digits = 2),
			HDIlo = round(summary.BESTout['muDiff','HDIlo'], digits = 2),
			HDIup = round(summary.BESTout['muDiff','HDIup'], digits = 2),
			effect.size = round(summary.BESTout['effSz','median'], digits = 2),
			number.cell.lines = length(drug.response.variant),
			HDI.greater.ROPE = BESTPower[1,'mean'],
			HDI.less.ROPE = BESTPower[2,'mean'],
			HDI.in.ROPE = BESTPower[2,'mean']
			);
	} else {
		# return data frame of hypothesis testing results
		data.frame(
			compound = compound,
			variant = variant,
			p.value = round(p.values, digits=4),
			median.diff = round(summary.BESTout['muDiff','median'], digits = 2),
			HDIlo = round(summary.BESTout['muDiff','HDIlo'], digits = 2),
			HDIup = round(summary.BESTout['muDiff','HDIup'], digits = 2),
			effect.size = round(summary.BESTout['effSz','median'], digits = 2),
			number.cell.lines = length(drug.response.variant)
			);
		}
	}

variant.drug.response <- foreach(compound = colnames(drug.response), .combine='rbind', .packages = 'BEST') %dopar%
	calculate.drug.response.differences(
		compound = compound,
		variant = opt$gene,
		mutations = mutations,
		drug.response = drug.response,
		repetitions = opt$rep,
		power = opt$power
		);

write.table(
	x =variant.drug.response,
	file = paste0(
		opt$output.dir,
		date,
		'_',
		opt$gene,
		'VariantDrugCombinations.txt'
		),
	sep = '\t',
	row.names = FALSE,
	quote = FALSE
	);