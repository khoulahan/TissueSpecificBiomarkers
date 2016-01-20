### test_neighbour_variant_drug_combinations.R ####################################################
suppressWarnings(suppressMessages(library(rjags)));
suppressWarnings(suppressMessages(library(BEST)));
suppressWarnings(suppressMessages(library(doParallel)));
suppressWarnings(suppressMessages(library(getopt)));
suppressWarnings(suppressMessages(library(org.Hs.eg.db)));
suppressWarnings(suppressMessages(library(annotate)));

date <- Sys.Date();

### OBTAIN COMMAND LINE ARGUMENTS #################################################################
spec = matrix(c(
	# required parameters
	'gene',				'g',	1,	"character",
	'mutations',		'm',	1,	"character",
	'drug.response',	'd',	1,	"character",
	'network',			'n',	1,	"character",
	'compound',			'c',	2,	"character",
	'output.dir',		'o',	2,	"character",
	'bayesian',			'b',	2,	"logical",
	'parallel',			'p',	2,	"logicial",
	'power',			'w',	2,	"logical",
	'rep',				'r',	2,	"integer",
	'cores',			's',	2,	"integer"	
	), byrow = TRUE, ncol = 4);
opt = getopt(spec);,

if (is.null(opt$compound)) {opt$compound = 'all'}
if (is.null(opt$bayesian)) {opt$bayesian = FALSE}
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
network <- read.delim(
	opt$network,
	sep = '\t',
	header = FALSE
	);

#### TEST ASSOCIATIONS ############################################################################
# set parallel background if parallel is set to true
if (opt$parallel) {
	cl <- makeCluster(opt$cores);
	registerDoParallel(cl, cores = opt$cores);
	}
# create function to be run in parallel
calculate.drug.response.differences <- function(compound, variants, mutations, drug.response, 
	bayesian = FALSE, repetitions = 200, power = FALSE) {

		cat(paste("Comparing response to", compound, "..."));

		# check that all variants are in mutation matrix
		in.mutations <- variants %in% colnames(mutations);
		if (any(!in.mutations)) {
			print(paste("Variant(s) not in mutation matrix:", variants[!in.mutations], collapse = ' '));
			print("Skipping variant and continuing with analysis ...")
			variants <- variants[in.mutations];
		}

		# find drug response of cell lines with variant and without
		if (length(variants) == 1) {
			is.variant <- mutations[,variants] == 1;
			names(is.variant) <- rownames(mutations);
		} else {
			is.variant <- apply(mutations[,variants], 1, function(x) {any(x == 1)});
		}

		drug.response.variant <- drug.response[rownames(drug.response) %in% names(which(is.variant)),compound];
		drug.response.no.variant <- drug.response[rownames(drug.response) %in% names(which(!is.variant)),compound];

		# ensure no NAs in drug response
		drug.response.variant <- drug.response.variant[!is.na(drug.response.variant)];
		drug.response.no.variant <- drug.response.no.variant[!is.na(drug.response.no.variant)];

		# ensure there are more than 5 cell lines with mutation 
		if (length(drug.response.variant) < 5 | length(drug.response.no.variant) < 5) {
			print("Not enough cell lines with/without variant to test ...")
			return(NULL);
			}

		# Wilcox rank sum test
		# not sure the direction therefore testing two sided
		p.values <- wilcox.test(
			drug.response.variant,
			drug.response.no.variant,
			alternative = 'two.sided'
			)$p.value;

		# set results dataframe
		results <- data.frame(
			compound = compound,
			variant = variants[1],
			p.value = round(p.values, digits=4),
			median.diff = NA,
			HDIlo = NA,
			HDIup = NA,
			effect.size = NA,
			number.cell.lines = length(drug.response.variant)
			);

		if (!bayesian) {
			# return wilcox rank sum test results
			return(results);

		} else {
			# bayesian estimation of posterior distributions
			BESTout 		<- BESTmcmc(drug.response.variant, drug.response.no.variant);
			summary.BESTout <- summary(BESTout);

			# add bayesian estimation results
			results$median.diff <- round(summary.BESTout['muDiff','median'], digits = 2);
			results$HDIlo 		<- round(summary.BESTout['muDiff','HDIlo'], digits = 2);
			results$HDIup 		<- round(summary.BESTout['muDiff','HDIup'], digits = 2);
			results$effect.size <- round(summary.BESTout['effSz','median'], digits = 2);

			if (!power) {
				# return bayesian estimation results
				return(results);
			} else {
				# calculate power of difference in median being outside ROPE of (-0.01, 0.01)
				BESTPower <- BESTpower(
					BESTout,
					N1 = length(drug.response.variant),
					N2 = length(drug.response.no.variant),
					ROPEm = c(-0.05,0.05),
					nRep = repetitions
					);
				# return data frame with power analysis results
				results <- data.frame(
					results,
					HDI.greater.ROPE = BESTPower[1,'mean'],
					HDI.less.ROPE = BESTPower[2,'mean'],
					HDI.in.ROPE = BESTPower[2,'mean']
					);
				return(results);
				}
			} 
		}

### FIND NEIGHBOURS IN NETWORK ####################################################################
# convert entrez IDs to Hugo symbols
network[,1] <- getSYMBOL(as.character(network[,1]), data='org.Hs.eg');
network[,2] <- getSYMBOL(as.character(network[,2]), data='org.Hs.eg');
# remove any genes that do not have a Hugo symbol 
network <- network[!is.na(network[,1]) & !is.na(network[,2]),];

neighbour.variants <- c(network[network[,1] == opt$gene, 2], network[network[,2] == opt$gene,1]);
# only keep variants present in mutation matrix
neighbour.variants <- neighbour.variants[neighbour.variants %in% colnames(mutations)];

# set compounds to test
if (opt$compound == 'all') {
	compounds <- colnames(drug.response)
} else {
	compounds <- opt$compound;
}

variant.drug.response <- foreach(compound = compounds, .combine='rbind', .packages = 'BEST') %dopar%
	calculate.drug.response.differences(
		compound = compound,
		variants = c(opt$gene, neighbour.variants),
		mutations = mutations,
		drug.response = drug.response,
		bayesian = opt$bayesian,
		repetitions = opt$rep,
		power = opt$power
		);

write.table(
	x = variant.drug.response,
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