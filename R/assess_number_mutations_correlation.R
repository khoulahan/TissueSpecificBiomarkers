### assess_number_mutations_correlation.R #########################################################
# Generates scatterplot to assess the correlation between number of cell lines mutated and p-value 
# for each variant for the specified compound. Tests for significant differences in mutational 
# numbers using Wilcox Rank Sum Test - significant difference in number of cell lines mutated in 
# significant variants versus non-significant variants
## Arguments:
## mutations			binary matrix of variants mutated in each cell line
## drug.response		matrix of drug response for each compound tested on each cell line
## hypothesis.results	results of systematic assessment of variant/drug combinations
## compound				compound to analysis
## filename				path to write scatterplot to
###################################################################################################
assess.number.mutations.correlation <- function(mutations, drug.response, hypothesis.results, compound, filename) {

	# filter hypothesis results down to only results for compound specified
	hypothesis.results <- hypothesis.results[hypothesis.results$compound == compound,];
	# filter drug response and mutations to only cell lines with drug response for compound specified 
	drug.response <- drug.response[!is.na(drug.response[,compound]),];
	mutations <- mutations[which(rownames(mutations) %in% rownames(drug.response)),];

	# find number of cell lines mutated for each variant in hypothesis results 
	mutations <- mutations[,which(colnames(mutations) %in% hypothesis.results$variant)];
	mutations <- mutations[,order(colnames(mutations))];
	hypothesis.results <- hypothesis.results[order(hypothesis.results$variant),];

	# error check
	if (!all(colnames(mutations) == hypothesis.results$variant)) {
		stop("Check that all variants listed in hypothesis results file are present in mutations matrix ...");
	}

	# generate plot data
	plot.data <- data.frame(
		variant = hypothesis.results$variant,
		p.value = hypothesis.results$p.value,
		num.mut = colSums(mutations),
		sig = NA
		);
	plot.data[plot.data$p.value < 0.05,'sig']	<- TRUE;
	plot.data[plot.data$p.value >= 0.05, 'sig'] <- FALSE;

	# generate scatterplot 
	create.scatterplot(
		p.value ~ num.mut,
		data = plot.data,
		filename = filename,
		main = compound,
		main.cex = 2,
		xlab.label = 'Number of Cell Lines Mutated',
		xlab.cex = 1.5,
		xaxis.cex = 1.5,
		ylab.label = 'P-value',
		ylab.cex = 1.5,
		yaxis.cex = 1.5,
		ylimits = c(0,1),
		yat = seq(0,1,0.2),
		yaxis.lab = c('0.0','0.2','0.4','0.6','0.8','1.0'),
		key = BoutrosLab.plotting.general::get.corr.key(
			plot.data$p.value,
			plot.data$num.mut,
			label.items = 'spearman',
			key.cex = 1.5,
			x.pos = 0.97,
			y.pos = 0.97
			)
		)

	# wilcox test 
	wilcox <- wilcox.test(
		plot.data[plot.data$sig == TRUE,'num.mut'],
		plot.data[plot.data$sig == FALSE,'num.mut'],
		alternative = 'greater'
		)
	print("Wilcox Rank Sum Test testing for greater number of mutated cell lines in significant variants.")
	print(paste("P-value:", wilcox$p.value))

	return(plot.data);
}