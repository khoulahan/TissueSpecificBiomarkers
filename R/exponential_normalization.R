### exponential_normalization.R ###################################################################
# normalizes gene expression data based on underlying exponential distribution 
###################################################################################################
exponential.normalization <- function(gene.expression) {
	# sample-wise normalization 
	normalized.gene.expression <- apply(
		gene.expression,
		1,
		function(sample) {
			sample <- rank(sample);
			sample <- max(sample) - sample;
			sample <- (sample-min(sample))/(max(sample)-min(sample));
			dexp(sample);
			}
		);
	return(t(normalized.gene.expression));
}