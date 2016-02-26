### FUNCTION TO RUN GLM ###########################################################################
apply.glm.gamma <- function(plot.data) {
	# remove any cell lines without drug response
	plot.data <- plot.data[!is.na(plot.data$drug.response),];
	# ensure ordered by variant and then cell line 
	plot.data <- plot.data[order(plot.data$gene, plot.data$cell.line),];

	# convert plot data into sample by feature matrix where features are tissue and binary variants
	all.variants <- as.character(unique(plot.data$gene));
	sample.features <- data.frame(
		matrix(
			nrow = sum(plot.data$gene == all.variants[1]),
			ncol = 2+length(all.variants)
			)
		);
	rownames(sample.features) <- plot.data[plot.data$gene == all.variants[1],'cell.line'];
	colnames(sample.features) <- c('drug.response','tissue',all.variants);

	sample.features$drug.response 	<- plot.data[plot.data$gene == all.variants[1],'drug.response'];
	sample.features$tissue 			<- unlist(
		sapply(
			as.character(plot.data[plot.data$gene == all.variants[1],'cell.line']), 
			function(x) {
				tmp <- unlist(strsplit(x, split = '_'));
				paste(tmp[2:length(tmp)], collapse = '_');
				}
			)
		);
	# change unknown tissues to unknown 
	sample.features[sample.features$tissue %in% c('NA_M059J','NA_SF8657','NA_SNUC2B'),'tissue'] <- 'UNKNOWN';

	for (i in 1:length(all.variants)) {
		cell.lines.with.variants 	<- as.character(
			plot.data[plot.data$gene == all.variants[i] & plot.data$variant == 'Variant','cell.line']
			);
		sample.features[cell.lines.with.variants,all.variants[i]] <- 1;
		}
	sample.features[is.na(sample.features)] <- 0;

	# remove cell lines with drug response equal to 0 
	sample.features <- sample.features[sample.features$drug.response != 0,];

	# test linear model for tissues and variants
	# using Gamma function due to the high skewness of drug response data
	formula.call <- formula(paste0(
		'drug.response~',
		paste(
			paste(colnames(sample.features)[2:ncol(sample.features)], collapse = '+'),
			paste(paste('tissue', colnames(sample.features)[3:ncol(sample.features)], sep = '*'), collapse = '+'),
			sep = '+'
			)
		));
	glm.results <- glm(formula.call, sample.features, family = Gamma(link='log'));
	return(glm.results);
	}