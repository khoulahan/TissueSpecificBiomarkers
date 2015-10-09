### create dotmap of results of glm ##############################################################
generate.glm.dotmap <- function(
	therapeutic.terms,
	biomarkers,
	mutations,
	drug.response,
	map.tissue,
	compounds,
	filename,
	dataset = 'CCLE',
	width = 4,
	height = 8,
	NA.spot.size = 3.5,
	spot.size.function = NULL
	) {

	# generate plot data for each compound specified 
	plot.data <- sapply(
		compounds,
		function(compound) {
			extract.biomarkers.drug.response(
				therapeutic.terms = therapeutic.terms, 
				biomarkers = biomarkers, 
				mutations = mutations, 
				drug.response = drug.response,
				compound = compound, 
				map.tissue = map.tissue,
				dataset = dataset
				);
			},
		simplify = FALSE
		);

	# run glm on each compound 
	glm.results <- lapply(
		plot.data,
		apply.glm.gamma
		);

	# create dotmap plot data
	coefficients.plot.data <- lapply(
		glm.results,
		function(result) {
			summary.glm <- summary(result);
			# extract coefficients
			data.frame(
				Features = rownames(summary.glm$coefficients), 
				Coefficients = summary.glm$coefficients[,'Estimate']
				);
			})
	all.coefficients.plot.data <- suppressWarnings(Reduce(function(...) merge(..., all=T, by = 'Features'), coefficients.plot.data));
	colnames(all.coefficients.plot.data) <- c('Features', compounds);
	rownames(all.coefficients.plot.data) <- all.coefficients.plot.data$Features;
	# remove column of feature names and intercept results
	all.coefficients.plot.data <- all.coefficients.plot.data[-1,-1];

	pvalues.plot.data <- lapply(
		glm.results,
		function(result) {
			summary.glm <- summary(result);
			# extract pvalues
			data.frame(
				Features = rownames(summary.glm$coefficients),
				Compound = summary.glm$coefficients[,'Pr(>|t|)']
				);
			});
	all.pvalues.plot.data <- suppressWarnings(Reduce(function(...) merge(..., all=T, by = 'Features'), pvalues.plot.data));
	colnames(all.pvalues.plot.data) <- c('Features', compounds);
	rownames(all.pvalues.plot.data) <- all.pvalues.plot.data$Features;
	# remove column of feature names and intercept results
	all.pvalues.plot.data <- all.pvalues.plot.data[-1,-1];

	# order data and bg data by feature coefficients
	all.pvalues.plot.data		<- all.pvalues.plot.data[order(rowSums(abs(all.coefficients.plot.data), na.rm=TRUE), decreasing = TRUE),]; 
	all.coefficients.plot.data	<- all.coefficients.plot.data[order(rowSums(abs(all.coefficients.plot.data), na.rm=TRUE), decreasing = TRUE),];

	# set yaxis labels 
	yaxis.lab <- rownames(all.coefficients.plot.data);
	yaxis.lab[grep('tissue',yaxis.lab)] <- unlist(
		sapply(
			yaxis.lab[grep('tissue',yaxis.lab)],
			function(x) {
				tmp <- unlist(strsplit(x, split = '_'));
				substr(tmp[1], start = 7, stop = nchar(tmp[1]));
				}
			)
		);
	# make tissue labels red
	tissue.colours <- rownames(all.coefficients.plot.data);
	tissue.colours[grep('tissue', tissue.colours)] <- 'firebrick2';
	tissue.colours[tissue.colours != 'firebrick2'] <- 'black';

	# adjust spot size and colour for dotmap 
	if (is.null(spot.size.function)) {
		spot.size.function <- function(x) { 0.1 + (2 * abs(x)); }
		}
	spot.colour.function <- function(x) {
			colours <- rep('white', length(x));
			colours[sign(x) == -1] <- 'dodgerblue2';
			colours[sign(x) ==  1] <- 'firebrick2';
			return(colours);
			}
	spot.outline.function <- function(x) {
			colours <- rep('black', length(x));
			colours[x == 0] <- 'transparent';
			return(colours);
			}

	create.dotmap(
		x = all.coefficients.plot.data,
		bg.data = -log10(all.pvalues.plot.data),
		filename = filename,
		spot.colour.function = spot.colour.function,
		spot.size.function = spot.size.function,
		pch = 21,
		xaxis.cex = 1,
		xaxis.rot = 90,
		yaxis.cex = 1,
		yaxis.lab = yaxis.lab,
		yaxis.col = rev(tissue.colours),
		height = height,
		width = width,
		colourkey = TRUE,
		colour.scheme = c('white', 'white', 'black'),
		colour.centering.value = -log10(0.5),
		at = seq(0,5,0.5),
		colourkey.labels.at = seq(0,5,1),
		colourkey.labels = c(
				expression(10^0),
				expression(10^-1),
				expression(10^-2),
				expression(10^-3),
				expression(10^-4),
				substitute(p.value<=10^-5, list(p.value = ''))
				),
		bg.alpha = 1,
		NA.spot.size = NA.spot.size,
		resolution = 500,
		key.top = 0.6,
		bottom.padding = 3
		);
}