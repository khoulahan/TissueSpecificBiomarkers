### create dotmap of results of glm ##############################################################
generate.glm.dotmap <- function(
	therapeutic.terms,
	biomarkers,
	mutations,
	drug.response,
	map.tissue,
	compounds,
	filename,
	key.sizes = seq(-1,1,0.4),
	dataset = 'CCLE',
	width = 4,
	height = 8,
	NA.spot.size = 3.5,
	spot.size.function = NULL,
	interactions.only = FALSE
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

	# run glm on each compound and generate new vector of compound names that ran through glm 
	glm.compounds <- c();
	glm.results <- lapply(
		seq_along(plot.data),
		function(i) {
			# first ensure that compound has at least one variant associated with at least one tissue
			if (length(levels(plot.data[[i]]$tissue)) == 1) {
				print(paste(names(plot.data)[i], "has only one variants associated with", as.character(levels(plot.data[[i]]$tissue))));
				print(paste("Cannot run glm. Skipping", names(plot.data)[i], "..."));
				NULL;
			} else {
				glm.compounds <<- c(glm.compounds, names(plot.data)[i]);
				print(paste("Modeling", names(plot.data)[i], "response ..."));
				apply.glm.gamma(plot.data[[i]])
				}
			}
		);
	# remove the NULL from the list 
	glm.results <- glm.results[!vapply(glm.results, is.null, logical(1))]

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
	colnames(all.coefficients.plot.data) <- c('Features', glm.compounds);
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
	colnames(all.pvalues.plot.data) <- c('Features', glm.compounds);
	rownames(all.pvalues.plot.data) <- all.pvalues.plot.data$Features;
	# remove column of feature names and intercept results
	all.pvalues.plot.data <- all.pvalues.plot.data[-1,-1];

	# order data and bg data by feature coefficients
	all.pvalues.plot.data		<- all.pvalues.plot.data[order(rowSums(abs(all.coefficients.plot.data), na.rm=TRUE), decreasing = TRUE),]; 
	all.coefficients.plot.data	<- all.coefficients.plot.data[order(rowSums(abs(all.coefficients.plot.data), na.rm=TRUE), decreasing = TRUE),];

	# only keep features that are significant (p-value < 0.05) for at least one gene
	features.to.keep 			<- apply(all.pvalues.plot.data, 1, function(x) {
		any(x <= 0.05, na.rm = TRUE)}) 
	all.pvalues.plot.data 		<- all.pvalues.plot.data[features.to.keep,];
	all.coefficients.plot.data 	<- all.coefficients.plot.data[features.to.keep,];

	if (interactions.only) {
		# if interactions only specified only plot interaction terms
		all.coefficients.plot.data 	<- all.coefficients.plot.data[grep(":", rownames(all.coefficients.plot.data)),];
		all.pvalues.plot.data 		<- all.pvalues.plot.data[grep(":", rownames(all.pvalues.plot.data)),];
	}

	# set yaxis labels 
	yaxis.lab <- rownames(all.coefficients.plot.data);
	yaxis.lab[grep('tissue',yaxis.lab)] <- unlist(
		sapply(
			yaxis.lab[grep('tissue',yaxis.lab)],
			function(x) {
				substr(x, start = 7, stop = nchar(x));
				}
			)
		);

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

	# generate key 
	#max.coefficient <- round(max(all.coefficients.plot.data, na.rm = TRUE));
	#key.sizes <- rev(seq(-max.coefficient,max.coefficient,round(max.coefficient/3)));

	# find mode of action of each drug 
	find.mode.of.action <- function(compounds, column) {
		unlist(
			sapply(
				compounds, 
				function(x) {
					tmp <- unique(therapeutic.terms[grep(x, therapeutic.terms[,column], ignore.case = TRUE),'Function']);
					# if more than one mode of action, choose first and tell user which one chosen
					if (length(tmp) > 1) {
						print(paste("Due to multiple modes of action, considering", x, "as", tmp[1]));
						tmp <- tmp[1]
						}
					return(as.character(tmp));
					}
				)
			);
		}

	if (dataset == 'CCLE') {
		mode.of.action <- find.mode.of.action(glm.compounds, 'AssociatedCompoundsCCLE');
	} else if (dataset == 'CTDD') {
		mode.of.action <- find.mode.of.action(glm.compounds, 'AssociatedCompoundsCTDD');
	} else if (dataset == 'Sanger') {
		mode.of.action <- find.mode.of.action(glm.compounds, 'AssociatedCompoundsSanger');
	}

	# add top covariate indicating the mode of action of the drugs
	mode.of.action.colours <- unlist(sapply(mode.of.action, function(x) {
		as.character(unique(therapeutic.terms[therapeutic.terms$Function == x,'AssociatedColour']))
		}));

	# create right hand covariate to specify number of cell lines associated with feature 
	number.cell.lines <- yaxis.lab;
	tissues.lab <- sapply(yaxis.lab, function(x) { sum(grepl(x, rownames(mutations)))});
	interactions <- yaxis.lab[grep(':', yaxis.lab)];
	number.cell.lines[yaxis.lab %in% colnames(mutations)] 	<- colSums(mutations[,number.cell.lines[yaxis.lab %in% colnames(mutations)]])
	number.cell.lines[tissues.lab != 0]						<- tissues.lab[tissues.lab != 0]; 
	number.cell.lines[yaxis.lab %in% interactions]			<- sapply(interactions, function(x) {
		tmp.tissue 	<- unlist(strsplit(x, split = ':'))[1];
		tmp.variant <- unlist(strsplit(x, split = ':'))[2];
		sum(mutations[grep(tmp.tissue, rownames(mutations)),tmp.variant]);
		});
	number.cell.lines[number.cell.lines == 'UNKNOWN'] 		<- 3;
	number.cell.lines										<- as.numeric(number.cell.lines);

	# set colours 
	cell.lines.colours <- number.cell.lines;
	cell.lines.colours[number.cell.lines < 5] 								<- 'lightcyan';
	cell.lines.colours[number.cell.lines >= 5 & number.cell.lines < 10] 	<- 'lightblue2';
	cell.lines.colours[number.cell.lines >= 10 & number.cell.lines < 20] 	<- 'dodgerblue3';
	cell.lines.colours[number.cell.lines >= 20] 							<- 'dodgerblue4';

	# top covariate
	top.covariate <- list(
		rect = list(
			col = 'black',
			fill = mode.of.action.colours,
			lwd = 1.5
			)
		);

	top.cov.grob <- covariates.grob(
		covariates = top.covariate,
		ord = c(1:length(mode.of.action.colours)),
		side = 'top'
		);

	right.covariate <- list(
		rect = list(
			col = 'black',
			fill = cell.lines.colours,
			lwd = 1.5
			)
		);
	right.cov.grob <- covariates.grob(
		covariates = right.covariate,
		ord = c(length(cell.lines.colours):1),
		side = 'right'
		);

	cov.legend <- list(
		legend = list(
			labels = as.character(unique(mode.of.action)),
			colours = unique(mode.of.action.colours),
			title = 'Mode of Action'
			),
		legend = list(
			labels = c('< 5','< 10', '< 20', '> 20'),
			colours = c('lightcyan','lightblue2','dodgerblue3','dodgerblue4'),
			title = 'Number of Cell Lines'
			)
		);

	cov.legend <- legend.grob(
		legends = cov.legend
		);

	create.dotmap(
		x = all.coefficients.plot.data,
		bg.data = -log10(all.pvalues.plot.data),
		filename = filename,
		main = NULL,
		spot.colour.function = spot.colour.function,
		spot.size.function = spot.size.function,
		pch = 21,
		xaxis.cex = 1,
		xaxis.rot = 90,
		yaxis.cex = 1,
		yaxis.lab = yaxis.lab,
		#yaxis.col = rev(tissue.colours),
		height = height,
		width = width,
		colourkey = TRUE,
		key.top = 1,
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
		bottom.padding = 3,
		key = NULL,
		# key = list(
		# 	space = 'right',
		# 	points = list(
		# 		cex = spot.size.function(key.sizes),
		# 		col = spot.colour.function(key.sizes),
		# 		pch = 19
		# 		),
		# 	text = list(
		# 		lab = as.character(key.sizes),
		# 		cex = 1,
		# 		adj = 1
		# 		),
		# 	padding.text = 2,
		# 	background = 'white'
		# 	),
		legend = list(
			top = list(
				fun = top.cov.grob
				),
			right = list(
				fun = right.cov.grob
				),
			left = list( 
				fun = cov.legend
				)
			)
		);
}