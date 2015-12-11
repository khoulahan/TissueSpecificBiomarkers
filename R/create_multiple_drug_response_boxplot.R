### create_multiple_drug_response_boxplot.R #######################################################
# creates boxplot subsetted into associated tissue and non-associated tissue for all compounds 
# present in the plot data matrix
## Arguments:
## plot.data		plot data matrix as outputted by extract_multiple_drug_response_boxplot.R
## filename			name of file to write plot to
## xaxis.labels		labels of xaxis - drug names in alphabetical order
## width			width of plot
###################################################################################################
create.multiple.drug.response.boxplot <- function(plot.data, filename, xaxis.labels, width) {
	# create group for plotting by combining tissue and variant information 
	plot.data$group <- apply(
		plot.data,
		1,
		function(x) {
			if (x['tissue'] != 'OTHER') {
				paste(x['compound'],'Associated',x['variant'],sep='_');
			} else {
				paste(x['compound'],'NonAssociated',x['variant'],sep='_');
				}
			}
		);

	# check that all groups exist
	all.groups <- expand.grid(unique(plot.data$compound), c('Associated','NonAssociated'), unique(plot.data$variant));
	all.groups <- apply(all.groups, 1, paste, collapse = '_');

	if (!all(all.groups %in% plot.data$group)) {
		missing.groups <- all.groups[which(!all.groups %in% plot.data$group)]
		for (missing.group in missing.groups) {
			# add tmp row with value above plot limits 
			tmp <- data.frame(matrix(nrow = 1, ncol = ncol(plot.data), dimnames = list(NULL, colnames(plot.data))))
			tmp$group <- missing.group;
			plot.data <- rbind(plot.data, tmp);
		}
	}
	# order plot data by groups
	plot.data <- plot.data[order(plot.data$group),];

	# tissue covariates colours
	tissue.covariates <- unique(plot.data$group);
	tissue.covariates[grep("NonAssociated", tissue.covariates)] <- default.colours(2, palette.type = 'spiral.sunrise')[1];
	tissue.covariates[grep("Associated", tissue.covariates)] 	<- default.colours(2, palette.type = "spiral.sunrise")[2];
	# variant covariates colours
	variant.covariates <- unique(plot.data$group);
	variant.covariates[grep("NoVariant", variant.covariates)] 	<- default.colours(2, palette.type = 'spiral.morning')[1];
	variant.covariates[grep("Variant", variant.covariates)] 	<- default.colours(2, palette.type = 'spiral.morning')[2];

	# create covariates
	covariate <- list(
		rect = list(
			col = 'black',
			fill = tissue.covariates,
			lwd = 1.5
			),
		rect = list(
			col = 'black',
			fill = variant.covariates,
			lwd = 1.5
			)
		);
	cov.grob <- covariates.grob(
		covariates = covariate,
		ord = c(1:length(unique(plot.data$group))),
		side = 'top'
		);
	# create covariate legend
	cov.legend <- list(
		legend = list(
			colours = default.colours(2, palette.type = 'spiral.sunrise'),
			labels = c('Non-Associated','Associated'),
			cex = 1,
			title = 'Tissue'
			),
		legend = list(
			colours = default.colours(2, palette.type = 'spiral.morning'),
			labels = c('No Variant','Variant'),
			cex = 1,
			title = 'Mutation Status'
			)
		);
	cov.legend <- legend.grob(
		legends = cov.legend,
		title.cex = 1.2,
		label.cex = 1,
		size = 2
		);

	# calculate rectangle break points
	groups 	<- unique(plot.data$group);
	if (length(xaxis.labels) > 1) {
		add.rectangle <- TRUE
		if (length(xaxis.labels) == 2) {
			xleft.rectangle <- 0;
			xright.rectangle <- 4.5;
		} else {
			xleft.rectangle <- c(0,seq(8.5,length(groups)-3,8));
			xright.rectangle <- seq(4.5,length(groups)+1,8);
		}
	} else {
		add.rectangle <- FALSE
	}

	# create boxplot
	create.boxplot(
		drug.response ~ group,
		data = plot.data,
		filename = filename,
		add.rectangle = add.rectangle,
		add.stripplot = TRUE,
		xright.rectangle = xright.rectangle,
		xleft.rectangle = xleft.rectangle,
		ybottom.rectangle = 0,
		col.rectangle = 'grey',
		alpha.rectangle = 0.5,
		xaxis.lab = xaxis.labels,
		xat = seq(2.5,length(groups),4),
		ytop.rectangle = 1.05,
		ylimits = c(0,1.05),
		yat = seq(0,1,0.2),
		yaxis.lab = c("0.0","0.2","0.4","0.6","0.8","1.0"),
		xaxis.cex = 1.2,
		yaxis.cex = 1.2,
		symbol.cex = 0.3,
		points.cex = 0.3,
		xlab.label = 'Compound',
		ylab.label = 'Activity Area',
		xlab.cex = 1.5,
		ylab.cex = 1.5,
		legend = list(
			top = list(
				fun = cov.grob
				),
			right = list(
				fun = cov.legend
				)
			),
		top.padding = 10,
		width = width,
		resolution = 500
		);
}