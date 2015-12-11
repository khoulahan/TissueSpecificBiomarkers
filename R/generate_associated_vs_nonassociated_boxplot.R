### Generate Boxplot of Associated vs Non-Associated Tissues ######################################
# groups all associated tissues into one category in an effort to increase power
# creates one boxplot for each variant associated with compound in plot data matrix
## Arguments:
## plot.data		Plot data matrix as outputted by extract_biomarkers_drug_response.R
## output.dir		Directory plots will be outputted to
## dataset			Dataset to analyze
## width			Width of plot
###################################################################################################
generate.associated.vs.nonassociated.boxplot <- function(
	plot.data,
	output.dir,
	dataset = 'CCLE',
	width = 6
	) {
		# load library
		library(BEST);

		# set list
		hypothesis.test.results <- list();

		# set date
		date <- Sys.Date();

		# create group for plotting by combining tissue and variant information 
		plot.data$group <- apply(
			plot.data,
			1,
			function(x) {
				if (x['tissue'] != 'OTHER') {
					paste('Associated',x['variant'],sep='_');
					} else {
						paste('NonAssociated',x['variant'],sep='_');
					}
				}
			);

		# remove any cell lines with no drug response
		plot.data <- plot.data[!is.na(plot.data$drug.response),];

		for (variant in unique(plot.data$gene)) {
			# subset plot data to only variant of interest
			plot.data.subset <- plot.data[plot.data$gene == variant,];
			
			# order plot data by groups
			plot.data.subset <- plot.data.subset[order(plot.data.subset$group),];

			# calculate posterior distributions of mu's for both population for each tissue
			legend.results <- c('Population Comparison:');
			for (tissue in c('Associated','NonAssociated')) {
				# subset groups 
				group1 <- plot.data.subset[plot.data.subset$group == paste(tissue, 'Variant', sep = '_') ,];
				group2 <- plot.data.subset[plot.data.subset$group == paste(tissue, 'NoVariant', sep = '_'),];

				# check that there are more than 5 samples
				if (nrow(group1) <= 5 | nrow(group2) <= 5) {
					next;
				}

				# bayesian estimation of posterior distributions
				BESTout <- BESTmcmc(group1$drug.response, group2$drug.response);
				summary.BESTout <- summary(BESTout);
				HDI.interval <- paste(
					round(summary.BESTout['muDiff','HDIlo'], digits = 2), 
					round(summary.BESTout['muDiff','HDIup'], digits = 2),
					sep = ':'
					);

				if (dataset == 'CCLE') {
					# CCLE response is area above curve 
					# therefore increased area = increased sensitivity and vice versa
					if (any(grepl('sensitivity|^response', group1$association))) {
						hypothesis <- 'greater';
					} else if (any(grepl('none', group1$association))) {
						hypothesis <- 'two.sided';
					} else {
						hypothesis <- 'less';
					}
				} else if (dataset %in% c('CTDD','CGP')) {
					# CTDD response is area under curve
					# therefore increased area = increased resistance and vice versa 
					if (any(grepl('sensitivity|^response', group1$association))) {
						hypothesis <- 'less';
					} else if (any(grepl('none', group1$association))) {
						hypothesis <- 'two.sided';
					} else {
						hypothesis <- 'greater';
					}
				}
				p.values <- wilcox.test(
					group1$drug.response,
					group2$drug.response,
					alternative = hypothesis
					)$p.value;

				legend.results <- c(
					legend.results,
					paste0(tissue, ':'),
					paste('P-value:', round(p.values, digits = 4)),
					paste('u1-u2:', round(summary.BESTout['muDiff','median'], digits = 2)),
					paste('HDI:', HDI.interval),
					''
					);

				#store hypothesis testing results in list 
				hypothesis.test.results[[paste(variant,tissue,sep='_')]] <- data.frame(
					variant = variant,
					tissue = tissue,
					p.value = round(p.values, digits=4),
					median.diff = round(summary.BESTout['muDiff','median'], digits = 2),
					HDIlo = round(summary.BESTout['muDiff','HDIlo'], digits = 2),
					HDIup = round(summary.BESTout['muDiff','HDIup'], digits = 2),
					effect.size = round(summary.BESTout['effSz','median'], digits = 2)
					);
			}

			# set point colours
			point.colours <- as.character(plot.data.subset$association);
			point.colours[grep('^sensitivity|response',point.colours)] <- 'firebrick3';
			point.colours[grep('resistance', point.colours)] <- 'dodgerblue3';
			point.colours[grep('none', point.colours)] <- 'black';

			# create boxplot
			create.boxplot(
				drug.response ~ group,
				data = plot.data.subset,
				filename = paste0(
					output.dir,
					date,
					'_',
					variant,
					'AllTissueBoxplot.tiff'
					),
				add.stripplot = TRUE,
				points.col = point.colours, 
				symbol.cex = 0.25,
				add.rectangle = TRUE,
				xright.rectangle = 2.5,
				xleft.rectangle = 0,
				ybottom.rectangle = 0,
				col.rectangle = 'grey',
				alpha.rectangle = 0.5,
				xaxis.lab = c('Associated Tissue','Non-Associated Tissue'),
				xat = seq(1.5, 4, 2),
				ytop.rectangle = round(max(plot.data.subset$drug.response)),
				ylimits = c(0,round(max(plot.data.subset$drug.response))),
				yat = seq(0,round(max(plot.data.subset$drug.response)),1),
				yaxis.lab = seq(0,round(max(plot.data.subset$drug.response)),1),
				xaxis.cex = 0.7,
				xaxis.rot = 90,
				legend = list(
					right = list(
						fun = draw.key,
						args = list(
							key = list(
								text = list(
									lab = legend.results,
									cex = 0.8,
									font = 2
									)
								)
							)
						)
					),
				yaxis.cex = 1,
				xlab.label = 'Tissue',
				ylab.label = 'Activity Area',
				xlab.cex = 1.5,
				ylab.cex = 1.5,
				width = width,
				resolution = 500
				);
		}
	hypothesis.test.results <- do.call(rbind,hypothesis.test.results);
	return(hypothesis.test.results);
	}
