### Generate Boxplot of All Associated Tissues vs Other Tissues ###########################################
# generate boxplot of all individual tissues associated with given drug in literature file
# generates one plot for each variant associated with the specified compound in the plot data matrix
## Arguments:
## plot.data		Plot data matrix as outputted by extract_biomarkers_drug_response.R
## output.dir		Directory plots will be outputted to
## dataset			Dataset to analyze
## width			Width of plot
###################################################################################################
generate.associated.tissue.boxplot <- function(
	plot.data,
	output.dir,
	dataset = 'CCLE',
	width = 6
	) {
		# load library
		library(BEST);
		# set date
		date <- Sys.Date();
		# set list
		hypothesis.test.results <- list();

		# create group for plotting by combining tissue and variant information 
		plot.data$group <- apply(
			plot.data,
			1,
			function(x) {
				paste(x['tissue'],x['variant'],sep='_');
				}
			);

		# remove any cell lines with no drug response
		plot.data <- plot.data[!is.na(plot.data$drug.response),];

		for (variant in unique(plot.data$gene)) {
			# subset plot data to only variant of interest
			plot.data.subset <- plot.data[plot.data$gene == variant,];
			# re factor group to drop unused levels
			# but ensure that both variant and no variant remain for each tissue type
			new.levels <- apply(
				expand.grid(unique(plot.data.subset$tissue), c('Variant','NoVariant')),
				1,
				function(x) {
					paste(x[1],x[2],sep = '_')
					}
				);
			new.levels <- new.levels[order(new.levels)];
			plot.data.subset$group <- factor(plot.data.subset$group, levels = new.levels, ordered = TRUE);

			# if not all levels are present all 0 entries for missing levels
			if (!all(new.levels %in% plot.data.subset$group)) {
				missing.levels <- new.levels[which(!new.levels %in% plot.data.subset$group)];

				for (i in missing.levels) {
					tissue <- unlist(strsplit(i, split='_'));
					tissue <- paste(tissue[1:(length(tissue)-1)], collapse = '_');

					tmp.row <- data.frame(
						drug.response = -1,
						cell.line = NA,
						tissue = tissue,
						gene = variant,
						variant = NA,
						association = 'none',
						group = i
						);
					plot.data.subset <- rbind(plot.data.subset, tmp.row);
				}
			}
			# order plot data by groups
			plot.data.subset <- plot.data.subset[order(plot.data.subset$group),];

			# calculate posterior distributions of mu's for both population for each tissue
			legend.results <- c('Population Comparison:');
			for (tissue in unique(plot.data.subset$tissue)) {
				# subset groups 
				group1 <- plot.data.subset[plot.data.subset$tissue == tissue & plot.data.subset$variant == 'Variant',];
				group2 <- plot.data.subset[plot.data.subset$tissue == tissue & plot.data.subset$variant == 'NoVariant',];

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

			# find number of groups
			num.groups <- length(levels(plot.data.subset$group));

			# set rectangle parameters
			if (num.groups == 2) {
				add.rectangle <- FALSE;
			} else if (num.groups == 4) {
				add.rectangle <- TRUE;
				xright.rectangle <- 0;
				xleft.rectangle <- 2.5;
			} else if (num.groups == 6) {
				add.rectangle <- TRUE;
				xright.rectangle <- c(0,4.5);
				xleft.rectangle <- c(2.5,6.5)
			} else {
				add.rectangle <- TRUE;
				xright.rectangle <- c(0,seq(4.5,num.groups,4));
				xleft.rectangle <- seq(2.5,num.groups+1,4);
			}

			# set point colours
			point.colours <- as.character(plot.data.subset$association);
			point.colours[grep('^sensitivity|response',point.colours)] <- 'firebrick3';
			point.colours[grep('resistance', point.colours)] <- 'dodgerblue3';
			point.colours[grep('none', point.colours)] <- 'black';

			# for sake of labeling change HAEMATOPOIETIC_AND_LYMPHOID_TISSUE to BLOOD
			plot.data.subset$tissue <- as.character(plot.data.subset$tissue);
			plot.data.subset[plot.data.subset$tissue == 'HAEMATOPOIETIC_AND_LYMPHOID_TISSUE','tissue'] <- 'BLOOD';

			# create boxplot
			create.boxplot(
				drug.response ~ group,
				data = plot.data.subset,
				filename = paste0(
					output.dir,
					date,
					'_',
					variant,
					'SubsetTissueBoxplot.tiff'
					),
				add.stripplot = TRUE,
				points.col = point.colours, 
				symbol.cex = 0.25,
				add.rectangle = add.rectangle,
				xright.rectangle = xright.rectangle,
				xleft.rectangle = xleft.rectangle,
				ybottom.rectangle = 0,
				col.rectangle = 'grey',
				alpha.rectangle = 0.5,
				xaxis.lab = unique(plot.data.subset$tissue),
				xat = seq(1.5, num.groups, 2),
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
