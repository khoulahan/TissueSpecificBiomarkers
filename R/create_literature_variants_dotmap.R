### create_literature_variants_dotmap.R ###########################################################
# pulls out all significant variant-compound combinations and plots as a dotmap
## dot size and colour indicates effect size and direction 
## background shading indicates FDR
## covariates indicating if result is in tissue reported in literature or non-reported tissues
## covariates also show compound analyzed
### Arguments:
### hypothesis.results		Hypothesis results file formed as aggregate from generate_associated_vs_nonassociated_boxplot.R
### map.drug.names			File mapping drug names in each dataset as well specifies associated colour
### filename				Filename to output plot to
### dataset					Dataset to plot
### width					Width of plot
### height 					Height of plot
### x.legend				X-coordinate of "FDR" label
### y.legend				Y-coordinate of "FDR" label
### spot.size.function		Function to determine spot size for dotmap
###################################################################################################
create.literature.variants.dotmap <- function(
	hypothesis.results, 
	map.drug.names, 
	filename, 
	dataset = 'CCLE', 
	width = 6, 
	height = 6,
	spot.size.function = function(x) { 0.5 + (1.1 * abs(x)) }
	) {

		# add dataset column 
		hypothesis.results$dataset <- dataset;
		
		### MAKE DATASETS CONSISTENT - MATCH EFFECT SIZE SIGN AND COMPOUND NAMES ###
		if (dataset == 'CCLE') {
			# change CCLE effect sizes sign to match CTDD and CGP
			hypothesis.results[hypothesis.results$dataset == 'CCLE','effect.size'] <- hypothesis.results[hypothesis.results$dataset == 'CCLE','effect.size']*(-1)
		} else if (dataset == 'CTDD') {
			# change all drug names to match CCLE
			hypothesis.results$compound <- unlist(sapply(as.character(hypothesis.results$compound), function(x) {
				map.drug.names[map.drug.names$CTDD == x & !is.na(map.drug.names$CTDD),'CCLE'];
				}))
		} else if (dataset == 'CGP') {
			# change all drug names to match CCLE
			hypothesis.results$compound <- unlist(sapply(as.character(hypothesis.results$compound), function(x) {
				map.drug.names[map.drug.names$CGP == x & !is.na(map.drug.names$CGP),'CGP']
				}))
		}

		### CORRECT P-VALUES & PULL OUT SIGNIFICANT VARIANTS ###
		# adjust pvalues
		hypothesis.results$p.value <- p.adjust(hypothesis.results$p.value, method = 'fdr')
		# only variants with significant results
		sig.variant.drug <- by(
			hypothesis.results,
			hypothesis.results[,c('variant','compound')],
			function(subset) {
				if (any(subset$p.value < 0.05)) {
					return(data.frame(variant = unique(subset$variant), drug = unique(subset$compound)))
				} else {
					return(NULL)
				}
			});
		sig.variant.drug <- do.call(rbind, sig.variant.drug)
		# filter hypothesis results to keep only significant compounds and variants
		hypothesis.results <- hypothesis.results[hypothesis.results$variant %in% sig.variant.drug$variant & hypothesis.results$compound %in% sig.variant.drug$drug,];

		### FORMAT DATA FOR PLOTTING ###
		# create identifying group column
		hypothesis.results$group <- apply(
			hypothesis.results,
			1,
			function(x) {
				paste(x['compound'],x['dataset'],x['tissue'], sep = '_')
				})
		# factor variant and group columns
		hypothesis.results$variant 	<- factor(hypothesis.results$variant);
		hypothesis.results$group 	<- factor(hypothesis.results$group);
		# reformat p-value, num cell lines and effect size
		reformat.data.frames <- function(hypothesis.results, variable) {
			tmp <- matrix(
				nrow=nlevels(hypothesis.results$group), 
				ncol=nlevels(hypothesis.results$variant),
				dimnames=list(levels(hypothesis.results$group), levels(hypothesis.results$variant))
				)
			tmp[cbind(hypothesis.results$group, hypothesis.results$variant)] <- hypothesis.results[,variable]
			tmp
		}
		pvalues <- reformat.data.frames(hypothesis.results,'p.value');
		literature <- reformat.data.frames(hypothesis.results, 'literature');
		effect.size <- reformat.data.frames(hypothesis.results,'effect.size');

		### ORDER PLOTTING DATA ###
		# order by pvalues
		max.effect.size <- apply(effect.size,2, function(x) {max(abs(x), na.rm = TRUE)});
		effect.size <- effect.size[,order(
			colSums(-log10(pvalues), na.rm = TRUE),
			max.effect.size,
			decreasing = TRUE)];
		literature <- literature[,order(
			colSums(-log10(pvalues), na.rm = TRUE),
			max.effect.size,
			decreasing = TRUE)];
		pvalues <- pvalues[,order(
			colSums(-log10(pvalues), na.rm = TRUE), 
			max.effect.size,
			decreasing = TRUE)];
		# order compounds as in map drug names file so that similar mode of actions are adjacent
		compounds.order <- sapply(as.character(map.drug.names$CCLE), function(x) { grep(x, rownames(effect.size))});
		compounds.order <- unlist(compounds.order)

		effect.size <- effect.size[compounds.order,];
		literature	<- literature[compounds.order,];
		pvalues 	<- pvalues[compounds.order,];

		### SPOT COLOUR FUNCTION ###
		spot.colour.function <- function(x) {
			colours <- rep('white', length(x));
			colours[sign(x) ==  1] <- 'dodgerblue2';
			colours[sign(x) == -1] <- 'firebrick2';
			return(colours);
			}

		### GENERATE COVARIATES ###
		# tissue covariate colours
		tissue.colours <- rownames(pvalues);
		tissue.colours[grep('NonAssociated', tissue.colours)] <- default.colours(2, palette.type = 'spiral.sunrise')[1];
		tissue.colours[grep('Associated', tissue.colours)] <- default.colours(2, palette.type = 'spiral.sunrise')[2];
		# compound covariate colours
		compound.colours <- rownames(pvalues);
		compounds.used <- rep(NA, nrow(map.drug.names));
		for (i in 1:nrow(map.drug.names)) {
			compounds.used[i] <- any(grepl(map.drug.names[i,'CCLE'], compound.colours));
			compound.colours[grep(map.drug.names[i,'CCLE'], compound.colours)] <- as.character(map.drug.names[i,'Colours']);
		}
		# create covariates
		covariate <- list(
			rect = list(
				col = 'black',
				fill = tissue.colours,
				lwd = 1.5
				),
			rect = list(
				col = 'black',
				fill = compound.colours,
				lwd = 1.5
				)
			);
		cov.grob <- covariates.grob(
			covariates = covariate,
			ord = c(length(tissue.colours):1),
			side = 'right'
			);

		### GENERATE COVARIATE LEGEND ###
		cov.legend <- list(
			legend = list(
				colours = default.colours(2, palette.type = 'spiral.sunrise'),
				labels = c('Non-Associated','Associated'),
				cex = 1,
				title = 'Tissue'
				),
			legend = list(
				colours = as.character(map.drug.names[compounds.used,'Colours']),
				labels = as.character(map.drug.names[compounds.used,'CCLE']),
				cex = 1,
				title = 'Compounds'
				),
			legend = list(
				colours = c('black','white','grey'),
				labels = c('Literature','Non-Literature','N/A'),
				cex = 1,
				title = 'Literature Reported'
				)
			);

		cov.legend <- legend.grob(
			legends = cov.legend,
			title.cex = 0.8,
			label.cex = 0.75,
			size = 1
			);

		### CREATE GRID LAYOUT FOR COVARIATES AND LEGEND ###
		right.layout <- grid.layout(
			nrow = 1,
			ncol = 2,
			width = unit(
				x = c(0,1),
				units = rep('lines',2)
				),
			heights = unit(
				x = c(1,1),
				units = rep('npc', 1)
				)
			);
		right.grob <- frameGrob(layout = right.layout);
		right.grob <- packGrob(
			frame = right.grob,
			grob = cov.grob,
			row = 1,
			col = 1
			);
		right.grob <- packGrob(
			frame = right.grob,
			grob = cov.legend,
			row = 1,
			col = 2
		 	);

		### CREATE DOT LEGEND ###
		key.sizes <- c(1.5, 1,0.5,-0.5,-1,-1.5)
		dot.grob <- draw.key(
			list(
				space = 'right',
				points = list(
					cex = spot.size.function(key.sizes),
					col = 'black',
					pch = 21,
					fill = spot.colour.function(key.sizes)
					),
				text = list(
					lab = c("+1.5","+1.0","+0.5","-0.5","-1.0","-1.5"),
					cex = 0.75,
					adj = 1
					),
				title = 'Effect Size',
				cex.title = 0.8,
				background = 'white',
				padding.text = 3
				)
			);

		### CREATE HEATMAP ###
		heatmap <- create.heatmap(
			literature,
			clustering.method = 'none',
			colour.scheme = c('black','white'),
			yaxis.lab = rep('', nrow(pvalues)),
			grid.col = TRUE,
			same.as.matrix = TRUE
			)

		### CREATE DOTMAP ###
		dotmap <- create.dotmap(
			x = effect.size,
			bg.data = -log10(pvalues),
			spot.colour.function = spot.colour.function,
			spot.size.function = spot.size.function,
			pch = 21,
			bg.alpha = 1,
			yaxis.lab = rep('', nrow(pvalues)),
			colour.scheme = c('white','black'),
			colourkey = FALSE,
			yaxis.cex = 0.5,
			xaxis.cex = 1,
			xaxis.rot = 90,
			NA.spot.size = 1.5,
			resolution = 500,
			at = seq(0,5,1),
			axis.bottom = 1.05,
			bottom.padding = 4
			)

		### CREATE MULTIPLOT ###
		create.multiplot(
			list(
				dotmap,
				heatmap
				),
			filename = filename,
			plot.layout = c(2,1),
			panel.widths = c(2.5,1.5),
			xaxis.cex = 0.8,
			xaxis.rot = 90,
			retrieve.plot.labels = TRUE,
			print.new.legend = TRUE,
			bottom.padding = 5,
			height = height,
			width = width,
			legend = list(
				right = list(
					fun = right.grob
					),
				left = list(
					fun = dot.grob
					),
				inside = list(
					fun = BoutrosLab.plotting.general::create.colourkey(
						x = effect.size,
						colour.scheme = c('white','black'),
						colourkey.labels.cex = 0.9,
						at = seq(0,5,1), 
						colourkey.labels.at = seq(0,5,1),
						colourkey.labels = c(
							expression(10^0),
							expression(10^-1),
							expression(10^-2),
							expression(10^-3),
							expression(10^-4),
							substitute(p.value<=10^-5, list(p.value = ''))
							),
						placement = viewport(just='left',x=0.5,y=-0.7,width=0.55)
						)
					)
				) 
			)
	}