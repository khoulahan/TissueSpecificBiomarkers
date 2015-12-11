#### generate variant drug association heatmap ####################################################
# generates heatmap showing variants each drug is associated with in the biomarkers file
## Arguments:
## therapeutic.terms	File mapping therapeutic terms in biomarkers file to drug names in each dataset
## biomarkers			Literature cultivated biomarkers file
## mutations			Mutations features matrix; binary matrix sample by gene
## drug response		Drug response feature matrix; sample by compound
## map.tissue			File mapping disease in biomarkers file to tissue specification in each dataset
## compounds			Compounds to include in plot
## filename				Filename to write plot to
## plot					Logical indicating whether to plot heatmap or just return plot data
## dataset				The dataset to analyze
## width				Width of plot
## height				Height of plot
###################################################################################################
generate.variant.drug.heatmap <- function(
	therapeutic.terms,
	biomarkers,
	mutations,
	drug.response,
	map.tissue,
	compounds,
	filename = NULL,
	plot = TRUE,
	dataset = 'CCLE',
	width = 5,
	height = 7
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

	# find the genes associated with each drug 
	associated.variants <- lapply(
		plot.data,
		function(mat) {
			as.character(unique(mat$gene))
			});
	unique.variants <- unique(unlist(associated.variants));

	# generate heatmap plot data
	heatmap.plot.data <- data.frame(matrix(nrow = length(compounds), ncol = length(unique.variants)));
	colnames(heatmap.plot.data) <- unique.variants;
	rownames(heatmap.plot.data) <- compounds;

	for (compound in compounds) {
		heatmap.plot.data[compound,associated.variants[[compound]]] <- 1;
	}
	heatmap.plot.data[is.na(heatmap.plot.data)] <- 0;

	# order heatmap data by most present genes
	heatmap.plot.data <- heatmap.plot.data[,order(colSums(heatmap.plot.data))]

	if (plot) {
		create.heatmap(
			heatmap.plot.data,
			filename,
			clustering.method = 'none',
			colour.scheme = c('white','seagreen'),
			xaxis.lab = NA,
			yaxis.lab = NA,
			print.colour.key = FALSE,
			xaxis.cex = 1.5,
			grid.row = TRUE,
			grid.col = TRUE,
			width = width,
			height = height,
			resolution = 500
			);
	}
	# return matrix
	return(heatmap.plot.data);
}