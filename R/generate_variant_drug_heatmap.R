#### generate variant drug association heatmap ####################################################
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