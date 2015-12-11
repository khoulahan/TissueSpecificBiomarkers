### extract_variant_multiple_datasets.R ###########################################################
# extract literature variants drug response for multiple datasets 
## Arguments:
## compound.data.grid 		a matrix of drugs, first column, and datasets, second column 
## additional.variant		vector of gene names to test on every drug regardless if associated in 
##							literature
## neighbour.variants		logical indicating if samples with variant in gene present in the same 
##							pathway should be included 
## network					pathway network to use if neighbour.variants set as TRUE
###################################################################################################
extract.multiple.datasets <- function(compound.dataset.grid, additional.variant, neighbour.variants = FALSE, network = NULL) {
	plot.data <- list();
	for (i in 1:nrow(compound.dataset.grid)) {
		# set mutation profile and drug response profile
		if (compound.dataset.grid[i,2] == 'CCLE') {
			mutations <- mutations.ccle;
			drug.response <- drug.response.ccle;
		} else if (compound.dataset.grid[i,2] == 'CTDD') {
			mutations <- mutations.ccle;
			drug.response <- drug.response.ctdd;
		} else if (compound.dataset.grid[i,2] == 'CGP') {
			mutations <- mutations.cgp;
			drug.response <- drug.response.cgp;
		}
		# ensure compound name is specific to dataset
		compound.name <- map.drug.names[map.drug.names$CCLE == as.character(compound.dataset.grid[i,1]),
			as.character(compound.dataset.grid[i,2])];
		# check if drug is in dataset
		if (is.na(compound.name)) {
			print(paste(compound.dataset.grid[i,1], "not found in", compound.dataset.grid[i,2]));
			next;
		}
		if (!neighbour.variants) {
			source('extract_biomarkers_drug_response.R')
			# extract biomarker drug response for all drugs associated with BRAF and all three datasets
			tmp <- extract.biomarkers.drug.response(
				therapeutic.terms = therapeutic.terms, 
				biomarkers = biomarkers, 
				mutations = mutations, 
				drug.response = drug.response,
				compound = compound.name, 
				map.tissue = map.tissue,
				dataset = compound.dataset.grid[i,2],
				additional.variant = additional.variant
				);
		} else {
			source('extract_biomarkers_neighbour_drug_response.R')
			# extract biomarker drug response for all drugs associated with BRAF and all three datasets
			tmp <- extract.biomarkers.neighbour.drug.response(
				therapeutic.terms = therapeutic.terms, 
				biomarkers = biomarkers, 
				mutations = mutations, 
				drug.response = drug.response,
				compound = compound.name, 
				map.tissue = map.tissue,
				neighbour = network,
				dataset = compound.dataset.grid[i,2],
				additional.variant = additional.variant
				);
			}
		# add compound and dataset
		tmp$compound <- compound.dataset.grid[i,1];
		tmp$dataset <- compound.dataset.grid[i,2];
		plot.data[[paste(compound.dataset.grid[i,1], compound.dataset.grid[i,2], sep = '_')]] <- tmp;
		}
	do.call(rbind, plot.data);
	}