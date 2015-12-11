### extract_biomarkers_neighbour_drug_response.R ##################################################
# extract drug response of cell lines with variant in associated biomarker plus its neighbours
## Arguments:
## therapeutic.terms	File mapping therapeutic terms in biomarkers file to drug names in each dataset
## biomarkers			Literature cultivated biomarkers file
## mutations			Mutations features matrix; binary matrix sample by gene
## drug response		Drug response feature matrix; sample by compound
## network				Pathway network used to determine genes in same pathway
## map.tissue			File mapping disease in biomarkers file to tissue specification in each dataset
## dataset				The dataset to analyze
###################################################################################################
extract.biomarkers.neighbour.drug.response <- function(
	therapeutic.terms,
	biomarkers,
	mutations,
	drug.response,
	compound,
	network,
	map.tissue,
	dataset = 'CCLE'
	) {
		# find all variants associated with specified compound
		source('find_associated_variants.R');
		associated.biomarkers.variants <- find.associated.variants(
			compound,
			biomarkers = biomarkers,
			therapeutic.terms = therapeutic.terms,
			dataset = dataset
			);

		drug.response.variants <- by(
			associated.biomarkers.variants,
			associated.biomarkers.variants$Gene,
			function(x) {
				tissues <- unique(map.tissue[map.tissue$disease %in% x[,'Disease'],grep(dataset, colnames(map.tissue))]);
				gene <- as.character(unique(x$Gene));
				# find neigbours of gene 
				neighbour.variants <- c(gene, network[network[,1] == gene, 2], network[network[,2] == gene,1]);
				# keep only neighbouring variants present in mutation matrix
				neighbour.variants <- neighbour.variants[neighbour.variants %in% colnames(mutations)];
				# find cell lines with variant
				if (length(neighbour.variants) == 1) {
					is.variant <- mutations[,neighbour.variants] == 1;
					names(is.variant) <- rownames(mutations);
				} else {
					is.variant <- apply(mutations[,neighbour.variants], 1, function(x) {any(x == 1)});
					}
				# get drug reponse of all tissues associated with gene of interest
				associated.drug.response <- sapply(
					tissues,
					function(tissue) {
						# find cell lines associated with tissue that have variant
						cell.lines.variant <- names(is.variant[is.variant & grepl(tissue, 
							names(is.variant), ignore.case = TRUE)]);

						# only keep cell lines that have drug response
						cell.lines.variant <- cell.lines.variant[cell.lines.variant %in% rownames(drug.response)];
						# find cell lines associated with tissue that do not have variant
						cell.lines.no.variant <- rownames(drug.response)[grep(tissue, rownames(drug.response), ignore.case = TRUE)];
						cell.lines.no.variant <- cell.lines.no.variant[!cell.lines.no.variant %in% cell.lines.variant];

						# find the most predominant therapeutic association 
						therapeutic.association <- summary(
							x[x$Disease %in% map.tissue[map.tissue[,grep(dataset, colnames(map.tissue))] == tissue,'disease'],'Association']
							);
						therapeutic.association <- names(therapeutic.association)[therapeutic.association == max(therapeutic.association)][1];

						data.frame(
							drug.response = c(
								drug.response[cell.lines.variant,compound],
								drug.response[cell.lines.no.variant,compound]
								),
							cell.line = c(
								cell.lines.variant,
								cell.lines.no.variant
								),
							tissue = rep(tissue, length(c(cell.lines.variant, cell.lines.no.variant))),
							gene = rep(gene, length(c(cell.lines.variant, cell.lines.no.variant))),
							variant = c(
								rep('Variant', length(cell.lines.variant)),
								rep('NoVariant', length(cell.lines.no.variant))
								),
							association = c(
								rep(therapeutic.association, length(cell.lines.variant)),
								rep('none', length(cell.lines.no.variant))
								)
							);
						},
					simplify = FALSE
					);
				associated.drug.response <- do.call(rbind, associated.drug.response);
				# add the remaining tissues binned as "other"
				associated.tissues <- paste(tissues, collapse = '|');
				# find all tissues not associated with variants
				non.associated.tissues <- rownames(drug.response)[grep(associated.tissues, rownames(drug.response), 
					invert = TRUE, ignore.case = TRUE)];
				# find cell lines not associated with tissue that have variant
				cell.lines.variant <- names(is.variant)[is.variant & names(is.variant) %in% non.associated.tissues];

				# find non associated cell lines without variant
				cell.lines.no.variant <- non.associated.tissues[!non.associated.tissues %in% cell.lines.variant];

				# generate drug response matrix of non associated cell lines
				non.associated.drug.response <- data.frame(
					drug.response = c(
						drug.response[cell.lines.variant,compound],
						drug.response[cell.lines.no.variant,compound]
						),
					cell.line = c(
						cell.lines.variant,
						cell.lines.no.variant
						),
					tissue = rep('OTHER', length(c(cell.lines.variant, cell.lines.no.variant))),
					gene = rep(gene, length(c(cell.lines.variant, cell.lines.no.variant))),
					variant = c(
						rep('Variant', length(cell.lines.variant)),
						rep('NoVariant', length(cell.lines.no.variant))
						),
					association = rep('none', length(c(cell.lines.variant, cell.lines.no.variant)))
					);

				return(rbind(associated.drug.response, non.associated.drug.response));
				}
			);
		drug.response.variants <- do.call(rbind, drug.response.variants);

		return(drug.response.variants);
	}