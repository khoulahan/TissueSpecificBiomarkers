extract.biomarkers.drug.response <- function(
	therapeutic.terms,
	biomarkers,
	mutations = NULL,
	cnv = NULL,
	drug.response,
	compound,
	map.tissue,
	dataset = 'CCLE',
	variant.type = 'mutation'
	) {

		# find all variants associated with specified compound
		source('~/Documents/TissueSpecificBiomarkers/find_associated_variants.R');
		associated.biomarkers.variants <- find.associated.variants(
			compound,
			biomarkers = biomarkers,
			therapeutic.terms = therapeutic.terms,
			dataset = dataset,
			variant.type = variant.type
			);

		# for each gene find tissues that associated with listed diseases and pull out drug response for all those tissues	
		drug.response.variants <- by(
			associated.biomarkers.variants,
			associated.biomarkers.variants$Gene,
			function(x) {
				tissues <- unique(map.tissue[map.tissue$disease %in% x[,'Disease'],grep(dataset, colnames(map.tissue))]);
				gene <- as.character(unique(x$Gene));
				variant.cat <- max(as.character(x$Variant));
				# check that gene is in mutations if not move onto next variant
				if (variant.type == 'mutation' & !gene %in% colnames(mutations) | variant.type == 'cnv' & !gene %in% colnames(cnv)) {
					return(NULL);
				}
				# get drug reponse of all tissues associated with gene of interest
				associated.drug.response <- sapply(
					tissues,
					function(tissue) {
						if (variant.type == 'mutation') {
							# find cell lines associated with tissue that have variant
							cell.lines.variant <- rownames(mutations[mutations[,gene] == 1,])[grep(tissue, 
								rownames(mutations[mutations[,gene] == 1,]), ignore.case = TRUE)];
						} else if (variant.type == 'cnv') {
							# find cell lines with amplification or deletion 
							copy.number.variant <- round(2*(2^cnv[,gene]));
							if (variant.cat == 'amplification') {
								cell.lines.variant <- rownames(cnv)[which(copy.number.variant > 2)];
							} else if (variant.cat == 'deletion') {
								cell.lines.variant <- rownames(cnv)[which(copy.number.variant < 2)];
							}
							cell.lines.variant <- cell.lines.variant[grep(tissue, cell.lines.variant, ignore.case = TRUE)];
						}
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
				# find non associated cell lines with variant
				if (variant.type == 'mutation') {
					# find cell lines not associated with tissue that have variant
					cell.lines.variant <- rownames(mutations[mutations[,gene] == 1,])[rownames(
						mutations[mutations[,gene] == 1,]) %in% non.associated.tissues];
				} else if (variant.type == 'cnv') {
					# find cell lines with amplification or deletion that are not associated with tissue
					copy.number.variant <- round(2*(2^cnv[,gene]));
					if (variant.cat == 'amplification') {
						cell.lines.variant <- rownames(cnv)[which(copy.number.variant > 2)];
					} else if (variant.cat == 'deletion') {
						cell.lines.variant <- rownames(cnv)[which(copy.number.variant < 2)];
					}
					cell.lines.variant <- cell.lines.variant[cell.lines.variant %in% non.associated.tissues];
				}

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
