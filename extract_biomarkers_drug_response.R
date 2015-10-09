extract.biomarkers.drug.response <- function(
	therapeutic.terms,
	biomarkers,
	mutations,
	drug.response,
	compound,
	map.tissue,
	dataset = 'CCLE'
	) {

		# define non silent mutations 
		non.silent.mutations <- c(
			'Frame_Shift_Del',
			'Frame_Shift_Ins',
			'Indel',
			'In_Frame_Del',
			'In_Frame_Ins',
			'Missense_Mutation',
			'Nonsense_Mutation',
			'Splice_Site',
			'Splice_Site_Del',
			'Splice_Site_DNP',
			'Splice_Site_Ins',
			'Splice_Site_SNP',
			'De_novo_Start_InFrame',
			'De_novo_Start_OutOfFrame',
			'Nonstop_Mutation',
			'Start_Codon_Del',
			'Stop_Codon_DNP',
			'Stop_Codon_Ins'
			);

		# find terms in biomarkers that match compound of interest
		if (dataset == 'CCLE') {
			associated.biomarkers.terms <- therapeutic.terms[grep(compound, therapeutic.terms$AssociatedCompoundsCCLE),];
			} else if (dataset == 'CTDD') {
				associated.biomarkers.terms <- therapeutic.terms[grep(compound, therapeutic.terms$AssociatedCompoundsCTDD),];
			} else {
				stop("Please specify appropriate dataset. Either CCLE or CTDD ...")
			}

		# if no associated terms, stop the function and report to user
		if (nrow(associated.biomarkers.terms) == 0) {
			stop("No therapeutic term in biomarkers associated with specified drug ...");
		}

		# pull out variants associated with compound
		# need to iterate over all the therapeutic.response columns
		associated.biomarkers.variants <- apply(
			associated.biomarkers.terms,
			1,
			function(term) {
				tmp.term <- term['BiomarkersTherapeuticTerm'];
				tmp.list <- list();
				for (i in 1:8) {
					colnames.keep <- c(
						'Disease',
						'Gene',
						'Variant',
						paste('Association',i,sep='_'),
						paste('Therapeutic.context',i,sep='_')
						);
					tmp.list[[i]] <- biomarkers[grep(tmp.term, biomarkers[,paste('Therapeutic.context', i, sep = '_')]), colnames.keep];
					colnames(tmp.list[[i]]) <- c('Disease','Gene','Variant','Association','Therapeutic.context');
					}
				do.call(rbind, tmp.list);
				}
			);
		associated.biomarkers.variants <- do.call(rbind, associated.biomarkers.variants);
		# remove duplicated rows
		associated.biomarkers.variants <- associated.biomarkers.variants[!duplicated(associated.biomarkers.variants),];
		# keep only mutations, not rearrangements or amplifications
		associated.biomarkers.variants <- associated.biomarkers.variants[!associated.biomarkers.variants$Variant %in% c('rearrangement','amplification'),];
		# re factor variants to drop levels
		associated.biomarkers.variants$Gene <- factor(associated.biomarkers.variants$Gene);

		if (nrow(associated.biomarkers.variants) == 0) {
			stop("No mutation variants associated with compound specified ...")
		}

		# for each gene find tissues that associated with listed diseases and pull out drug response for all those tissues	
		drug.response.variants <- by(
			associated.biomarkers.variants,
			associated.biomarkers.variants$Gene,
			function(x) {
				tissues <- unique(map.tissue[map.tissue$disease %in% x[,'Disease'],'tissue']);
				gene <- as.character(unique(x$Gene));
				# find which cell lines have variant
				mutations.variants <- mutations[mutations$Hugo_Symbol == gene & mutations$Variant_Classification %in% non.silent.mutations,];
				# get drug reponse of all tissues associated with gene of interest
				associated.drug.response <- sapply(
					tissues,
					function(tissue) {
						# find cell lines associated with tissue that have variant
						cell.lines.variant <- as.character(
							unique(mutations.variants[grep(tissue, mutations.variants$Tumor_Sample_Barcode),'Tumor_Sample_Barcode'])
							);
						# only keep cell lines that have drug response
						cell.lines.variant <- cell.lines.variant[cell.lines.variant %in% rownames(drug.response)];
						# find cell lines associated with tissue that do not have variant
						cell.lines.no.variant <- rownames(drug.response)[grep(tissue, rownames(drug.response), ignore.case = TRUE)];
						cell.lines.no.variant <- cell.lines.no.variant[!cell.lines.no.variant %in% cell.lines.variant];

						# find the most predominant therapeutic association 
						therapeutic.association <- summary(
							x[x$Disease %in% map.tissue[map.tissue$tissue == tissue,'disease'],'Association']
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
				cell.lines.variant <- as.character(
					unique(
						mutations.variants$Tumor_Sample_Barcode[mutations.variants$Tumor_Sample_Barcode %in% non.associated.tissues]
						)
					);
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

				rbind(associated.drug.response, non.associated.drug.response);
				}
			);
		drug.response.variants <- do.call(rbind, drug.response.variants);

		return(drug.response.variants);
	}
