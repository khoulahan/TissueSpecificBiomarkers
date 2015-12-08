### find_associated_variants.R ####################################################################
# find all variants associated with specified drug and return in a more condensed format 
###################################################################################################
find.associated.variants <- function(
	compound, 
	biomarkers, 
	therapeutic.terms, 
	dataset = 'CCLE', 
	variant.type = 'mutation'
	) {

		# find terms in biomarkers that match compound of interest
		if (dataset == 'CCLE') {
			associated.biomarkers.terms <- therapeutic.terms[grep(compound, therapeutic.terms$AssociatedCompoundsCCLE),];
			} else if (dataset == 'CTDD') {
				associated.biomarkers.terms <- therapeutic.terms[grep(compound, therapeutic.terms$AssociatedCompoundsCTDD),];
			} else if (dataset == 'CGP') {
				associated.biomarkers.terms <- therapeutic.terms[grep(compound, therapeutic.terms$AssociatedCompoundsCGP),];
			} else {
				stop("Please specify appropriate dataset. Either CCLE, CTDD or CGP ...")
			}

		# if no associated terms, stop the function and report to user
		if (nrow(associated.biomarkers.terms) == 0) {
			stop(paste("No therapeutic term in biomarkers associated with", compound,"..."));
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
		# separate mutations from amplifications and rearrangments
		if (variant.type == 'cnv') {
			associated.biomarkers.variants <- associated.biomarkers.variants[associated.biomarkers.variants$Variant %in% c('amplification','deletion'),];
		} else if (variant.type == 'mutation') {
			associated.biomarkers.variants <- associated.biomarkers.variants[!associated.biomarkers.variants$Variant %in% c('rearrangement','amplification','splice variant mRNA','deletion'),];
		} else {
			stop("Please specify valid variant type. Options are 'mutation' or 'cnv'");
		}
		# re factor variants to drop levels
		associated.biomarkers.variants$Gene <- factor(associated.biomarkers.variants$Gene);

		if (nrow(associated.biomarkers.variants) == 0) {
			stop(paste("No", variant.type, "variants associated with", compound, "..."));
		}

		return(associated.biomarkers.variants);
	}
