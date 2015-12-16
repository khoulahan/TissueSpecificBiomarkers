### generate_drug_response_curve.R ################################################################
# generate drug response curve for area under the curve 
## Arguments:
## drug.response		drug response feature matrix
## compound				name of compound to plot; NOTE: needs to match colnames of drug response 
##						matrix
## filename				filename to write plot to 
###################################################################################################
generate.drug.response.curve <- function(drug.response, compound, filename) {
	# get drug response for compound
	drug.response.compound <- drug.response[,compound];
	# remove NA
	drug.response.compound <- drug.response.compound[!is.na(drug.response.compound)];
	# order response
	drug.response.compound <- drug.response.compound[order(drug.response.compound)];
	# find quantiles
	summary.drug.response <- summary(drug.response.compound);

	if (length(drug.response.compound) == 0) {
		stop(paste("No drug response values for", compound, '...'));
	}

	# create plot data
	plot.data <- data.frame(
		drug.response = drug.response.compound,
		index = 1:length(drug.response.compound)
		);

	# plot curve
	create.scatterplot(
		drug.response ~ index,
		data = plot.data,
		filename = filename,
		type = 'l',
		lwd = 3,
		main = compound,
		abline.v = c(
			plot.data[plot.data$drug.response == max(plot.data[plot.data$drug.response <= summary.drug.response['1st Qu.'],'drug.response']),'index'], 
			plot.data[plot.data$drug.response == min(plot.data[plot.data$drug.response >= summary.drug.response['3rd Qu.'],'drug.response']),'index']
			),
		abline.col = 'firebrick3',
		main.cex = 2,
		xaxis.cex = 1.5,
		yaxis.cex = 1.5,
		xlab.cex = 1.5,
		xlab.label = 'Index',
		ylab.label = 'Area Under Curve',
		ylab.cex = 1.5
		)
	}