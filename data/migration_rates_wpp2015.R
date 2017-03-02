# Construct historical migration rates from wpp2015 datasets

migration_rates_wpp2015 <- local({
	library(wpp2015)
	data(UNlocations, package="wpp2015")
	countries <- subset(UNlocations, location_type==4)$country_code

	# population & migration
	data(pop, package="wpp2015")
	popu <- subset(pop, country_code %in% countries)

	data(migration, package="wpp2015")
	# put into same order as popu
	mig <- merge(popu[,'country_code', drop=FALSE], migration[,1:which(colnames(migration)=="2010-2015")], sort=FALSE)

	# put popu and mig into same shape
	popu.matrix <- popu[,4:ncol(popu)]
	rownames(popu.matrix) <- popu$country_code
	mig.matrix <- mig[,3:ncol(mig)]
	rownames(mig.matrix) <- mig$country_code

	rates <- mig.matrix / (popu.matrix - mig.matrix)
	cbind(data.frame(country_code=rownames(rates), name=popu$name), rates)
})
