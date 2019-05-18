# Construct historical migration rates from wpp2017 datasets

migration_rates_wpp2019 <- local({
	library(wpp2019)
	data(UNlocations, package="wpp2019")
	countries <- subset(UNlocations, location_type==4)$country_code

	# population & migration
	data(pop, package="wpp2019")
	popu <- subset(pop, country_code %in% countries)

	data(migration, package="wpp2019")
	# put into same order as popu
	mig <- merge(popu[,'country_code', drop=FALSE], migration[,1:which(colnames(migration)=="2015-2020")], sort=FALSE)

	# put popu and mig into same shape
	popu.matrix <- popu[,4:ncol(popu)]
	rownames(popu.matrix) <- popu$country_code
	mig.matrix <- mig[,3:ncol(mig)]
	rownames(mig.matrix) <- mig$country_code

	rates <- mig.matrix / (popu.matrix - mig.matrix)
	rates <- cbind(data.frame(country_code=rownames(rates), name=popu$name), rates)
	rates
})
