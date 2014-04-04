library(bayesPop)
start.test <- function(name) cat('\n<=== Starting test of', name,'====\n')
test.ok <- function(name) cat('\n==== Test of', name, 'OK.===>\n')

test.prediction <- function() {
	test.name <- 'Running prediction'
	start.test(test.name)
	set.seed(1)
	sim.dir <- tempfile()
	pred <- pop.predict(countries=c(528,218,450), present.year=2010,
				nr.traj = 3, verbose=FALSE, output.dir=sim.dir)
	s <- summary(pred)
	stopifnot(s$nr.traj == 3)
	stopifnot(s$nr.countries == 3)
	stopifnot(length(s$projection.years) == 18)
	test.ok(test.name)

	# aggregate
	test.name <- 'Running aggregation'
	start.test(test.name)
	aggr <- pop.aggregate(pred, c(900,904))
	stopifnot(nrow(aggr$countries) == 2)
	test.ok(test.name)
	
	# aggregate with user-defined groupings
	test.name <- 'Running aggregation with user-defined groupings'
	start.test(test.name)	
	UNlocs <- cbind(UNlocations, agcode_10=99)
	UNlocs[UNlocs$country_code %in% c(218,528), 'agcode_10'] <- 9000
	UNlocs <- rbind(UNlocs, 
				data.frame(name='my_aggregation', country_code=9000, reg_code=-1, reg_name='', area_code=-1, area_name='', 
							location_type=10, agcode_10=-1))
	locfile <- tempfile()
	write.table(UNlocs, file=locfile, sep='\t')
	aggr1 <- pop.aggregate(pred, 9000, my.location.file=locfile)
	unlink(locfile)
	stopifnot(length(aggr1$aggregated.countries[['9000']]) == 2)
	stopifnot(all(is.element(c(218, 528), aggr1$aggregated.countries[['9000']]))) 
	test.ok(test.name)
	unlink(sim.dir, recursive=TRUE)
}

test.expressions <- function() {
	test.name <- 'Population expressions'
	start.test(test.name)
	sim.dir <- tempfile()
	pred <- pop.predict(countries=c(528,218,450, 242, 458), present.year=2010,
				nr.traj = 3, verbose=FALSE, output.dir=sim.dir)
	filename <- tempfile()
	png(filename=filename)
	pop.trajectories.plot(pred, expression='P528_F[1]')
	pop.byage.plot(pred, expression='P528_F{} / PNL_M{}')
	dev.off()
	size <- file.info(filename)['size']
	unlink(filename)
	stopifnot(size > 0)
	pop.trajectories.table(pred, expression='P242 / (P528 + P218 + P450 + P242 + P458)')
	
	write.pop.projection.summary(pred, expression="PXXX[1] / PXXX", output.dir=sim.dir)
	t <- read.table(file.path(sim.dir, 'projection_summary_expression.csv'), sep=',', header=TRUE)
	stopifnot(all(dim(t) == c(25,22)))
	
	write.pop.projection.summary(pred, expression="GXXX[1:10]", output.dir=sim.dir) # migration
	t <- read.table(file.path(sim.dir, 'projection_summary_expression.csv'), sep=',', header=TRUE)
	stopifnot(all(dim(t) == c(25,22)))
	
	aggr <- pop.aggregate(pred, 900)
	pop.trajectories.table(pred, expression='P528_M / P900')
	write.pop.projection.summary(pred, expression="PXXX_M / P900_M", output.dir=sim.dir)
	t <- read.table(file.path(sim.dir, 'projection_summary_expression.csv'), sep=',', header=TRUE)
	stopifnot(all(dim(t) == c(25,22)))
	
	test.ok(test.name)
	unlink(sim.dir, recursive=TRUE)
}

test.expressions.with.VE <- function(map=TRUE) {
	test.name <- 'Expressions with vital events'
	start.test(test.name)
	sim.dir <- tempfile()
	pred <- pop.predict(countries=c(528, 218), present.year=2010,
				nr.traj = 3, verbose=FALSE, output.dir=sim.dir, keep.vital.events=TRUE)
	filename <- tempfile()
	png(filename=filename)
	pop.trajectories.plot(pred, expression='F528_F[10]')
	dev.off()
	size <- file.info(filename)['size']
	unlink(filename)
	stopifnot(size > 0)
	pop.trajectories.table(pred, expression='D528 / (DNLD + D218)')
	
	write.pop.projection.summary(pred, expression="BXXX[5] / BXXX", output.dir=sim.dir)
	t <- read.table(file.path(sim.dir, 'projection_summary_expression.csv'), sep=',', header=TRUE)
	stopifnot(all(dim(t) == c(10,22))) # 2 countries 5 rows each
	
	write.pop.projection.summary(pred, expression="pop.combine(BXXX[5], BXXX, '/')", output.dir=sim.dir)
	t <- read.table(file.path(sim.dir, 'projection_summary_expression.csv'), sep=',', header=TRUE)
	stopifnot(all(dim(t) == c(10,22))) # 2 countries 5 rows each
	
	t <- pop.byage.table(pred, expression='M528_M{}')
	stopifnot(all(dim(t) == c(27,5)))
	write.pop.projection.summary(pred, expression="SXXX_M{0}", output.dir=sim.dir)
	t <- read.table(file.path(sim.dir, 'projection_summary_expression.csv'), sep=',', header=TRUE)
	stopifnot(all(dim(t) == c(10,22)))
		
	filename <- tempfile()
	png(filename=filename)
	pop.byage.plot(pred, expression='log(QEC_M{age.index01(27)})', year=2050)
	pop.byage.plot(pred, expression='log(QECU_M{age.index01(21)})', year=2008)
	pop.byage.plot(pred, expression='M218_F{age.index05(27)}', year=2050)
	pop.trajectories.plot(pred, expression="pop.apply(P528_F{4:10}, gmedian, cats=seq(15, by=5, length=8))")
	pop.byage.plot(pred, expression="pop.combine(M218_F{age.index05(27)}, P218, '/')", year=2050)
	pop.byage.plot(pred, expression="pop.combine(M218_F{age.index05(27)}, P218, '/')", year=1970)
	pop.trajectories.plot(pred, expression="pop.combine(B218 - D218, G218, '+', split.along='traj')")
	if(map) pop.map(pred, expression="pop.combine(PXXX_M, P528, '/', split.along='country')", year=1980)
	dev.off()
	size <- file.info(filename)['size']
	unlink(filename)
	stopifnot(size > 0)
	write.pop.projection.summary(pred, expression="QXXX_F[0]", output.dir=sim.dir)
	t <- read.table(file.path(sim.dir, 'projection_summary_expression.csv'), sep=',', header=TRUE)
	stopifnot(all(dim(t) == c(10,22)))
	
	filename <- tempfile()
	png(filename=filename)
	pop.pyramid(pred, 218, indicator='D')
	pop.pyramid(pred, 218, indicator='B')
	dev.off()
	size <- file.info(filename)['size']
	unlink(filename)
	stopifnot(size > 0)
	test.ok(test.name)
	unlink(sim.dir, recursive=TRUE)
}
