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
	test.name <- 'Running aggregation (country type)'
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
	
	test.name <- 'Running prediction with 1 trajectory'
	start.test(test.name)
	pred <- pop.predict(countries=528, present.year=2010, keep.vital.events=TRUE,
				nr.traj = 1, verbose=FALSE, output.dir=sim.dir, replace.output=TRUE, end.year=2040)
	# check that it took the median TFR and not high or low
	tfr <- get.pop("F528", pred)
	stopifnot(all(round(tfr[1,1,,1],2) == c(1.75, 1.77, 1.79, 1.81, 1.82, 1.84, 1.85))) # WPP 2012 data
	
	pred <- pop.predict(countries=528, present.year=2010, keep.vital.events=TRUE,
				nr.traj = 3, verbose=FALSE, output.dir=sim.dir, replace.output=TRUE, end.year=2040,
				inputs=list(tfr.file='median_', e0M.file='median_'))
	tfr <- get.pop("F528", pred)
	stopifnot(all(round(tfr[1,1,,1],2) == c(1.75, 1.77, 1.79, 1.81, 1.82, 1.84, 1.85))) # WPP 2012 data
	stopifnot(pred$nr.traj==1) # even though we want 3 trajectories, only one available, because we take TFR median
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
	pop.trajectories.plot(pred, expression="pop.combine(G218, P218, '/', split.along='traj')")
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
	
	write.pop.projection.summary(pred, output.dir=sim.dir)
	t <- read.table(file.path(sim.dir, 'projection_summary_tpop.csv'), sep=',', header=TRUE)
	stopifnot(all(dim(t) == c(10,22)))
	t <- read.table(file.path(sim.dir, 'projection_summary_asfrage.csv'), sep=',', header=TRUE)
	stopifnot(all(dim(t) == c(70,23)))
	
	test.ok(test.name)
	unlink(sim.dir, recursive=TRUE)
}

test.prediction.with.prob.migration <- function() {
	test.name <- 'Running prediction with probabilstic migration'
	start.test(test.name)
	set.seed(1)
	# create migration files with two countries and two trajectories
	migMfile <- tempfile()
	migFfile <- tempfile()
	sim.dir <- tempfile()
	nr.traj <- 1
	time <- 5
	ncountries <- 2
	
	write.migration <- function(nr.traj) {
		nrows.country <- nr.traj*21*time
		mig <- data.frame(LocID=rep(c(528,218), each=nrows.country), Year=rep(rep(seq(2013, by=5, length=time), each=21), times=nr.traj*ncountries),
						Trajectory=rep(rep(1:nr.traj, each=21*time), times=ncountries), 
						Age=c(paste(seq(0,95,by=5), seq(4,99,by=5), sep='-'), '100+'), Migration=0)
		migM <- migF <- mig
		migM$Migration <- rnorm(nrow(mig), mean=rep(c(3,0), each=nrows.country), sd=rep(c(2, 1), each=nrows.country))
		migF$Migration <- rnorm(nrow(mig), mean=rep(c(3,0), each=nrows.country), sd=rep(c(2, 1), each=nrows.country))
		write.csv(migM, file=migMfile, row.names=FALSE)
		write.csv(migF, file=migFfile, row.names=FALSE)
	}
	write.migration(nr.traj=2)
	pred <- pop.predict(countries=c(528,218), present.year=2010, end.year=2033,
				verbose=FALSE, output.dir=sim.dir, keep.vital.events=TRUE, replace.output=TRUE,
				inputs=list(migMtraj=migMfile, migFtraj=migFfile))
	s <- summary(pred)
	# should have 3 trajectories because TFR has 3
	stopifnot(s$nr.traj == 3)
	stopifnot(s$nr.countries == 2)
	stopifnot(length(s$projection.years) == 5)
	mgr <- get.pop("G528", pred)
	stopifnot(dim(mgr)[4] == 3) # migration is re-sampled to 3 trajs
	
	write.migration(nr.traj=5)
	pred <- pop.predict(countries=c(528,218), present.year=2010, end.year=2033,
				verbose=FALSE, output.dir=sim.dir, keep.vital.events=TRUE, replace.output=TRUE,
				inputs=list(migMtraj=migMfile, migFtraj=migFfile))
	stopifnot(pred$nr.traj == 5)
	stopifnot(dim(get.pop("G218", pred))[4] == 5)
	
	pred <- pop.predict(countries=c(528,218), present.year=2010, end.year=2033,
				verbose=FALSE, output.dir=sim.dir, keep.vital.events=TRUE, replace.output=TRUE,
				inputs=list(migMtraj=migMfile, migFtraj=migFfile), nr.traj=1)
	stopifnot(pred$nr.traj == 1)
	stopifnot(dim(get.pop("G218", pred))[4] == 1)
	
	write.migration(nr.traj=1)
	pred <- pop.predict(countries=c(528,218), present.year=2010, end.year=2033,
				verbose=FALSE, output.dir=sim.dir, keep.vital.events=TRUE, replace.output=TRUE,
				inputs=list(migMtraj=migMfile)) # female is taken the default one (only works if male has 1 trajectory)
	stopifnot(pred$nr.traj == 3)
	stopifnot(dim(get.pop("G218_M", pred))[4] == 3)
	
	pred <- pop.predict(countries=c(528,218), present.year=2010, end.year=2033,
				verbose=FALSE, output.dir=sim.dir, keep.vital.events=TRUE, replace.output=TRUE,
				inputs=list(migFtraj=migFfile)) # male is taken the default one (only works if male has 1 trajectory)
	stopifnot(pred$nr.traj == 3)
	stopifnot(dim(get.pop("G218_F", pred))[4] == 3)

	
	test.ok(test.name)
	unlink(sim.dir, recursive=TRUE)
	unlink(migMfile)
	unlink(migFfile)
}

test.regional.aggregation <-function() {
	test.name <- 'Regional aggregation'
	start.test(test.name)
	regions <- c(900, 908, 904)
	sim.dir.tfr <- tempfile()
	sim.dir.e0 <- tempfile()
	sim.dir.pop <- tempfile()
	# Estimate TFR parameters
	run.tfr.mcmc(iter=25, nr.chains=1, output.dir=sim.dir.tfr, seed=1)
	run.tfr.mcmc.extra(sim.dir=sim.dir.tfr, countries=regions, burnin=0)
	# Predict TFR 
	tfr.predict(sim.dir=sim.dir.tfr, burnin=5, save.as.ascii=0, use.tfr3=FALSE)
	# Estimate e0 parameters 
	run.e0.mcmc(sex='F', iter=25, nr.chains=1, thin=1, output.dir=sim.dir.e0, seed=1)
	run.e0.mcmc.extra(sim.dir=sim.dir.e0, countries=regions, burnin=0)
	# Predict female and male e0
	warn <- options('warn'); options(warn=-1) # the joined estimation and pop projection has some warnings which can be ignored
	e0.predict(sim.dir=sim.dir.e0, burnin=5, save.as.ascii=0)
	# Population prediction
	pred <- pop.predict(output.dir=sim.dir.pop, verbose=TRUE, 
    			inputs = list(tfr.sim.dir=sim.dir.tfr, e0F.sim.dir=sim.dir.e0, e0M.sim.dir='joint_'))
    options(warn=warn$warn)
	stopifnot(pred$nr.traj==20)
	aggr <- pop.aggregate(pred, regions=regions, input.type="region", verbose=TRUE)
	stopifnot(setequal(aggr$countries$code, regions))
	unlink(sim.dir.tfr, recursive=TRUE)
	unlink(sim.dir.e0, recursive=TRUE)
	unlink(sim.dir.pop, recursive=TRUE)
}