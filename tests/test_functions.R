library(bayesPop)
library(data.table)

start.test <- function(name) cat('\n<=== Starting test of', name,'====\n')
test.ok <- function(name) cat('\n==== Test of', name, 'OK.===>\n')

test.prediction <- function(parallel = FALSE) {
	test.name <- paste('Running prediction', if(parallel) 'in parallel' else '')
	start.test(test.name)
	set.seed(1)
	sim.dir <- tempfile()	
	pred <- pop.predict(countries=c(528,218,450), 
				nr.traj = 3, verbose=FALSE, output.dir=sim.dir, parallel = parallel)
	s <- summary(pred)
	stopifnot(s$nr.traj == 3)
	stopifnot(s$nr.countries == 3)
	stopifnot(length(s$projection.years) == 16)
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
	UNlocs <- cbind(UNlocations[,1:7], agcode_10=99)
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
	
	test.name <- paste('Running prediction with 1 trajectory', if(parallel) 'in parallel' else '')
	start.test(test.name)
	pred <- pop.predict(countries=528, keep.vital.events=TRUE,
				nr.traj = 1, verbose=FALSE, output.dir=sim.dir, replace.output=TRUE, end.year=2040,
				parallel = parallel)
	# check that it took the median TFR and not high or low
	tfr <- get.pop("F528", pred)
	# extract data from wpp2019
	data(tfrprojMed, package = "wpp2019")
	tfr.should.be <- round(subset(tfrprojMed, country_code == 528)[3:6], 2)
	stopifnot(all(round(tfr[1,1,2:dim(tfr)[3],1],2) == tfr.should.be))
	
	# check that writing summary is OK
	write.pop.projection.summary(pred, what=c("popsexage", "popage"), output.dir=sim.dir, include.means = TRUE)
	t <- read.table(file.path(sim.dir, 'projection_summary_popsexage.csv'), sep=',', header=TRUE)
	s <- summary(pred, country = 528)
	stopifnot(round(s$projections[1,1]) == sum(t[t$variant == "median", "X2020"]))
	stopifnot('mean' %in% t$variant)
	
	t <- read.table(file.path(sim.dir, 'projection_summary_popage.csv'), sep=',', header=TRUE)
	s <- pop.byage.table(pred, country = 528, year = 2040)
	stopifnot(round(s[2,1]) == t[t$variant == "median" & t$age == "5-9", "X2040"])
	stopifnot('mean' %in% t$variant)
	
	pred <- pop.predict(countries=528, keep.vital.events=TRUE,
				nr.traj = 3, verbose=FALSE, output.dir=sim.dir, replace.output=TRUE, end.year=2040,
				inputs=list(tfr.file='median_', e0M.file='median_'), parallel = parallel)
	tfr <- get.pop("F528", pred)
	stopifnot(all(round(tfr[1,1,2:dim(tfr)[3],1],2) == tfr.should.be))
	stopifnot(pred$nr.traj==1) # even though we want 3 trajectories, only one is available, because we take TFR median

	# check that writing summary is OK
	write.pop.projection.summary(pred, what=c("popsexage"), output.dir=sim.dir)
	t <- read.table(file.path(sim.dir, 'projection_summary_popsexage.csv'), sep=',', header=TRUE)
	s <- summary(pred, country = 528)
	stopifnot(round(s$projections[1,1]) == sum(t[t$variant == "median", "X2020"]))
	
	test.ok(test.name)
	unlink(sim.dir, recursive=TRUE)
}

test.expressions <- function(parallel = FALSE) {
	test.name <- paste('Population expressions', if(parallel) 'in parallel' else '')
	start.test(test.name)
	sim.dir <- tempfile()
	pred <- pop.predict(countries=c(528,218,450, 242, 458), nr.traj = 3, verbose=FALSE, 
	                    output.dir=sim.dir, parallel = parallel)
	filename <- tempfile()
	png(filename=filename)
	pop.trajectories.plot(pred, expression='P528_F[1]')
	pop.byage.plot(pred, expression='P528_F{} / PNL_M{}')
	dev.off()
	size <- file.info(filename)['size']
	unlink(filename)
	stopifnot(size > 0)
	pop.trajectories.table(pred, expression='P242 / (P528 + P218 + P450 + P242 + P458)')
	
	write.pop.projection.summary(pred, expression="PXXX[1] / PXXX", output.dir=sim.dir, include.observed=TRUE)
	t <- read.table(file.path(sim.dir, 'projection_summary_expression.csv'), sep=',', header=TRUE)
	stopifnot(all(dim(t) == c(25,34)))
	
	write.pop.projection.summary(pred, expression="GXXX[1:10]", output.dir=sim.dir, include.observed=TRUE) # migration
	t <- read.table(file.path(sim.dir, 'projection_summary_expression.csv'), sep=',', header=TRUE)
	stopifnot(all(dim(t) == c(25,34)))
	
	aggr <- pop.aggregate(pred, 900)
	pop.trajectories.table(pred, expression='P528_M / P900')
	write.pop.projection.summary(pred, expression="PXXX_M / P900_M", output.dir=sim.dir)
	t <- read.table(file.path(sim.dir, 'projection_summary_expression.csv'), sep=',', header=TRUE)
	stopifnot(all(dim(t) == c(25,20)))
	
	test.ok(test.name)
	unlink(sim.dir, recursive=TRUE)
}

test.expressions.with.VE <- function(map=TRUE, parallel = FALSE) {
	test.name <- paste('Expressions with vital events', if(parallel) 'in parallel' else '')
	start.test(test.name)
	sim.dir <- tempfile()
	pred <- pop.predict(countries=c(528, 218), nr.traj = 3, verbose=FALSE, output.dir=sim.dir, 
	                    keep.vital.events=TRUE, parallel = parallel)

	filename <- tempfile()
	png(filename=filename)
	pop.trajectories.plot(pred, expression='F528_F[10]')
	dev.off()
	size <- file.info(filename)['size']
	unlink(filename)
	stopifnot(size > 0)
	pop.trajectories.table(pred, expression='D528 / (DNLD + D218)')
	pop.trajectories.table(pred, expression='F528_F[4]/(R528_F[4]/100)') # gives TFR
	pop.trajectories.table(pred, expression=mac.expression("ECU")) # MAC
	
	write.pop.projection.summary(pred, expression="FXXX", output.dir=sim.dir, include.observed=TRUE) # TFR
	tb <- read.table(file.path(sim.dir, 'projection_summary_expression.csv'), sep=',', header=TRUE)
	stopifnot(all(tb[, 5:ncol(tb)] < 7))
	
	write.pop.projection.summary(pred, expression="BXXX[5] / BXXX", output.dir=sim.dir)
	t1 <- read.table(file.path(sim.dir, 'projection_summary_expression.csv'), sep=',', header=TRUE)
	stopifnot(all(dim(t1) == c(10,20))) # 2 countries 5 rows each
	
	write.pop.projection.summary(pred, expression="pop.combine(BXXX[5], BXXX, '/')", output.dir=sim.dir)
	t2 <- read.table(file.path(sim.dir, 'projection_summary_expression.csv'), sep=',', header=TRUE)
	stopifnot(identical(t1, t2))
	
	t <- pop.byage.table(pred, expression='M528_M{}')
	stopifnot(all(dim(t) == c(27,5)))
	write.pop.projection.summary(pred, expression="SXXX_M{0}", output.dir=sim.dir)
	t <- read.table(file.path(sim.dir, 'projection_summary_expression.csv'), sep=',', header=TRUE)
	stopifnot(all(dim(t) == c(10,20)))
		
	filename <- tempfile()
	png(filename=filename)
	pop.byage.plot(pred, expression='log(QEC_M{age.index01(27)})', year=2050)
	pop.byage.plot(pred, expression='log(QECU_M{age.index01(21)})', year=2018)
	pop.byage.plot(pred, expression='M218_F{age.index05(27)}', year=2050)
	pop.trajectories.plot(pred, expression="pop.apply(P528_F{4:10}, gmedian, cats=seq(15, by=5, length=8))")
	pop.byage.plot(pred, expression="pop.combine(M218_F{age.index05(27)}, P218, '/')", year=2050)
	pop.byage.plot(pred, expression="pop.combine(M218_F{age.index05(27)}, P218, '/')", year=1990)
	pop.trajectories.plot(pred, expression="pop.combine(B218 - D218, G218, '+', split.along='traj')")
	pop.trajectories.plot(pred, expression="pop.combine(G218, P218, '/', split.along='traj')")
	if(map) {
	    pop.map(pred, expression="pop.combine(PXXX_M, P528, '/', split.along='country')", year=1980)
	    pop.ggmap(pred, expression="pop.combine(PXXX_M, P528, '/', split.along='country')", year=2050)
	}
	dev.off()
	size <- file.info(filename)['size']
	unlink(filename)
	stopifnot(size > 0)
	write.pop.projection.summary(pred, expression="QXXX_F[0]", output.dir=sim.dir)
	t <- read.table(file.path(sim.dir, 'projection_summary_expression.csv'), sep=',', header=TRUE)
	stopifnot(all(dim(t) == c(10,20)))
	
	filename <- tempfile()
	png(filename=filename)
	pop.pyramid(pred, 218)
	pop.pyramid(pred, 528, year=2052, proportion=TRUE)
	pop.pyramid(pred, 218, indicator='D')
	pop.pyramid(pred, 218, indicator='B', year=2100)
	pop.pyramid(pred, 218, indicator='D', proportion=TRUE, year=2100)
	pop.trajectories.pyramid(pred, 218, proportion=TRUE)
	pop.trajectories.pyramid(pred, 528, year=2052, proportion=TRUE)
	pop.trajectories.pyramid(pred, 218, indicator='B')
	pop.trajectories.pyramid(pred, 218, indicator='D', year=2100)
	pop.trajectories.pyramid(pred, 218, indicator='B', year=2100)
	pop.trajectories.pyramid(pred, 218, indicator='B', proportion=TRUE, year=2100)
	dev.off()
	size <- file.info(filename)['size']
	unlink(filename)
	stopifnot(size > 0)
	
	write.pop.projection.summary(pred, output.dir=sim.dir)
	t <- read.table(file.path(sim.dir, 'projection_summary_tpop.csv'), sep=',', header=TRUE)
	stopifnot(all(dim(t) == c(10,20)))
	t <- read.table(file.path(sim.dir, 'projection_summary_asfrage.csv'), sep=',', header=TRUE)
	stopifnot(all(dim(t) == c(70,21)))
	
	write.pop.projection.summary(pred, output.dir=sim.dir, include.observed=TRUE)
	t <- read.table(file.path(sim.dir, 'projection_summary_tpop.csv'), sep=',', header=TRUE)
	stopifnot(all(dim(t) == c(10,34)))
	
	test.ok(test.name)
	unlink(sim.dir, recursive=TRUE)
}

test.prediction.with.prob.migration <- function(parallel = FALSE) {
	test.name <- paste('Running prediction with probabilistic migration', if(parallel) 'in parallel' else '')
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
	pred <- pop.predict(countries=c(528,218), end.year=2033,
				verbose=FALSE, output.dir=sim.dir, keep.vital.events=TRUE, replace.output=TRUE,
				inputs=list(migMtraj=migMfile, migFtraj=migFfile), parallel = parallel)
	s <- summary(pred)
	# should have 3 trajectories because TFR has 3
	stopifnot(s$nr.traj == 3)
	stopifnot(s$nr.countries == 2)
	stopifnot(length(s$projection.years) == 3)
	mgr <- get.pop("G528", pred)
	stopifnot(dim(mgr)[4] == 3) # migration is re-sampled to 3 trajs
	
	write.migration(nr.traj=5)
	pred <- pop.predict(countries=c(528,218), end.year=2033,
				verbose=FALSE, output.dir=sim.dir, keep.vital.events=TRUE, replace.output=TRUE,
				inputs=list(migMtraj=migMfile, migFtraj=migFfile), parallel = parallel)
	stopifnot(pred$nr.traj == 5)
	stopifnot(dim(get.pop("G218", pred))[4] == 5)
	
	pred <- pop.predict(countries=c(528,218), end.year=2033,
				verbose=FALSE, output.dir=sim.dir, keep.vital.events=TRUE, replace.output=TRUE,
				inputs=list(migMtraj=migMfile, migFtraj=migFfile), nr.traj=1, parallel = parallel)
	stopifnot(pred$nr.traj == 1)
	stopifnot(dim(get.pop("G218", pred))[4] == 1)
	
	write.migration(nr.traj=1)
	pred <- pop.predict(countries=c(528,218), end.year=2033,
				verbose=FALSE, output.dir=sim.dir, keep.vital.events=TRUE, replace.output=TRUE, parallel = parallel,
				inputs=list(migMtraj=migMfile)) # female is taken the default one (only works if male has 1 trajectory)
	stopifnot(pred$nr.traj == 3)
	stopifnot(dim(get.pop("G218_M", pred))[4] == 3)
	
	pred <- pop.predict(countries=c(528,218), end.year=2033,
				verbose=FALSE, output.dir=sim.dir, keep.vital.events=TRUE, replace.output=TRUE, parallel = parallel,
				inputs=list(migFtraj=migFfile)) # male is taken the default one (only works if male has 1 trajectory)
	stopifnot(pred$nr.traj == 3)
	stopifnot(dim(get.pop("G218_F", pred))[4] == 3)

	
	test.ok(test.name)
	unlink(sim.dir, recursive=TRUE)
	unlink(migMfile)
	unlink(migFfile)
}

test.regional.aggregation <-function(parallel = FALSE) {
	test.name <- paste('Regional aggregation', if(parallel) 'in parallel' else '')
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
	pred <- pop.predict(output.dir=sim.dir.pop, verbose=TRUE, parallel = parallel,
    			inputs = list(tfr.sim.dir=sim.dir.tfr, 
    			              e0F.sim.dir=sim.dir.e0, e0M.sim.dir='joint_'))
    options(warn=warn$warn)
	stopifnot(pred$nr.traj==20)
	aggr <- pop.aggregate(pred, regions=regions, input.type="region", verbose=TRUE)
	stopifnot(setequal(aggr$countries$code, regions))
	test.ok(test.name)
	unlink(sim.dir.tfr, recursive=TRUE)
	unlink(sim.dir.e0, recursive=TRUE)
	unlink(sim.dir.pop, recursive=TRUE)
}

test.life.table <- function(parallel = FALSE){
	test.name <- paste('Life Tables', if(parallel) 'in parallel' else '')
	start.test(test.name)
	sim.dir <- tempfile()
	# this is the Example from LifeTableMx
	pred <- pop.predict(countries="Ecuador", output.dir=sim.dir, wpp.year=2015, parallel = parallel,
    			present.year=2015, keep.vital.events=TRUE, fixed.mx=TRUE, fixed.pasfr=TRUE)
	# get male mortality rates from current year for age groups 0-1, 1-4, 5-9, ...
	mx <- pop.byage.table(pred, expression="MEC_M{c(-1,0,2:27)}")[,1]
	LT <- LifeTableMx(mx)
	stopifnot(all(dim(LT) == c(28,10)))
	stopifnot(!any(is.na(LT)))
	mxf <- pop.byage.table(pred, expression="MEC_F{age.index01(27)}", year=2020)[,1]
	LT <- LifeTableMx(mxf, sex="Female", include01=FALSE)
	stopifnot(all(dim(LT) == c(27,10)))
	stopifnot(!any(is.na(LT)))
	sx1 <- as.double(LifeTableMxCol(mx, 'sx', age05=c(FALSE, FALSE, TRUE)))
	sx2 <- get.pop.exba("SEC_M{1:27}", pred, observed=TRUE)
	sx2 <- as.double(sx2[,ncol(sx2)])
	sxpred <- get.pop.exba("SEC_M{1:27}", pred, observed=FALSE)
	sx3 <- as.double(sxpred[,1,1])
	stopifnot(all.equal(sx1[1:19], sx2[1:19]) && all.equal(sx2[1:19], sx3[1:19]))
	sx4 <- as.double(LifeTableMxCol(pop.byage.table(pred, expression="MEC_M{age.index01(27)}", year=2053)[,1], 'sx', age05=c(FALSE, FALSE, TRUE)))
	sx5 <- as.double(sxpred[,"2053",1])
	stopifnot(all.equal(sx4, sx5))
	test.ok(test.name)
	unlink(sim.dir, recursive=TRUE)
}

test.adjustment <- function() {
    test.name <- 'Adjustments'
    start.test(test.name)
    sim.dir <- file.path(find.package("bayesPop"), "ex-data", "Pop")
    pred <- get.pop.prediction(sim.dir)
    med <- pop.trajectories.table(pred, "Ecuador")[,"median"]
    adj.med <- pop.trajectories.table(pred, "Ecuador", adjust=TRUE)[,"median"]
    # from wpp2017 popproj dataset:
    #data(popproj, package="wpp2017")
    #should.be <- unlist(subset(popproj, name=="Ecuador")[,c("2080", "2090", "2100")])
    should.be <- c(24876.80, 24725.84, 24320.58) 
    stopifnot(all.equal(adj.med[c("2080", "2090", "2100")], should.be, 
                        tolerance = 0.01, check.attributes = FALSE))
    test.ok(test.name)
}

test.subnat <- function(mig.age.method = "rc") {
  test.name <- "Subnational projections with national TFR"
  start.test(test.name)
  data.dir <- file.path(find.package("bayesPop"), "extdata")
  # Use national data for tfr and e0
  sim.dir <- tempfile()
  pred <- pop.predict.subnat(output.dir = sim.dir,
                             locations = file.path(data.dir, "CANlocations.txt"),
                             inputs = list(popM = file.path(data.dir, "CANpopM.txt"),
                                           popF = file.path(data.dir, "CANpopF.txt"),
                                           patterns = file.path(data.dir, "CANpatterns.txt")
                                           ),
                             mig.age.method = mig.age.method)
  ct <- get.countries.table(pred)
  stopifnot(nrow(ct) == 13) # 13 sub-regions of Canada
  stopifnot(dim(get.pop("P658", pred))[3] == 9) # projection until 2060
  stopifnot(pred$inputs$mig.age.method == if(is.null(mig.age.method)) "rc" else mig.age.method)
  
  aggr <- pop.aggregate.subnat(pred, regions = 124, 
                locations = file.path(data.dir, "CANlocations.txt"))
  ct <- get.countries.table(aggr)
  stopifnot(nrow(ct) == 1)
  stopifnot(dim(get.pop("P124", aggr))[3] == 9) # projection until 2060
  unlink(sim.dir, recursive=TRUE)
  test.ok(test.name)
}

test.subnat.with.subnat.tfr.e0 <- function() {
  test.name <- "Subnational projections with subnational TFR and e0"
  start.test(test.name)
  # TFR projections
  data.dir <- file.path(find.package("bayesPop"), "extdata")
  my.subtfr.file <- file.path(find.package("bayesTFR"), 'extdata', 'subnational_tfr_template.txt')
  tfr.nat.dir <- file.path(find.package("bayesTFR"), "ex-data", "bayesTFR.output")
  tfr.reg.dir <- tempfile()
  tfr.preds <- tfr.predict.subnat(124, my.tfr.file = my.subtfr.file,
                                  sim.dir = tfr.nat.dir, output.dir = tfr.reg.dir, start.year = 2013)
  
  # e0 projections
  my.sube0.file <- file.path(find.package("bayesLife"), 'extdata', 'subnational_e0_template.txt')
  e0.nat.dir <- file.path(find.package("bayesLife"), "ex-data", "bayesLife.output")
  e0.reg.dir <- tempfile()
  e0.preds <- e0.predict.subnat(124, my.e0.file=my.sube0.file,
                             sim.dir=e0.nat.dir, output.dir=e0.reg.dir, 
                             predict.jmale = TRUE, my.e0M.file = my.sube0.file)
  
  # Pop projections
  sim.dir <- tempfile()
  pred <- pop.predict.subnat(output.dir = sim.dir,
                             locations = file.path(data.dir, "CANlocations.txt"),
                             inputs = list(popM = file.path(data.dir, "CANpopM.txt"),
                                           popF = file.path(data.dir, "CANpopF.txt"),
                                           patterns = file.path(data.dir, "CANpatterns.txt"),
                                           tfr.sim.dir = file.path(tfr.reg.dir, "subnat", "c124"),
                                           e0F.sim.dir = file.path(e0.reg.dir, "subnat_ar1", "c124"),
                                           e0M.sim.dir = "joint_"
                             ), verbose = FALSE)
  stopifnot(all(dim(get.pop("P658", pred)) == c(1,1,9,30))) # 30 trajectories because TFR example has 30 national trajs. 
  aggr <- pop.aggregate.subnat(pred, regions = 124, 
                               locations = file.path(data.dir, "CANlocations.txt"))
  tbl <- pop.trajectories.table(aggr, "Canada")
  stopifnot(nrow(tbl) == 28)
  stopifnot(summary(aggr)$nr.traj == 30)

  unlink(sim.dir, recursive = TRUE)
  unlink(tfr.reg.dir, recursive = TRUE)
  unlink(e0.reg.dir, recursive = TRUE)
  test.ok(test.name)
}

test.prediction.with.patterns <- function() {
    test.name <- paste('Running prediction with different patterns')
    start.test(test.name)
    data("vwBaseYear2019")
    patterns <- vwBaseYear2019
    cntries <- c(528,218,450)
    patterns <- subset(patterns, country_code %in% cntries)
    #patterns$AgeMortProjMethod1[] <- "LC"
    patterns$AgeMortProjMethod1[] <- "modPMD"
    patterns$AgeMortProjMethod2[] <- ""
    patterns$SmoothDFLatestAgeMortalityPattern <- c(4, 0, 19)
    patterns$SmoothLatestAgeMortalityPattern[] <- 1
    patterns$LatestAgeMortalityPattern <- c("c(-2, 3)", -2, "c(1, 1)") # the last one should give a warning
    pattern.file <- tempfile()
    write.table(patterns, file = pattern.file, sep = "\t", row.names = FALSE)
    set.seed(1)
    sim.dir <- tempfile()	
    pred <- pop.predict(countries = cntries, 
                        nr.traj = 3, verbose=TRUE, output.dir=sim.dir,
                        inputs = list(patterns = pattern.file),
                        keep.vital.events = TRUE, lc.for.all = FALSE, replace.output = TRUE
                        )
    # there isn't really a way to test that the patterns were used, 
    # apart from exploring the resulting mx, e.g.
    # pop.byage.plot(pred, expression = "log(M218_M{age.index01(27)})", year = 2023)
    # pop.byage.plot(pred, expression = "log(M218_M{age.index01(27)})", add = TRUE)
    s <- summary(pred)
    stopifnot(s$nr.traj == 3)
    stopifnot(s$nr.countries == 3)
    stopifnot(length(s$projection.years) == 16)
    unlink(sim.dir, recursive = TRUE)
    unlink(pattern.file)
    test.ok(test.name)
}    

test.age.specific.migration <- function(){
    test.name <- paste('Generating age-specific migration')
    start.test(test.name)
    
    # Taken from Examples which are "dontrun"
    asmig <- age.specific.migration()
    stopifnot(all(dim(asmig$male) == c(4221, 33)))
    stopifnot(all(dim(asmig$female) == c(4221, 33)))
    
    # disaggregate WPP 2019 migration for all countries, one sex
    data(migration, package = "wpp2019")
    # assuming equal sex migration ratio
    asmig.all <- migration.totals2age(migration, scale = 0.5, method = "rc") 
    # result for the US in 2095-2100
    mig1sex.us <- subset(asmig.all, country_code == 840)[["2095-2100"]]
    stopifnot(length(mig1sex.us) == 21) 
    stopifnot(mig1sex.us[1] > 2*mig1sex.us[3] && 8*mig1sex.us[3] < mig1sex.us[5] && mig1sex.us[21] == 0)
    # check that the sum is half of the original total
    stopifnot(sum(mig1sex.us) == subset(migration, country_code == 840)[["2095-2100"]]/2)
    
    test.ok(test.name)
}

test.different.migration.methods <- function(wpp.year = 2024, annual = TRUE) {
    step <- if(annual) 1 else 5
    test.name <- paste0('Running different migration methods ', 'for a ', step, '-year simulation with wpp', wpp.year)
    start.test(test.name)
    set.seed(1)
    # create migration files with two countries and two trajectories
    migfile <- tempfile()
    sim.dir <- tempfile()
    time <- 5
    countries <- c(528,218)
    ncountries <- length(countries)
    start.years <- if(annual) list("2024" = 2023, "2022" = 2021, "2019" = 2015) else list("2024" = 2020, "2022" = 2020, "2019" = 2015)
    present.year <- start.years[[as.character(wpp.year)]]
    midyr <- if(annual) 1 else 3
    write.migration <- function(nr.traj, rates = FALSE) {
        nrows.country <- nr.traj*time
        mig <- data.frame(LocID=rep(c(528,218), each=nrows.country), 
                          Year=rep(seq(present.year + midyr, by=step, length=time), times=nr.traj*ncountries),
                          Trajectory=rep(rep(1:nr.traj, each=time), times=ncountries), 
                          Migration=0)
        if(!rates)
            mig$Migration <- rnorm(nrow(mig), mean=rep(c(3,0), each=nrows.country), sd=rep(c(2, 1), each=nrows.country))
        else mig$Migration <- rnorm(nrow(mig), mean=rep(c(0.005,0.0005), each=nrows.country), sd=rep(c(0.01, 0.001), each=nrows.country))
        write.csv(mig, file=migfile, row.names=FALSE)
    }
    # Test "fdm" with 1 trajectory
    nr.traj <- 1
    write.migration(nr.traj = nr.traj)
    pred <- pop.predict(countries=countries, end.year= present.year + time*step, present.year = present.year,
                        annual = annual, wpp.year = wpp.year, nr.traj = nr.traj,
                        verbose=FALSE, output.dir=sim.dir, keep.vital.events=TRUE, replace.output=TRUE,
                        inputs=list(migtraj=migfile), mig.age.method = "fdmnop")
    s <- summary(pred)
    stopifnot(s$nr.traj == nr.traj)
    stopifnot(s$nr.countries == 2)
    stopifnot(length(s$projection.years) == time)
    mgr <- get.pop("G528", pred)
    stopifnot(dim(mgr)[4] == nr.traj)

    # Test "auto" with 5 trajectories
    nr.traj <- 5
    write.migration(nr.traj = nr.traj)
    pred <- pop.predict(countries = countries, end.year = present.year + time*step, present.year = present.year,
                        annual = annual, wpp.year = wpp.year, nr.traj = nr.traj,
                        verbose=FALSE, output.dir=sim.dir, keep.vital.events=TRUE, replace.output=TRUE,
                        inputs=list(migtraj=migfile), mig.age.method = "auto")
    stopifnot(pred$nr.traj == nr.traj)
    stopifnot(dim(get.pop("G218", pred))[4] == nr.traj)
    stopifnot(pred$inputs$mig.age.method == if(wpp.year == 2019 & !annual) "residual" else "fdmp")

    # Test "fdmnop" counts with 5 trajectories
    nr.traj <- 5
    write.migration(nr.traj = nr.traj)
    pred <- pop.predict(countries = countries, end.year = present.year + time*step, present.year = present.year,
                        annual = annual, wpp.year = wpp.year, nr.traj = nr.traj,
                        verbose=FALSE, output.dir=sim.dir, keep.vital.events=TRUE, replace.output=TRUE,
                        inputs=list(migtraj=migfile), mig.age.method = "fdmnop")
    stopifnot(pred$nr.traj == nr.traj)
    stopifnot(dim(get.pop("G218", pred))[4] == nr.traj)
    stopifnot(pred$inputs$mig.age.method == "fdmnop")

    # Test "fdmp" rates with 5 trajectories
    nr.traj <- 5
    write.migration(nr.traj = nr.traj, rates = TRUE)
    pred <- pop.predict(countries = countries, end.year = present.year + time*step, present.year = present.year,
                        annual = annual, wpp.year = wpp.year, nr.traj = nr.traj,
                        verbose=FALSE, output.dir=sim.dir, keep.vital.events=TRUE, replace.output=TRUE,
                        inputs=list(migtraj=migfile), mig.age.method = "fdmp",
                        mig.is.rate = c(FALSE, TRUE))
    stopifnot(pred$nr.traj == nr.traj)
    stopifnot(dim(get.pop("G218", pred))[4] == nr.traj)
    stopifnot(pred$inputs$mig.age.method == "fdmp")

    #options(error=quote(dump.frames("last.dump", TRUE)))
    #load("last.dump.rda"); debugger()
    # Test "rc" rates with 1 trajectories with external RC proportions
    nr.traj <- 1
    write.migration(nr.traj = nr.traj, rates = TRUE)
    rc.schedules <- subset(DemoTools::mig_un_families, family == "Male Labor")
    pred <- pop.predict(countries = countries, end.year = present.year + time*step, present.year = present.year,
                        annual = annual, wpp.year = wpp.year, nr.traj = nr.traj,
                        verbose=FALSE, output.dir=sim.dir, keep.vital.events=TRUE, replace.output=TRUE,
                        inputs=list(migtraj=migfile), mig.age.method = "rc", 
                        mig.rc.fam = rc.schedules,
                        mig.is.rate = c(FALSE, TRUE))
    stopifnot(pred$nr.traj == nr.traj)
    stopifnot(dim(get.pop("G218", pred))[4] == nr.traj)
    stopifnot(pred$inputs$mig.age.method == "rc")
    # check the migration sex ratio, if it is equal to the one we passed via the mig.rc.fam argument
    sex.ratio.in <- sum(subset(rc.schedules, mig_sign == "Inmigration" & sex == "Male" & age < 101)$prop)/
        sum(subset(rc.schedules, mig_sign == "Inmigration" & age < 101)$prop)
    sex.ratio.out <- sum(subset(rc.schedules, mig_sign == "Emigration" & sex == "Male" & age < 101)$prop)/
        sum(subset(rc.schedules, mig_sign == "Emigration" & age < 101)$prop)
    should.be.sr <- ifelse(get.pop.ex("G218", pred)[-1] >= 0, sex.ratio.in, sex.ratio.out)
    stopifnot(all.equal(get.pop.ex("G218_M/G218", pred)[-1], should.be.sr, tolerance = if(annual) 1e-6 else 1e-3)) # for 5-year sim, the should be ratios are a little off
    test.ok(test.name)
    unlink(sim.dir, recursive=TRUE)
    unlink(migfile)
}

test.probabilistic.fdmp <- function(wpp.year = 2024, annual = TRUE) {
    step <- if(annual) 1 else 5
    test.name <- paste0('Running probabilistic FDMp ', 'for a ', step, '-year simulation with wpp', wpp.year)
    start.test(test.name)
    set.seed(1)
    # create migration files with two countries and two trajectories
    migfile <- tempfile()
    migfdmtraj <- tempfile()
    sim.dir <- tempfile()
    time <- 5
    country <- 250
    ages <- bayesPop:::ages.all(annual, observed = TRUE)
    nages <- length(ages)
    age.labels <- bayesPop:::get.age.labels(ages, single.year = annual, last.open = TRUE)
    start.years <- if(annual) list("2024" = 2023, "2022" = 2021, "2019" = 2015) else list("2024" = 2020, "2022" = 2020, "2019" = 2015)
    present.year <- start.years[[as.character(wpp.year)]]
    midyr <- if(annual) 1 else 3
    write.migration <- function(nr.traj, nr.traj.fdm, rates = TRUE) {
        # assemble trajectories of total migration
        nrows.country <- nr.traj*time
        mig <- data.frame(LocID=rep(country, each=nrows.country), 
                          Year=rep(seq(present.year + midyr, by=step, length=time), times=nr.traj),
                          Trajectory=rep(1:nr.traj, each=time), 
                          Migration=0)
        if(!rates)
            mig$Migration <- rnorm(nrow(mig), mean=rep(10, each=nrows.country), sd=rep(5, each=nrows.country))
        else mig$Migration <- rnorm(nrow(mig), mean=rep(0.01, each=nrows.country), sd=rep(0.01, each=nrows.country))
        # assemble FDM trajectories of age-specific RC curves
        nrows.country <- nr.traj.fdm*nages*2
        fdmtraj <- data.table(LocID=rep(country, each=nrows.country), 
                            Trajectory=rep(rep(1:nr.traj.fdm, each=nages), times = 2), 
                            Age=rep(age.labels, times = 2*nr.traj.fdm),
                            Parameter = rep(c("in", "out"), each = nr.traj.fdm*nages),
                            Value=rep(rcastro.schedule(annual), 2))
        fdmtraj$Value <- rnorm(nrow(fdmtraj), mean=fdmtraj$Value, sd=rep(0.001, each=nrows.country))
        fdmtraj[, Value := Value/sum(Value), by = c("Trajectory", "Parameter")]
        write.csv(mig, file=migfile, row.names=FALSE)
        fwrite(fdmtraj, file=migfdmtraj)
    }
    # Test "fdmp" with 5 mig trajectories and 4 RC-in/out trajectories
    nr.traj <- 5
    nr.traj.fdm <- 4
    write.migration(nr.traj = nr.traj, nr.traj.fdm = nr.traj.fdm, rates = FALSE)
    pred <- pop.predict(countries=country, end.year= present.year + time*step, present.year = present.year,
                        annual = annual, wpp.year = wpp.year, nr.traj = nr.traj,
                        verbose=FALSE, output.dir=sim.dir, keep.vital.events=TRUE, replace.output=TRUE,
                        inputs=list(migtraj=migfile, migFDMtraj = migfdmtraj), mig.age.method = "fdmp")

    s <- summary(pred)
    stopifnot(s$nr.traj == nr.traj)
    stopifnot(s$nr.countries == 1)
    stopifnot(length(s$projection.years) == time)
    mgr <- get.pop("G250", pred)
    stopifnot(dim(mgr)[4] == nr.traj) 
    stopifnot(nrow(pred$inputs$migFDMpred) == nr.traj.fdm * 2 * nages)
    
    nr.traj.fdm <- 6
    write.migration(nr.traj = nr.traj, nr.traj.fdm = nr.traj.fdm, rates = TRUE)
    pred <- pop.predict(countries=country, end.year= present.year + (time+1)*step, present.year = present.year,
                        annual = annual, wpp.year = wpp.year, nr.traj = nr.traj,
                        verbose=FALSE, output.dir=sim.dir, keep.vital.events=TRUE, replace.output=TRUE,
                        inputs=list(migtraj=migfile, migFDMtraj = migfdmtraj), mig.age.method = "fdmp",
                        mig.is.rate = c(FALSE, TRUE)
                        )
    mgr <- get.pop("G250", pred)
    stopifnot(dim(mgr)[4] == nr.traj) 
    stopifnot(ncol(pred$inputs$migFDMpred) == 5)
    stopifnot(nrow(pred$inputs$migFDMpred) == nr.traj.fdm * 2 * nages)
    
    test.ok(test.name)
    unlink(sim.dir, recursive=TRUE)
    unlink(migfile)
    unlink(migfdmtraj)
}


#TODO: test project.pasfr function