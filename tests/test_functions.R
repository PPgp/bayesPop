
test.prediction <- function() {
	set.seed(1)
	sim.dir <- tempfile()
	pred <- pop.predict(countries=c(528,218,450), wpp.year=2010, present.year=2010,
				nr.traj = 10, verbose=FALSE, output.dir=sim.dir)
	s <- summary(pred)
	stopifnot(s$nr.traj == 10)
	stopifnot(s$nr.countries == 3)
	stopifnot(length(s$projection.years) == 18)
	
	# aggregate
	aggr <- pop.aggregate(pred, c(900,904))
	stopifnot(nrow(aggr$countries) == 2)
	unlink(sim.dir, recursive=TRUE)
}

test.expressions <- function() {
	sim.dir <- tempfile()
	pred <- pop.predict(countries=c(528,218,450, 242, 458), wpp.year=2010, present.year=2010,
				nr.traj = 10, verbose=FALSE, output.dir=sim.dir)
	filename <- tempfile()
	png(filename=filename)
	pop.trajectories.plot(pred, expression='P528_F[1]')
	dev.off()
	size <- file.info(filename)['size']
	unlink(filename)
	stopifnot(size > 0)
	pop.trajectories.table(pred, expression='P242 / (P528 + P218 + P450 + P242 + P458)')
	
	write.pop.projection.summary(pred, expression="PXXX[1] / PXXX", output.dir=sim.dir)
	t <- read.table(file.path(sim.dir, 'projection_summary_expression.csv'), sep=',', header=TRUE)
	stopifnot(all(dim(t) == c(25,22)))
	
	aggr <- pop.aggregate(pred, 900)
	pop.trajectories.table(pred, expression='P528_M / P900')
	write.pop.projection.summary(pred, expression="PXXX_M / P900_M", output.dir=sim.dir)
	t <- read.table(file.path(sim.dir, 'projection_summary_expression.csv'), sep=',', header=TRUE)
	stopifnot(all(dim(t) == c(25,22)))

	unlink(sim.dir, recursive=TRUE)
}