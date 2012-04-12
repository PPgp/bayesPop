
test.prediction <- function() {
	set.seed(1)
	sim.dir <- tempfile()
	pred <- pop.predict(countries=c(528,218,450), wpp.year=2010, present.year=2010,
				nr.traj = 10, verbose=FALSE, output.dir=sim.dir)
	s <- summary(pred)
	stopifnot(s$nr.traj == 10)
	stopifnot(s$nr.countries == 3)
	stopifnot(length(s$projection.years) == 18)
	unlink(sim.dir, recursive=TRUE)
}