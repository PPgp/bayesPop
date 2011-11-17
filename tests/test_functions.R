
test.prediction <- function() {
	set.seed(1)
	sim.dir <- tempfile()
	pred <- pop.predict(countries=c(4,8,12), wpp.year=2010, present.year=2010,
				#inputs=list(#e0M.sim.dir='/Users/hana/bayespop/R/LE/3x100TM',
									#e0F.sim.dir='/Users/hana/bayespop/R/LE/3x100TF',
									#tfr.sim.dir='/Users/hana/bayespop/R/TFR/5x8000'
									#e0M.sim.dir='/Users/hana/bayespop/R/LE/2x30M',
									#e0F.sim.dir='/Users/hana/bayespop/R/LE/2x30F',
									#tfr.sim.dir='/Users/hana/bayespop/R/TFR/bayesTFR.output'
				#					),
				nr.traj = 10, verbose=TRUE, output.dir=sim.dir)
	summary(pred)
	unlink(sim.dir, recursive=TRUE)
}