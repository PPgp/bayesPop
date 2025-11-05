library(bayesPop)
source('test_functions.R')

CRAN <- TRUE

test.expressions()

if(!CRAN) {
# longer tests 
    warn <- options('warn')
    options(warn=2)
    test.expressions(parallel = TRUE)
	test.prediction()
	test.prediction(parallel = TRUE)
	test.prediction.with.prob.migration()
	test.prediction.with.prob.migration(parallel = TRUE)
	test.expressions.with.VE(map=FALSE)
    test.regional.aggregation()
    test.regional.aggregation(parallel = TRUE)
    test.life.table()
    test.life.table(parallel = TRUE)
    test.adjustment()
    test.subnat()
    options(warn=warn$warn)
    test.subnat(mig.age.method = NULL) # default is fdmp (will generate a warning)
    test.subnat.with.subnat.tfr.e0()
    test.prediction.with.patterns()
	test.expressions.with.VE(map=TRUE, parallel = TRUE) # generates warnings
	test.age.specific.migration()
	for(yr in c(2024, 2022, 2019))
	    test.different.migration.methods(wpp.year = yr, annual = FALSE)
	for(yr in c(2024, 2022))
	    test.different.migration.methods(wpp.year = yr)
	test.probabilistic.fdmp()
	test.probabilistic.fdmp(annual = FALSE)
}
