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
	#test.expressions.with.VE(map=TRUE)
	test.expressions.with.VE(map=FALSE)
	test.expressions.with.VE(map=FALSE, parallel = TRUE)
	test.regional.aggregation()
	test.regional.aggregation(parallel = TRUE)
	test.life.table()
	test.life.table(parallel = TRUE)
	test.adjustment()
	test.subnat()
	options(warn=warn$warn)
	test.subnat.with.subnat.tfr.e0()
}
