library(bayesPop)
source('test_functions.R')

CRAN <- TRUE

test.expressions()
test.adjustment()
if(!CRAN) {
# longer tests 
    warn <- options('warn')
    options(warn=2)
	test.prediction()
	test.prediction.with.prob.migration()
	test.expressions.with.VE(map=FALSE)
	test.regional.aggregation()
	test.life.table()
	test.adjustment()
	test.subnat()
	options(warn=warn$warn)
}
