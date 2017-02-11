library(bayesPop)
source('test_functions.R')

CRAN <- TRUE
warn <- options('warn')
options(warn=2)
test.expressions()

if(!CRAN) {
# longer tests 
	test.prediction()
	test.prediction.with.prob.migration()
	test.expressions.with.VE(map=FALSE)
	test.regional.aggregation()
	test.life.table()
}
options(warn=warn$warn)