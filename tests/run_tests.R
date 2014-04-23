library(bayesPop)
source('test_functions.R')

warn <- options('warn')
options(warn=2)
test.prediction()
test.prediction.with.prob.migration()
test.expressions()

test.expressions.with.VE(map=FALSE)
try(options(warn=warn), silent=TRUE)