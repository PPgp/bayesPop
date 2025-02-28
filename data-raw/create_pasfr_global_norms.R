library(bayesPop)

# Create PASFR global norms dataset from two simulations, 5x5 and 1x1
# these can be toy simulation as it only needs the input object
sim.dir5 <- "~/bayespop/R/Pop/wpp2021/public_inputs/sim5x5_mig"
sim.dir1 <- "~/bayespop/R/Pop/wpp2021/public_inputs/sim1x1"

pred5 <- get.pop.prediction(sim.dir5)
pred1 <- get.pop.prediction(sim.dir1)

pasfr.glob.norms5 <- pred5$inputs$PASFRnorms['PasfrGlobalNorm']
pasfr.glob.norms1 <- pred1$inputs$PASFRnorms['PasfrGlobalNorm']

save(pasfr.glob.norms5, pasfr.glob.norms1, file = "pasfr_global_norms.rda")