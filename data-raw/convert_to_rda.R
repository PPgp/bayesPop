# convert ASCII into rda to be used in the bayesPop package

library(data.table)
e <- new.env()

dir <- "."
datasets <- c("MLTbx", paste0("vwBaseYear", c(2010, 2012, 2015, 2017, 2019, 2022, 2024)))
datasets <- "vwBaseYear2024"
#datasets <- "rcFDM"
#datasets <- "MLTbx"

for(f in datasets){
    e[[f]] <- as.data.frame(fread(file.path(dir, paste0(f, ".txt"))))
    if(f == "MLTbx") {
        rown <- e[[f]][,1]
        e[[f]] <- e[[f]][, -1]
        rownames(e[[f]]) <- rown
        colnames(e[[f]]) <- NULL
    }
    save(list = f, envir = e, file = paste0(f, ".rda"), compress = "bzip2")
}

