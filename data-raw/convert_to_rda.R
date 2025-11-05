# convert ASCII into rda to be used in the bayesPop package

library(data.table)
e <- new.env()

dir <- "."
datasets <- c("MLTbx", paste0("vwBaseYear", c(2010, 2012, 2015, 2017, 2019, 2022, 2024)))
datasets <- "vwBaseYear2024"
#datasets <- c("rc1FDM", "rc5FDM")
#datasets <- "MLTbx"

add5year.settings <- function(df){
    # convert some settings into 5-year equivalents
    vw <- data.table(df)
    cols <- colnames(vw)
    LAMP <- which(cols == "LatestAgeMortalityPattern")
    DF <- which(cols == "SmoothDFLatestAgeMortalityPattern")
    
    vw[, LatestAgeMortalityPattern1 := LatestAgeMortalityPattern]
    age.pat.val <- sapply(vw[,LatestAgeMortalityPattern], function(x) eval(parse(text = x)))
    
    age.pat.val5 <- lapply(age.pat.val, 
                           function(x) if(length(x) == 1) {if(x < 0) floor(x) else ceiling(x)} 
                           else c(floor(x[1]/5), ceiling(x[2]/5)))
    vw[, LatestAgeMortalityPattern := sapply(age.pat.val5, 
                                             function(x) if(length(x) == 1) paste(x) else paste0("c(", x[1], ",", x[2], ")"))]
    vw[, SmoothDFLatestAgeMortalityPattern1 := SmoothDFLatestAgeMortalityPattern]
    vw[, SmoothDFLatestAgeMortalityPattern := floor(SmoothDFLatestAgeMortalityPattern1/5)]
    
    # reorder columns
    setcolorder(vw, c(1:LAMP, ncol(vw)-1, (LAMP+1):DF, ncol(vw)))
    
    # convert non-ascii
    replace.names <- list("384"="Cote d'Ivoire",
                          "638"="Reunion",
                          "531"="Curacao",
                          "652"="Saint Bartelemy",
                          "792"="Turkiye"
    )
    for(code in names(replace.names)) {
        vw[country_code == code, country := replace.names[[code]]]
    }
    return(data.frame(vw))
}

for(f in datasets){
    e[[f]] <- as.data.frame(fread(file.path(dir, paste0(f, ".txt"))))
    if(f == "MLTbx") {
        rown <- e[[f]][,1]
        e[[f]] <- e[[f]][, -1]
        rownames(e[[f]]) <- rown
        colnames(e[[f]]) <- NULL
    }
    if( f %in% c("vwBaseYear2024", "vwBaseYear2022"))
        e[[f]] <- add5year.settings(e[[f]])
    
    if(f == "vwBaseYear2024"){
        fdmpat <- as.data.frame(fread(file.path(dir, "mig_fdm_patterns.csv")))
        e[[f]] <- merge(e[[f]][, !colnames(e[[f]]) %in% setdiff(colnames(fdmpat), "country_code")], fdmpat, 
                        by = "country_code", all.x = TRUE, sort = FALSE)
    }
    save(list = f, envir = e, file = paste0(f, ".rda"), compress = "bzip2")
}

