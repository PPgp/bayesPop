\name{adjust.to.aggregation}
\alias{adjust.to.aggregation}

\title{
Population Adjustments To Match Aggregated Targets
}
\description{
Creates a file with adjusted sex- and age-specific population counts for each location and year, where the sum of the adjusted counts over all locations matches the given target for each year, sex and age. 
}
\usage{
adjust.to.aggregation(pop.pred, target.file, output.file = "adjusted_population.txt", 
    target.code = NULL, variant.name = "median", target.id.column = "country_code", 
    output.id.column = "reg_code", stat = "mean", exclude.codes = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pop.pred}{
%%     ~~Describe \code{pop.pred} here~~
}
  \item{target.file}{
%%     ~~Describe \code{agfile} here~~
}
  \item{output.file}{
%%     ~~Describe \code{output.file} here~~
}
  \item{target.code}{
%%     ~~Describe \code{adjust.to.code} here~~
}
  \item{variant.name}{
%%     ~~Describe \code{variant.name} here~~
}
  \item{target.id.column}{
%%     ~~Describe \code{ag.id.column} here~~
}
  \item{output.id.column}{
%%     ~~Describe \code{output.id.column} here~~
}
  \item{stat}{
%%     ~~Describe \code{stat} here~~
}
  \item{exclude.codes}{
%%     ~~Describe \code{exclude.codes} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or standard data sets, see data().

## The function is currently defined as
function (pop.pred, agfile, output.file = "adjusted_population.txt", 
    adjust.to.code = NULL, variant.name = "median", ag.id.column = "country_code", 
    output.id.column = "reg_code", stat = "mean", exclude.codes = NULL) 
{
    if (!has.pop.aggregation(pop.pred = pop.pred)) 
        stop("The pop.pred object does not contain any aggregation. Consider running pop.aggregate().")
    if (!stat \%in\% c("median", "mean")) 
        stop("Argument 'stat' must be either 'mean' or 'median'.")
    agdat <- data.table::fread(agfile)
    if ("variant" \%in\% colnames(agdat)) {
        if (!tolower(variant.name) \%in\% tolower(agdat$variant)) 
            stop("Value '", variant.name, "' not found in the column 'variant' of the agfile. Consider changing the 'variant.name' argument.")
        agdat <- agdat[variant == variant.name][, `:=`(variant, 
            NULL)]
    }
    if (ag.id.column \%in\% colnames(agdat)) {
        if (!is.null(adjust.to.code)) {
            agdat <- agdat[get(ag.id.column) == adjust.to.code][, 
                `:=`((ag.id.column), NULL)]
            if (nrow(agdat) == 0) 
                stop("Value ", adjust.to.code, " not found in agfile.")
        }
        else {
            if (length((adj.to.code <- unique(agdat[[ag.id.column]]))) > 
                1) {
                stop("agfile contains more than one location. Set 'adjust.to.code' to specifify which location to use.")
            }
        }
    }
    else {
        if (!is.null(adjust.to.code) && any(duplicated(agdat[, 
            c("sex", "age"), with = FALSE]))) 
            stop("Column ", ag.id.column, " not found in agfile. Consider changing the 'ag.id.column' argument.")
    }
    if (any(duplicated(agdat[, c("sex", "age"), with = FALSE]))) 
        stop("The agfile contains duplicates of sex and age combinations.")
    if (pop.pred$annual) 
        agdat[, `:=`(age, as.integer(gsub("+", "", age, fixed = TRUE)))]
    year.cols <- grep("^[0-9]{4}", colnames(agdat), value = TRUE)
    agdat.long <- data.table::melt(agdat, id.vars = c("sex", 
        "age"), measure.vars = year.cols, variable.name = "year", 
        variable.factor = FALSE, value.name = "target")
    agdat.long[, `:=`(year, as.integer(year))]
    pop.aggr <- get.pop.aggregation(pop.pred = pop.pred)
    agsimlocs <- pop.aggr$countries
    if (nrow(agsimlocs) > 1) {
        if (is.null(adj.to.code)) 
            stop("The aggregation of pop.pred contain more than one aggregated locations. Consider setting 'adj.to.code'.")
        agsim.id <- adj.to.code
    }
    else agsim.id <- agsimlocs$code
    popaggr.traj <- rbind(get.pop.exba(paste0("P", agsim.id, 
        "_M{}"), pop.aggr, as.dt = TRUE)[, `:=`(sex, "male")], 
        get.pop.exba(paste0("P", agsim.id, "_F{}"), pop.aggr, 
            as.dt = TRUE)[, `:=`(sex, "female")])
    popaggr.stat <- popaggr.traj[, list(sim = do.call(stat, list(indicator))), 
        by = c("year", "sex", "age")]
    targets <- merge(agdat.long, popaggr.stat, by = c("year", 
        "sex", "age"))[, `:=`(shift, sim - target)]
    pop.traj <- rbind(get.pop.exba(paste0("PXXX_M{}"), pop.pred, 
        as.dt = TRUE)[, `:=`(sex, "male")], get.pop.exba(paste0("PXXX_F{}"), 
        pop.pred, as.dt = TRUE)[, `:=`(sex, "female")])
    pop.stat <- pop.traj[, list(sim = do.call(stat, list(indicator))), 
        by = c("country_code", "year", "sex", "age")]
    pop.stat[, `:=`(share, sim/sum(sim)), by = c("year", "sex", 
        "age")]
    totpop <- pop.stat[, list(totpop = sum(sim)), by = c("country_code", 
        "year")][, `:=`(totshare, totpop/sum(totpop)), by = c("year")]
    pop.stat[totpop, `:=`(totshare, i.totshare), on = c("country_code", 
        "year")]
    pop.stat[is.na(share), `:=`(share, totshare)]
    if (!is.null(exclude.codes)) {
        pop.stat[country_id \%in\% exclude.codes, `:=`(share, 0)]
    }
    pop.stat[, `:=`(share, share/sum(share)), by = c("year", 
        "sex", "age")]
    pop.stat[targets, `:=`(totshift, i.shift), on = c("year", 
        "sex", "age")]
    pop.stat[totshift > 0 & sim == 0, `:=`(share, 0)]
    pop.stat[, `:=`(share, share/sum(share)), by = c("year", 
        "sex", "age")][is.na(share), `:=`(share, 0)]
    pop.stat[, `:=`(shift, totshift * share)][, `:=`(simadj, 
        pmax(0, sim - shift))]
    respop <- data.table::dcast(pop.stat, country_code + sex + 
        age ~ year, value.var = "simadj")
    data.table::setnames(respop, "country_code", output.id.column)
    data.table::fwrite(respop, file = output.file, sep = "\t")
    cat("\nAdjustment file written into", output.file, "\n")
    return(output.file)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
