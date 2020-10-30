\name{UpdateSoil}
\alias{UpdateSoil}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  UpdateSoil
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
UpdateSoil(state, parms, site, general.info, weather)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{state}{
%%     ~~Describe \code{state} here~~
}
  \item{parms}{
%%     ~~Describe \code{parms} here~~
}
  \item{site}{
%%     ~~Describe \code{site} here~~
}
  \item{general.info}{
%%     ~~Describe \code{general.info} here~~
}
  \item{weather}{
%%     ~~Describe \code{weather} here~~
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
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (state, parms, site, general.info, weather) 
{
    GPP <- state[["GPP"]]
    NPP <- state[["NPP"]]
    Wdl <- state[["Wdl"]]
    Wds <- state[["Wds"]]
    Wdr <- state[["Wdr"]]
    YlC <- state[["YlC"]]
    YrC <- state[["YrC"]]
    OC <- state[["OC"]]
    YlN <- state[["YlN"]]
    YrN <- state[["YrN"]]
    ON <- state[["ON"]]
    Nav <- state[["Nav"]]
    Rs <- state[["Rs"]]
    difWl <- state[["difWl"]]
    difLitter <- state[["difLitter"]]
    difRoots <- state[["difRoots"]]
    dg <- state[["dg"]]
    fT <- state[["fT"]]
    fSW <- state[["fSW"]]
    klmax <- parms[["klmax"]]
    krmax <- parms[["krmax"]]
    komax <- parms[["komax"]]
    hc <- parms[["hc"]]
    qir <- parms[["qir"]]
    qil <- parms[["qil"]]
    qbc <- parms[["qbc"]]
    qh <- parms[["qh"]]
    el <- parms[["el"]]
    er <- parms[["er"]]
    Nf <- parms[["Nf"]]
    Navm <- parms[["Navm"]]
    Navx <- parms[["Navx"]]
    if (state[["t"]] == 0) {
        parms.soil <- general.info$parms.soil
        parms.sw.site <- parms.soil[which(parms.soil$soilclass == 
            site[["soilclass"]]), ]
        SWconst <- parms.sw.site[["SWconst"]]
        SWpower <- parms.sw.site[["SWpower"]]
        ASW <- state[["ASW"]]
        MaxASW <- site[["MaxASW"]]
        MoistRatio <- ASW/MaxASW
        Tmin <- parms[["Tmin"]]
        Tmax <- parms[["Tmax"]]
        Topt <- parms[["Topt"]]
        Tav <- weather[1, "Tmean"]
        dg <- ((state[["Wsbr"]] * 1000/state[["N"]])/parms[["aS"]])^(1/parms[["nS"]])
        fSW <- 1/(1 + ((1 - MoistRatio)/SWconst)^SWpower)
        if (Tav < Tmin | Tav > Tmax) {
            fT <- 0
        }
        else {
            fT <- ((Tav - Tmin)/(Topt - Tmin)) * ((Tmax - Tav)/(Tmax - 
                Topt))^((Tmax - Topt)/(Topt - Tmin))
        }
        Wds <- 0
        Wdl <- 0
        Wdr <- 0
        difLitter <- 0
        difRoots <- 0
        coarseDifRoots <- 0
        fineDifRoots <- 0
        difWl <- 0
    }
    fineCoarseRatio <- 0.1276 + 1462.2671 * exp(-1.7958 * dg)
    coarseDifRoots <- difRoots/(1 + fineCoarseRatio)
    fineDifRoots <- difRoots - coarseDifRoots
    kr <- krmax * fSW * fT
    kl <- klmax * fSW * fT
    ko <- komax * fSW * fT
    YlCflx <- kl * (1 - hc) * YlC
    YrCflx <- kr * (1 - hc) * YrC
    OCflx <- ko * OC
    hl <- kl * hc * YlC
    hr <- kr * hc * YrC
    hNl <- kl * hc * (YlN/qh)
    hNr <- kr * hc * (YrN/qh)
    YlNflx <- kl * ((1 - hc)/(1 - el)) * (YlN - el * (YlC/qbc))
    YlNflx <- max(0, YlNflx)
    YrNflx <- kr * ((1 - hc)/(1 - er)) * (YrN - er * (YrC/qbc))
    YrNflx <- max(0, YrNflx)
    ONflx <- ko * ON
    ONflx <- max(0, ONflx)
    Rs <- YlCflx + YrCflx + OCflx
    YrC <- YrC + ((Wds + coarseDifRoots)/2) - YrCflx - hr
    YlC <- YlC + ((difLitter + fineDifRoots + Wdl + Wdr)/2) - 
        YlCflx - hl
    OC <- OC + hl + hr - OCflx
    YrN <- YrN + ((Wds + coarseDifRoots)/(2 * qir)) - YrNflx - 
        hNr
    YlN <- YlN + ((difLitter + fineDifRoots + Wdl + Wdr)/(2 * 
        qil)) - YlNflx - hNl
    ON <- ON + hNr + hNl - ONflx
    totC <- YrC + YlC + OC
    totN <- YrN + YlN + ON
    Navflx <- YrNflx + YlNflx + ONflx
    Un <- difWl * Nf
    Nav <- Nav + Navflx - Un
    Nav <- max(Navm, Nav)
    if (Nav > Navx) {
        Nleach <- Nav - Navx
        Nav <- Navx
    }
    else {
        Nleach <- 0
    }
    fN <- (Nav - Navm)/(Navx - Navm)
    if (state[["t"]] != 0) {
        NEE <- Rs - NPP
        Ra <- GPP - NPP
        Reco <- Ra + Rs
    }
    else {
        NEE <- 0
        Ra <- 0
        Reco <- 0
    }
    state[c("YrC", "YlC", "OC", "YrN", "YlN", "ON", "kl", "kr", 
        "ko", "hl", "hr", "hNl", "hNr", "YrCflx", "YlCflx", "OCflx", 
        "YrNflx", "YlNflx", "ONflx", "Navflx", "totC", "totN", 
        "Un", "Nav", "Nleach", "fN", "NEE", "Reco", "Ra", "Rs")] <- c(YrC, 
        YlC, OC, YrN, YlN, ON, kl, kr, ko, hl, hr, hNl, hNr, 
        YrCflx, YlCflx, OCflx, YrNflx, YlNflx, ONflx, Navflx, 
        totC, totN, Un, Nav, Nleach, fN, NEE, Reco, Ra, Rs)
    return(state)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
