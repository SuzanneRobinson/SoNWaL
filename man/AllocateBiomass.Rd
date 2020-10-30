\name{AllocateBiomass}
\alias{AllocateBiomass}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  AllocateBiomass
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
AllocateBiomass(state, site, parms, weather)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{state}{
%%     ~~Describe \code{state} here~~
}
  \item{site}{
%%     ~~Describe \code{site} here~~
}
  \item{parms}{
%%     ~~Describe \code{parms} here~~
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
function (state, site, parms, weather) 
{
    current.month <- weather[["Month"]]
    leaf.grow <- parms[["leaf.grow"]]
    leaf.fall <- parms[["leaf.fall"]]
    m0 <- parms[["m0"]]
    pRx <- parms[["pRx"]]
    pRn <- parms[["pRn"]]
    fN <- state[["fN"]]
    PhysMod <- state[["PhysMod"]]
    m <- m0 + (1 - m0) * fN
    pR <- pRx * pRn/(pRn + (pRx - pRn) * m * PhysMod)
    pFS2 <- parms[["pFS2"]]
    pFS20 <- parms[["pFS20"]]
    pfsPower <- log(pFS20/pFS2)/log(10)
    pfsConst <- pFS2/(2^pfsPower)
    dg <- state[["dg"]]
    pFS <- pfsConst * dg^pfsPower
    pS <- (1 - pR)/(pFS + 1)
    pF <- 1 - pR - pS
    NPP <- state[["NPP"]]
    difWl <- NPP * pF
    difWr <- NPP * pR
    difWsbr <- NPP * pS
    gammaFx <- parms[["gammaFx"]]
    gammaF0 <- parms[["gammaF0"]]
    tgammaF <- parms[["tgammaF"]]
    if (isTRUE(isDormantSeason(current.month, leaf.grow, leaf.fall))) {
        if (current.month - 1 == 0) {
            if (isFALSE(isDormantSeason(12, leaf.grow, leaf.fall))) {
                Littfall <- state[["WlDormant"]]
            }
            else {
                Littfall <- 0
            }
        }
        else {
            if (isFALSE(isDormantSeason(current.month - 1, leaf.grow, 
                leaf.fall))) {
                Littfall <- state[["WlDormant"]]
            }
            else {
                Littfall <- 0
            }
        }
        state[["Wl"]] <- 0
        difRoots <- 0
        state[["Wr"]] <- state[["Wr"]] + difRoots
        state[["Wsbr"]] <- state[["Wsbr"]]
        difLitter <- ifelse(Littfall == 0, state[["Wlitt"]] * 
            Littfall, Littfall)
        TotalLitter <- state[["Wlitt"]] + difLitter
    }
    else if ((current.month - 1 == 0 && isFALSE(isDormantSeason(current.month, 
        leaf.grow, leaf.fall)) && isTRUE(isDormantSeason(12, 
        leaf.grow, leaf.fall))) || (isFALSE(isDormantSeason(current.month, 
        leaf.grow, leaf.fall)) && isTRUE(isDormantSeason(current.month - 
        1, leaf.grow, leaf.fall)))) {
        Wl <- state[["WlDormant"]]
        NPPDebt <- NPP - Wl
        t <- state[["t"]]
        Littfall <- gammaFx * gammaF0/(gammaF0 + (gammaFx - gammaF0) * 
            exp(-12 * log(1 + gammaFx/gammaF0) * t/tgammaF))
        difLitter <- Littfall * Wl
        Wr <- state[["Wr"]]
        Rttover <- parms[["Rttover"]]
        difRoots <- Rttover * Wr
        if (NPPDebt > 0) {
            difWl <- NPPDebt * pF
            difWr <- NPPDebt * pR
            difWsbr <- NPPDebt * pS
        }
        else {
            difWl <- 0
            difWr <- 0
            difWsbr <- 0
        }
        state[["Wl"]] <- Wl
        state[["Wr"]] <- Wr + difWr - difRoots
        state[["Wsbr"]] <- state[["Wsbr"]] + difWsbr
        TotalLitter <- state[["Wlitt"]] + difLitter
    }
    else if (leaf.grow != 0) {
        Wl <- state[["Wl"]]
        t <- state[["t"]]
        Littfall <- gammaFx * gammaF0/(gammaF0 + (gammaFx - gammaF0) * 
            exp(-12 * log(1 + gammaFx/gammaF0) * t/tgammaF))
        difLitter <- Littfall * Wl
        Wr <- state[["Wr"]]
        Rttover <- parms[["Rttover"]]
        difRoots <- Rttover * Wr
        NPPDebt <- NPP - Wl
        if (NPPDebt < 0) {
            NPPDebt <- NPPDebt + NPP
            if (NPPDebt > 0) {
                difWl <- NPPDebt * pF
                difWr <- NPPDebt * pR
                difWsbr <- NPPDebt * pS
            }
            else {
                difWl <- 0
                difWr <- 0
                difWsbr <- 0
            }
        }
        else if (NPPDebt >= 0) {
            difWl <- NPP * pF
            difWr <- NPP * pR
            difWsbr <- NPP * pS
        }
        state[["Wl"]] <- Wl + difWl - difLitter
        state[["Wr"]] <- Wr + difWr - difRoots
        state[["Wsbr"]] <- state[["Wsbr"]] + difWsbr
        TotalLitter <- state[["Wlitt"]] + difLitter
    }
    else {
        Wl <- state[["Wl"]]
        t <- state[["t"]]
        Littfall <- gammaFx * gammaF0/(gammaF0 + (gammaFx - gammaF0) * 
            exp(-12 * log(1 + gammaFx/gammaF0) * t/tgammaF))
        difLitter <- Littfall * Wl
        state[["Wlitt"]] <- state[["Wlitt"]] + difLitter
        Wr <- state[["Wr"]]
        Rttover <- parms[["Rttover"]]
        difRoots <- Rttover * Wr
        Wsbr <- state[["Wsbr"]]
        TotalLitter <- state[["Wlitt"]]
        state[["Wl"]] <- Wl + difWl - difLitter
        state[["Wr"]] <- Wr + difWr - difRoots
        state[["Wsbr"]] <- state[["Wsbr"]] + difWsbr
        TotalLitter <- state[["Wlitt"]] + difLitter
    }
    if (current.month + 1 == 13) {
        if (isFALSE(isDormantSeason(current.month, leaf.grow, 
            leaf.fall)) && isTRUE(isDormantSeason(1, leaf.grow, 
            leaf.fall))) {
            state[["WlDormant"]] <- state[["Wl"]]
        }
    }
    else {
        if (isFALSE(isDormantSeason(current.month, leaf.grow, 
            leaf.fall)) && isTRUE(isDormantSeason(current.month + 
            1, leaf.grow, leaf.fall))) {
            state[["WlDormant"]] <- state[["Wl"]]
        }
    }
    state[c("pR", "pFS", "pS", "pF", "difWl", "difWr", "difWsbr", 
        "Littfall", "difLitter", "difRoots", "TotalLitter")] <- c(pR, 
        pFS, pS, pF, difWl, difWr, difWsbr, Littfall, difLitter, 
        difRoots, TotalLitter)
    return(state)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
