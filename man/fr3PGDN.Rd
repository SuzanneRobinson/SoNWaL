\name{fr3PGDN}
\alias{fr3PGDN}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  fr3PGDN
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
fr3PGDN(weather, presc, t = 0, N = 2500, Wl = 0.01, WlDormant = 0, Wr = 0.01, Wsbr = 0.1, Wlitt = 0, YrC = 18.2, YlC = 6.5, OC = 122.77, YrN = 0, YlN = 1.52, ON = 6.9, Nav = 3, rotation = 1, cycle = 1, rm.sprouts = F, nyears = 35, initial.month = 1, latitude = 57.06, soilclass = -1, ASW = 200, MaxASW = 300, MinASW = 0, CO2 = 400, pFS2 = 0.76, pFS20 = 0.46, aS = 0.0572, nS = 2.4568, pRx = 0.48, pRn = 0.21, Tmin = -5, Topt = 15, Tmax = 35, kF = 1, SWconst0 = 0.7, SWpower0 = 9, m0 = 0.018, fN0 = 0.296, fNn = 1.27, MaxAge = 265.6, nAge = 3.545, rAge = 0.796, gammaFx = 0.025, gammaF0 = 0.0022, tgammaF = 60, Rttover = 0.07, MaxCond = 0.02, LAIgcx = 3.33, BLcond = 0.2, wSx1000 = 255, thinPower = 1.5, mF = 0.2, mR = 0.2, mS = 0.2, SLA0 = 6, SLA1 = 4, tSLA = 6, k = 0.52, fullCanAge = 18, MaxIntcptn = 0.15, LAImaxIntcptn = 5, alpha = 0.05, Y = 0.49, poolFractn = 0, e20 = 2.2, rhoAir = 1.2, lambda = 2460000, VPDconv = 0.000622, fracBB0 = 0.3, fracBB1 = 0.1, tBB = 10, rhoMin = 0.39, rhoMax = 0.39, tRho = 5, Qa = -90, Qb = 0.8, gDM_mol = 24, molPAR_MJ = 2.3, CoeffCond = 0.05, fCalpha700 = 1.433, fCg700 = 0.451, fCalphax = 2.33333333333333, fCg0 = 1.75, MinCond = 0.015, klmax = 0.031532358, krmax = 0.0042394358, komax = 0.000458861, hc = 0.23, qir = 334.29, qil = 49.08, qh = 23.63, qbc = 2.21, el = 0.25, er = 0.56, Nf = 0.0684, Navm = 0.01, Navx = 20, leaf.grow = 4, leaf.fall = 11, Wl.s = 0.01, Wsbr.s = 0.1, Wr.s = 0.01, pWl.sprouts = 0.5, pWsbr.sprouts = 0.9, cod.pred = "3PG", cod.clim = "Average")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{weather}{
%%     ~~Describe \code{weather} here~~
}
  \item{presc}{
%%     ~~Describe \code{presc} here~~
}
  \item{t}{
%%     ~~Describe \code{t} here~~
}
  \item{N}{
%%     ~~Describe \code{N} here~~
}
  \item{Wl}{
%%     ~~Describe \code{Wl} here~~
}
  \item{WlDormant}{
%%     ~~Describe \code{WlDormant} here~~
}
  \item{Wr}{
%%     ~~Describe \code{Wr} here~~
}
  \item{Wsbr}{
%%     ~~Describe \code{Wsbr} here~~
}
  \item{Wlitt}{
%%     ~~Describe \code{Wlitt} here~~
}
  \item{YrC}{
%%     ~~Describe \code{YrC} here~~
}
  \item{YlC}{
%%     ~~Describe \code{YlC} here~~
}
  \item{OC}{
%%     ~~Describe \code{OC} here~~
}
  \item{YrN}{
%%     ~~Describe \code{YrN} here~~
}
  \item{YlN}{
%%     ~~Describe \code{YlN} here~~
}
  \item{ON}{
%%     ~~Describe \code{ON} here~~
}
  \item{Nav}{
%%     ~~Describe \code{Nav} here~~
}
  \item{rotation}{
%%     ~~Describe \code{rotation} here~~
}
  \item{cycle}{
%%     ~~Describe \code{cycle} here~~
}
  \item{rm.sprouts}{
%%     ~~Describe \code{rm.sprouts} here~~
}
  \item{nyears}{
%%     ~~Describe \code{nyears} here~~
}
  \item{initial.month}{
%%     ~~Describe \code{initial.month} here~~
}
  \item{latitude}{
%%     ~~Describe \code{latitude} here~~
}
  \item{soilclass}{
%%     ~~Describe \code{soilclass} here~~
}
  \item{ASW}{
%%     ~~Describe \code{ASW} here~~
}
  \item{MaxASW}{
%%     ~~Describe \code{MaxASW} here~~
}
  \item{MinASW}{
%%     ~~Describe \code{MinASW} here~~
}
  \item{CO2}{
%%     ~~Describe \code{CO2} here~~
}
  \item{pFS2}{
%%     ~~Describe \code{pFS2} here~~
}
  \item{pFS20}{
%%     ~~Describe \code{pFS20} here~~
}
  \item{aS}{
%%     ~~Describe \code{aS} here~~
}
  \item{nS}{
%%     ~~Describe \code{nS} here~~
}
  \item{pRx}{
%%     ~~Describe \code{pRx} here~~
}
  \item{pRn}{
%%     ~~Describe \code{pRn} here~~
}
  \item{Tmin}{
%%     ~~Describe \code{Tmin} here~~
}
  \item{Topt}{
%%     ~~Describe \code{Topt} here~~
}
  \item{Tmax}{
%%     ~~Describe \code{Tmax} here~~
}
  \item{kF}{
%%     ~~Describe \code{kF} here~~
}
  \item{SWconst0}{
%%     ~~Describe \code{SWconst0} here~~
}
  \item{SWpower0}{
%%     ~~Describe \code{SWpower0} here~~
}
  \item{m0}{
%%     ~~Describe \code{m0} here~~
}
  \item{fN0}{
%%     ~~Describe \code{fN0} here~~
}
  \item{fNn}{
%%     ~~Describe \code{fNn} here~~
}
  \item{MaxAge}{
%%     ~~Describe \code{MaxAge} here~~
}
  \item{nAge}{
%%     ~~Describe \code{nAge} here~~
}
  \item{rAge}{
%%     ~~Describe \code{rAge} here~~
}
  \item{gammaFx}{
%%     ~~Describe \code{gammaFx} here~~
}
  \item{gammaF0}{
%%     ~~Describe \code{gammaF0} here~~
}
  \item{tgammaF}{
%%     ~~Describe \code{tgammaF} here~~
}
  \item{Rttover}{
%%     ~~Describe \code{Rttover} here~~
}
  \item{MaxCond}{
%%     ~~Describe \code{MaxCond} here~~
}
  \item{LAIgcx}{
%%     ~~Describe \code{LAIgcx} here~~
}
  \item{BLcond}{
%%     ~~Describe \code{BLcond} here~~
}
  \item{wSx1000}{
%%     ~~Describe \code{wSx1000} here~~
}
  \item{thinPower}{
%%     ~~Describe \code{thinPower} here~~
}
  \item{mF}{
%%     ~~Describe \code{mF} here~~
}
  \item{mR}{
%%     ~~Describe \code{mR} here~~
}
  \item{mS}{
%%     ~~Describe \code{mS} here~~
}
  \item{SLA0}{
%%     ~~Describe \code{SLA0} here~~
}
  \item{SLA1}{
%%     ~~Describe \code{SLA1} here~~
}
  \item{tSLA}{
%%     ~~Describe \code{tSLA} here~~
}
  \item{k}{
%%     ~~Describe \code{k} here~~
}
  \item{fullCanAge}{
%%     ~~Describe \code{fullCanAge} here~~
}
  \item{MaxIntcptn}{
%%     ~~Describe \code{MaxIntcptn} here~~
}
  \item{LAImaxIntcptn}{
%%     ~~Describe \code{LAImaxIntcptn} here~~
}
  \item{alpha}{
%%     ~~Describe \code{alpha} here~~
}
  \item{Y}{
%%     ~~Describe \code{Y} here~~
}
  \item{poolFractn}{
%%     ~~Describe \code{poolFractn} here~~
}
  \item{e20}{
%%     ~~Describe \code{e20} here~~
}
  \item{rhoAir}{
%%     ~~Describe \code{rhoAir} here~~
}
  \item{lambda}{
%%     ~~Describe \code{lambda} here~~
}
  \item{VPDconv}{
%%     ~~Describe \code{VPDconv} here~~
}
  \item{fracBB0}{
%%     ~~Describe \code{fracBB0} here~~
}
  \item{fracBB1}{
%%     ~~Describe \code{fracBB1} here~~
}
  \item{tBB}{
%%     ~~Describe \code{tBB} here~~
}
  \item{rhoMin}{
%%     ~~Describe \code{rhoMin} here~~
}
  \item{rhoMax}{
%%     ~~Describe \code{rhoMax} here~~
}
  \item{tRho}{
%%     ~~Describe \code{tRho} here~~
}
  \item{Qa}{
%%     ~~Describe \code{Qa} here~~
}
  \item{Qb}{
%%     ~~Describe \code{Qb} here~~
}
  \item{gDM_mol}{
%%     ~~Describe \code{gDM_mol} here~~
}
  \item{molPAR_MJ}{
%%     ~~Describe \code{molPAR_MJ} here~~
}
  \item{CoeffCond}{
%%     ~~Describe \code{CoeffCond} here~~
}
  \item{fCalpha700}{
%%     ~~Describe \code{fCalpha700} here~~
}
  \item{fCg700}{
%%     ~~Describe \code{fCg700} here~~
}
  \item{fCalphax}{
%%     ~~Describe \code{fCalphax} here~~
}
  \item{fCg0}{
%%     ~~Describe \code{fCg0} here~~
}
  \item{MinCond}{
%%     ~~Describe \code{MinCond} here~~
}
  \item{klmax}{
%%     ~~Describe \code{klmax} here~~
}
  \item{krmax}{
%%     ~~Describe \code{krmax} here~~
}
  \item{komax}{
%%     ~~Describe \code{komax} here~~
}
  \item{hc}{
%%     ~~Describe \code{hc} here~~
}
  \item{qir}{
%%     ~~Describe \code{qir} here~~
}
  \item{qil}{
%%     ~~Describe \code{qil} here~~
}
  \item{qh}{
%%     ~~Describe \code{qh} here~~
}
  \item{qbc}{
%%     ~~Describe \code{qbc} here~~
}
  \item{el}{
%%     ~~Describe \code{el} here~~
}
  \item{er}{
%%     ~~Describe \code{er} here~~
}
  \item{Nf}{
%%     ~~Describe \code{Nf} here~~
}
  \item{Navm}{
%%     ~~Describe \code{Navm} here~~
}
  \item{Navx}{
%%     ~~Describe \code{Navx} here~~
}
  \item{leaf.grow}{
%%     ~~Describe \code{leaf.grow} here~~
}
  \item{leaf.fall}{
%%     ~~Describe \code{leaf.fall} here~~
}
  \item{Wl.s}{
%%     ~~Describe \code{Wl.s} here~~
}
  \item{Wsbr.s}{
%%     ~~Describe \code{Wsbr.s} here~~
}
  \item{Wr.s}{
%%     ~~Describe \code{Wr.s} here~~
}
  \item{pWl.sprouts}{
%%     ~~Describe \code{pWl.sprouts} here~~
}
  \item{pWsbr.sprouts}{
%%     ~~Describe \code{pWsbr.sprouts} here~~
}
  \item{cod.pred}{
%%     ~~Describe \code{cod.pred} here~~
}
  \item{cod.clim}{
%%     ~~Describe \code{cod.clim} here~~
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
function (weather, presc, t = 0, N = 2500, Wl = 0.01, WlDormant = 0, 
    Wr = 0.01, Wsbr = 0.1, Wlitt = 0, YrC = 18.2, YlC = 6.5, 
    OC = 122.77, YrN = 0, YlN = 1.52, ON = 6.9, Nav = 3, rotation = 1, 
    cycle = 1, rm.sprouts = F, nyears = 35, initial.month = 1, 
    latitude = 57.06, soilclass = -1, ASW = 200, MaxASW = 300, 
    MinASW = 0, CO2 = 400, pFS2 = 0.76, pFS20 = 0.46, aS = 0.0572, 
    nS = 2.4568, pRx = 0.48, pRn = 0.21, Tmin = -5, Topt = 15, 
    Tmax = 35, kF = 1, SWconst0 = 0.7, SWpower0 = 9, m0 = 0.018, 
    fN0 = 0.296, fNn = 1.27, MaxAge = 265.6, nAge = 3.545, rAge = 0.796, 
    gammaFx = 0.025, gammaF0 = 0.0022, tgammaF = 60, Rttover = 0.07, 
    MaxCond = 0.02, LAIgcx = 3.33, BLcond = 0.2, wSx1000 = 255, 
    thinPower = 1.5, mF = 0.2, mR = 0.2, mS = 0.2, SLA0 = 6, 
    SLA1 = 4, tSLA = 6, k = 0.52, fullCanAge = 18, MaxIntcptn = 0.15, 
    LAImaxIntcptn = 5, alpha = 0.05, Y = 0.49, poolFractn = 0, 
    e20 = 2.2, rhoAir = 1.2, lambda = 2460000, VPDconv = 0.000622, 
    fracBB0 = 0.3, fracBB1 = 0.1, tBB = 10, rhoMin = 0.39, rhoMax = 0.39, 
    tRho = 5, Qa = -90, Qb = 0.8, gDM_mol = 24, molPAR_MJ = 2.3, 
    CoeffCond = 0.05, fCalpha700 = 1.433, fCg700 = 0.451, fCalphax = 2.33333333333333, 
    fCg0 = 1.75, MinCond = 0.015, klmax = 0.031532358, krmax = 0.0042394358, 
    komax = 0.000458861, hc = 0.23, qir = 334.29, qil = 49.08, 
    qh = 23.63, qbc = 2.21, el = 0.25, er = 0.56, Nf = 0.0684, 
    Navm = 0.01, Navx = 20, leaf.grow = 4, leaf.fall = 11, Wl.s = 0.01, 
    Wsbr.s = 0.1, Wr.s = 0.01, pWl.sprouts = 0.5, pWsbr.sprouts = 0.9, 
    cod.pred = "3PG", cod.clim = "Average") 
{
    parms <- c(pFS2, pFS20, aS, nS, pRx, pRn, Tmin, Topt, Tmax, 
        kF, SWconst0, SWpower0, m0, fN0, fNn, MaxAge, nAge, rAge, 
        gammaFx, gammaF0, tgammaF, Rttover, MaxCond, LAIgcx, 
        BLcond, wSx1000, thinPower, mF, mR, mS, SLA0, SLA1, tSLA, 
        k, fullCanAge, MaxIntcptn, LAImaxIntcptn, alpha, Y, poolFractn, 
        e20, rhoAir, lambda, VPDconv, fracBB0, fracBB1, tBB, 
        rhoMin, rhoMax, tRho, Qa, Qb, gDM_mol, molPAR_MJ, CoeffCond, 
        fCalpha700, fCg700, fCalphax, fCg0, klmax, krmax, komax, 
        hc, qir, qil, qh, qbc, el, er, Nf, Navm, Navx, MinCond, 
        Wl.s, Wsbr.s, Wr.s, pWl.sprouts, pWsbr.sprouts, leaf.grow, 
        leaf.fall)
    names(parms) <- c("pFS2", "pFS20", "aS", "nS", "pRx", "pRn", 
        "Tmin", "Topt", "Tmax", "kF", "SWconst0", "SWpower0", 
        "m0", "fN0", "fNn", "MaxAge", "nAge", "rAge", "gammaFx", 
        "gammaF0", "tgammaF", "Rttover", "MaxCond", "LAIgcx", 
        "BLcond", "wSx1000", "thinPower", "mF", "mR", "mS", "SLA0", 
        "SLA1", "tSLA", "k", "fullCanAge", "MaxIntcptn", "LAImaxIntcptn", 
        "alpha", "Y", "poolFractn", "e20", "rhoAir", "lambda", 
        "VPDconv", "fracBB0", "fracBB1", "tBB", "rhoMin", "rhoMax", 
        "tRho", "Qa", "Qb", "gDM_mol", "molPAR_MJ", "CoeffCond", 
        "fCalpha700", "fCg700", "fCalphax", "fCg0", "klmax", 
        "krmax", "komax", "hc", "qir", "qil", "qh", "qbc", "el", 
        "er", "Nf", "Navm", "Navx", "MinCond", "Wl.s", "Wsbr.s", 
        "Wr.s", "pWl.sprouts", "pWsbr.sprouts", "leaf.grow", 
        "leaf.fall")
    vars.ini <- c(t, N, Wl, WlDormant, Wr, Wsbr, Wlitt, YrC, 
        YlC, OC, YrN, YlN, ON, Nf, Nav, rotation, cycle, rm.sprouts, 
        nyears, initial.month)
    names(vars.ini) <- c("t", "N", "Wl", "WlDormant", "Wr", "Wsbr", 
        "Wlitt", "YrC", "YlC", "OC", "YrN", "YlN", "ON", "Nf", 
        "Nav", "rotation", "cycle", "rm.sprouts", "nyears", "initial.month")
    site.info <- c(latitude, soilclass, ASW, MaxASW, MinASW, 
        CO2)
    names(site.info) <- c("latitude", "soilclass", "ASW", 
        "MaxASW", "MinASW", "CO2")
    parms.general <- list(daysinmonth = c(Jan = 31, Feb = 28, 
        Mar = 31, Apr = 30, May = 31, Jun = 30, Jul = 31, Aug = 31, 
        Sep = 30, Oct = 31, Nov = 30, Dec = 31), parms.soil = data.frame(soilclass.name = c("Sandy", 
        "Sandy loam", "Clay loam", "Clay", "Non standard", "No effect of ASW"), 
        soilclass = c(1, 2, 3, 4, -1, 0), SWconst = c(0.7, 0.6, 
            0.5, 0.4, parms[["SWconst0"]], 1), SWpower = c(9, 
            7, 5, 3, parms[["SWpower0"]], 1)))
    proj <- RunModel(stand.init = vars.ini, weather = weather, 
        site = site.info, general.info = parms.general, presc = presc, 
        parms = parms, cod.pred = cod.pred, cod.clim = cod.clim)
    return(as.data.frame(proj))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
