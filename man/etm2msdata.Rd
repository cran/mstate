% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/etm2msdata.R
\name{etm2msdata}
\alias{etm2msdata}
\alias{tra2trans}
\title{Converts between etm and msdata format}
\usage{
etm2msdata(etmdata, id, tra, covs)
}
\arguments{
\item{etmdata}{Multi-state data in \code{etm} format}

\item{id}{Column name identifying the subject id}

\item{tra}{Transition matrix in \code{etm} format}

\item{covs}{Vector of column names containing covariates to be included}
}
\description{
Converts multi-state data back and forth between etm and msdata formats.
Covariates have to be dealt with separately.
}
\details{
\code{msdata2etm} will convert from \code{msdata} format to \code{etm}
format; \code{etm2msdata} will convert from \code{etm} format to
\code{msdata} format. Both \code{msdata2etm} and \code{etm2msdata} work with
basic time-fixed covariates. Time-dependent covariates are not supported.
The function \code{msdata2etm} will work for transition-specific covariates,
but the result does not really make much sense when used in etm.
}
\examples{

# Transition matrix for illness-death model
tmat <- trans.illdeath()
# Data in wide format, for transition 1 this is dataset E1 of
# Therneau & Grambsch (T&G)
tg <- data.frame(id=1:6,illt=c(1,1,6,6,8,9),ills=c(1,0,1,1,0,1),
                 dt=c(5,1,9,7,8,12),ds=c(1,1,1,1,1,1),
                 x1=c(1,1,1,0,0,0),x2=c(6:1))
# Data in long format using msprep
tglong <- msprep(time=c(NA,"illt","dt"),status=c(NA,"ills","ds"),
                 data=tg,keep=c("x1","x2"),trans=tmat, id="id")
# Same thing in etm format
tra <- trans2tra(tmat)
tgetm <- msdata2etm(tglong, id="id")
tgetm <- msdata2etm(tglong, id="id", covs=c("x1", "x2")) # with covariates
# And back
etm2msdata(tgetm, id="id", tra=tra)
etm2msdata(tgetm, id="id", tra=tra, covs=c("x1", "x2")) # with covariates
}
\author{
Hein Putter \email{H.Putter@lumc.nl}
}
\keyword{datagen}
