% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vis.multiple.pt.R
\name{vis.multiple.pt}
\alias{vis.multiple.pt}
\title{Visualise multiple probtrans objects}
\usage{
vis.multiple.pt(
  x,
  from = 1,
  to,
  xlab = "Time",
  ylab = "Probability",
  xlim = NULL,
  ylim = NULL,
  cols,
  lwd,
  labels,
  conf.int = 0.95,
  conf.type = c("log", "plain", "none"),
  legend.title
)
}
\arguments{
\item{x}{A list of \code{"probtrans"} objects}

\item{from}{The starting state from which the probabilities are used to plot
Numeric, as in \code{plot.probtrans}}

\item{to}{(Numeric) destination state}

\item{xlab}{A title for the x-axis; default is \code{"Time"}}

\item{ylab}{A title for the y-axis; default is \code{"Probability"}}

\item{xlim}{The x limits of the plot(s), default is range of time}

\item{ylim}{The y limits of the plot(s); if ylim is specified for
type="separate", then all plots use the same ylim for y limits}

\item{cols}{A vector specifying colors for the different transitions;
default is a palette from green to red, when type=\code{"filled"} (reordered
according to \code{ord}, and 1 (black), otherwise}

\item{lwd}{The line width, see \code{\link{par}}; default is 1}

\item{labels}{Character vector labelling each element of x (e.g. label
for a reference patient) - so labels = c("Patient 1", "Patient 2")}

\item{conf.int}{Confidence level (\%) from 0-1 for probabilities, 
default is 0.95 (95\% CI). Setting to 0 removes the CIs.}

\item{conf.type}{Type of confidence interval - either "log" or "plain" . See
function details for details.}

\item{legend.title}{Character - title of legend}
}
\value{
A ggplot object.
}
\description{
Helper function allowing to visualise state probabilities for 
different reference patients/covariates. Multiple \code{"probtrans"} objects
are thus needed.
}
\examples{

library(ggplot2)

data("aidssi")
head(aidssi)
si <- aidssi

# Prepare transition matrix
tmat <- trans.comprisk(2, names = c("event-free", "AIDS", "SI"))

# Run msprep
si$stat1 <- as.numeric(si$status == 1)
si$stat2 <- as.numeric(si$status == 2)

silong <- msprep(
time = c(NA, "time", "time"), 
status = c(NA, "stat1", "stat2"), 
data = si, keep = "ccr5", trans = tmat
)

# Run cox model
silong <- expand.covs(silong, "ccr5")
c1 <- coxph(Surv(time, status) ~ ccr5WM.1 + ccr5WM.2 + strata(trans),
            data = silong)
            
# 1. Prepare patient data - both CCR5 genotypes
WW <- data.frame(
ccr5WM.1 = c(0, 0), 
ccr5WM.2 = c(0, 0), 
trans = c(1, 2), 
strata = c(1, 2)
)

WM <- data.frame(
ccr5WM.1 = c(1, 0), 
ccr5WM.2 = c(0, 1),
trans = c(1, 2), 
strata = c(1, 2)
)

# 2. Make msfit objects
msf.WW <- msfit(c1, WW, trans = tmat)
msf.WM <- msfit(c1, WM, trans = tmat)

# 3. Make probtrans objects
pt.WW <- probtrans(msf.WW, predt = 0)
pt.WM <- probtrans(msf.WM, predt = 0)           

# Plot - see vignette for more details
vis.multiple.pt(
x = list(pt.WW, pt.WM), 
from = 1,
to = 2, 
conf.type = "log",
cols = c(1, 2),
labels = c("Pat WW", "Pat WM"),
legend.title = "Ref patients"
)

}
\author{
Edouard F. Bonneville \email{e.f.bonneville@lumc.nl}
}
