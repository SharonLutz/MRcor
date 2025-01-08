\name{MRcorDiff}
\alias{MRcorDiff}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A Capitalized Title (ideally limited to 65 characters)
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
MRcorDiff(n = 1e+05, nSNPx = 10, MAFx = c(rep(0.5, 10)), nSNPy = 10, MAFy = c(rep(0.5, 10)), nCov = 3, meanC = rep(0, 3), sdC = c(rep(1, 3)), deltaC = c(rep(0.005, 3)), betaC = c(rep(0.005, 3)), betaGY = c(rep(0.2, 10)), betaGX = c(rep(0, 10)), deltaGX = c(rep(0.2, 10)), betaX = 0.2, sdX = 1, sdY = 1, SEED = 1, sig.level = 0.05, nSims = 5, table.name = "MRcor", pFilterMRCD = 5e-08, pFilterBDC = 0.05/20, approaches = c("All"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
%%     ~~Describe \code{n} here~~
}
  \item{nSNPx}{
%%     ~~Describe \code{nSNPx} here~~
}
  \item{MAFx}{
%%     ~~Describe \code{MAFx} here~~
}
  \item{nSNPy}{
%%     ~~Describe \code{nSNPy} here~~
}
  \item{MAFy}{
%%     ~~Describe \code{MAFy} here~~
}
  \item{nCov}{
%%     ~~Describe \code{nCov} here~~
}
  \item{meanC}{
%%     ~~Describe \code{meanC} here~~
}
  \item{sdC}{
%%     ~~Describe \code{sdC} here~~
}
  \item{deltaC}{
%%     ~~Describe \code{deltaC} here~~
}
  \item{betaC}{
%%     ~~Describe \code{betaC} here~~
}
  \item{betaGY}{
%%     ~~Describe \code{betaGY} here~~
}
  \item{betaGX}{
%%     ~~Describe \code{betaGX} here~~
}
  \item{deltaGX}{
%%     ~~Describe \code{deltaGX} here~~
}
  \item{betaX}{
%%     ~~Describe \code{betaX} here~~
}
  \item{sdX}{
%%     ~~Describe \code{sdX} here~~
}
  \item{sdY}{
%%     ~~Describe \code{sdY} here~~
}
  \item{SEED}{
%%     ~~Describe \code{SEED} here~~
}
  \item{sig.level}{
%%     ~~Describe \code{sig.level} here~~
}
  \item{nSims}{
%%     ~~Describe \code{nSims} here~~
}
  \item{table.name}{
%%     ~~Describe \code{table.name} here~~
}
  \item{pFilterMRCD}{
%%     ~~Describe \code{pFilterMRCD} here~~
}
  \item{pFilterBDC}{
%%     ~~Describe \code{pFilterBDC} here~~
}
  \item{approaches}{
%%     ~~Describe \code{approaches} here~~
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
MRcorDiff(n = 1000, nSNPx = 10, MAFx = c(rep(0.5, 10)), nSNPy = 10, MAFy = c(rep(0.5, 10)), nCov = 3, meanC = rep(0, 3), sdC = c(rep(1, 3)), deltaC = c(rep(0.005, 3)), betaC = c(rep(0.005, 3)), betaGY = c(rep(0.2, 10)), betaGX = c(rep(0, 10)), deltaGX = c(rep(0.2, 10)), betaX = 0.2, sdX = 1, sdY = 1, SEED = 1, sig.level = 0.05, nSims = 5, table.name = "MRcor", pFilterMRCD = 5e-08, pFilterBDC = 0.05/20, approaches = c("All"))

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
