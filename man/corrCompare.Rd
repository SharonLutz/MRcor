\name{corrCompare}
\alias{corrCompare}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A Capitalized Title (ideally limited to 65 characters)
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
corrCompare(n = 1e+05, nSNP1 = 10, MAF1 = c(rep(0.5, 10)), nSNP2 = 10, MAF2 = c(rep(0.5, 10)), nCov = 3, meanC = rep(0, 3), sdC = c(rep(1, 3)), deltaC = c(rep(0.005, 3)), betaC = c(rep(0.005, 3)), betaG2 = c(rep(0.2, 10)), betaG1 = c(rep(0, 10)), deltaG1 = c(rep(0.2, 10)), betaX = 0.2, sdX = 1, sdY = 1, SEED = 1, sig.level = 0.05, nSims = 5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
%%     ~~Describe \code{n} here~~
}
  \item{nSNP1}{
%%     ~~Describe \code{nSNP1} here~~
}
  \item{MAF1}{
%%     ~~Describe \code{MAF1} here~~
}
  \item{nSNP2}{
%%     ~~Describe \code{nSNP2} here~~
}
  \item{MAF2}{
%%     ~~Describe \code{MAF2} here~~
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
  \item{betaG2}{
%%     ~~Describe \code{betaG2} here~~
}
  \item{betaG1}{
%%     ~~Describe \code{betaG1} here~~
}
  \item{deltaG1}{
%%     ~~Describe \code{deltaG1} here~~
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
corrCompare(n = 1000, nSNP1 = 10, MAF1 = c(rep(0.5, 10)), nSNP2 = 10, MAF2 = c(rep(0.5, 10)), nCov = 3, meanC = rep(0, 3), sdC = c(rep(1, 3)), deltaC = c(rep(0.005, 3)), betaC = c(rep(0.005, 3)), betaG2 = c(rep(0.2, 10)), betaG1 = c(rep(0, 10)), deltaG1 = c(rep(0.2, 10)), betaX = 0.2, sdX = 1, sdY = 1, SEED = 1, sig.level = 0.05, nSims = 5)

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
