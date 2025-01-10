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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or standard data sets, see data().

## The function is currently defined as
function (n = 1e+05, nSNP1 = 10, MAF1 = c(rep(0.5, 10)), nSNP2 = 10, 
    MAF2 = c(rep(0.5, 10)), nCov = 3, meanC = rep(0, 3), sdC = c(rep(1, 
        3)), deltaC = c(rep(0.005, 3)), betaC = c(rep(0.005, 
        3)), betaG2 = c(rep(0.2, 10)), betaG1 = c(rep(0, 10)), 
    deltaG1 = c(rep(0.2, 10)), betaX = 0.2, sdX = 1, sdY = 1, 
    SEED = 1, sig.level = 0.05, nSims = 5) 
{
    if (length(MAF1) != nSNP1) {
        stop("length(MAF1) must equal nSNP1.")
    }
    if (length(MAF2) != nSNP2) {
        stop("length(MAF2) must equal nSNP2")
    }
    if (length(deltaG1) != nSNP1) {
        stop("length(deltaG1) must equal nSNP1")
    }
    if (length(betaG2) != nSNP2) {
        stop("length(betaG2) must equal nSNP2")
    }
    if (length(betaG1) != nSNP2) {
        stop("length(betaG1) must equal nSNP2")
    }
    if (length(meanC) != nCov) {
        stop("length(meanC) must equal nCov")
    }
    if (length(sdC) != nCov) {
        stop("length(sdC) must equal nCov")
    }
    if (length(deltaC) != nCov) {
        stop("ncol(deltaC) must equal nCov")
    }
    if (length(betaC) != nCov) {
        stop("length(betaC) must equal nCov")
    }
    colnames.matR <- c("corGX", "corGXformulaMLR", "corGXresid", 
        "corGY", "corGYformulaMLR", "corGYresid")
    matR <- matrix(0, ncol = length(colnames.matR), nrow = nSims)
    colnames(matR) <- colnames.matR
    colnames.matBeta <- c("corGX", "corGXformulaMLR", "corGXresid", 
        "corGY", "corGYformulaMLR", "corGYresid")
    matBeta <- matrix(0, ncol = length(colnames.matBeta), nrow = nSims)
    colnames(matBeta) <- colnames.matBeta
    colnames.matSE <- c("corGX", "corGXformulaMLR", "corGXresid", 
        "corGY", "corGYformulaMLR", "corGYresid")
    matSE <- matrix(0, ncol = length(colnames.matSE), nrow = nSims)
    colnames(matSE) <- colnames.matSE
    for (i in 1:nSims) {
        printCut <- 10
        if (floor(i/printCut) == ceiling(i/printCut)) {
            print(paste(i, "of", nSims, "simulations"))
        }
        set.seed(SEED + i - 1)
        deltaG1 <- matrix(deltaG1, ncol = 1, nrow = nSNP1)
        betaC <- matrix(betaC, ncol = 1, nrow = nCov)
        betaG1 <- matrix(betaG1, ncol = 1, nrow = nSNP2)
        betaG2 <- matrix(betaG2, ncol = 1, nrow = nSNP2)
        deltaC <- matrix(deltaC, ncol = 1, nrow = nCov)
        matG1 <- matrix(0, ncol = nSNP1, nrow = n)
        for (ss in 1:nSNP1) {
            matG1[, ss] <- rbinom(n, 2, MAF1[ss])
        }
        matG2 <- matrix(0, ncol = nSNP2, nrow = n)
        for (ss in 1:nSNP2) {
            matG2[, ss] <- rbinom(n, 2, MAF2[ss])
        }
        matC <- matrix(0, ncol = nCov, nrow = n)
        for (cc in 1:nCov) {
            matC[, cc] <- rnorm(n, meanC[cc], sdC[cc])
        }
        x <- rnorm(n, (matG1 \%*\% deltaG1 + matC \%*\% deltaC), 
            sdX)
        y <- rnorm(n, (matG2 \%*\% betaG2 + matG1 \%*\% betaG1 + 
            matC \%*\% betaC + betaX * x), sdY)
        modelRx <- lm(x ~ matC)
        xr <- modelRx$res
        modelRy <- lm(y ~ matC)
        yr <- modelRy$res
        modelXS <- summary(lm(x ~ matG1[, 1]))
        modelXS_Coef <- modelXS$coef[2, c(1, 2)]
        matBeta[i, "corGX"] <- modelXS_Coef[1]
        matSE[i, "corGX"] <- modelXS_Coef[2]
        modelYS <- summary(lm(y ~ matG1[, 1]))
        modelYS_Coef <- modelYS$coef[2, c(1, 2)]
        matBeta[i, "corGY"] <- modelYS_Coef[1]
        matSE[i, "corGY"] <- modelYS_Coef[2]
        modelxr <- summary(lm(xr ~ matG1[, 1]))
        modelxr_Coef <- modelxr$coef[2, c(1, 2)]
        matBeta[i, "corGXresid"] <- modelxr_Coef[1]
        matSE[i, "corGXresid"] <- modelxr_Coef[2]
        modelyr <- summary(lm(yr ~ matG1[, 1]))
        modelyr_Coef <- modelyr$coef[2, c(1, 2)]
        matBeta[i, "corGYresid"] <- modelyr_Coef[1]
        matSE[i, "corGYresid"] <- modelyr_Coef[2]
        modelXM <- summary(lm(x ~ matG1[, 1] + matC))
        modelXM_Coef <- modelXM$coef[2, c(1, 2)]
        matBeta[i, "corGXformulaMLR"] <- modelXM_Coef[1]
        matSE[i, "corGXformulaMLR"] <- modelXM_Coef[2]
        modelYM <- summary(lm(y ~ matG1[, 1] + matC))
        modelYM_Coef <- modelYM$coef[2, c(1, 2)]
        matBeta[i, "corGYformulaMLR"] <- modelYM_Coef[1]
        matSE[i, "corGYformulaMLR"] <- modelYM_Coef[2]
        corGXslr <- modelXS_Coef[1]/(sqrt(modelXS_Coef[1]^2 + 
            (n - 2) * modelXS_Coef[2]^2))
        matR[i, "corGX"] <- corGXslr
        corGXmlr <- modelXM_Coef[1]/(sqrt(modelXM_Coef[1]^2 + 
            (n - 2) * modelXM_Coef[2]^2))
        matR[i, "corGXformulaMLR"] <- corGXmlr
        corGYslr <- modelYS_Coef[1]/(sqrt(modelYS_Coef[1]^2 + 
            (n - 2) * modelYS_Coef[2]^2))
        matR[i, "corGY"] <- corGYslr
        corGYmlr <- modelYM_Coef[1]/(sqrt(modelYM_Coef[1]^2 + 
            (n - 2) * modelYM_Coef[2]^2))
        matR[i, "corGYformulaMLR"] <- corGYmlr
        corGXslrRes <- modelxr_Coef[1]/(sqrt(modelxr_Coef[1]^2 + 
            (n - 2) * modelxr_Coef[2]^2))
        matR[i, "corGXresid"] <- corGXslrRes
        corGYslrRes <- modelyr_Coef[1]/(sqrt(modelyr_Coef[1]^2 + 
            (n - 2) * modelyr_Coef[2]^2))
        matR[i, "corGYresid"] <- corGYslrRes
    }
    matCorDiff <- matrix(NA, nrow = 2, ncol = 3)
    colnames(matCorDiff) <- c("min", "mean", "max")
    rownames(matCorDiff) <- c("corGXformulaMLR-corGXresid", "corGYformulaMLR-corGYresid")
    matCorDiff[1, "min"] <- min(matR[, "corGXformulaMLR"] - matR[, 
        "corGXresid"])
    matCorDiff[1, "max"] <- max(matR[, "corGXformulaMLR"] - matR[, 
        "corGXresid"])
    matCorDiff[1, "mean"] <- mean(matR[, "corGXformulaMLR"] - 
        matR[, "corGXresid"])
    matCorDiff[2, "min"] <- min(matR[, "corGYformulaMLR"] - matR[, 
        "corGYresid"])
    matCorDiff[2, "max"] <- max(matR[, "corGYformulaMLR"] - matR[, 
        "corGYresid"])
    matCorDiff[2, "mean"] <- mean(matR[, "corGYformulaMLR"] - 
        matR[, "corGYresid"])
    write.table(matCorDiff, paste0("corrCompareMatDiffCor_nSims", 
        nSims, "n", n, "nSNP1_", nSNP1, "nSNP2_", nSNP2, "maf1_", 
        MAF1[1] * 100, "maf2_", MAF2[1] * 100, "nCov", nCov, 
        "betaX", betaX * 100, "deltaC1_", deltaC[1] * 100, "betaC1_", 
        betaC[1] * 100, "betaG1_", betaG1[1] * 100, "betaG2_", 
        betaG2[1] * 100, "deltaG1_", deltaG1[1] * 100, ".txt"), 
        quote = F, row.names = T)
    colnames.matMeanCor <- c("min", "mean", "max")
    rownames.matMeanCor <- c("corGX", "corGXresid", "corGXformulaMLR", 
        "corGX-corGXresid", "corGX-corGXformulaMLR", "corGXformulaMLR-corGXresid", 
        "corGY", "corGYresid", "corGYformulaMLR", "corGY-corGYresid", 
        "corGY-corGYformulaMLR", "corGYformulaMLR-corGYresid")
    matMeanCor <- matrix(NA, nrow = length(rownames.matMeanCor), 
        ncol = length(colnames.matMeanCor))
    colnames(matMeanCor) <- colnames.matMeanCor
    rownames(matMeanCor) <- rownames.matMeanCor
    matMeanCor["corGX", "min"] <- min(matR[, "corGX"])
    matMeanCor["corGX", "max"] <- max(matR[, "corGX"])
    matMeanCor["corGX", "mean"] <- mean(matR[, "corGX"])
    matMeanCor["corGXresid", "min"] <- min(matR[, "corGXresid"])
    matMeanCor["corGXresid", "max"] <- max(matR[, "corGXresid"])
    matMeanCor["corGXresid", "mean"] <- mean(matR[, "corGXresid"])
    matMeanCor["corGXformulaMLR", "min"] <- min(matR[, "corGXformulaMLR"])
    matMeanCor["corGXformulaMLR", "max"] <- max(matR[, "corGXformulaMLR"])
    matMeanCor["corGXformulaMLR", "mean"] <- mean(matR[, "corGXformulaMLR"])
    matMeanCor["corGX-corGXresid", "min"] <- min(matR[, "corGX"] - 
        matR[, "corGXresid"])
    matMeanCor["corGX-corGXresid", "max"] <- max(matR[, "corGX"] - 
        matR[, "corGXresid"])
    matMeanCor["corGX-corGXresid", "mean"] <- mean(matR[, "corGX"] - 
        matR[, "corGXresid"])
    matMeanCor["corGX-corGXformulaMLR", "min"] <- min(matR[, 
        "corGX"] - matR[, "corGXformulaMLR"])
    matMeanCor["corGX-corGXformulaMLR", "max"] <- max(matR[, 
        "corGX"] - matR[, "corGXformulaMLR"])
    matMeanCor["corGX-corGXformulaMLR", "mean"] <- mean(matR[, 
        "corGX"] - matR[, "corGXformulaMLR"])
    matMeanCor["corGXformulaMLR-corGXresid", "min"] <- min(matR[, 
        "corGXformulaMLR"] - matR[, "corGXresid"])
    matMeanCor["corGXformulaMLR-corGXresid", "max"] <- max(matR[, 
        "corGXformulaMLR"] - matR[, "corGXresid"])
    matMeanCor["corGXformulaMLR-corGXresid", "mean"] <- mean(matR[, 
        "corGXformulaMLR"] - matR[, "corGXresid"])
    matMeanCor["corGY", "min"] <- min(matR[, "corGY"])
    matMeanCor["corGY", "max"] <- max(matR[, "corGY"])
    matMeanCor["corGY", "mean"] <- mean(matR[, "corGY"])
    matMeanCor["corGYresid", "min"] <- min(matR[, "corGYresid"])
    matMeanCor["corGYresid", "max"] <- max(matR[, "corGYresid"])
    matMeanCor["corGYresid", "mean"] <- mean(matR[, "corGYresid"])
    matMeanCor["corGYformulaMLR", "min"] <- min(matR[, "corGYformulaMLR"])
    matMeanCor["corGYformulaMLR", "max"] <- max(matR[, "corGYformulaMLR"])
    matMeanCor["corGYformulaMLR", "mean"] <- mean(matR[, "corGYformulaMLR"])
    matMeanCor["corGY-corGYresid", "min"] <- min(matR[, "corGY"] - 
        matR[, "corGYresid"])
    matMeanCor["corGY-corGYresid", "max"] <- max(matR[, "corGY"] - 
        matR[, "corGYresid"])
    matMeanCor["corGY-corGYresid", "mean"] <- mean(matR[, "corGY"] - 
        matR[, "corGYresid"])
    matMeanCor["corGY-corGYformulaMLR", "min"] <- min(matR[, 
        "corGY"] - matR[, "corGYformulaMLR"])
    matMeanCor["corGY-corGYformulaMLR", "max"] <- max(matR[, 
        "corGY"] - matR[, "corGYformulaMLR"])
    matMeanCor["corGY-corGYformulaMLR", "mean"] <- mean(matR[, 
        "corGY"] - matR[, "corGYformulaMLR"])
    matMeanCor["corGYformulaMLR-corGYresid", "min"] <- min(matR[, 
        "corGYformulaMLR"] - matR[, "corGYresid"])
    matMeanCor["corGYformulaMLR-corGYresid", "max"] <- max(matR[, 
        "corGYformulaMLR"] - matR[, "corGYresid"])
    matMeanCor["corGYformulaMLR-corGYresid", "mean"] <- mean(matR[, 
        "corGYformulaMLR"] - matR[, "corGYresid"])
    write.table(matMeanCor, file = paste0("corrCompareMatCor_nSims", 
        nSims, "n", n, "nSNP1_", nSNP1, "nSNP2_", nSNP2, "maf1_", 
        MAF1[1] * 100, "maf2_", MAF2[1] * 100, "nCov", nCov, 
        "betaX", betaX * 100, "deltaC1_", deltaC[1] * 100, "betaC1_", 
        betaC[1] * 100, "betaG1_", betaG1[1] * 100, "betaG2_", 
        betaG2[1] * 100, "deltaG1_", deltaG1[1] * 100, ".txt"), 
        quote = F, row.names = T)
    colnames.matMeanBeta <- c("betaMin", "betaMean", "betaMax", 
        "seMin", "seMean", "seMax")
    rownames.matMeanBeta <- c("corGX", "corGXresid", "corGXformulaMLR", 
        "corGY", "corGYresid", "corGYformulaMLR")
    matMeanBeta <- matrix(NA, nrow = length(rownames.matMeanBeta), 
        ncol = length(colnames.matMeanBeta))
    colnames(matMeanBeta) <- colnames.matMeanBeta
    rownames(matMeanBeta) <- rownames.matMeanBeta
    matMeanBeta["corGX", "betaMin"] <- min(matBeta[, "corGX"])
    matMeanBeta["corGX", "betaMax"] <- max(matBeta[, "corGX"])
    matMeanBeta["corGX", "betaMean"] <- mean(matBeta[, "corGX"])
    matMeanBeta["corGX", "seMin"] <- min(matSE[, "corGX"])
    matMeanBeta["corGX", "seMax"] <- max(matSE[, "corGX"])
    matMeanBeta["corGX", "seMean"] <- mean(matSE[, "corGX"])
    matMeanBeta["corGXresid", "betaMin"] <- min(matBeta[, "corGXresid"])
    matMeanBeta["corGXresid", "betaMax"] <- max(matBeta[, "corGXresid"])
    matMeanBeta["corGXresid", "betaMean"] <- mean(matBeta[, "corGXresid"])
    matMeanBeta["corGXresid", "seMin"] <- min(matSE[, "corGXresid"])
    matMeanBeta["corGXresid", "seMax"] <- max(matSE[, "corGXresid"])
    matMeanBeta["corGXresid", "seMean"] <- mean(matSE[, "corGXresid"])
    matMeanBeta["corGXformulaMLR", "betaMin"] <- min(matBeta[, 
        "corGXformulaMLR"])
    matMeanBeta["corGXformulaMLR", "betaMax"] <- max(matBeta[, 
        "corGXformulaMLR"])
    matMeanBeta["corGXformulaMLR", "betaMean"] <- mean(matBeta[, 
        "corGXformulaMLR"])
    matMeanBeta["corGXformulaMLR", "seMin"] <- min(matSE[, "corGXformulaMLR"])
    matMeanBeta["corGXformulaMLR", "seMax"] <- max(matSE[, "corGXformulaMLR"])
    matMeanBeta["corGXformulaMLR", "seMean"] <- mean(matSE[, 
        "corGXformulaMLR"])
    matMeanBeta["corGY", "betaMin"] <- min(matBeta[, "corGY"])
    matMeanBeta["corGY", "betaMax"] <- max(matBeta[, "corGY"])
    matMeanBeta["corGY", "betaMean"] <- mean(matBeta[, "corGY"])
    matMeanBeta["corGY", "seMin"] <- min(matSE[, "corGY"])
    matMeanBeta["corGY", "seMax"] <- max(matSE[, "corGY"])
    matMeanBeta["corGY", "seMean"] <- mean(matSE[, "corGY"])
    matMeanBeta["corGYresid", "betaMin"] <- min(matBeta[, "corGYresid"])
    matMeanBeta["corGYresid", "betaMax"] <- max(matBeta[, "corGYresid"])
    matMeanBeta["corGYresid", "betaMean"] <- mean(matBeta[, "corGYresid"])
    matMeanBeta["corGYresid", "seMin"] <- min(matSE[, "corGYresid"])
    matMeanBeta["corGYresid", "seMax"] <- max(matSE[, "corGYresid"])
    matMeanBeta["corGYresid", "seMean"] <- mean(matSE[, "corGYresid"])
    matMeanBeta["corGYformulaMLR", "betaMin"] <- min(matBeta[, 
        "corGYformulaMLR"])
    matMeanBeta["corGYformulaMLR", "betaMax"] <- max(matBeta[, 
        "corGYformulaMLR"])
    matMeanBeta["corGYformulaMLR", "betaMean"] <- mean(matBeta[, 
        "corGYformulaMLR"])
    matMeanBeta["corGYformulaMLR", "seMin"] <- min(matSE[, "corGYformulaMLR"])
    matMeanBeta["corGYformulaMLR", "seMax"] <- max(matSE[, "corGYformulaMLR"])
    matMeanBeta["corGYformulaMLR", "seMean"] <- mean(matSE[, 
        "corGYformulaMLR"])
    write.table(matMeanBeta, file = paste0("corrCompareMatBetaSE_nSims", 
        nSims, "n", n, "nSNP1_", nSNP1, "nSNP2_", nSNP2, "maf1_", 
        MAF1[1] * 100, "maf2_", MAF2[1] * 100, "nCov", nCov, 
        "betaX", betaX * 100, "deltaC1_", deltaC[1] * 100, "betaC1_", 
        betaC[1] * 100, "betaG1_", betaG1[1] * 100, "betaG2_", 
        betaG2[1] * 100, "deltaG1_", deltaG1[1] * 100, ".txt"), 
        quote = F, row.names = T)
    pdf(paste("corrCompare_GX_CorPlot_nSims", nSims, "n", n, 
        "nSNP1_", nSNP1, "nSNP2_", nSNP2, "maf1_", MAF1[1] * 
            100, "maf2_", MAF2[1] * 100, "nCov", nCov, "betaX", 
        betaX * 100, "deltaC1_", deltaC[1] * 100, "betaC1_", 
        betaC[1] * 100, "betaG1_", betaG1[1] * 100, "betaG2_", 
        betaG2[1] * 100, "deltaG1_", deltaG1[1] * 100, ".pdf", 
        sep = ""))
    boxplot(matR[, c("corGX", "corGXresid", "corGXformulaMLR")], 
        names = c("cor(G,X)", "cor(G,Xresid)", "Formula:MLR"), 
        cex.axis = 0.75)
    dev.off()
    pdf(paste("corrCompare_GY_CorPlot_nSims", nSims, "n", n, 
        "nSNP1_", nSNP1, "nSNP2_", nSNP2, "maf1_", MAF1[1] * 
            100, "maf2_", MAF2[1] * 100, "nCov", nCov, "betaX", 
        betaX * 100, "deltaC1_", deltaC[1] * 100, "betaC1_", 
        betaC[1] * 100, "betaG1_", betaG1[1] * 100, "betaG2_", 
        betaG2[1] * 100, "deltaG1_", deltaG1[1] * 100, ".pdf", 
        sep = ""))
    boxplot(matR[, c("corGY", "corGYresid", "corGYformulaMLR")], 
        names = c("cor(G,Y)", "cor(G,Yresid)", "Formula:MLR"), 
        cex.axis = 0.75)
    dev.off()
    pdf(paste("corrCompare_GX_DiffPlot_nSims", nSims, "n", n, 
        "nSNP1_", nSNP1, "nSNP2_", nSNP2, "maf1_", MAF1[1] * 
            100, "maf2_", MAF2[1] * 100, "nCov", nCov, "betaX", 
        betaX * 100, "deltaC1_", deltaC[1] * 100, "betaC1_", 
        betaC[1] * 100, "betaG1_", betaG1[1] * 100, "betaG2_", 
        betaG2[1] * 100, "deltaG1_", deltaG1[1] * 100, ".pdf", 
        sep = ""))
    boxplot(matR[, "corGX"] - matR[, "corGXresid"], matR[, "corGX"] - 
        matR[, "corGXformulaMLR"], matR[, "corGXformulaMLR"] - 
        matR[, "corGXresid"], names = c("cor(G,X)-cor(G,Xresid)", 
        "cor(G,X)-Formula:MLR", "Formula:MLR-cor(G,Xresid)"), 
        cex.axis = 0.65)
    dev.off()
    pdf(paste("corrCompare_GY_DiffPlot_nSims", nSims, "n", n, 
        "nSNP1_", nSNP1, "nSNP2_", nSNP2, "maf1_", MAF1[1] * 
            100, "maf2_", MAF2[1] * 100, "nCov", nCov, "betaX", 
        betaX * 100, "deltaC1_", deltaC[1] * 100, "betaC1_", 
        betaC[1] * 100, "betaG1_", betaG1[1] * 100, "betaG2_", 
        betaG2[1] * 100, "deltaG1_", deltaG1[1] * 100, ".pdf", 
        sep = ""))
    boxplot(matR[, "corGY"] - matR[, "corGYresid"], matR[, "corGY"] - 
        matR[, "corGYformulaMLR"], matR[, "corGYformulaMLR"] - 
        matR[, "corGYresid"], names = c("cor(G,Y)-cor(G,Yresid)", 
        "cor(G,Y)-Formula:MLR", "Formula:MLR-cor(G,Yresid)"), 
        cex.axis = 0.65)
    dev.off()
    return(list(matMeanCor = matMeanCor, matMeanBeta = matMeanBeta, 
        matCorDiff = matCorDiff))
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
