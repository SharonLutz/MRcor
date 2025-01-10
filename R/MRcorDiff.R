MRcorDiff <-
function(n=100000, nSNPx = 10, MAFx = c(rep(0.5,  10)), nSNPy = 10, MAFy = c(rep(0.5, 10)),
           nCov = 3, meanC = rep(0, 3), sdC = c(rep(1, 3)), deltaC =  c(rep(0.005,3)), betaC=c(rep(0.005,3)),
           betaGY = c(rep(0.2, 10)), betaGX = c(rep(0, 10)), deltaGX = c(rep(0.2, 10)), betaX=0.2,
           sdX=1,sdY=1, SEED=1, sig.level=0.05, nSims=5, table.name="MRcor", pFilterMRCD=5E-8, pFilterBDC=0.05/20,
           approaches=c("All"), libPath=NULL, plot.pdf=F){
    
    if(length(MAFx)!=nSNPx){stop("length(MAFx) must equal nSNPx.")}
    if(length(MAFy)!=nSNPy){stop("length(MAFy) must equal nSNPy")}
    if(length(deltaGX)!=nSNPx){stop("length(deltaGX) must equal nSNPx")}
    if(length(betaGY)!=nSNPy){stop("length(betaGY) must equal nSNPy")}
    if(length(betaGX)!=nSNPy){stop("length(betaGX) must equal nSNPy")}
    
    if(length(meanC)!=nCov){stop("length(meanC) must equal nCov")}
    if(length(sdC)!=nCov){stop("length(sdC) must equal nCov")}
    if(length(deltaC)!=nCov){stop("ncol(deltaC) must equal nCov")}
    if(length(betaC)!=nCov){stop("length(betaC) must equal nCov")}
    
    
    if(is.null(libPath)){
      library(MRCD)
      library(MRcML)
      library(BiDirectCausal)
      library(TwoSampleMR)
      library(pals)
    }
    if(!is.null(libPath)){
      library(MRCD, lib.loc=libPath)
      library(MRcML, lib.loc=libPath)
      library(BiDirectCausal, lib.loc=libPath)
      library(TwoSampleMR, lib.loc=libPath)
      library(pals, lib.loc=libPath)
    }
    
    
    ################################################################################
    # Matrix to save Results
    ################################################################################
    #save results for type 1 error rate betaX=0 and power betaX>0
    All <- c("CDcML", "MRcML", "CDRatio", "CDEgger", "CDGLS",
             "CDRatio.slrF", "CDEgger.slrF", "CDGLS.slrF", 
             "CDRatio.mlrF", "CDEgger.mlrF", "CDGLS.mlrF",
             "cdRatio.S", "cdEgger.S", "cdRatio.NoS", "cdEgger.NoS",
             "MRS.IVW", "MRS.wMedian", "MRS.Egger")
    
    colnames.1 <- NULL
    if("All" %in% approaches){colnames.1 <- c(colnames.1, All)}
    if("CDcML" %in% approaches){colnames.1 <- c(colnames.1, "CDcML")}
    if("MRcML" %in% approaches){colnames.1 <- c(colnames.1, "MRcML")}
    
    if("CDRatio" %in% approaches){colnames.1 <- c(colnames.1, "CDRatio")}
    if("CDEgger" %in% approaches){colnames.1 <- c(colnames.1, "CDEgger")}
    if("CDGLS" %in% approaches){colnames.1 <- c(colnames.1, "CDGLS")}
    
    if("CDRatio.slrF" %in% approaches){colnames.1 <- c(colnames.1, "CDRatio.slrF")}
    if("CDEgger.slrF" %in% approaches){colnames.1 <- c(colnames.1, "CDEgger.slrF")}
    if("CDGLS.slrF" %in% approaches){colnames.1 <- c(colnames.1, "CDGLS.slrF")}
    
    if("CDRatio.mlrF" %in% approaches){colnames.1 <- c(colnames.1, "CDRatio.mlrF")}
    if("CDEgger.mlrF" %in% approaches){colnames.1 <- c(colnames.1, "CDEgger.mlrF")}
    if("CDGLS.mlrF" %in% approaches){colnames.1 <- c(colnames.1, "CDGLS.mlrF")}
    
    if("cdRatio.S" %in% approaches){colnames.1 <- c(colnames.1, "cdRatio.S")}
    if("cdEgger.S" %in% approaches){colnames.1 <- c(colnames.1, "cdEgger.S")}
    
    if("cdRatio.NoS" %in% approaches){colnames.1 <- c(colnames.1, "cdRatio.NoS")}
    if("cdEgger.NoS" %in% approaches){colnames.1 <- c(colnames.1, "cdEgger.NoS")}
    
    if("MRS.IVW" %in% approaches){colnames.1 <- c(colnames.1, "MRS.IVW")}
    if("MRS.wMedian" %in% approaches){colnames.1 <- c(colnames.1, "MRS.wMedian")}
    if("MRS.Egger" %in% approaches){colnames.1 <- c(colnames.1, "MRS.Egger")}
    
    # make sure each method is only one column - no duplicate column names
    colnames.mat1 <- unique(colnames.1)
    
    matR <- matrix(0, ncol=length(colnames.mat1), nrow=29)
    
    # expand this out to include res1,res2,res3
    colnames(matR) <- colnames.mat1
    rownames(matR) <- c("slr1_res1_mlr1", "slr1_res1_mlr2", "slr1_res1_mlr3",
                        "slr1_res2_mlr1", "slr1_res2_mlr2", "slr1_res2_mlr3",
                        "slr1_res3_mlr1", "slr1_res3_mlr2", "slr1_res3_mlr3",
                        
                        "slr2_res1_mlr1", "slr2_res1_mlr2", "slr2_res1_mlr3",
                        "slr2_res2_mlr1", "slr2_res2_mlr2", "slr2_res2_mlr3",
                        "slr2_res3_mlr1", "slr2_res3_mlr2", "slr2_res3_mlr3",
                        
                        "slr3_res1_mlr1", "slr3_res1_mlr2", "slr3_res1_mlr3",
                        "slr3_res2_mlr1", "slr3_res2_mlr2", "slr3_res2_mlr3",
                        "slr3_res3_mlr1", "slr3_res3_mlr2", "slr3_res3_mlr3",
                        "allFiltered.MLR", "allFiltered.SLR")
    
    #storing matrices/vectors
    corGXms<-rep(0,nSims)
    corGYms<-rep(0,nSims)
    corGXsr<-rep(0,nSims)
    corGYsr<-rep(0,nSims)
    corGXmr<-rep(0,nSims)
    corGYmr<-rep(0,nSims)
    
    
    for(i in 1:nSims){
      printCut<-1
      if(floor(i/printCut)==ceiling(i/printCut)){print(paste(i,"of",nSims, "simulations"))}
      
      set.seed(SEED+i-1)
      
      ###################################
      # matrix for deltaGX, deltaC, betaGX, betaGY, betaC
      deltaGX<-matrix(deltaGX, ncol=1, nrow=nSNPx)
      betaC<-matrix(betaC, ncol=1, nrow=nCov)
      betaGX<- matrix(betaGX,ncol=1,nrow=nSNPy)
      betaGY<- matrix(betaGY,ncol=1,nrow=nSNPy)
      deltaC<-matrix(deltaC, ncol=1, nrow=nCov)
      
      # create SNP matrix for X
      matGX<-matrix(0,ncol=nSNPx,nrow=n)
      for(ss in 1:nSNPx){
        matGX[,ss]<-rbinom(n,2, MAFx[ss])
      }
      
      # create SNP matrix for Y
      matGY<-matrix(0,ncol=nSNPy,nrow=n)
      for(ss in 1:nSNPy){
        matGY[,ss]<-rbinom(n,2, MAFy[ss])
      }
      
      # create covariate matrix
      matC<-matrix(0, ncol=nCov,nrow=n)
      for(cc in 1:nCov){
        matC[,cc]<-rnorm(n, meanC[cc], sdC[cc])
      }
      
      x <- rnorm(n, (matGX%*%deltaGX + matC%*%deltaC), sdX) # trait X
      y <- rnorm(n, (matGY%*%betaGY + matGX%*%betaGX+ matC%*%betaC+betaX*x), sdY)
      
      # residuals
      modelRx <- lm(x~matC)
      xr <- modelRx$res
      
      modelRy <- lm(y~matC)
      yr <- modelRy$res
      
      ###################################
      # correlations
      ###################################
      # simple linear regression (SLR)
      modelXS<-summary(lm(x~matGX[,1]))$coef[2,c(1,2)]
      modelYS<-summary(lm(y~matGX[,1]))$coef[2,c(1,2)]
      
      modelxr<-summary(lm(xr~matGX[,1]))$coef[2,c(1,2)]
      modelyr<-summary(lm(yr~matGX[,1]))$coef[2,c(1,2)]
      
      #multiple linear regression (MLR)
      modelXM<-summary(lm(x~matGX[,1]+matC))$coef[2,c(1,2)]
      modelYM<-summary(lm(y~matGX[,1]+matC))$coef[2,c(1,2)]
      
      # correlation
      corGXslr<-modelXS[1]/(sqrt( modelXS[1]^2+(n-2)* modelXS[2]^2 ))
      corGXmlr<-modelXM[1]/(sqrt( modelXM[1]^2+(n-2)* modelXM[2]^2 ))
      
      corGYslr<-modelYS[1]/(sqrt( modelYS[1]^2+(n-2)* modelYS[2]^2 ))
      corGYmlr<-modelYM[1]/(sqrt( modelYM[1]^2+(n-2)* modelYM[2]^2 ))
      
      corGXslrRes<-modelxr[1]/(sqrt( modelxr[1]^2+(n-2)* modelxr[2]^2 ))
      corGYslrRes<-modelyr[1]/(sqrt( modelyr[1]^2+(n-2)* modelyr[2]^2 ))
      
      # save difference
      corGXms[i]<-corGXmlr - corGXslr
      corGYms[i]<-corGYmlr - corGYslr
      
      corGXsr[i] <- corGXslr-corGXslrRes
      corGYsr[i] <- corGYslr-corGYslrRes
      
      corGXmr[i] <- corGXmlr-corGXslrRes
      corGYmr[i] <- corGYmlr-corGYslrRes
      
      ###################################
      # GX --> X --> Y
      bxA <- NULL
      sxA <- NULL
      pxA <- NULL
      
      byA <- NULL
      syA <- NULL
      pyA <- NULL
      
      for(ll in 1:ncol(matGX)){
        bxA <- c(bxA, summary(lm(x~matGX[,ll]+matC))$coef[2,1])
        sxA <- c(sxA, summary(lm(x~matGX[,ll]+matC))$coef[2,2])
        pxA <- c(pxA, summary(lm(x~matGX[,ll]+matC))$coef[2,4])
        
        byA <- c(byA, summary(lm(y~matGX[,ll]+matC))$coef[2,1])
        syA <- c(syA, summary(lm(y~matGX[,ll]+matC))$coef[2,2])
        pyA <- c(pyA, summary(lm(y~matGX[,ll]+matC))$coef[2,4])
      }
      
      # save for cd-ratio, cd-gls, cd-egger
     # a.betaGX_X <- bxA[pxA < pFilterMRCD & pyA < pFilterMRCD]
      a.betaGX_X <- bxA
      a.seGX_X <- sxA
      
      a.betaGX_Y <- byA
      a.seGX_Y <- syA
      
      for(ll in 1:ncol(matGY)){
        bxA <- c(bxA, summary(lm(x~matGY[,ll]+matC))$coef[2,1])
        sxA <- c(sxA, summary(lm(x~matGY[,ll]+matC))$coef[2,2])
        
        byA <- c(byA, summary(lm(y~matGY[,ll]+matC))$coef[2,1])
        syA <- c(syA, summary(lm(y~matGY[,ll]+matC))$coef[2,2])
      }
      
      bxU <- NULL
      sxU <- NULL
      pxU <- NULL
      
      byU <- NULL
      syU <- NULL
      pyU <- NULL
      
      for(ll in 1:ncol(matGX)){
        bxU <- c(bxU, summary(lm(x~matGX[,ll]))$coef[2,1])
        sxU <- c(sxU, summary(lm(x~matGX[,ll]))$coef[2,2])
        pxU <- c(pxU, summary(lm(x~matGX[,ll]))$coef[2,4])
        
        byU <- c(byU, summary(lm(y~matGX[,ll]))$coef[2,1])
        syU <- c(syU, summary(lm(y~matGX[,ll]))$coef[2,2])
        pyU <- c(pyU, summary(lm(y~matGX[,ll]))$coef[2,4])
      }
      
      # save for cd-ratio, cd-gls, cd-egger
      u.betaGX_X <- bxU
      u.seGX_X <- sxU
      
      u.betaGX_Y <- byU
      u.seGX_Y <- syU
      
      for(ll in 1:ncol(matGY)){
        bxU <- c(bxU, summary(lm(x~matGY[,ll]))$coef[2,1])
        sxU <- c(sxU, summary(lm(x~matGY[,ll]))$coef[2,2])

        byU <- c(byU, summary(lm(y~matGY[,ll]))$coef[2,1])
        syU <- c(syU, summary(lm(y~matGY[,ll]))$coef[2,2])
      }
      
      # save betas for residuals and GX, GY
      bxR <- NULL
      sxR <- NULL
      pxR <- NULL

      byR <- NULL
      syR <- NULL
      pyR <- NULL

      for(ll in 1:ncol(matGX)){
        bxR <- c(bxR, summary(lm(xr~matGX[,ll]))$coef[2,1])
        sxR <- c(sxR, summary(lm(xr~matGX[,ll]))$coef[2,2])
        pxR <- c(pxR, summary(lm(xr~matGX[,ll]))$coef[2,4])
        
        byR <- c(byR, summary(lm(yr~matGX[,ll]))$coef[2,1])
        syR <- c(syR, summary(lm(yr~matGX[,ll]))$coef[2,2])
        pyR <- c(pyR, summary(lm(yr~matGX[,ll]))$coef[2,4])
      }
      
      # save for cd-ratio, cd-gls, cd-egger
      r.betaGX_X <- bxR
      r.seGX_X <- sxR
      
      r.betaGX_Y <- byR
      r.seGX_Y <- syR
      
      for(ll in 1:ncol(matGY)){
        bxR <- c(bxR, summary(lm(xr~matGY[,ll]))$coef[2,1])
        sxR <- c(sxR, summary(lm(xr~matGY[,ll]))$coef[2,2])
        
        byR <- c(byR, summary(lm(yr~matGY[,ll]))$coef[2,1])
        syR <- c(syR, summary(lm(yr~matGY[,ll]))$coef[2,2])
      }
      
      # save MLR filtered and SLR filtered SNPs for MRCD methods
      a.betaGX_X.MLR <- a.betaGX_X[pxA < pFilterMRCD & pyA < pFilterMRCD]
      a.seGX_X.MLR <- a.seGX_X[pxA < pFilterMRCD & pyA < pFilterMRCD]
      a.betaGX_Y.MLR <- a.betaGX_Y[pxA < pFilterMRCD & pyA < pFilterMRCD]
      a.seGX_Y.MLR <- a.seGX_Y[pxA < pFilterMRCD & pyA < pFilterMRCD]
      
      a.betaGX_X.SLR <- a.betaGX_X[pxU < pFilterMRCD & pyU < pFilterMRCD]
      a.seGX_X.SLR <- a.seGX_X[pxU < pFilterMRCD & pyU < pFilterMRCD]
      a.betaGX_Y.SLR <- a.betaGX_Y[pxU < pFilterMRCD & pyU < pFilterMRCD]
      a.seGX_Y.SLR <- a.seGX_Y[pxU < pFilterMRCD & pyU < pFilterMRCD]
      
      u.betaGX_X.MLR <- u.betaGX_X[pxA < pFilterMRCD & pyA < pFilterMRCD]
      u.seGX_X.MLR <- u.seGX_X[pxA < pFilterMRCD & pyA < pFilterMRCD]
      u.betaGX_Y.MLR <- u.betaGX_Y[pxA < pFilterMRCD & pyA < pFilterMRCD]
      u.seGX_Y.MLR <- u.seGX_Y[pxA < pFilterMRCD & pyA < pFilterMRCD]
      
      u.betaGX_X.SLR <- u.betaGX_X[pxU < pFilterMRCD & pyU < pFilterMRCD]
      u.seGX_X.SLR <- u.seGX_X[pxU < pFilterMRCD & pyU < pFilterMRCD]
      u.betaGX_Y.SLR <- u.betaGX_Y[pxU < pFilterMRCD & pyU < pFilterMRCD]
      u.seGX_Y.SLR <- u.seGX_Y[pxU < pFilterMRCD & pyU < pFilterMRCD]
      
      r.betaGX_X.MLR <- r.betaGX_X[pxA < pFilterMRCD & pyA < pFilterMRCD]
      r.seGX_X.MLR <- r.seGX_X[pxA < pFilterMRCD & pyA < pFilterMRCD]
      r.betaGX_Y.MLR <- r.betaGX_Y[pxA < pFilterMRCD & pyA < pFilterMRCD]
      r.seGX_Y.MLR <- r.seGX_Y[pxA < pFilterMRCD & pyA < pFilterMRCD]
      
      r.betaGX_X.SLR <- r.betaGX_X[pxU < pFilterMRCD & pyU < pFilterMRCD]
      r.seGX_X.SLR <- r.seGX_X[pxU < pFilterMRCD & pyU < pFilterMRCD]
      r.betaGX_Y.SLR <- r.betaGX_Y[pxU < pFilterMRCD & pyU < pFilterMRCD]
      r.seGX_Y.SLR <- r.seGX_Y[pxU < pFilterMRCD & pyU < pFilterMRCD]
      

      if(length(u.betaGX_X.MLR) <= 1 & "CDRatio.mlrF" %in% colnames.mat1){matR["allFiltered.MLR", "CDRatio.mlrF"] <- matR["allFiltered.MLR", "CDRatio.mlrF"] + 1}
      if(length(u.betaGX_X.MLR) <= 1 & "CDEgger.mlrF" %in% colnames.mat1){matR["allFiltered.MLR", "CDEgger.mlrF"] <- matR["allFiltered.MLR", "CDEgger.mlrF"] + 1}
      if(length(u.betaGX_X.MLR) <= 1 & "CDGLS.mlrF" %in% colnames.mat1){matR["allFiltered.MLR", "CDGLS.mlrF"] <- matR["allFiltered.MLR", "CDGLS.mlrF"] + 1}
 
     if(length(u.betaGX_X.SLR) <= 1 & "CDRatio.slrF" %in% colnames.mat1){matR["allFiltered.SLR", "CDRatio.slrF"] <- matR["allFiltered.SLR", "CDRatio.slrF"] + 1}
     if(length(u.betaGX_X.SLR) <= 1 & "CDEgger.slrF" %in% colnames.mat1){matR["allFiltered.SLR", "CDEgger.slrF"] <- matR["allFiltered.SLR", "CDEgger.slrF"] + 1}
     if(length(u.betaGX_X.SLR) <= 1 & "CDGLS.slrF" %in% colnames.mat1){matR["allFiltered.SLR", "CDGLS.slrF"] <- matR["allFiltered.SLR", "CDGLS.slrF"] + 1}
     
          
      if("CDRatio" %in% colnames.mat1 | "CDEgger" %in% colnames.mat1 | "CDGLS" %in% colnames.mat1){
        ######################################################################################################
        # Unadjusted -- CD-Ratio, CD-Egger, CD-GLS (1-Sample, multiple SNPs):
        ######################################################################################################
        uR <- NA
        uG <- NA
        uE <- NA
        
        u.tmp.mat <- matrix(NA,ncol=12,nrow=length(u.betaGX_X))
        u.tmp.df <- as.data.frame(u.tmp.mat)
        colnames(u.tmp.df) <- c("chr", "pos", "rsid", "A1", "A2", "beta_T1",
                                "se_T1", "N_T1", "beta_T2", "se_T2", "N_T2", "loci")
        u.tmp.df$rsid<-paste0("SNP", c(1:length(u.betaGX_X)))
        u.tmp.df$beta_T1 <- u.betaGX_X
        u.tmp.df$se_T1 <- u.seGX_X
        u.tmp.df$N_T1 <- n
        u.tmp.df$beta_T2 <- u.betaGX_Y
        u.tmp.df$se_T2 <- u.seGX_Y
        u.tmp.df$N_T2 <- n
        
        u.pruned1 <- list(sig_part=u.tmp.df)
        u.cd3 <- CD_3_methods_Independent(u.pruned1$sig_part)
        
        u.ratio.YX <- u.cd3$CD_Ratio_result$T1toT2
        u.ratio.XY <- u.cd3$CD_Ratio_result$T2toT1
        u.egger.YX<- u.cd3$CD_Egger_result$T1toT2
        u.egger.XY <- u.cd3$CD_Egger_result$T2toT1
        u.gls.YX <- u.cd3$CD_GLS_result$T1toT2
        u.gls.XY <- u.cd3$CD_GLS_result$T2toT1
        
        # Confidence intervals for K - needed for use in decision rules
        # CD-ratio CIs
        u.lowerCIyx.cdRatio <- u.ratio.YX["K"] - qnorm(1-sig.level/2)*u.ratio.YX["se(K)"]
        u.upperCIyx.cdRatio <- u.ratio.YX["K"] + qnorm(1-sig.level/2)*u.ratio.YX["se(K)"]
        u.lowerCIxy.cdRatio <- u.ratio.XY["K"] - qnorm(1-sig.level/2)*u.ratio.XY["se(K)"]
        u.upperCIxy.cdRatio <- u.ratio.XY["K"] + qnorm(1-sig.level/2)*u.ratio.XY["se(K)"]
        
        # CD-Egger CIs
        u.lowerCIyx.cdEgger <- u.egger.YX["K"] - qnorm(1-sig.level/2)*u.egger.YX["se(K)"]
        u.upperCIyx.cdEgger <- u.egger.YX["K"] + qnorm(1-sig.level/2)*u.egger.YX["se(K)"]
        u.lowerCIxy.cdEgger <- u.egger.XY["K"] - qnorm(1-sig.level/2)*u.egger.XY["se(K)"]
        u.upperCIxy.cdEgger <- u.egger.XY["K"] + qnorm(1-sig.level/2)*u.egger.XY["se(K)"]
        
        # CD-GLS CIs
        u.lowerCIyx.cdGLS <- u.gls.YX["K"] - qnorm(1-sig.level/2)*u.gls.YX["se(K)"]
        u.upperCIyx.cdGLS <- u.gls.YX["K"] + qnorm(1-sig.level/2)*u.gls.YX["se(K)"]
        u.lowerCIxy.cdGLS <- u.gls.XY["K"] - qnorm(1-sig.level/2)*u.gls.XY["se(K)"]
        u.upperCIxy.cdGLS <- u.gls.XY["K"] + qnorm(1-sig.level/2)*u.gls.XY["se(K)"]
        
        # case 1: X->Y if CIxy not inside [-1,1] and CIyx inside [-1,0) OR (0,1]
        # case 2: X<-Y if CIyx not inside [-1,1] and CIxy inside [-1,0) OR (0,1]
        # case 3: neither otherwise
        # CD-Ratio decisions
        if(u.upperCIxy.cdRatio < u.lowerCIxy.cdRatio){
          stop("LowerCIxy > UpperCIxy")}
        if(u.upperCIyx.cdRatio < u.lowerCIyx.cdRatio){
          stop("LowerCIyx > UpperCIyx")}
        
        if(u.upperCIxy.cdEgger < u.lowerCIxy.cdEgger){
          stop("LowerCIxyEgger > UpperCIxyEgger")}
        if(u.upperCIyx.cdEgger < u.lowerCIyx.cdEgger){
          stop("LowerCIyxEgger > UpperCIyxEgger")}
        
        if(u.upperCIxy.cdGLS < u.lowerCIxy.cdGLS){
          stop("LowerCIxy > UpperCIxy")}
        if(u.upperCIyx.cdGLS < u.lowerCIyx.cdGLS){
          stop("LowerCIyx > UpperCIyx")}
        
        # CD-Ratio decisions
        if("CDRatio" %in% colnames.mat1){
          if((u.upperCIxy.cdRatio < (-1) | u.lowerCIxy.cdRatio>1) &
             ((u.lowerCIyx.cdRatio>=(-1) & u.upperCIyx.cdRatio<0) | (u.lowerCIyx.cdRatio>(0) & u.upperCIyx.cdRatio<=1)) ){
            uR <- 1
          }else if((u.upperCIyx.cdRatio < (-1) | u.lowerCIyx.cdRatio>1) &
                   ((u.lowerCIxy.cdRatio>=(-1) & u.upperCIxy.cdRatio< 0) | (u.lowerCIxy.cdRatio>(0) & u.upperCIxy.cdRatio<= 1)) ){
            uR <- 2
          }else{
            uR <- 3
          }
        }
        
        # CD-Egger decisions
        if("CDEgger" %in% colnames.mat1){
          if((u.upperCIxy.cdEgger < (-1)  | u.lowerCIxy.cdEgger>1) &
             ((u.lowerCIyx.cdEgger>=(-1) & u.upperCIyx.cdEgger<0) | (u.lowerCIyx.cdEgger>(0) & u.upperCIyx.cdEgger<=1))){
            uE <- 1
          }else if((u.upperCIyx.cdEgger < (-1) | u.lowerCIyx.cdEgger>1) &
                   ((u.lowerCIxy.cdEgger>=(-1)  & u.upperCIxy.cdEgger< 0) | (u.lowerCIxy.cdEgger>(0)  & u.upperCIxy.cdEgger<=1))){
            uE <- 2
          }else{
            uE <- 3
          }
        }
        
        # CD-GLS decisions
        if("CDGLS" %in% colnames.mat1){
          if((u.upperCIxy.cdGLS < (-1) | u.lowerCIxy.cdGLS>1) &
             ( (u.lowerCIyx.cdGLS>= (-1) & u.upperCIyx.cdGLS<0) | (u.lowerCIyx.cdGLS> (0) & u.upperCIyx.cdGLS<=1)) ){
            uG <- 1
          }else if((u.upperCIyx.cdGLS < (-1) | u.lowerCIyx.cdGLS>1) &
                   ((u.lowerCIxy.cdGLS>= (-1) & u.upperCIxy.cdGLS<0) | (u.lowerCIxy.cdGLS> (0) & u.upperCIxy.cdGLS<=1)) ){
            uG <- 2
          }else{
            uG <- 3
          }
        }
        
        ######################################################################################################
        # Adjusted -- CD-Ratio, CD-Egger, CD-GLS (1-Sample, multiple SNPs):
        ######################################################################################################
        a.tmp.mat <- matrix(NA,ncol=12,nrow=length(a.betaGX_X))
        a.tmp.df <- as.data.frame(a.tmp.mat)
        colnames(a.tmp.df) <- c("chr", "pos", "rsid", "A1", "A2", "beta_T1",
                                "se_T1", "N_T1", "beta_T2", "se_T2", "N_T2", "loci")
        a.tmp.df$rsid<-paste0("SNP", c(1:length(a.betaGX_X)))
        a.tmp.df$beta_T1 <- a.betaGX_X
        a.tmp.df$se_T1 <- a.seGX_X
        a.tmp.df$N_T1 <- n
        a.tmp.df$beta_T2 <- a.betaGX_Y
        a.tmp.df$se_T2 <- a.seGX_Y
        a.tmp.df$N_T2 <- n
        a.pruned1 <- list(sig_part=a.tmp.df)
        
        a.cd3 <- CD_3_methods_Independent(a.pruned1$sig_part)
        
        # CD-RATIO
        if("CDRatio" %in% colnames.mat1){
          aR <- NA
          a.ratio.YX <- a.cd3$CD_Ratio_result$T1toT2
          a.ratio.XY <- a.cd3$CD_Ratio_result$T2toT1
          
          # Confidence intervals for K - needed for use in decision rules
          # CD-ratio CIs
          a.lowerCIyx.cdRatio <- a.ratio.YX["K"] - qnorm(1-sig.level/2)*a.ratio.YX["se(K)"]
          a.upperCIyx.cdRatio <- a.ratio.YX["K"] + qnorm(1-sig.level/2)*a.ratio.YX["se(K)"]
          a.lowerCIxy.cdRatio <- a.ratio.XY["K"] - qnorm(1-sig.level/2)*a.ratio.XY["se(K)"]
          a.upperCIxy.cdRatio <- a.ratio.XY["K"] + qnorm(1-sig.level/2)*a.ratio.XY["se(K)"]
          
          # case 1: X->Y if CIxy not inside [-1,1] and CIyx inside [-1,0) OR (0,1]
          # case 2: X<-Y if CIyx not inside [-1,1] and CIxy inside [-1,0) OR (0,1]
          # case 3: neither otherwise
          # CD-Ratio decisions
          if(a.upperCIxy.cdRatio < a.lowerCIxy.cdRatio){
            stop("LowerCIxy > UpperCIxy")}
          if(a.upperCIyx.cdRatio < a.lowerCIyx.cdRatio){
            stop("LowerCIyx > UpperCIyx")}
          
          # CD-Ratio decisions
          if((a.upperCIxy.cdRatio < (-1) | a.lowerCIxy.cdRatio>1) &
             ((a.lowerCIyx.cdRatio>=(-1) & a.upperCIyx.cdRatio<0) | (a.lowerCIyx.cdRatio>(0) & a.upperCIyx.cdRatio<=1)) ){
            aR <- 1
          }else if((a.upperCIyx.cdRatio < (-1) | a.lowerCIyx.cdRatio>1) &
                   ((a.lowerCIxy.cdRatio>=(-1) & a.upperCIxy.cdRatio< 0) | (a.lowerCIxy.cdRatio>(0) & a.upperCIxy.cdRatio<= 1)) ){
            aR <- 2
          }else{
            aR <- 3
          }
        }
        
        
        # START CD-GLS
        if("CDGLS" %in% colnames.mat1){
          aG <- NA
          
          a.gls.YX <- a.cd3$CD_GLS_result$T1toT2
          a.gls.XY <- a.cd3$CD_GLS_result$T2toT1
          
          # CD-GLS CIs - Confidence intervals for K - needed for use in decision rules
          a.lowerCIyx.cdGLS <- a.gls.YX["K"] - qnorm(1-sig.level/2)*a.gls.YX["se(K)"]
          a.upperCIyx.cdGLS <- a.gls.YX["K"] + qnorm(1-sig.level/2)*a.gls.YX["se(K)"]
          a.lowerCIxy.cdGLS <- a.gls.XY["K"] - qnorm(1-sig.level/2)*a.gls.XY["se(K)"]
          a.upperCIxy.cdGLS <- a.gls.XY["K"] + qnorm(1-sig.level/2)*a.gls.XY["se(K)"]
          
          # case 1: X->Y if CIxy not inside [-1,1] and CIyx inside [-1,0) OR (0,1]
          # case 2: X<-Y if CIyx not inside [-1,1] and CIxy inside [-1,0) OR (0,1]
          # case 3: neither otherwise
          if(a.upperCIxy.cdGLS < a.lowerCIxy.cdGLS){
            stop("LowerCIxy > UpperCIxy")}
          if(a.upperCIyx.cdGLS < a.lowerCIyx.cdGLS){
            stop("LowerCIyx > UpperCIyx")}
          
          # CD-GLS decisions
          if((a.upperCIxy.cdGLS < (-1) | a.lowerCIxy.cdGLS>1) &
             ((a.lowerCIyx.cdGLS>= (-1) & a.upperCIyx.cdGLS<0) | (a.lowerCIyx.cdGLS> (0) & a.upperCIyx.cdGLS<=1)) ){
            aG <- 1
          }else if((a.upperCIyx.cdGLS < (-1) | a.lowerCIyx.cdGLS>1) &
                   ((a.lowerCIxy.cdGLS>= (-1) & a.upperCIxy.cdGLS<0) | (a.lowerCIxy.cdGLS> (0) & a.upperCIxy.cdGLS<=1)) ){
            aG <- 2
          }else{
            aG <- 3
          }
        }# end CD-GLS
        
        # START CDEGGER
        if("CDEgger" %in% colnames.mat1){
          aE <- NA
          a.egger.YX<- a.cd3$CD_Egger_result$T1toT2
          a.egger.XY <- a.cd3$CD_Egger_result$T2toT1
          
          # CD-Egger CIs - Confidence intervals for K - needed for use in decision rules
          a.lowerCIyx.cdEgger <- a.egger.YX["K"] - qnorm(1-sig.level/2)*a.egger.YX["se(K)"]
          a.upperCIyx.cdEgger <- a.egger.YX["K"] + qnorm(1-sig.level/2)*a.egger.YX["se(K)"]
          a.lowerCIxy.cdEgger <- a.egger.XY["K"] - qnorm(1-sig.level/2)*a.egger.XY["se(K)"]
          a.upperCIxy.cdEgger <- a.egger.XY["K"] + qnorm(1-sig.level/2)*a.egger.XY["se(K)"]
          
          # case 1: X->Y if CIxy not inside [-1,1] and CIyx inside [-1,0) OR (0,1]
          # case 2: X<-Y if CIyx not inside [-1,1] and CIxy inside [-1,0) OR (0,1]
          # case 3: neither otherwise
          if(a.upperCIxy.cdEgger < a.lowerCIxy.cdEgger){
            stop("LowerCIxyEgger > UpperCIxyEgger")}
          if(a.upperCIyx.cdEgger < a.lowerCIyx.cdEgger){
            stop("LowerCIyxEgger > UpperCIyxEgger")}
          
          # CD-Egger decisions
          if((a.upperCIxy.cdEgger < (-1)  | a.lowerCIxy.cdEgger>1) &
             ((a.lowerCIyx.cdEgger>=(-1) & a.upperCIyx.cdEgger<0) | (a.lowerCIyx.cdEgger>(0) & a.upperCIyx.cdEgger<=1))){
            aE <- 1
          }else if((a.upperCIyx.cdEgger < (-1) | a.lowerCIyx.cdEgger>1) &
                   ((a.lowerCIxy.cdEgger>=(-1)  & a.upperCIxy.cdEgger< 0) | (a.lowerCIxy.cdEgger>(0)  & a.upperCIxy.cdEgger<=1))){
            aE <- 2
          }else{
            aE <- 3
          }
          
        }
        
        ######################################################################################################
        # Residual -- CD-Ratio, CD-Egger, CD-GLS (1-Sample, multiple SNPs):
        ######################################################################################################
        rR <- NA
        rG <- NA
        rE <- NA
        
        r.tmp.mat <- matrix(NA,ncol=12,nrow=length(r.betaGX_X))
        r.tmp.df <- as.data.frame(r.tmp.mat)
        colnames(r.tmp.df) <- c("chr", "pos", "rsid", "A1", "A2", "beta_T1",
                                "se_T1", "N_T1", "beta_T2", "se_T2", "N_T2", "loci")
        r.tmp.df$rsid<-paste0("SNP", c(1:length(r.betaGX_X)))
        r.tmp.df$beta_T1 <- r.betaGX_X
        r.tmp.df$se_T1 <- r.seGX_X
        r.tmp.df$N_T1 <- n
        r.tmp.df$beta_T2 <- r.betaGX_Y
        r.tmp.df$se_T2 <- r.seGX_Y
        r.tmp.df$N_T2 <- n
        
        r.pruned1 <- list(sig_part=r.tmp.df)
        r.cd3 <- CD_3_methods_Independent(r.pruned1$sig_part)
        
        r.ratio.YX <- r.cd3$CD_Ratio_result$T1toT2
        r.ratio.XY <- r.cd3$CD_Ratio_result$T2toT1
        r.egger.YX<- r.cd3$CD_Egger_result$T1toT2
        r.egger.XY <- r.cd3$CD_Egger_result$T2toT1
        r.gls.YX <- r.cd3$CD_GLS_result$T1toT2
        r.gls.XY <- r.cd3$CD_GLS_result$T2toT1
        
        # Confidence intervals for K - needed for use in decision rules
        # CD-ratio CIs
        r.lowerCIyx.cdRatio <- r.ratio.YX["K"] - qnorm(1-sig.level/2)*r.ratio.YX["se(K)"]
        r.upperCIyx.cdRatio <- r.ratio.YX["K"] + qnorm(1-sig.level/2)*r.ratio.YX["se(K)"]
        r.lowerCIxy.cdRatio <- r.ratio.XY["K"] - qnorm(1-sig.level/2)*r.ratio.XY["se(K)"]
        r.upperCIxy.cdRatio <- r.ratio.XY["K"] + qnorm(1-sig.level/2)*r.ratio.XY["se(K)"]
        
        # CD-Egger CIs
        r.lowerCIyx.cdEgger <- r.egger.YX["K"] - qnorm(1-sig.level/2)*r.egger.YX["se(K)"]
        r.upperCIyx.cdEgger <- r.egger.YX["K"] + qnorm(1-sig.level/2)*r.egger.YX["se(K)"]
        r.lowerCIxy.cdEgger <- r.egger.XY["K"] - qnorm(1-sig.level/2)*r.egger.XY["se(K)"]
        r.upperCIxy.cdEgger <- r.egger.XY["K"] + qnorm(1-sig.level/2)*r.egger.XY["se(K)"]
        
        # CD-GLS CIs
        r.lowerCIyx.cdGLS <- r.gls.YX["K"] - qnorm(1-sig.level/2)*r.gls.YX["se(K)"]
        r.upperCIyx.cdGLS <- r.gls.YX["K"] + qnorm(1-sig.level/2)*r.gls.YX["se(K)"]
        r.lowerCIxy.cdGLS <- r.gls.XY["K"] - qnorm(1-sig.level/2)*r.gls.XY["se(K)"]
        r.upperCIxy.cdGLS <- r.gls.XY["K"] + qnorm(1-sig.level/2)*r.gls.XY["se(K)"]
        
        
        # case 1: X->Y if CIxy not inside [-1,1] and CIyx inside [-1,0) OR (0,1]
        # case 2: X<-Y if CIyx not inside [-1,1] and CIxy inside [-1,0) OR (0,1]
        # case 3: neither otherwise
        # CD-Ratio decisions
        
        if(r.upperCIxy.cdRatio < r.lowerCIxy.cdRatio){
          stop("LowerCIxy > UpperCIxy")}
        if(r.upperCIyx.cdRatio < r.lowerCIyx.cdRatio){
          stop("LowerCIyx > UpperCIyx")}
        
        if(r.upperCIxy.cdEgger < r.lowerCIxy.cdEgger){
          stop("LowerCIxyEgger > UpperCIxyEgger")}
        if(r.upperCIyx.cdEgger < r.lowerCIyx.cdEgger){
          stop("LowerCIyxEgger > UpperCIyxEgger")}
        
        if(r.upperCIxy.cdGLS < r.lowerCIxy.cdGLS){
          stop("LowerCIxy > UpperCIxy")}
        if(r.upperCIyx.cdGLS < r.lowerCIyx.cdGLS){
          stop("LowerCIyx > UpperCIyx")}
        
        
        # CD-Ratio decisions
        if("CDRatio" %in% colnames.mat1){
          if((r.upperCIxy.cdRatio < (-1) | r.lowerCIxy.cdRatio>1) &
             ((r.lowerCIyx.cdRatio>=(-1) & r.upperCIyx.cdRatio<0) | (r.lowerCIyx.cdRatio>(0) & r.upperCIyx.cdRatio<=1)) ){
            rR <- 1
          }else if((r.upperCIyx.cdRatio < (-1) | r.lowerCIyx.cdRatio>1) &
                   ((r.lowerCIxy.cdRatio>=(-1) & r.upperCIxy.cdRatio< 0) | (r.lowerCIxy.cdRatio>(0) & r.upperCIxy.cdRatio<= 1)) ){
            rR <- 2
          }else{
            rR <- 3
          }
        }
        
        # CD-Egger decisions
        if("CDEgger" %in% colnames.mat1){
          if((r.upperCIxy.cdEgger < (-1)  | r.lowerCIxy.cdEgger>1) &
             ((r.lowerCIyx.cdEgger>=(-1) & r.upperCIyx.cdEgger<0) | (r.lowerCIyx.cdEgger>(0) & r.upperCIyx.cdEgger<=1))){
            rE <- 1
          }else if((r.upperCIyx.cdEgger < (-1) | r.lowerCIyx.cdEgger>1) &
                   ((r.lowerCIxy.cdEgger>=(-1)  & r.upperCIxy.cdEgger< 0) | (r.lowerCIxy.cdEgger>(0)  & r.upperCIxy.cdEgger<=1))){
            rE <- 2
          }else{
            rE <- 3
          }
        }
        
        # CD-GLS decisions
        if("CDGLS" %in% colnames.mat1){
          if((r.upperCIxy.cdGLS < (-1) | r.lowerCIxy.cdGLS>1) &
             ((r.lowerCIyx.cdGLS>= (-1) & r.upperCIyx.cdGLS<0) | (r.lowerCIyx.cdGLS> (0) & r.upperCIyx.cdGLS<=1)) ){
            rG <- 1
          }else if((r.upperCIyx.cdGLS < (-1) | r.lowerCIyx.cdGLS>1) &
                   ((r.lowerCIxy.cdGLS>= (-1) & r.upperCIxy.cdGLS<0) | (r.lowerCIxy.cdGLS> (0) & r.upperCIxy.cdGLS<=1)) ){
            rG <- 2
          }else{
            rG <- 3
          }
        }
        
        
        # SAVE CD-RATIO RESULTS
        if("CDRatio" %in% colnames.mat1){
          if(uR==1 & aR==1 & rR==1){matR["slr1_res1_mlr1", "CDRatio"] <- matR["slr1_res1_mlr1", "CDRatio"]+1}
          if(uR==1 & aR==1 & rR==2){matR["slr1_res2_mlr1", "CDRatio"] <- matR["slr1_res2_mlr1", "CDRatio"]+1}
          if(uR==1 & aR==1 & rR==3){matR["slr1_res3_mlr1", "CDRatio"] <- matR["slr1_res3_mlr1", "CDRatio"]+1}
          
          if(uR==1 & aR==2 & rR==1){matR["slr1_res1_mlr2", "CDRatio"] <- matR["slr1_res1_mlr2", "CDRatio"]+1}
          if(uR==1 & aR==2 & rR==2){matR["slr1_res2_mlr2", "CDRatio"] <- matR["slr1_res2_mlr2", "CDRatio"]+1}
          if(uR==1 & aR==2 & rR==3){matR["slr1_res3_mlr2", "CDRatio"] <- matR["slr1_res3_mlr2", "CDRatio"]+1}
          
          if(uR==1 & aR==3 & rR==1){matR["slr1_res1_mlr3", "CDRatio"] <- matR["slr1_res1_mlr3", "CDRatio"]+1}
          if(uR==1 & aR==3 & rR==2){matR["slr1_res2_mlr3", "CDRatio"] <- matR["slr1_res2_mlr3", "CDRatio"]+1}
          if(uR==1 & aR==3 & rR==3){matR["slr1_res3_mlr3", "CDRatio"] <- matR["slr1_res3_mlr3", "CDRatio"]+1}
          
          if(uR==2 & aR==1 & rR==1){matR["slr2_res1_mlr1", "CDRatio"] <- matR["slr2_res1_mlr1", "CDRatio"]+1}
          if(uR==2 & aR==1 & rR==2){matR["slr2_res2_mlr1", "CDRatio"] <- matR["slr2_res2_mlr1", "CDRatio"]+1}
          if(uR==2 & aR==1 & rR==3){matR["slr2_res3_mlr1", "CDRatio"] <- matR["slr2_res3_mlr1", "CDRatio"]+1}
          
          if(uR==2 & aR==2 & rR==1){matR["slr2_res1_mlr2", "CDRatio"] <- matR["slr2_res1_mlr2", "CDRatio"]+1}
          if(uR==2 & aR==2 & rR==2){matR["slr2_res2_mlr2", "CDRatio"] <- matR["slr2_res2_mlr2", "CDRatio"]+1}
          if(uR==2 & aR==2 & rR==3){matR["slr2_res3_mlr2", "CDRatio"] <- matR["slr2_res3_mlr2", "CDRatio"]+1}
          
          if(uR==2 & aR==3 & rR==1){matR["slr2_res1_mlr3", "CDRatio"] <- matR["slr2_res1_mlr3", "CDRatio"]+1}
          if(uR==2 & aR==3 & rR==2){matR["slr2_res2_mlr3", "CDRatio"] <- matR["slr2_res2_mlr3", "CDRatio"]+1}
          if(uR==2 & aR==3 & rR==3){matR["slr2_res3_mlr3", "CDRatio"] <- matR["slr2_res3_mlr3", "CDRatio"]+1}
          
          if(uR==3 & aR==1 & rR==1){matR["slr3_res1_mlr1", "CDRatio"] <- matR["slr3_res1_mlr1", "CDRatio"]+1}
          if(uR==3 & aR==1 & rR==2){matR["slr3_res2_mlr1", "CDRatio"] <- matR["slr3_res2_mlr1", "CDRatio"]+1}
          if(uR==3 & aR==1 & rR==3){matR["slr3_res3_mlr1", "CDRatio"] <- matR["slr3_res3_mlr1", "CDRatio"]+1}
          
          if(uR==3 & aR==2 & rR==1){matR["slr3_res1_mlr2", "CDRatio"] <- matR["slr3_res1_mlr2", "CDRatio"]+1}
          if(uR==3 & aR==2 & rR==2){matR["slr3_res2_mlr2", "CDRatio"] <- matR["slr3_res2_mlr2", "CDRatio"]+1}
          if(uR==3 & aR==2 & rR==3){matR["slr3_res3_mlr2", "CDRatio"] <- matR["slr3_res3_mlr2", "CDRatio"]+1}
          
          if(uR==3 & aR==3 & rR==1){matR["slr3_res1_mlr3", "CDRatio"] <- matR["slr3_res1_mlr3", "CDRatio"]+1}
          if(uR==3 & aR==3 & rR==2){matR["slr3_res2_mlr3", "CDRatio"] <- matR["slr3_res2_mlr3", "CDRatio"]+1}
          if(uR==3 & aR==3 & rR==3){matR["slr3_res3_mlr3", "CDRatio"] <- matR["slr3_res3_mlr3", "CDRatio"]+1}
        }
        
        # Save CD-GLS results
        if("CDGLS" %in% colnames.mat1){
          
          if(uG==1 & aG==1 & rG==1){matR["slr1_res1_mlr1", "CDGLS"] <- matR["slr1_res1_mlr1", "CDGLS"]+1}
          if(uG==1 & aG==1 & rG==2){matR["slr1_res2_mlr1", "CDGLS"] <- matR["slr1_res2_mlr1", "CDGLS"]+1}
          if(uG==1 & aG==1 & rG==3){matR["slr1_res3_mlr1", "CDGLS"] <- matR["slr1_res3_mlr1", "CDGLS"]+1}
          
          if(uG==1 & aG==2 & rG==1){matR["slr1_res1_mlr2", "CDGLS"] <- matR["slr1_res1_mlr2", "CDGLS"]+1}
          if(uG==1 & aG==2 & rG==2){matR["slr1_res2_mlr2", "CDGLS"] <- matR["slr1_res2_mlr2", "CDGLS"]+1}
          if(uG==1 & aG==2 & rG==3){matR["slr1_res3_mlr2", "CDGLS"] <- matR["slr1_res3_mlr2", "CDGLS"]+1}
          
          if(uG==1 & aG==3 & rG==1){matR["slr1_res1_mlr3", "CDGLS"] <- matR["slr1_res1_mlr3", "CDGLS"]+1}
          if(uG==1 & aG==3 & rG==2){matR["slr1_res2_mlr3", "CDGLS"] <- matR["slr1_res2_mlr3", "CDGLS"]+1}
          if(uG==1 & aG==3 & rG==3){matR["slr1_res3_mlr3", "CDGLS"] <- matR["slr1_res3_mlr3", "CDGLS"]+1}
          
          if(uG==2 & aG==1 & rG==1){matR["slr2_res1_mlr1", "CDGLS"] <- matR["slr2_res1_mlr1", "CDGLS"]+1}
          if(uG==2 & aG==1 & rG==2){matR["slr2_res2_mlr1", "CDGLS"] <- matR["slr2_res2_mlr1", "CDGLS"]+1}
          if(uG==2 & aG==1 & rG==3){matR["slr2_res3_mlr1", "CDGLS"] <- matR["slr2_res3_mlr1", "CDGLS"]+1}
          
          if(uG==2 & aG==2 & rG==1){matR["slr2_res1_mlr2", "CDGLS"] <- matR["slr2_res1_mlr2", "CDGLS"]+1}
          if(uG==2 & aG==2 & rG==2){matR["slr2_res2_mlr2", "CDGLS"] <- matR["slr2_res2_mlr2", "CDGLS"]+1}
          if(uG==2 & aG==2 & rG==3){matR["slr2_res3_mlr2", "CDGLS"] <- matR["slr2_res3_mlr2", "CDGLS"]+1}
          
          if(uG==2 & aG==3 & rG==1){matR["slr2_res1_mlr3", "CDGLS"] <- matR["slr2_res1_mlr3", "CDGLS"]+1}
          if(uG==2 & aG==3 & rG==2){matR["slr2_res2_mlr3", "CDGLS"] <- matR["slr2_res2_mlr3", "CDGLS"]+1}
          if(uG==2 & aG==3 & rG==3){matR["slr2_res3_mlr3", "CDGLS"] <- matR["slr2_res3_mlr3", "CDGLS"]+1}
          
          if(uG==3 & aG==1 & rG==1){matR["slr3_res1_mlr1", "CDGLS"] <- matR["slr3_res1_mlr1", "CDGLS"]+1}
          if(uG==3 & aG==1 & rG==2){matR["slr3_res2_mlr1", "CDGLS"] <- matR["slr3_res2_mlr1", "CDGLS"]+1}
          if(uG==3 & aG==1 & rG==3){matR["slr3_res3_mlr1", "CDGLS"] <- matR["slr3_res3_mlr1", "CDGLS"]+1}
          
          if(uG==3 & aG==2 & rG==1){matR["slr3_res1_mlr2", "CDGLS"] <- matR["slr3_res1_mlr2", "CDGLS"]+1}
          if(uG==3 & aG==2 & rG==2){matR["slr3_res2_mlr2", "CDGLS"] <- matR["slr3_res2_mlr2", "CDGLS"]+1}
          if(uG==3 & aG==2 & rG==3){matR["slr3_res3_mlr2", "CDGLS"] <- matR["slr3_res3_mlr2", "CDGLS"]+1}
          
          if(uG==3 & aG==3 & rG==1){matR["slr3_res1_mlr3", "CDGLS"] <- matR["slr3_res1_mlr3", "CDGLS"]+1}
          if(uG==3 & aG==3 & rG==2){matR["slr3_res2_mlr3", "CDGLS"] <- matR["slr3_res2_mlr3", "CDGLS"]+1}
          if(uG==3 & aG==3 & rG==3){matR["slr3_res3_mlr3", "CDGLS"] <- matR["slr3_res3_mlr3", "CDGLS"]+1}
        }
        
        # save CD-Egger results
        if("CDEgger" %in% colnames.mat1){
          
          if(uE==1 & aE==1 & rE==1){matR["slr1_res1_mlr1", "CDEgger"] <- matR["slr1_res1_mlr1", "CDEgger"]+1}
          if(uE==1 & aE==1 & rE==2){matR["slr1_res2_mlr1", "CDEgger"] <- matR["slr1_res2_mlr1", "CDEgger"]+1}
          if(uE==1 & aE==1 & rE==3){matR["slr1_res3_mlr1", "CDEgger"] <- matR["slr1_res3_mlr1", "CDEgger"]+1}
          
          if(uE==1 & aE==2 & rE==1){matR["slr1_res1_mlr2", "CDEgger"] <- matR["slr1_res1_mlr2", "CDEgger"]+1}
          if(uE==1 & aE==2 & rE==2){matR["slr1_res2_mlr2", "CDEgger"] <- matR["slr1_res2_mlr2", "CDEgger"]+1}
          if(uE==1 & aE==2 & rE==3){matR["slr1_res3_mlr2", "CDEgger"] <- matR["slr1_res3_mlr2", "CDEgger"]+1}
          
          if(uE==1 & aE==3 & rE==1){matR["slr1_res1_mlr3", "CDEgger"] <- matR["slr1_res1_mlr3", "CDEgger"]+1}
          if(uE==1 & aE==3 & rE==2){matR["slr1_res2_mlr3", "CDEgger"] <- matR["slr1_res2_mlr3", "CDEgger"]+1}
          if(uE==1 & aE==3 & rE==3){matR["slr1_res3_mlr3", "CDEgger"] <- matR["slr1_res3_mlr3", "CDEgger"]+1}
          
          
          if(uE==2 & aE==1 & rE==1){matR["slr2_res1_mlr1", "CDEgger"] <- matR["slr2_res1_mlr1", "CDEgger"]+1}
          if(uE==2 & aE==1 & rE==2){matR["slr2_res2_mlr1", "CDEgger"] <- matR["slr2_res2_mlr1", "CDEgger"]+1}
          if(uE==2 & aE==1 & rE==3){matR["slr2_res3_mlr1", "CDEgger"] <- matR["slr2_res3_mlr1", "CDEgger"]+1}
          
          if(uE==2 & aE==2 & rE==1){matR["slr2_res1_mlr2", "CDEgger"] <- matR["slr2_res1_mlr2", "CDEgger"]+1}
          if(uE==2 & aE==2 & rE==2){matR["slr2_res2_mlr2", "CDEgger"] <- matR["slr2_res2_mlr2", "CDEgger"]+1}
          if(uE==2 & aE==2 & rE==3){matR["slr2_res3_mlr2", "CDEgger"] <- matR["slr2_res3_mlr2", "CDEgger"]+1}
          
          if(uE==2 & aE==3 & rE==1){matR["slr2_res1_mlr3", "CDEgger"] <- matR["slr2_res1_mlr3", "CDEgger"]+1}
          if(uE==2 & aE==3 & rE==2){matR["slr2_res2_mlr3", "CDEgger"] <- matR["slr2_res2_mlr3", "CDEgger"]+1}
          if(uE==2 & aE==3 & rE==3){matR["slr2_res3_mlr3", "CDEgger"] <- matR["slr2_res3_mlr3", "CDEgger"]+1}
          
          if(uE==3 & aE==1 & rE==1){matR["slr3_res1_mlr1", "CDEgger"] <- matR["slr3_res1_mlr1", "CDEgger"]+1}
          if(uE==3 & aE==1 & rE==2){matR["slr3_res2_mlr1", "CDEgger"] <- matR["slr3_res2_mlr1", "CDEgger"]+1}
          if(uE==3 & aE==1 & rE==3){matR["slr3_res3_mlr1", "CDEgger"] <- matR["slr3_res3_mlr1", "CDEgger"]+1}
          
          if(uE==3 & aE==2 & rE==1){matR["slr3_res1_mlr2", "CDEgger"] <- matR["slr3_res1_mlr2", "CDEgger"]+1}
          if(uE==3 & aE==2 & rE==2){matR["slr3_res2_mlr2", "CDEgger"] <- matR["slr3_res2_mlr2", "CDEgger"]+1}
          if(uE==3 & aE==2 & rE==3){matR["slr3_res3_mlr2", "CDEgger"] <- matR["slr3_res3_mlr2", "CDEgger"]+1}
          
          if(uE==3 & aE==3 & rE==1){matR["slr3_res1_mlr3", "CDEgger"] <- matR["slr3_res1_mlr3", "CDEgger"]+1}
          if(uE==3 & aE==3 & rE==2){matR["slr3_res2_mlr3", "CDEgger"] <- matR["slr3_res2_mlr3", "CDEgger"]+1}
          if(uE==3 & aE==3 & rE==3){matR["slr3_res3_mlr3", "CDEgger"] <- matR["slr3_res3_mlr3", "CDEgger"]+1}
        } # End CD-Egger
        
        
      } # END CD-RATIO / CD-EGGER / CD-GLS no filtering
      
      
      if(("CDRatio.slrF" %in% colnames.mat1 | "CDEgger.slrF" %in% colnames.mat1 | "CDGLS.slrF" %in% colnames.mat1) & (length(u.betaGX_X.SLR) > 1)){
        ######################################################################################################
        # Unadjusted -- CD-Ratio, CD-Egger, CD-GLS (1-Sample, multiple SNPs):
        ######################################################################################################
        uR <- NA
        uG <- NA
        uE <- NA
        
        u.tmp.mat <- matrix(NA,ncol=12,nrow=length(u.betaGX_X.SLR))
        u.tmp.df <- as.data.frame(u.tmp.mat)
        colnames(u.tmp.df) <- c("chr", "pos", "rsid", "A1", "A2", "beta_T1",
                                "se_T1", "N_T1", "beta_T2", "se_T2", "N_T2", "loci")
        u.tmp.df$rsid<-paste0("SNP", c(1:length(u.betaGX_X.SLR)))
        u.tmp.df$beta_T1 <- u.betaGX_X.SLR
        u.tmp.df$se_T1 <- u.seGX_X.SLR
        u.tmp.df$N_T1 <- n
        u.tmp.df$beta_T2 <- u.betaGX_Y.SLR
        u.tmp.df$se_T2 <- u.seGX_Y.SLR
        u.tmp.df$N_T2 <- n
        
        u.pruned1 <- list(sig_part=u.tmp.df)
        u.cd3 <- CD_3_methods_Independent(u.pruned1$sig_part)
        
        cdDecision.Fun <- function(cd3){
          # Confidence intervals for K - needed for use in decision rules
          # CD-ratio CIs
          
          ratio.YX <- cd3$CD_Ratio_result$T1toT2
          ratio.XY <- cd3$CD_Ratio_result$T2toT1
          egger.YX<- cd3$CD_Egger_result$T1toT2
          egger.XY <- cd3$CD_Egger_result$T2toT1
          gls.YX <- cd3$CD_GLS_result$T1toT2
          gls.XY <- cd3$CD_GLS_result$T2toT1
          
          # Confidence intervals for K - needed for use in decision rules
          # CD-ratio CIs
          lowerCIyx.cdRatio <- ratio.YX["K"] - qnorm(1-sig.level/2)*ratio.YX["se(K)"]
          upperCIyx.cdRatio <- ratio.YX["K"] + qnorm(1-sig.level/2)*ratio.YX["se(K)"]
          lowerCIxy.cdRatio <- ratio.XY["K"] - qnorm(1-sig.level/2)*ratio.XY["se(K)"]
          upperCIxy.cdRatio <- ratio.XY["K"] + qnorm(1-sig.level/2)*ratio.XY["se(K)"]
          
          # CD-Egger CIs
          lowerCIyx.cdEgger <- egger.YX["K"] - qnorm(1-sig.level/2)*egger.YX["se(K)"]
          upperCIyx.cdEgger <- egger.YX["K"] + qnorm(1-sig.level/2)*egger.YX["se(K)"]
          lowerCIxy.cdEgger <- egger.XY["K"] - qnorm(1-sig.level/2)*egger.XY["se(K)"]
          upperCIxy.cdEgger <- egger.XY["K"] + qnorm(1-sig.level/2)*egger.XY["se(K)"]
          
          # CD-GLS CIs
          lowerCIyx.cdGLS <- gls.YX["K"] - qnorm(1-sig.level/2)*gls.YX["se(K)"]
          upperCIyx.cdGLS <- gls.YX["K"] + qnorm(1-sig.level/2)*gls.YX["se(K)"]
          lowerCIxy.cdGLS <- gls.XY["K"] - qnorm(1-sig.level/2)*gls.XY["se(K)"]
          upperCIxy.cdGLS <- gls.XY["K"] + qnorm(1-sig.level/2)*gls.XY["se(K)"]
          
          # case 1: X->Y if CIxy not inside [-1,1] and CIyx inside [-1,0) OR (0,1]
          # case 2: X<-Y if CIyx not inside [-1,1] and CIxy inside [-1,0) OR (0,1]
          # case 3: neither otherwise
          # CD-Ratio decisions
          if(upperCIxy.cdRatio < lowerCIxy.cdRatio){
            stop("LowerCIxy > UpperCIxy")}
          if(upperCIyx.cdRatio < lowerCIyx.cdRatio){
            stop("LowerCIyx > UpperCIyx")}
          
          if(upperCIxy.cdEgger < lowerCIxy.cdEgger){
            stop("LowerCIxyEgger > UpperCIxyEgger")}
          if(upperCIyx.cdEgger < lowerCIyx.cdEgger){
            stop("LowerCIyxEgger > UpperCIyxEgger")}
          
          if(upperCIxy.cdGLS < lowerCIxy.cdGLS){
            stop("LowerCIxy > UpperCIxy")}
          if(upperCIyx.cdGLS < lowerCIyx.cdGLS){
            stop("LowerCIyx > UpperCIyx")}
         
          tmpR <- NA
          # CD-Ratio decisions
            if((upperCIxy.cdRatio < (-1) | lowerCIxy.cdRatio>1) &
               ((lowerCIyx.cdRatio>=(-1) & upperCIyx.cdRatio<0) | (lowerCIyx.cdRatio>(0) & upperCIyx.cdRatio<=1)) ){
              tmpR <- 1
            }else if((upperCIyx.cdRatio < (-1) | lowerCIyx.cdRatio>1) &
                     ((lowerCIxy.cdRatio>=(-1) & upperCIxy.cdRatio< 0) | (lowerCIxy.cdRatio>(0) & upperCIxy.cdRatio<= 1)) ){
              tmpR <- 2
            }else{
              tmpR <- 3
            }
          
          # CD-Egger decisions
          tmpE <- NA
            if((upperCIxy.cdEgger < (-1)  | lowerCIxy.cdEgger>1) &
               ((lowerCIyx.cdEgger>=(-1) & upperCIyx.cdEgger<0) | (lowerCIyx.cdEgger>(0) & upperCIyx.cdEgger<=1))){
              tmpE <- 1
            }else if((upperCIyx.cdEgger < (-1) | lowerCIyx.cdEgger>1) &
                     ((lowerCIxy.cdEgger>=(-1)  & upperCIxy.cdEgger< 0) | (lowerCIxy.cdEgger>(0)  & upperCIxy.cdEgger<=1))){
              tmpE <- 2
            }else{
              tmpE <- 3
            }
          
          # CD-GLS decisions
          tmpG <- NA
            if((upperCIxy.cdGLS < (-1) | lowerCIxy.cdGLS>1) &
               ( (lowerCIyx.cdGLS>= (-1) & upperCIyx.cdGLS<0) | (lowerCIyx.cdGLS> (0) & upperCIyx.cdGLS<=1)) ){
              tmpG <- 1
            }else if((upperCIyx.cdGLS < (-1) | lowerCIyx.cdGLS>1) &
                     ((lowerCIxy.cdGLS>= (-1) & upperCIxy.cdGLS<0) | (lowerCIxy.cdGLS> (0) & upperCIxy.cdGLS<=1)) ){
              tmpG <- 2
            }else{
              tmpG <- 3
            }
          
        return(list("ratio.Case" = tmpR, "egger.Case"=tmpE, "gls.Case" = tmpG))  
        }
        
        
        # unadjusted cases
        u.Out <- cdDecision.Fun(cd3=u.cd3)
        uR <- u.Out$ratio.Case
        uG <- u.Out$gls.Case
        uE <- u.Out$egger.Case
        
        ######################################################################################################
        # Adjusted -- CD-Ratio, CD-Egger, CD-GLS (1-Sample, multiple SNPs):
        ######################################################################################################
        a.tmp.mat <- matrix(NA,ncol=12,nrow=length(a.betaGX_X.SLR))
        a.tmp.df <- as.data.frame(a.tmp.mat)
        colnames(a.tmp.df) <- c("chr", "pos", "rsid", "A1", "A2", "beta_T1",
                                "se_T1", "N_T1", "beta_T2", "se_T2", "N_T2", "loci")
        a.tmp.df$rsid<-paste0("SNP", c(1:length(a.betaGX_X.SLR)))
        a.tmp.df$beta_T1 <- a.betaGX_X.SLR
        a.tmp.df$se_T1 <- a.seGX_X.SLR
        a.tmp.df$N_T1 <- n
        a.tmp.df$beta_T2 <- a.betaGX_Y.SLR
        a.tmp.df$se_T2 <- a.seGX_Y.SLR
        a.tmp.df$N_T2 <- n
        a.pruned1 <- list(sig_part=a.tmp.df)
        
        a.cd3 <- CD_3_methods_Independent(a.pruned1$sig_part)
        aR <- NA
        aG <- NA
        aE <- NA

        a.Out <- cdDecision.Fun(cd3=a.cd3)
        aR <- a.Out$ratio.Case
        aG <- a.Out$gls.Case
        aE <- a.Out$egger.Case        
        
        ######################################################################################################
        # Residual -- CD-Ratio, CD-Egger, CD-GLS (1-Sample, multiple SNPs):
        ######################################################################################################
        rR <- NA
        rG <- NA
        rE <- NA
        
        r.tmp.mat <- matrix(NA,ncol=12,nrow=length(r.betaGX_X.SLR))
        r.tmp.df <- as.data.frame(r.tmp.mat)
        colnames(r.tmp.df) <- c("chr", "pos", "rsid", "A1", "A2", "beta_T1",
                                "se_T1", "N_T1", "beta_T2", "se_T2", "N_T2", "loci")
        r.tmp.df$rsid<-paste0("SNP", c(1:length(r.betaGX_X.SLR)))
        r.tmp.df$beta_T1 <- r.betaGX_X.SLR
        r.tmp.df$se_T1 <- r.seGX_X.SLR
        r.tmp.df$N_T1 <- n
        r.tmp.df$beta_T2 <- r.betaGX_Y.SLR
        r.tmp.df$se_T2 <- r.seGX_Y.SLR
        r.tmp.df$N_T2 <- n
        
        r.pruned1 <- list(sig_part=r.tmp.df)
        r.cd3 <- CD_3_methods_Independent(r.pruned1$sig_part)
        
        r.Out <- cdDecision.Fun(cd3=r.cd3)
        rR <- r.Out$ratio.Case
        rG <- r.Out$gls.Case
        rE <- r.Out$egger.Case
        
        
        # SAVE CD-RATIO RESULTS
        if("CDRatio.slrF" %in% colnames.mat1){
          if(uR==1 & aR==1 & rR==1){matR["slr1_res1_mlr1", "CDRatio.slrF"] <- matR["slr1_res1_mlr1", "CDRatio.slrF"]+1}
          if(uR==1 & aR==1 & rR==2){matR["slr1_res2_mlr1", "CDRatio.slrF"] <- matR["slr1_res2_mlr1", "CDRatio.slrF"]+1}
          if(uR==1 & aR==1 & rR==3){matR["slr1_res3_mlr1", "CDRatio.slrF"] <- matR["slr1_res3_mlr1", "CDRatio.slrF"]+1}
          
          if(uR==1 & aR==2 & rR==1){matR["slr1_res1_mlr2", "CDRatio.slrF"] <- matR["slr1_res1_mlr2", "CDRatio.slrF"]+1}
          if(uR==1 & aR==2 & rR==2){matR["slr1_res2_mlr2", "CDRatio.slrF"] <- matR["slr1_res2_mlr2", "CDRatio.slrF"]+1}
          if(uR==1 & aR==2 & rR==3){matR["slr1_res3_mlr2", "CDRatio.slrF"] <- matR["slr1_res3_mlr2", "CDRatio.slrF"]+1}
          
          if(uR==1 & aR==3 & rR==1){matR["slr1_res1_mlr3", "CDRatio.slrF"] <- matR["slr1_res1_mlr3", "CDRatio.slrF"]+1}
          if(uR==1 & aR==3 & rR==2){matR["slr1_res2_mlr3", "CDRatio.slrF"] <- matR["slr1_res2_mlr3", "CDRatio.slrF"]+1}
          if(uR==1 & aR==3 & rR==3){matR["slr1_res3_mlr3", "CDRatio.slrF"] <- matR["slr1_res3_mlr3", "CDRatio.slrF"]+1}
          
          if(uR==2 & aR==1 & rR==1){matR["slr2_res1_mlr1", "CDRatio.slrF"] <- matR["slr2_res1_mlr1", "CDRatio.slrF"]+1}
          if(uR==2 & aR==1 & rR==2){matR["slr2_res2_mlr1", "CDRatio.slrF"] <- matR["slr2_res2_mlr1", "CDRatio.slrF"]+1}
          if(uR==2 & aR==1 & rR==3){matR["slr2_res3_mlr1", "CDRatio.slrF"] <- matR["slr2_res3_mlr1", "CDRatio.slrF"]+1}
          
          if(uR==2 & aR==2 & rR==1){matR["slr2_res1_mlr2", "CDRatio.slrF"] <- matR["slr2_res1_mlr2", "CDRatio.slrF"]+1}
          if(uR==2 & aR==2 & rR==2){matR["slr2_res2_mlr2", "CDRatio.slrF"] <- matR["slr2_res2_mlr2", "CDRatio.slrF"]+1}
          if(uR==2 & aR==2 & rR==3){matR["slr2_res3_mlr2", "CDRatio.slrF"] <- matR["slr2_res3_mlr2", "CDRatio.slrF"]+1}
          
          if(uR==2 & aR==3 & rR==1){matR["slr2_res1_mlr3", "CDRatio.slrF"] <- matR["slr2_res1_mlr3", "CDRatio.slrF"]+1}
          if(uR==2 & aR==3 & rR==2){matR["slr2_res2_mlr3", "CDRatio.slrF"] <- matR["slr2_res2_mlr3", "CDRatio.slrF"]+1}
          if(uR==2 & aR==3 & rR==3){matR["slr2_res3_mlr3", "CDRatio.slrF"] <- matR["slr2_res3_mlr3", "CDRatio.slrF"]+1}
          
          if(uR==3 & aR==1 & rR==1){matR["slr3_res1_mlr1", "CDRatio.slrF"] <- matR["slr3_res1_mlr1", "CDRatio.slrF"]+1}
          if(uR==3 & aR==1 & rR==2){matR["slr3_res2_mlr1", "CDRatio.slrF"] <- matR["slr3_res2_mlr1", "CDRatio.slrF"]+1}
          if(uR==3 & aR==1 & rR==3){matR["slr3_res3_mlr1", "CDRatio.slrF"] <- matR["slr3_res3_mlr1", "CDRatio.slrF"]+1}
          
          if(uR==3 & aR==2 & rR==1){matR["slr3_res1_mlr2", "CDRatio.slrF"] <- matR["slr3_res1_mlr2", "CDRatio.slrF"]+1}
          if(uR==3 & aR==2 & rR==2){matR["slr3_res2_mlr2", "CDRatio.slrF"] <- matR["slr3_res2_mlr2", "CDRatio.slrF"]+1}
          if(uR==3 & aR==2 & rR==3){matR["slr3_res3_mlr2", "CDRatio.slrF"] <- matR["slr3_res3_mlr2", "CDRatio.slrF"]+1}
          
          if(uR==3 & aR==3 & rR==1){matR["slr3_res1_mlr3", "CDRatio.slrF"] <- matR["slr3_res1_mlr3", "CDRatio.slrF"]+1}
          if(uR==3 & aR==3 & rR==2){matR["slr3_res2_mlr3", "CDRatio.slrF"] <- matR["slr3_res2_mlr3", "CDRatio.slrF"]+1}
          if(uR==3 & aR==3 & rR==3){matR["slr3_res3_mlr3", "CDRatio.slrF"] <- matR["slr3_res3_mlr3", "CDRatio.slrF"]+1}
        }
        
        # Save CD-GLS results
        if("CDGLS.slrF" %in% colnames.mat1){
          
          if(uG==1 & aG==1 & rG==1){matR["slr1_res1_mlr1", "CDGLS.slrF"] <- matR["slr1_res1_mlr1", "CDGLS.slrF"]+1}
          if(uG==1 & aG==1 & rG==2){matR["slr1_res2_mlr1", "CDGLS.slrF"] <- matR["slr1_res2_mlr1", "CDGLS.slrF"]+1}
          if(uG==1 & aG==1 & rG==3){matR["slr1_res3_mlr1", "CDGLS.slrF"] <- matR["slr1_res3_mlr1", "CDGLS.slrF"]+1}
          
          if(uG==1 & aG==2 & rG==1){matR["slr1_res1_mlr2", "CDGLS.slrF"] <- matR["slr1_res1_mlr2", "CDGLS.slrF"]+1}
          if(uG==1 & aG==2 & rG==2){matR["slr1_res2_mlr2", "CDGLS.slrF"] <- matR["slr1_res2_mlr2", "CDGLS.slrF"]+1}
          if(uG==1 & aG==2 & rG==3){matR["slr1_res3_mlr2", "CDGLS.slrF"] <- matR["slr1_res3_mlr2", "CDGLS.slrF"]+1}
          
          if(uG==1 & aG==3 & rG==1){matR["slr1_res1_mlr3", "CDGLS.slrF"] <- matR["slr1_res1_mlr3", "CDGLS.slrF"]+1}
          if(uG==1 & aG==3 & rG==2){matR["slr1_res2_mlr3", "CDGLS.slrF"] <- matR["slr1_res2_mlr3", "CDGLS.slrF"]+1}
          if(uG==1 & aG==3 & rG==3){matR["slr1_res3_mlr3", "CDGLS.slrF"] <- matR["slr1_res3_mlr3", "CDGLS.slrF"]+1}
          
          if(uG==2 & aG==1 & rG==1){matR["slr2_res1_mlr1", "CDGLS.slrF"] <- matR["slr2_res1_mlr1", "CDGLS.slrF"]+1}
          if(uG==2 & aG==1 & rG==2){matR["slr2_res2_mlr1", "CDGLS.slrF"] <- matR["slr2_res2_mlr1", "CDGLS.slrF"]+1}
          if(uG==2 & aG==1 & rG==3){matR["slr2_res3_mlr1", "CDGLS.slrF"] <- matR["slr2_res3_mlr1", "CDGLS.slrF"]+1}
          
          if(uG==2 & aG==2 & rG==1){matR["slr2_res1_mlr2", "CDGLS.slrF"] <- matR["slr2_res1_mlr2", "CDGLS.slrF"]+1}
          if(uG==2 & aG==2 & rG==2){matR["slr2_res2_mlr2", "CDGLS.slrF"] <- matR["slr2_res2_mlr2", "CDGLS.slrF"]+1}
          if(uG==2 & aG==2 & rG==3){matR["slr2_res3_mlr2", "CDGLS.slrF"] <- matR["slr2_res3_mlr2", "CDGLS.slrF"]+1}
          
          if(uG==2 & aG==3 & rG==1){matR["slr2_res1_mlr3", "CDGLS.slrF"] <- matR["slr2_res1_mlr3", "CDGLS.slrF"]+1}
          if(uG==2 & aG==3 & rG==2){matR["slr2_res2_mlr3", "CDGLS.slrF"] <- matR["slr2_res2_mlr3", "CDGLS.slrF"]+1}
          if(uG==2 & aG==3 & rG==3){matR["slr2_res3_mlr3", "CDGLS.slrF"] <- matR["slr2_res3_mlr3", "CDGLS.slrF"]+1}
          
          if(uG==3 & aG==1 & rG==1){matR["slr3_res1_mlr1", "CDGLS.slrF"] <- matR["slr3_res1_mlr1", "CDGLS.slrF"]+1}
          if(uG==3 & aG==1 & rG==2){matR["slr3_res2_mlr1", "CDGLS.slrF"] <- matR["slr3_res2_mlr1", "CDGLS.slrF"]+1}
          if(uG==3 & aG==1 & rG==3){matR["slr3_res3_mlr1", "CDGLS.slrF"] <- matR["slr3_res3_mlr1", "CDGLS.slrF"]+1}
          
          if(uG==3 & aG==2 & rG==1){matR["slr3_res1_mlr2", "CDGLS.slrF"] <- matR["slr3_res1_mlr2", "CDGLS.slrF"]+1}
          if(uG==3 & aG==2 & rG==2){matR["slr3_res2_mlr2", "CDGLS.slrF"] <- matR["slr3_res2_mlr2", "CDGLS.slrF"]+1}
          if(uG==3 & aG==2 & rG==3){matR["slr3_res3_mlr2", "CDGLS.slrF"] <- matR["slr3_res3_mlr2", "CDGLS.slrF"]+1}
          
          if(uG==3 & aG==3 & rG==1){matR["slr3_res1_mlr3", "CDGLS.slrF"] <- matR["slr3_res1_mlr3", "CDGLS.slrF"]+1}
          if(uG==3 & aG==3 & rG==2){matR["slr3_res2_mlr3", "CDGLS.slrF"] <- matR["slr3_res2_mlr3", "CDGLS.slrF"]+1}
          if(uG==3 & aG==3 & rG==3){matR["slr3_res3_mlr3", "CDGLS.slrF"] <- matR["slr3_res3_mlr3", "CDGLS.slrF"]+1}
        }
        
        # save CD-Egger results
        if("CDEgger.slrF" %in% colnames.mat1){
          
          if(uE==1 & aE==1 & rE==1){matR["slr1_res1_mlr1", "CDEgger.slrF"] <- matR["slr1_res1_mlr1", "CDEgger.slrF"]+1}
          if(uE==1 & aE==1 & rE==2){matR["slr1_res2_mlr1", "CDEgger.slrF"] <- matR["slr1_res2_mlr1", "CDEgger.slrF"]+1}
          if(uE==1 & aE==1 & rE==3){matR["slr1_res3_mlr1", "CDEgger.slrF"] <- matR["slr1_res3_mlr1", "CDEgger.slrF"]+1}
          
          if(uE==1 & aE==2 & rE==1){matR["slr1_res1_mlr2", "CDEgger.slrF"] <- matR["slr1_res1_mlr2", "CDEgger.slrF"]+1}
          if(uE==1 & aE==2 & rE==2){matR["slr1_res2_mlr2", "CDEgger.slrF"] <- matR["slr1_res2_mlr2", "CDEgger.slrF"]+1}
          if(uE==1 & aE==2 & rE==3){matR["slr1_res3_mlr2", "CDEgger.slrF"] <- matR["slr1_res3_mlr2", "CDEgger.slrF"]+1}
          
          if(uE==1 & aE==3 & rE==1){matR["slr1_res1_mlr3", "CDEgger.slrF"] <- matR["slr1_res1_mlr3", "CDEgger.slrF"]+1}
          if(uE==1 & aE==3 & rE==2){matR["slr1_res2_mlr3", "CDEgger.slrF"] <- matR["slr1_res2_mlr3", "CDEgger.slrF"]+1}
          if(uE==1 & aE==3 & rE==3){matR["slr1_res3_mlr3", "CDEgger.slrF"] <- matR["slr1_res3_mlr3", "CDEgger.slrF"]+1}
          
          
          if(uE==2 & aE==1 & rE==1){matR["slr2_res1_mlr1", "CDEgger.slrF"] <- matR["slr2_res1_mlr1", "CDEgger.slrF"]+1}
          if(uE==2 & aE==1 & rE==2){matR["slr2_res2_mlr1", "CDEgger.slrF"] <- matR["slr2_res2_mlr1", "CDEgger.slrF"]+1}
          if(uE==2 & aE==1 & rE==3){matR["slr2_res3_mlr1", "CDEgger.slrF"] <- matR["slr2_res3_mlr1", "CDEgger.slrF"]+1}
          
          if(uE==2 & aE==2 & rE==1){matR["slr2_res1_mlr2", "CDEgger.slrF"] <- matR["slr2_res1_mlr2", "CDEgger.slrF"]+1}
          if(uE==2 & aE==2 & rE==2){matR["slr2_res2_mlr2", "CDEgger.slrF"] <- matR["slr2_res2_mlr2", "CDEgger.slrF"]+1}
          if(uE==2 & aE==2 & rE==3){matR["slr2_res3_mlr2", "CDEgger.slrF"] <- matR["slr2_res3_mlr2", "CDEgger.slrF"]+1}
          
          if(uE==2 & aE==3 & rE==1){matR["slr2_res1_mlr3", "CDEgger.slrF"] <- matR["slr2_res1_mlr3", "CDEgger.slrF"]+1}
          if(uE==2 & aE==3 & rE==2){matR["slr2_res2_mlr3", "CDEgger.slrF"] <- matR["slr2_res2_mlr3", "CDEgger.slrF"]+1}
          if(uE==2 & aE==3 & rE==3){matR["slr2_res3_mlr3", "CDEgger.slrF"] <- matR["slr2_res3_mlr3", "CDEgger.slrF"]+1}
          
          if(uE==3 & aE==1 & rE==1){matR["slr3_res1_mlr1", "CDEgger.slrF"] <- matR["slr3_res1_mlr1", "CDEgger.slrF"]+1}
          if(uE==3 & aE==1 & rE==2){matR["slr3_res2_mlr1", "CDEgger.slrF"] <- matR["slr3_res2_mlr1", "CDEgger.slrF"]+1}
          if(uE==3 & aE==1 & rE==3){matR["slr3_res3_mlr1", "CDEgger.slrF"] <- matR["slr3_res3_mlr1", "CDEgger.slrF"]+1}
          
          if(uE==3 & aE==2 & rE==1){matR["slr3_res1_mlr2", "CDEgger.slrF"] <- matR["slr3_res1_mlr2", "CDEgger.slrF"]+1}
          if(uE==3 & aE==2 & rE==2){matR["slr3_res2_mlr2", "CDEgger.slrF"] <- matR["slr3_res2_mlr2", "CDEgger.slrF"]+1}
          if(uE==3 & aE==2 & rE==3){matR["slr3_res3_mlr2", "CDEgger.slrF"] <- matR["slr3_res3_mlr2", "CDEgger.slrF"]+1}
          
          if(uE==3 & aE==3 & rE==1){matR["slr3_res1_mlr3", "CDEgger.slrF"] <- matR["slr3_res1_mlr3", "CDEgger.slrF"]+1}
          if(uE==3 & aE==3 & rE==2){matR["slr3_res2_mlr3", "CDEgger.slrF"] <- matR["slr3_res2_mlr3", "CDEgger.slrF"]+1}
          if(uE==3 & aE==3 & rE==3){matR["slr3_res3_mlr3", "CDEgger.slrF"] <- matR["slr3_res3_mlr3", "CDEgger.slrF"]+1}
        } # End CD-Egger
        
        
      } # END CD-RATIO / CD-EGGER / CD-GLS SLR FILTERING
      
      # START CD-RATIO / CD-EGGER / CD-GLS MLR FILTERING
      
      
      if(("CDRatio.mlrF" %in% colnames.mat1 | "CDEgger.mlrF" %in% colnames.mat1 | "CDGLS.mlrF" %in% colnames.mat1) & (length(u.betaGX_X.MLR) > 1)){

        ######################################################################################################
        # Unadjusted -- CD-Ratio, CD-Egger, CD-GLS (1-Sample, multiple SNPs):
        ######################################################################################################
        uR <- NA
        uG <- NA
        uE <- NA
        
        u.tmp.mat <- matrix(NA,ncol=12,nrow=length(u.betaGX_X.MLR))
        u.tmp.df <- as.data.frame(u.tmp.mat)
        colnames(u.tmp.df) <- c("chr", "pos", "rsid", "A1", "A2", "beta_T1",
                                "se_T1", "N_T1", "beta_T2", "se_T2", "N_T2", "loci")
        u.tmp.df$rsid<-paste0("SNP", c(1:length(u.betaGX_X.MLR)))
        u.tmp.df$beta_T1 <- u.betaGX_X.MLR
        u.tmp.df$se_T1 <- u.seGX_X.MLR
        u.tmp.df$N_T1 <- n
        u.tmp.df$beta_T2 <- u.betaGX_Y.MLR
        u.tmp.df$se_T2 <- u.seGX_Y.MLR
        u.tmp.df$N_T2 <- n
        
        u.pruned1 <- list(sig_part=u.tmp.df)
        u.cd3 <- CD_3_methods_Independent(u.pruned1$sig_part)
        
        cdDecision.Fun <- function(cd3){
          # Confidence intervals for K - needed for use in decision rules
          # CD-ratio CIs
          
          ratio.YX <- cd3$CD_Ratio_result$T1toT2
          ratio.XY <- cd3$CD_Ratio_result$T2toT1
          egger.YX<- cd3$CD_Egger_result$T1toT2
          egger.XY <- cd3$CD_Egger_result$T2toT1
          gls.YX <- cd3$CD_GLS_result$T1toT2
          gls.XY <- cd3$CD_GLS_result$T2toT1
          
          # Confidence intervals for K - needed for use in decision rules
          # CD-ratio CIs
          lowerCIyx.cdRatio <- ratio.YX["K"] - qnorm(1-sig.level/2)*ratio.YX["se(K)"]
          upperCIyx.cdRatio <- ratio.YX["K"] + qnorm(1-sig.level/2)*ratio.YX["se(K)"]
          lowerCIxy.cdRatio <- ratio.XY["K"] - qnorm(1-sig.level/2)*ratio.XY["se(K)"]
          upperCIxy.cdRatio <- ratio.XY["K"] + qnorm(1-sig.level/2)*ratio.XY["se(K)"]
          
          # CD-Egger CIs
          lowerCIyx.cdEgger <- egger.YX["K"] - qnorm(1-sig.level/2)*egger.YX["se(K)"]
          upperCIyx.cdEgger <- egger.YX["K"] + qnorm(1-sig.level/2)*egger.YX["se(K)"]
          lowerCIxy.cdEgger <- egger.XY["K"] - qnorm(1-sig.level/2)*egger.XY["se(K)"]
          upperCIxy.cdEgger <- egger.XY["K"] + qnorm(1-sig.level/2)*egger.XY["se(K)"]
          
          # CD-GLS CIs
          lowerCIyx.cdGLS <- gls.YX["K"] - qnorm(1-sig.level/2)*gls.YX["se(K)"]
          upperCIyx.cdGLS <- gls.YX["K"] + qnorm(1-sig.level/2)*gls.YX["se(K)"]
          lowerCIxy.cdGLS <- gls.XY["K"] - qnorm(1-sig.level/2)*gls.XY["se(K)"]
          upperCIxy.cdGLS <- gls.XY["K"] + qnorm(1-sig.level/2)*gls.XY["se(K)"]
          
          # case 1: X->Y if CIxy not inside [-1,1] and CIyx inside [-1,0) OR (0,1]
          # case 2: X<-Y if CIyx not inside [-1,1] and CIxy inside [-1,0) OR (0,1]
          # case 3: neither otherwise
          # CD-Ratio decisions
          if(upperCIxy.cdRatio < lowerCIxy.cdRatio){
            stop("LowerCIxy > UpperCIxy")}
          if(upperCIyx.cdRatio < lowerCIyx.cdRatio){
            stop("LowerCIyx > UpperCIyx")}
          
          if(upperCIxy.cdEgger < lowerCIxy.cdEgger){
            stop("LowerCIxyEgger > UpperCIxyEgger")}
          if(upperCIyx.cdEgger < lowerCIyx.cdEgger){
            stop("LowerCIyxEgger > UpperCIyxEgger")}
          
          if(upperCIxy.cdGLS < lowerCIxy.cdGLS){
            stop("LowerCIxy > UpperCIxy")}
          if(upperCIyx.cdGLS < lowerCIyx.cdGLS){
            stop("LowerCIyx > UpperCIyx")}
          
          tmpR <- NA
          # CD-Ratio decisions
          if((upperCIxy.cdRatio < (-1) | lowerCIxy.cdRatio>1) &
             ((lowerCIyx.cdRatio>=(-1) & upperCIyx.cdRatio<0) | (lowerCIyx.cdRatio>(0) & upperCIyx.cdRatio<=1)) ){
            tmpR <- 1
          }else if((upperCIyx.cdRatio < (-1) | lowerCIyx.cdRatio>1) &
                   ((lowerCIxy.cdRatio>=(-1) & upperCIxy.cdRatio< 0) | (lowerCIxy.cdRatio>(0) & upperCIxy.cdRatio<= 1)) ){
            tmpR <- 2
          }else{
            tmpR <- 3
          }
          
          # CD-Egger decisions
          tmpE <- NA
          if((upperCIxy.cdEgger < (-1)  | lowerCIxy.cdEgger>1) &
             ((lowerCIyx.cdEgger>=(-1) & upperCIyx.cdEgger<0) | (lowerCIyx.cdEgger>(0) & upperCIyx.cdEgger<=1))){
            tmpE <- 1
          }else if((upperCIyx.cdEgger < (-1) | lowerCIyx.cdEgger>1) &
                   ((lowerCIxy.cdEgger>=(-1)  & upperCIxy.cdEgger< 0) | (lowerCIxy.cdEgger>(0)  & upperCIxy.cdEgger<=1))){
            tmpE <- 2
          }else{
            tmpE <- 3
          }
          
          # CD-GLS decisions
          tmpG <- NA
          if((upperCIxy.cdGLS < (-1) | lowerCIxy.cdGLS>1) &
             ( (lowerCIyx.cdGLS>= (-1) & upperCIyx.cdGLS<0) | (lowerCIyx.cdGLS> (0) & upperCIyx.cdGLS<=1)) ){
            tmpG <- 1
          }else if((upperCIyx.cdGLS < (-1) | lowerCIyx.cdGLS>1) &
                   ((lowerCIxy.cdGLS>= (-1) & upperCIxy.cdGLS<0) | (lowerCIxy.cdGLS> (0) & upperCIxy.cdGLS<=1)) ){
            tmpG <- 2
          }else{
            tmpG <- 3
          }
          
          return(list("ratio.Case" = tmpR, "egger.Case"=tmpE, "gls.Case" = tmpG))  
        }
        
        
        # unadjusted cases
        u.Out <- cdDecision.Fun(cd3=u.cd3)
        uR <- u.Out$ratio.Case
        uG <- u.Out$gls.Case
        uE <- u.Out$egger.Case
        
        ######################################################################################################
        # Adjusted -- CD-Ratio, CD-Egger, CD-GLS (1-Sample, multiple SNPs):
        ######################################################################################################
        a.tmp.mat <- matrix(NA,ncol=12,nrow=length(a.betaGX_X.MLR))
        a.tmp.df <- as.data.frame(a.tmp.mat)
        colnames(a.tmp.df) <- c("chr", "pos", "rsid", "A1", "A2", "beta_T1",
                                "se_T1", "N_T1", "beta_T2", "se_T2", "N_T2", "loci")
        a.tmp.df$rsid<-paste0("SNP", c(1:length(a.betaGX_X.MLR)))
        a.tmp.df$beta_T1 <- a.betaGX_X.MLR
        a.tmp.df$se_T1 <- a.seGX_X.MLR
        a.tmp.df$N_T1 <- n
        a.tmp.df$beta_T2 <- a.betaGX_Y.MLR
        a.tmp.df$se_T2 <- a.seGX_Y.MLR
        a.tmp.df$N_T2 <- n
        a.pruned1 <- list(sig_part=a.tmp.df)
        
        a.cd3 <- CD_3_methods_Independent(a.pruned1$sig_part)
        aR <- NA
        aG <- NA
        aE <- NA
        
        a.Out <- cdDecision.Fun(cd3=a.cd3)
        aR <- a.Out$ratio.Case
        aG <- a.Out$gls.Case
        aE <- a.Out$egger.Case        
        
        ######################################################################################################
        # Residual -- CD-Ratio, CD-Egger, CD-GLS (1-Sample, multiple SNPs):
        ######################################################################################################
        rR <- NA
        rG <- NA
        rE <- NA
        
        r.tmp.mat <- matrix(NA,ncol=12,nrow=length(r.betaGX_X.MLR))
        r.tmp.df <- as.data.frame(r.tmp.mat)
        colnames(r.tmp.df) <- c("chr", "pos", "rsid", "A1", "A2", "beta_T1",
                                "se_T1", "N_T1", "beta_T2", "se_T2", "N_T2", "loci")
        r.tmp.df$rsid<-paste0("SNP", c(1:length(r.betaGX_X.MLR)))
        r.tmp.df$beta_T1 <- r.betaGX_X.MLR
        r.tmp.df$se_T1 <- r.seGX_X.MLR
        r.tmp.df$N_T1 <- n
        r.tmp.df$beta_T2 <- r.betaGX_Y.MLR
        r.tmp.df$se_T2 <- r.seGX_Y.MLR
        r.tmp.df$N_T2 <- n
        
        r.pruned1 <- list(sig_part=r.tmp.df)
        r.cd3 <- CD_3_methods_Independent(r.pruned1$sig_part)
        
        r.Out <- cdDecision.Fun(cd3=r.cd3)
        rR <- r.Out$ratio.Case
        rG <- r.Out$gls.Case
        rE <- r.Out$egger.Case
        
        
        # SAVE CD-RATIO RESULTS
        if("CDRatio.mlrF" %in% colnames.mat1){
          if(uR==1 & aR==1 & rR==1){matR["slr1_res1_mlr1", "CDRatio.mlrF"] <- matR["slr1_res1_mlr1", "CDRatio.mlrF"]+1}
          if(uR==1 & aR==1 & rR==2){matR["slr1_res2_mlr1", "CDRatio.mlrF"] <- matR["slr1_res2_mlr1", "CDRatio.mlrF"]+1}
          if(uR==1 & aR==1 & rR==3){matR["slr1_res3_mlr1", "CDRatio.mlrF"] <- matR["slr1_res3_mlr1", "CDRatio.mlrF"]+1}
          
          if(uR==1 & aR==2 & rR==1){matR["slr1_res1_mlr2", "CDRatio.mlrF"] <- matR["slr1_res1_mlr2", "CDRatio.mlrF"]+1}
          if(uR==1 & aR==2 & rR==2){matR["slr1_res2_mlr2", "CDRatio.mlrF"] <- matR["slr1_res2_mlr2", "CDRatio.mlrF"]+1}
          if(uR==1 & aR==2 & rR==3){matR["slr1_res3_mlr2", "CDRatio.mlrF"] <- matR["slr1_res3_mlr2", "CDRatio.mlrF"]+1}
          
          if(uR==1 & aR==3 & rR==1){matR["slr1_res1_mlr3", "CDRatio.mlrF"] <- matR["slr1_res1_mlr3", "CDRatio.mlrF"]+1}
          if(uR==1 & aR==3 & rR==2){matR["slr1_res2_mlr3", "CDRatio.mlrF"] <- matR["slr1_res2_mlr3", "CDRatio.mlrF"]+1}
          if(uR==1 & aR==3 & rR==3){matR["slr1_res3_mlr3", "CDRatio.mlrF"] <- matR["slr1_res3_mlr3", "CDRatio.mlrF"]+1}
          
          if(uR==2 & aR==1 & rR==1){matR["slr2_res1_mlr1", "CDRatio.mlrF"] <- matR["slr2_res1_mlr1", "CDRatio.mlrF"]+1}
          if(uR==2 & aR==1 & rR==2){matR["slr2_res2_mlr1", "CDRatio.mlrF"] <- matR["slr2_res2_mlr1", "CDRatio.mlrF"]+1}
          if(uR==2 & aR==1 & rR==3){matR["slr2_res3_mlr1", "CDRatio.mlrF"] <- matR["slr2_res3_mlr1", "CDRatio.mlrF"]+1}
          
          if(uR==2 & aR==2 & rR==1){matR["slr2_res1_mlr2", "CDRatio.mlrF"] <- matR["slr2_res1_mlr2", "CDRatio.mlrF"]+1}
          if(uR==2 & aR==2 & rR==2){matR["slr2_res2_mlr2", "CDRatio.mlrF"] <- matR["slr2_res2_mlr2", "CDRatio.mlrF"]+1}
          if(uR==2 & aR==2 & rR==3){matR["slr2_res3_mlr2", "CDRatio.mlrF"] <- matR["slr2_res3_mlr2", "CDRatio.mlrF"]+1}
          
          if(uR==2 & aR==3 & rR==1){matR["slr2_res1_mlr3", "CDRatio.mlrF"] <- matR["slr2_res1_mlr3", "CDRatio.mlrF"]+1}
          if(uR==2 & aR==3 & rR==2){matR["slr2_res2_mlr3", "CDRatio.mlrF"] <- matR["slr2_res2_mlr3", "CDRatio.mlrF"]+1}
          if(uR==2 & aR==3 & rR==3){matR["slr2_res3_mlr3", "CDRatio.mlrF"] <- matR["slr2_res3_mlr3", "CDRatio.mlrF"]+1}
          
          if(uR==3 & aR==1 & rR==1){matR["slr3_res1_mlr1", "CDRatio.mlrF"] <- matR["slr3_res1_mlr1", "CDRatio.mlrF"]+1}
          if(uR==3 & aR==1 & rR==2){matR["slr3_res2_mlr1", "CDRatio.mlrF"] <- matR["slr3_res2_mlr1", "CDRatio.mlrF"]+1}
          if(uR==3 & aR==1 & rR==3){matR["slr3_res3_mlr1", "CDRatio.mlrF"] <- matR["slr3_res3_mlr1", "CDRatio.mlrF"]+1}
          
          if(uR==3 & aR==2 & rR==1){matR["slr3_res1_mlr2", "CDRatio.mlrF"] <- matR["slr3_res1_mlr2", "CDRatio.mlrF"]+1}
          if(uR==3 & aR==2 & rR==2){matR["slr3_res2_mlr2", "CDRatio.mlrF"] <- matR["slr3_res2_mlr2", "CDRatio.mlrF"]+1}
          if(uR==3 & aR==2 & rR==3){matR["slr3_res3_mlr2", "CDRatio.mlrF"] <- matR["slr3_res3_mlr2", "CDRatio.mlrF"]+1}
          
          if(uR==3 & aR==3 & rR==1){matR["slr3_res1_mlr3", "CDRatio.mlrF"] <- matR["slr3_res1_mlr3", "CDRatio.mlrF"]+1}
          if(uR==3 & aR==3 & rR==2){matR["slr3_res2_mlr3", "CDRatio.mlrF"] <- matR["slr3_res2_mlr3", "CDRatio.mlrF"]+1}
          if(uR==3 & aR==3 & rR==3){matR["slr3_res3_mlr3", "CDRatio.mlrF"] <- matR["slr3_res3_mlr3", "CDRatio.mlrF"]+1}
        }
        
        # Save CD-GLS results
        if("CDGLS.mlrF" %in% colnames.mat1){
          
          if(uG==1 & aG==1 & rG==1){matR["slr1_res1_mlr1", "CDGLS.mlrF"] <- matR["slr1_res1_mlr1", "CDGLS.mlrF"]+1}
          if(uG==1 & aG==1 & rG==2){matR["slr1_res2_mlr1", "CDGLS.mlrF"] <- matR["slr1_res2_mlr1", "CDGLS.mlrF"]+1}
          if(uG==1 & aG==1 & rG==3){matR["slr1_res3_mlr1", "CDGLS.mlrF"] <- matR["slr1_res3_mlr1", "CDGLS.mlrF"]+1}
          
          if(uG==1 & aG==2 & rG==1){matR["slr1_res1_mlr2", "CDGLS.mlrF"] <- matR["slr1_res1_mlr2", "CDGLS.mlrF"]+1}
          if(uG==1 & aG==2 & rG==2){matR["slr1_res2_mlr2", "CDGLS.mlrF"] <- matR["slr1_res2_mlr2", "CDGLS.mlrF"]+1}
          if(uG==1 & aG==2 & rG==3){matR["slr1_res3_mlr2", "CDGLS.mlrF"] <- matR["slr1_res3_mlr2", "CDGLS.mlrF"]+1}
          
          if(uG==1 & aG==3 & rG==1){matR["slr1_res1_mlr3", "CDGLS.mlrF"] <- matR["slr1_res1_mlr3", "CDGLS.mlrF"]+1}
          if(uG==1 & aG==3 & rG==2){matR["slr1_res2_mlr3", "CDGLS.mlrF"] <- matR["slr1_res2_mlr3", "CDGLS.mlrF"]+1}
          if(uG==1 & aG==3 & rG==3){matR["slr1_res3_mlr3", "CDGLS.mlrF"] <- matR["slr1_res3_mlr3", "CDGLS.mlrF"]+1}
          
          if(uG==2 & aG==1 & rG==1){matR["slr2_res1_mlr1", "CDGLS.mlrF"] <- matR["slr2_res1_mlr1", "CDGLS.mlrF"]+1}
          if(uG==2 & aG==1 & rG==2){matR["slr2_res2_mlr1", "CDGLS.mlrF"] <- matR["slr2_res2_mlr1", "CDGLS.mlrF"]+1}
          if(uG==2 & aG==1 & rG==3){matR["slr2_res3_mlr1", "CDGLS.mlrF"] <- matR["slr2_res3_mlr1", "CDGLS.mlrF"]+1}
          
          if(uG==2 & aG==2 & rG==1){matR["slr2_res1_mlr2", "CDGLS.mlrF"] <- matR["slr2_res1_mlr2", "CDGLS.mlrF"]+1}
          if(uG==2 & aG==2 & rG==2){matR["slr2_res2_mlr2", "CDGLS.mlrF"] <- matR["slr2_res2_mlr2", "CDGLS.mlrF"]+1}
          if(uG==2 & aG==2 & rG==3){matR["slr2_res3_mlr2", "CDGLS.mlrF"] <- matR["slr2_res3_mlr2", "CDGLS.mlrF"]+1}
          
          if(uG==2 & aG==3 & rG==1){matR["slr2_res1_mlr3", "CDGLS.mlrF"] <- matR["slr2_res1_mlr3", "CDGLS.mlrF"]+1}
          if(uG==2 & aG==3 & rG==2){matR["slr2_res2_mlr3", "CDGLS.mlrF"] <- matR["slr2_res2_mlr3", "CDGLS.mlrF"]+1}
          if(uG==2 & aG==3 & rG==3){matR["slr2_res3_mlr3", "CDGLS.mlrF"] <- matR["slr2_res3_mlr3", "CDGLS.mlrF"]+1}
          
          if(uG==3 & aG==1 & rG==1){matR["slr3_res1_mlr1", "CDGLS.mlrF"] <- matR["slr3_res1_mlr1", "CDGLS.mlrF"]+1}
          if(uG==3 & aG==1 & rG==2){matR["slr3_res2_mlr1", "CDGLS.mlrF"] <- matR["slr3_res2_mlr1", "CDGLS.mlrF"]+1}
          if(uG==3 & aG==1 & rG==3){matR["slr3_res3_mlr1", "CDGLS.mlrF"] <- matR["slr3_res3_mlr1", "CDGLS.mlrF"]+1}
          
          if(uG==3 & aG==2 & rG==1){matR["slr3_res1_mlr2", "CDGLS.mlrF"] <- matR["slr3_res1_mlr2", "CDGLS.mlrF"]+1}
          if(uG==3 & aG==2 & rG==2){matR["slr3_res2_mlr2", "CDGLS.mlrF"] <- matR["slr3_res2_mlr2", "CDGLS.mlrF"]+1}
          if(uG==3 & aG==2 & rG==3){matR["slr3_res3_mlr2", "CDGLS.mlrF"] <- matR["slr3_res3_mlr2", "CDGLS.mlrF"]+1}
          
          if(uG==3 & aG==3 & rG==1){matR["slr3_res1_mlr3", "CDGLS.mlrF"] <- matR["slr3_res1_mlr3", "CDGLS.mlrF"]+1}
          if(uG==3 & aG==3 & rG==2){matR["slr3_res2_mlr3", "CDGLS.mlrF"] <- matR["slr3_res2_mlr3", "CDGLS.mlrF"]+1}
          if(uG==3 & aG==3 & rG==3){matR["slr3_res3_mlr3", "CDGLS.mlrF"] <- matR["slr3_res3_mlr3", "CDGLS.mlrF"]+1}
        }
        
        # save CD-Egger results
        if("CDEgger.mlrF" %in% colnames.mat1){
          
          if(uE==1 & aE==1 & rE==1){matR["slr1_res1_mlr1", "CDEgger.mlrF"] <- matR["slr1_res1_mlr1", "CDEgger.mlrF"]+1}
          if(uE==1 & aE==1 & rE==2){matR["slr1_res2_mlr1", "CDEgger.mlrF"] <- matR["slr1_res2_mlr1", "CDEgger.mlrF"]+1}
          if(uE==1 & aE==1 & rE==3){matR["slr1_res3_mlr1", "CDEgger.mlrF"] <- matR["slr1_res3_mlr1", "CDEgger.mlrF"]+1}
          
          if(uE==1 & aE==2 & rE==1){matR["slr1_res1_mlr2", "CDEgger.mlrF"] <- matR["slr1_res1_mlr2", "CDEgger.mlrF"]+1}
          if(uE==1 & aE==2 & rE==2){matR["slr1_res2_mlr2", "CDEgger.mlrF"] <- matR["slr1_res2_mlr2", "CDEgger.mlrF"]+1}
          if(uE==1 & aE==2 & rE==3){matR["slr1_res3_mlr2", "CDEgger.mlrF"] <- matR["slr1_res3_mlr2", "CDEgger.mlrF"]+1}
          
          if(uE==1 & aE==3 & rE==1){matR["slr1_res1_mlr3", "CDEgger.mlrF"] <- matR["slr1_res1_mlr3", "CDEgger.mlrF"]+1}
          if(uE==1 & aE==3 & rE==2){matR["slr1_res2_mlr3", "CDEgger.mlrF"] <- matR["slr1_res2_mlr3", "CDEgger.mlrF"]+1}
          if(uE==1 & aE==3 & rE==3){matR["slr1_res3_mlr3", "CDEgger.mlrF"] <- matR["slr1_res3_mlr3", "CDEgger.mlrF"]+1}
          
          
          if(uE==2 & aE==1 & rE==1){matR["slr2_res1_mlr1", "CDEgger.mlrF"] <- matR["slr2_res1_mlr1", "CDEgger.mlrF"]+1}
          if(uE==2 & aE==1 & rE==2){matR["slr2_res2_mlr1", "CDEgger.mlrF"] <- matR["slr2_res2_mlr1", "CDEgger.mlrF"]+1}
          if(uE==2 & aE==1 & rE==3){matR["slr2_res3_mlr1", "CDEgger.mlrF"] <- matR["slr2_res3_mlr1", "CDEgger.mlrF"]+1}
          
          if(uE==2 & aE==2 & rE==1){matR["slr2_res1_mlr2", "CDEgger.mlrF"] <- matR["slr2_res1_mlr2", "CDEgger.mlrF"]+1}
          if(uE==2 & aE==2 & rE==2){matR["slr2_res2_mlr2", "CDEgger.mlrF"] <- matR["slr2_res2_mlr2", "CDEgger.mlrF"]+1}
          if(uE==2 & aE==2 & rE==3){matR["slr2_res3_mlr2", "CDEgger.mlrF"] <- matR["slr2_res3_mlr2", "CDEgger.mlrF"]+1}
          
          if(uE==2 & aE==3 & rE==1){matR["slr2_res1_mlr3", "CDEgger.mlrF"] <- matR["slr2_res1_mlr3", "CDEgger.mlrF"]+1}
          if(uE==2 & aE==3 & rE==2){matR["slr2_res2_mlr3", "CDEgger.mlrF"] <- matR["slr2_res2_mlr3", "CDEgger.mlrF"]+1}
          if(uE==2 & aE==3 & rE==3){matR["slr2_res3_mlr3", "CDEgger.mlrF"] <- matR["slr2_res3_mlr3", "CDEgger.mlrF"]+1}
          
          if(uE==3 & aE==1 & rE==1){matR["slr3_res1_mlr1", "CDEgger.mlrF"] <- matR["slr3_res1_mlr1", "CDEgger.mlrF"]+1}
          if(uE==3 & aE==1 & rE==2){matR["slr3_res2_mlr1", "CDEgger.mlrF"] <- matR["slr3_res2_mlr1", "CDEgger.mlrF"]+1}
          if(uE==3 & aE==1 & rE==3){matR["slr3_res3_mlr1", "CDEgger.mlrF"] <- matR["slr3_res3_mlr1", "CDEgger.mlrF"]+1}
          
          if(uE==3 & aE==2 & rE==1){matR["slr3_res1_mlr2", "CDEgger.mlrF"] <- matR["slr3_res1_mlr2", "CDEgger.mlrF"]+1}
          if(uE==3 & aE==2 & rE==2){matR["slr3_res2_mlr2", "CDEgger.mlrF"] <- matR["slr3_res2_mlr2", "CDEgger.mlrF"]+1}
          if(uE==3 & aE==2 & rE==3){matR["slr3_res3_mlr2", "CDEgger.mlrF"] <- matR["slr3_res3_mlr2", "CDEgger.mlrF"]+1}
          
          if(uE==3 & aE==3 & rE==1){matR["slr3_res1_mlr3", "CDEgger.mlrF"] <- matR["slr3_res1_mlr3", "CDEgger.mlrF"]+1}
          if(uE==3 & aE==3 & rE==2){matR["slr3_res2_mlr3", "CDEgger.mlrF"] <- matR["slr3_res2_mlr3", "CDEgger.mlrF"]+1}
          if(uE==3 & aE==3 & rE==3){matR["slr3_res3_mlr3", "CDEgger.mlrF"] <- matR["slr3_res3_mlr3", "CDEgger.mlrF"]+1}
        } # End CD-Egger
        
        
      }
      # END CD-RATIO / CD-EGGER / CD-GLS MLR FILTERING
      
      
      if("cdRatio.S" %in% colnames.mat1 | "cdEgger.S" %in% colnames.mat1 | "cdEgger.NoS" %in% colnames.mat1 | "cdRatio.NoS" %in% colnames.mat1){
        bdcd.ratio.U <- NA
        bdcd.ratio.A <- NA
        bdcd.ratio.R <- NA
        
        bdcd.egger.U <- NA
        bdcd.egger.A <- NA
        bdcd.egger.R <- NA
        
        
        CD.Out.U <- BiDirCDMethod(b_X = bxU,
                                  b_Y = byU,
                                  se_X = sxU,
                                  se_Y = syU,
                                  n_X = n,
                                  n_Y = n,
                                  sig.cutoff = pFilterBDC, random.seed = 0)
        
        CD.Out.A <- BiDirCDMethod(b_X = bxA,
                                  b_Y = byA,
                                  se_X = sxA,
                                  se_Y = syA,
                                  n_X = n,
                                  n_Y = n,
                                  sig.cutoff = pFilterBDC, random.seed = 0)
        
        CD.Out.R <- BiDirCDMethod(b_X = bxR,
                                  b_Y = byR,
                                  se_X = sxR,
                                  se_Y = syR,
                                  n_X = n,
                                  n_Y = n,
                                  sig.cutoff = pFilterBDC, random.seed = 0)
        
        cdCasesFun <- function(K.ratio.XY.S, seK.ratio.XY.S, K.egger.XY.S, seK.egger.XY.S, 
                               K.ratio.YX.S, seK.ratio.YX.S, K.egger.YX.S, seK.egger.YX.S){
          bdcd.ratio.tmp <- NA
          bdcd.egger.tmp <- NA
          
          ci.low.ratioXY.S <- K.ratio.XY.S - (qnorm(0.975)*seK.ratio.XY.S)
          ci.upp.ratioXY.S <- K.ratio.XY.S + (qnorm(0.975)*seK.ratio.XY.S)
          
          ci.low.ratioYX.S <- K.ratio.YX.S - (qnorm(0.975)*seK.ratio.YX.S)
          ci.upp.ratioYX.S <- K.ratio.YX.S + (qnorm(0.975)*seK.ratio.YX.S)
          
          ci.low.eggerXY.S <- K.egger.XY.S - (qnorm(0.975)*seK.egger.XY.S)
          ci.upp.eggerXY.S <- K.egger.XY.S + (qnorm(0.975)*seK.egger.XY.S)
          
          ci.low.eggerYX.S <- K.egger.YX.S - (qnorm(0.975)*seK.egger.YX.S)
          ci.upp.eggerYX.S <- K.egger.YX.S + (qnorm(0.975)*seK.egger.YX.S)
          
          #  DECISION RULES CD Ratio - .S
          if((ci.low.ratioXY.S>0 & ci.upp.ratioXY.S<1)| (ci.low.ratioXY.S>-1 & ci.upp.ratioXY.S<0)){
            if((ci.low.ratioYX.S>0 & ci.upp.ratioYX.S<1) | (ci.low.ratioYX.S>-1 & ci.upp.ratioYX.S<0)){ bdcd.ratio.tmp <- 3
            }else{bdcd.ratio.tmp <- 1}
          }else if((ci.low.ratioYX.S>0 & ci.upp.ratioYX.S<1) | (ci.low.ratioYX.S>-1 & ci.upp.ratioYX.S<0)){
            if((ci.low.ratioXY.S>0 & ci.upp.ratioXY.S<1)| (ci.low.ratioXY.S>-1 & ci.upp.ratioXY.S<0)){ bdcd.ratio.tmp <- 3
            }else{bdcd.ratio.tmp <- 2}
          }else{bdcd.ratio.tmp <- 3}
          
          #  DECISION RULES CD Egger
          if((ci.low.eggerXY.S>0 & ci.upp.eggerXY.S<1)| (ci.low.eggerXY.S>-1 & ci.upp.eggerXY.S<0)){
            if((ci.low.eggerYX.S>0 & ci.upp.eggerYX.S<1) | (ci.low.eggerYX.S>-1 & ci.upp.eggerYX.S<0)){bdcd.egger.tmp <- 3
            }else{bdcd.egger.tmp <- 1}
          }else if((ci.low.eggerYX.S>0 & ci.upp.eggerYX.S<1) | (ci.low.eggerYX.S>-1 & ci.upp.eggerYX.S<0)){
            if((ci.low.eggerXY.S>0 & ci.upp.eggerXY.S<1)| (ci.low.eggerXY.S>-1 & ci.upp.eggerXY.S<0)){bdcd.egger.tmp <- 3
            }else{bdcd.egger.tmp <- 2}
          }else{bdcd.egger.tmp <- 3}
          
          return(list("bdcd.ratio" = bdcd.ratio.tmp, "bdcd.egger" = bdcd.egger.tmp))
        }
        
        
        if("cdRatio.S" %in% colnames.mat1 | "cdEgger.S" %in% colnames.mat1){

            CD.Cases.U <- cdCasesFun(K.ratio.XY.S = CD.Out.U$CDRatio.XtoY.est.S,
                                     seK.ratio.XY.S = CD.Out.U$CDRatio.XtoY.se.S,
                                     K.egger.XY.S = CD.Out.U$CDEgger.XtoY.est.S,
                                     seK.egger.XY.S = CD.Out.U$CDEgger.XtoY.se.S,
                                     K.ratio.YX.S = CD.Out.U$CDRatio.YtoX.est.S,
                                     seK.ratio.YX.S = CD.Out.U$CDRatio.YtoX.se.S,
                                     K.egger.YX.S = CD.Out.U$CDEgger.YtoX.est.S,
                                     seK.egger.YX.S = CD.Out.U$CDEgger.YtoX.se.S)
            
            bdcd.ratio.S.U <- CD.Cases.U$bdcd.ratio
            bdcd.egger.S.U <- CD.Cases.U$bdcd.egger
            
            CD.Cases.A <- cdCasesFun(K.ratio.XY.S = CD.Out.A$CDRatio.XtoY.est.S,
                                     seK.ratio.XY.S = CD.Out.A$CDRatio.XtoY.se.S,
                                     K.egger.XY.S = CD.Out.A$CDEgger.XtoY.est.S,
                                     seK.egger.XY.S = CD.Out.A$CDEgger.XtoY.se.S,
                                     K.ratio.YX.S = CD.Out.A$CDRatio.YtoX.est.S,
                                     seK.ratio.YX.S = CD.Out.A$CDRatio.YtoX.se.S,
                                     K.egger.YX.S = CD.Out.A$CDEgger.YtoX.est.S,
                                     seK.egger.YX.S = CD.Out.A$CDEgger.YtoX.se.S)
            bdcd.ratio.S.A <- CD.Cases.A$bdcd.ratio
            bdcd.egger.S.A <- CD.Cases.A$bdcd.egger
            
            CD.Cases.R <- cdCasesFun(K.ratio.XY.S = CD.Out.R$CDRatio.XtoY.est.S,
                                     seK.ratio.XY.S = CD.Out.R$CDRatio.XtoY.se.S,
                                     K.egger.XY.S = CD.Out.R$CDEgger.XtoY.est.S,
                                     seK.egger.XY.S = CD.Out.R$CDEgger.XtoY.se.S,
                                     K.ratio.YX.S = CD.Out.R$CDRatio.YtoX.est.S,
                                     seK.ratio.YX.S = CD.Out.R$CDRatio.YtoX.se.S,
                                     K.egger.YX.S = CD.Out.R$CDEgger.YtoX.est.S,
                                     seK.egger.YX.S = CD.Out.R$CDEgger.YtoX.se.S)
            
            bdcd.ratio.S.R <- CD.Cases.R$bdcd.ratio
            bdcd.egger.S.R <- CD.Cases.R$bdcd.egger
            
            if("cdRatio.S" %in% colnames.mat1){
                # decision rules - cd ratio bidircausal screening
                if(bdcd.ratio.S.U==1 & bdcd.ratio.S.A==1 & bdcd.ratio.S.R==1){matR["slr1_res1_mlr1", "cdRatio.S"] <- matR["slr1_res1_mlr1", "cdRatio.S"]+1}
                if(bdcd.ratio.S.U==1 & bdcd.ratio.S.A==1 & bdcd.ratio.S.R==2){matR["slr1_res2_mlr1", "cdRatio.S"] <- matR["slr1_res2_mlr1", "cdRatio.S"]+1}
                if(bdcd.ratio.S.U==1 & bdcd.ratio.S.A==1 & bdcd.ratio.S.R==3){matR["slr1_res3_mlr1", "cdRatio.S"] <- matR["slr1_res3_mlr1", "cdRatio.S"]+1}
                
                if(bdcd.ratio.S.U==1 & bdcd.ratio.S.A==2 & bdcd.ratio.S.R==1){matR["slr1_res1_mlr2", "cdRatio.S"] <- matR["slr1_res1_mlr2", "cdRatio.S"]+1}
                if(bdcd.ratio.S.U==1 & bdcd.ratio.S.A==2 & bdcd.ratio.S.R==2){matR["slr1_res2_mlr2", "cdRatio.S"] <- matR["slr1_res2_mlr2", "cdRatio.S"]+1}
                if(bdcd.ratio.S.U==1 & bdcd.ratio.S.A==2 & bdcd.ratio.S.R==3){matR["slr1_res3_mlr2", "cdRatio.S"] <- matR["slr1_res3_mlr2", "cdRatio.S"]+1}
                
                if(bdcd.ratio.S.U==1 & bdcd.ratio.S.A==3 & bdcd.ratio.S.R==1){matR["slr1_res1_mlr3", "cdRatio.S"] <- matR["slr1_res1_mlr3", "cdRatio.S"]+1}
                if(bdcd.ratio.S.U==1 & bdcd.ratio.S.A==3 & bdcd.ratio.S.R==2){matR["slr1_res2_mlr3", "cdRatio.S"] <- matR["slr1_res2_mlr3", "cdRatio.S"]+1}
                if(bdcd.ratio.S.U==1 & bdcd.ratio.S.A==3 & bdcd.ratio.S.R==3){matR["slr1_res3_mlr3", "cdRatio.S"] <- matR["slr1_res3_mlr3", "cdRatio.S"]+1}
                
                if(bdcd.ratio.S.U==2 & bdcd.ratio.S.A==1 & bdcd.ratio.S.R==1){matR["slr2_res1_mlr1", "cdRatio.S"] <- matR["slr2_res1_mlr1", "cdRatio.S"]+1}
                if(bdcd.ratio.S.U==2 & bdcd.ratio.S.A==1 & bdcd.ratio.S.R==2){matR["slr2_res2_mlr1", "cdRatio.S"] <- matR["slr2_res2_mlr1", "cdRatio.S"]+1}
                if(bdcd.ratio.S.U==2 & bdcd.ratio.S.A==1 & bdcd.ratio.S.R==3){matR["slr2_res3_mlr1", "cdRatio.S"] <- matR["slr2_res3_mlr1", "cdRatio.S"]+1}
                
                if(bdcd.ratio.S.U==2 & bdcd.ratio.S.A==2 & bdcd.ratio.S.R==1){matR["slr2_res1_mlr2", "cdRatio.S"] <- matR["slr2_res1_mlr2", "cdRatio.S"]+1}
                if(bdcd.ratio.S.U==2 & bdcd.ratio.S.A==2 & bdcd.ratio.S.R==2){matR["slr2_res2_mlr2", "cdRatio.S"] <- matR["slr2_res2_mlr2", "cdRatio.S"]+1}
                if(bdcd.ratio.S.U==2 & bdcd.ratio.S.A==2 & bdcd.ratio.S.R==3){matR["slr2_res3_mlr2", "cdRatio.S"] <- matR["slr2_res3_mlr2", "cdRatio.S"]+1}
                
                if(bdcd.ratio.S.U==2 & bdcd.ratio.S.A==3 & bdcd.ratio.S.R==1){matR["slr2_res1_mlr3", "cdRatio.S"] <- matR["slr2_res1_mlr3", "cdRatio.S"]+1}
                if(bdcd.ratio.S.U==2 & bdcd.ratio.S.A==3 & bdcd.ratio.S.R==2){matR["slr2_res2_mlr3", "cdRatio.S"] <- matR["slr2_res2_mlr3", "cdRatio.S"]+1}
                if(bdcd.ratio.S.U==2 & bdcd.ratio.S.A==3 & bdcd.ratio.S.R==3){matR["slr2_res3_mlr3", "cdRatio.S"] <- matR["slr2_res3_mlr3", "cdRatio.S"]+1}
                
                if(bdcd.ratio.S.U==3 & bdcd.ratio.S.A==1 & bdcd.ratio.S.R==1){matR["slr3_res1_mlr1", "cdRatio.S"] <- matR["slr3_res1_mlr1", "cdRatio.S"]+1}
                if(bdcd.ratio.S.U==3 & bdcd.ratio.S.A==1 & bdcd.ratio.S.R==2){matR["slr3_res2_mlr1", "cdRatio.S"] <- matR["slr3_res2_mlr1", "cdRatio.S"]+1}
                if(bdcd.ratio.S.U==3 & bdcd.ratio.S.A==1 & bdcd.ratio.S.R==3){matR["slr3_res3_mlr1", "cdRatio.S"] <- matR["slr3_res3_mlr1", "cdRatio.S"]+1}
                
                if(bdcd.ratio.S.U==3 & bdcd.ratio.S.A==2 & bdcd.ratio.S.R==1){matR["slr3_res1_mlr2", "cdRatio.S"] <- matR["slr3_res1_mlr2", "cdRatio.S"]+1}
                if(bdcd.ratio.S.U==3 & bdcd.ratio.S.A==2 & bdcd.ratio.S.R==2){matR["slr3_res2_mlr2", "cdRatio.S"] <- matR["slr3_res2_mlr2", "cdRatio.S"]+1}
                if(bdcd.ratio.S.U==3 & bdcd.ratio.S.A==2 & bdcd.ratio.S.R==3){matR["slr3_res3_mlr2", "cdRatio.S"] <- matR["slr3_res3_mlr2", "cdRatio.S"]+1}
                
                if(bdcd.ratio.S.U==3 & bdcd.ratio.S.A==3 & bdcd.ratio.S.R==1){matR["slr3_res1_mlr3", "cdRatio.S"] <- matR["slr3_res1_mlr3", "cdRatio.S"]+1}
                if(bdcd.ratio.S.U==3 & bdcd.ratio.S.A==3 & bdcd.ratio.S.R==2){matR["slr3_res2_mlr3", "cdRatio.S"] <- matR["slr3_res2_mlr3", "cdRatio.S"]+1}
                if(bdcd.ratio.S.U==3 & bdcd.ratio.S.A==3 & bdcd.ratio.S.R==3){matR["slr3_res3_mlr3", "cdRatio.S"] <- matR["slr3_res3_mlr3", "cdRatio.S"]+1}
           }
            
            if("cdEgger.S" %in% colnames.mat1){

                # decision rules - cd egger bidircausal screening
                if(bdcd.egger.S.U==1 & bdcd.egger.S.A==1 & bdcd.egger.S.R==1){matR["slr1_res1_mlr1", "cdEgger.S"] <- matR["slr1_res1_mlr1", "cdEgger.S"]+1}
                if(bdcd.egger.S.U==1 & bdcd.egger.S.A==1 & bdcd.egger.S.R==2){matR["slr1_res2_mlr1", "cdEgger.S"] <- matR["slr1_res2_mlr1", "cdEgger.S"]+1}
                if(bdcd.egger.S.U==1 & bdcd.egger.S.A==1 & bdcd.egger.S.R==3){matR["slr1_res3_mlr1", "cdEgger.S"] <- matR["slr1_res3_mlr1", "cdEgger.S"]+1}
                
                if(bdcd.egger.S.U==1 & bdcd.egger.S.A==2 & bdcd.egger.S.R==1){matR["slr1_res1_mlr2", "cdEgger.S"] <- matR["slr1_res1_mlr2", "cdEgger.S"]+1}
                if(bdcd.egger.S.U==1 & bdcd.egger.S.A==2 & bdcd.egger.S.R==2){matR["slr1_res2_mlr2", "cdEgger.S"] <- matR["slr1_res2_mlr2", "cdEgger.S"]+1}
                if(bdcd.egger.S.U==1 & bdcd.egger.S.A==2 & bdcd.egger.S.R==3){matR["slr1_res3_mlr2", "cdEgger.S"] <- matR["slr1_res3_mlr2", "cdEgger.S"]+1}
                
                if(bdcd.egger.S.U==1 & bdcd.egger.S.A==3 & bdcd.egger.S.R==1){matR["slr1_res1_mlr3", "cdEgger.S"] <- matR["slr1_res1_mlr3", "cdEgger.S"]+1}
                if(bdcd.egger.S.U==1 & bdcd.egger.S.A==3 & bdcd.egger.S.R==2){matR["slr1_res2_mlr3", "cdEgger.S"] <- matR["slr1_res2_mlr3", "cdEgger.S"]+1}
                if(bdcd.egger.S.U==1 & bdcd.egger.S.A==3 & bdcd.egger.S.R==3){matR["slr1_res3_mlr3", "cdEgger.S"] <- matR["slr1_res3_mlr3", "cdEgger.S"]+1}
                
                if(bdcd.egger.S.U==2 & bdcd.egger.S.A==1 & bdcd.egger.S.R==1){matR["slr2_res1_mlr1", "cdEgger.S"] <- matR["slr2_res1_mlr1", "cdEgger.S"]+1}
                if(bdcd.egger.S.U==2 & bdcd.egger.S.A==1 & bdcd.egger.S.R==2){matR["slr2_res2_mlr1", "cdEgger.S"] <- matR["slr2_res2_mlr1", "cdEgger.S"]+1}
                if(bdcd.egger.S.U==2 & bdcd.egger.S.A==1 & bdcd.egger.S.R==3){matR["slr2_res3_mlr1", "cdEgger.S"] <- matR["slr2_res3_mlr1", "cdEgger.S"]+1}
                
                if(bdcd.egger.S.U==2 & bdcd.egger.S.A==2 & bdcd.egger.S.R==1){matR["slr2_res1_mlr2", "cdEgger.S"] <- matR["slr2_res1_mlr2", "cdEgger.S"]+1}
                if(bdcd.egger.S.U==2 & bdcd.egger.S.A==2 & bdcd.egger.S.R==2){matR["slr2_res2_mlr2", "cdEgger.S"] <- matR["slr2_res2_mlr2", "cdEgger.S"]+1}
                if(bdcd.egger.S.U==2 & bdcd.egger.S.A==2 & bdcd.egger.S.R==3){matR["slr2_res3_mlr2", "cdEgger.S"] <- matR["slr2_res3_mlr2", "cdEgger.S"]+1}
                
                if(bdcd.egger.S.U==2 & bdcd.egger.S.A==3 & bdcd.egger.S.R==1){matR["slr2_res1_mlr3", "cdEgger.S"] <- matR["slr2_res1_mlr3", "cdEgger.S"]+1}
                if(bdcd.egger.S.U==2 & bdcd.egger.S.A==3 & bdcd.egger.S.R==2){matR["slr2_res2_mlr3", "cdEgger.S"] <- matR["slr2_res2_mlr3", "cdEgger.S"]+1}
                if(bdcd.egger.S.U==2 & bdcd.egger.S.A==3 & bdcd.egger.S.R==3){matR["slr2_res3_mlr3", "cdEgger.S"] <- matR["slr2_res3_mlr3", "cdEgger.S"]+1}
                
                if(bdcd.egger.S.U==3 & bdcd.egger.S.A==1 & bdcd.egger.S.R==1){matR["slr3_res1_mlr1", "cdEgger.S"] <- matR["slr3_res1_mlr1", "cdEgger.S"]+1}
                if(bdcd.egger.S.U==3 & bdcd.egger.S.A==1 & bdcd.egger.S.R==2){matR["slr3_res2_mlr1", "cdEgger.S"] <- matR["slr3_res2_mlr1", "cdEgger.S"]+1}
                if(bdcd.egger.S.U==3 & bdcd.egger.S.A==1 & bdcd.egger.S.R==3){matR["slr3_res3_mlr1", "cdEgger.S"] <- matR["slr3_res3_mlr1", "cdEgger.S"]+1}
                
                if(bdcd.egger.S.U==3 & bdcd.egger.S.A==2 & bdcd.egger.S.R==1){matR["slr3_res1_mlr2", "cdEgger.S"] <- matR["slr3_res1_mlr2", "cdEgger.S"]+1}
                if(bdcd.egger.S.U==3 & bdcd.egger.S.A==2 & bdcd.egger.S.R==2){matR["slr3_res2_mlr2", "cdEgger.S"] <- matR["slr3_res2_mlr2", "cdEgger.S"]+1}
                if(bdcd.egger.S.U==3 & bdcd.egger.S.A==2 & bdcd.egger.S.R==3){matR["slr3_res3_mlr2", "cdEgger.S"] <- matR["slr3_res3_mlr2", "cdEgger.S"]+1}
                
                if(bdcd.egger.S.U==3 & bdcd.egger.S.A==3 & bdcd.egger.S.R==1){matR["slr3_res1_mlr3", "cdEgger.S"] <- matR["slr3_res1_mlr3", "cdEgger.S"]+1}
                if(bdcd.egger.S.U==3 & bdcd.egger.S.A==3 & bdcd.egger.S.R==2){matR["slr3_res2_mlr3", "cdEgger.S"] <- matR["slr3_res2_mlr3", "cdEgger.S"]+1}
                if(bdcd.egger.S.U==3 & bdcd.egger.S.A==3 & bdcd.egger.S.R==3){matR["slr3_res3_mlr3", "cdEgger.S"] <- matR["slr3_res3_mlr3", "cdEgger.S"]+1}
            }
            
        }
        
        ##################################
        # Save no screening results
        ##################################
        
        if("cdEgger.NoS" %in% colnames.mat1 | "cdRatio.NoS" %in% colnames.mat1){

        
            bdcd.ratio.NoS.U <- NA
            bdcd.ratio.NoS.A <- NA
            bdcd.ratio.NoS.R <- NA
            
            bdcd.egger.NoS.U <- NA
            bdcd.egger.NoS.A <- NA
            bdcd.egger.NoS.R <- NA
            
            CD.Cases.NoS.U <- cdCasesFun(K.ratio.XY.S = CD.Out.U$CDRatio.XtoY.est.NoS,
                                     seK.ratio.XY.S = CD.Out.U$CDRatio.XtoY.se.NoS,
                                     K.egger.XY.S = CD.Out.U$CDEgger.XtoY.est.NoS,
                                     seK.egger.XY.S = CD.Out.U$CDEgger.XtoY.se.NoS,
                                     K.ratio.YX.S = CD.Out.U$CDRatio.YtoX.est.NoS,
                                     seK.ratio.YX.S = CD.Out.U$CDRatio.YtoX.se.NoS,
                                     K.egger.YX.S = CD.Out.U$CDEgger.YtoX.est.NoS,
                                     seK.egger.YX.S = CD.Out.U$CDEgger.YtoX.se.NoS)
            bdcd.ratio.NoS.U <- CD.Cases.NoS.U$bdcd.ratio
            bdcd.egger.NoS.U <- CD.Cases.NoS.U$bdcd.egger
            
            CD.Cases.NoS.A <- cdCasesFun(K.ratio.XY.S = CD.Out.A$CDRatio.XtoY.est.NoS,
                                     seK.ratio.XY.S = CD.Out.A$CDRatio.XtoY.se.NoS,
                                     K.egger.XY.S = CD.Out.A$CDEgger.XtoY.est.NoS,
                                     seK.egger.XY.S = CD.Out.A$CDEgger.XtoY.se.NoS,
                                     K.ratio.YX.S = CD.Out.A$CDRatio.YtoX.est.NoS,
                                     seK.ratio.YX.S = CD.Out.A$CDRatio.YtoX.se.NoS,
                                     K.egger.YX.S = CD.Out.A$CDEgger.YtoX.est.NoS,
                                     seK.egger.YX.S = CD.Out.A$CDEgger.YtoX.se.NoS)
            bdcd.ratio.NoS.A <- CD.Cases.NoS.A$bdcd.ratio
            bdcd.egger.NoS.A <- CD.Cases.NoS.A$bdcd.egger
            
            CD.Cases.NoS.R <- cdCasesFun(K.ratio.XY.S = CD.Out.R$CDRatio.XtoY.est.NoS,
                                     seK.ratio.XY.S = CD.Out.R$CDRatio.XtoY.se.NoS,
                                     K.egger.XY.S = CD.Out.R$CDEgger.XtoY.est.NoS,
                                     seK.egger.XY.S = CD.Out.R$CDEgger.XtoY.se.NoS,
                                     K.ratio.YX.S = CD.Out.R$CDRatio.YtoX.est.NoS,
                                     seK.ratio.YX.S = CD.Out.R$CDRatio.YtoX.se.NoS,
                                     K.egger.YX.S = CD.Out.R$CDEgger.YtoX.est.NoS,
                                     seK.egger.YX.S = CD.Out.R$CDEgger.YtoX.se.NoS)
            bdcd.ratio.NoS.R <- CD.Cases.NoS.R$bdcd.ratio
            bdcd.egger.NoS.R <- CD.Cases.NoS.R$bdcd.egger
            
            
            if("cdRatio.NoS" %in% colnames.mat1){
                # decision rules - cd ratio bidircausal no screening
                if(bdcd.ratio.NoS.U==1 & bdcd.ratio.NoS.A==1 & bdcd.ratio.NoS.R==1){matR["slr1_res1_mlr1", "cdRatio.NoS"] <- matR["slr1_res1_mlr1", "cdRatio.NoS"]+1}
                if(bdcd.ratio.NoS.U==1 & bdcd.ratio.NoS.A==1 & bdcd.ratio.NoS.R==2){matR["slr1_res2_mlr1", "cdRatio.NoS"] <- matR["slr1_res2_mlr1", "cdRatio.NoS"]+1}
                if(bdcd.ratio.NoS.U==1 & bdcd.ratio.NoS.A==1 & bdcd.ratio.NoS.R==3){matR["slr1_res3_mlr1", "cdRatio.NoS"] <- matR["slr1_res3_mlr1", "cdRatio.NoS"]+1}
                
                if(bdcd.ratio.NoS.U==1 & bdcd.ratio.NoS.A==2 & bdcd.ratio.NoS.R==1){matR["slr1_res1_mlr2", "cdRatio.NoS"] <- matR["slr1_res1_mlr2", "cdRatio.NoS"]+1}
                if(bdcd.ratio.NoS.U==1 & bdcd.ratio.NoS.A==2 & bdcd.ratio.NoS.R==2){matR["slr1_res2_mlr2", "cdRatio.NoS"] <- matR["slr1_res2_mlr2", "cdRatio.NoS"]+1}
                if(bdcd.ratio.NoS.U==1 & bdcd.ratio.NoS.A==2 & bdcd.ratio.NoS.R==3){matR["slr1_res3_mlr2", "cdRatio.NoS"] <- matR["slr1_res3_mlr2", "cdRatio.NoS"]+1}
                
                if(bdcd.ratio.NoS.U==1 & bdcd.ratio.NoS.A==3 & bdcd.ratio.NoS.R==1){matR["slr1_res1_mlr3", "cdRatio.NoS"] <- matR["slr1_res1_mlr3", "cdRatio.NoS"]+1}
                if(bdcd.ratio.NoS.U==1 & bdcd.ratio.NoS.A==3 & bdcd.ratio.NoS.R==2){matR["slr1_res2_mlr3", "cdRatio.NoS"] <- matR["slr1_res2_mlr3", "cdRatio.NoS"]+1}
                if(bdcd.ratio.NoS.U==1 & bdcd.ratio.NoS.A==3 & bdcd.ratio.NoS.R==3){matR["slr1_res3_mlr3", "cdRatio.NoS"] <- matR["slr1_res3_mlr3", "cdRatio.NoS"]+1}
                
                if(bdcd.ratio.NoS.U==2 & bdcd.ratio.NoS.A==1 & bdcd.ratio.NoS.R==1){matR["slr2_res1_mlr1", "cdRatio.NoS"] <- matR["slr2_res1_mlr1", "cdRatio.NoS"]+1}
                if(bdcd.ratio.NoS.U==2 & bdcd.ratio.NoS.A==1 & bdcd.ratio.NoS.R==2){matR["slr2_res2_mlr1", "cdRatio.NoS"] <- matR["slr2_res2_mlr1", "cdRatio.NoS"]+1}
                if(bdcd.ratio.NoS.U==2 & bdcd.ratio.NoS.A==1 & bdcd.ratio.NoS.R==3){matR["slr2_res3_mlr1", "cdRatio.NoS"] <- matR["slr2_res3_mlr1", "cdRatio.NoS"]+1}
                
                if(bdcd.ratio.NoS.U==2 & bdcd.ratio.NoS.A==2 & bdcd.ratio.NoS.R==1){matR["slr2_res1_mlr2", "cdRatio.NoS"] <- matR["slr2_res1_mlr2", "cdRatio.NoS"]+1}
                if(bdcd.ratio.NoS.U==2 & bdcd.ratio.NoS.A==2 & bdcd.ratio.NoS.R==2){matR["slr2_res2_mlr2", "cdRatio.NoS"] <- matR["slr2_res2_mlr2", "cdRatio.NoS"]+1}
                if(bdcd.ratio.NoS.U==2 & bdcd.ratio.NoS.A==2 & bdcd.ratio.NoS.R==3){matR["slr2_res3_mlr2", "cdRatio.NoS"] <- matR["slr2_res3_mlr2", "cdRatio.NoS"]+1}
                
                if(bdcd.ratio.NoS.U==2 & bdcd.ratio.NoS.A==3 & bdcd.ratio.NoS.R==1){matR["slr2_res1_mlr3", "cdRatio.NoS"] <- matR["slr2_res1_mlr3", "cdRatio.NoS"]+1}
                if(bdcd.ratio.NoS.U==2 & bdcd.ratio.NoS.A==3 & bdcd.ratio.NoS.R==2){matR["slr2_res2_mlr3", "cdRatio.NoS"] <- matR["slr2_res2_mlr3", "cdRatio.NoS"]+1}
                if(bdcd.ratio.NoS.U==2 & bdcd.ratio.NoS.A==3 & bdcd.ratio.NoS.R==3){matR["slr2_res3_mlr3", "cdRatio.NoS"] <- matR["slr2_res3_mlr3", "cdRatio.NoS"]+1}
                
                if(bdcd.ratio.NoS.U==3 & bdcd.ratio.NoS.A==1 & bdcd.ratio.NoS.R==1){matR["slr3_res1_mlr1", "cdRatio.NoS"] <- matR["slr3_res1_mlr1", "cdRatio.NoS"]+1}
                if(bdcd.ratio.NoS.U==3 & bdcd.ratio.NoS.A==1 & bdcd.ratio.NoS.R==2){matR["slr3_res2_mlr1", "cdRatio.NoS"] <- matR["slr3_res2_mlr1", "cdRatio.NoS"]+1}
                if(bdcd.ratio.NoS.U==3 & bdcd.ratio.NoS.A==1 & bdcd.ratio.NoS.R==3){matR["slr3_res3_mlr1", "cdRatio.NoS"] <- matR["slr3_res3_mlr1", "cdRatio.NoS"]+1}
                
                if(bdcd.ratio.NoS.U==3 & bdcd.ratio.NoS.A==2 & bdcd.ratio.NoS.R==1){matR["slr3_res1_mlr2", "cdRatio.NoS"] <- matR["slr3_res1_mlr2", "cdRatio.NoS"]+1}
                if(bdcd.ratio.NoS.U==3 & bdcd.ratio.NoS.A==2 & bdcd.ratio.NoS.R==2){matR["slr3_res2_mlr2", "cdRatio.NoS"] <- matR["slr3_res2_mlr2", "cdRatio.NoS"]+1}
                if(bdcd.ratio.NoS.U==3 & bdcd.ratio.NoS.A==2 & bdcd.ratio.NoS.R==3){matR["slr3_res3_mlr2", "cdRatio.NoS"] <- matR["slr3_res3_mlr2", "cdRatio.NoS"]+1}
                
                if(bdcd.ratio.NoS.U==3 & bdcd.ratio.NoS.A==3 & bdcd.ratio.NoS.R==1){matR["slr3_res1_mlr3", "cdRatio.NoS"] <- matR["slr3_res1_mlr3", "cdRatio.NoS"]+1}
                if(bdcd.ratio.NoS.U==3 & bdcd.ratio.NoS.A==3 & bdcd.ratio.NoS.R==2){matR["slr3_res2_mlr3", "cdRatio.NoS"] <- matR["slr3_res2_mlr3", "cdRatio.NoS"]+1}
                if(bdcd.ratio.NoS.U==3 & bdcd.ratio.NoS.A==3 & bdcd.ratio.NoS.R==3){matR["slr3_res3_mlr3", "cdRatio.NoS"] <- matR["slr3_res3_mlr3", "cdRatio.NoS"]+1}
            }
            
            if("cdEgger.NoS" %in% colnames.mat1){
                # decision rules - cd egger bidircausal no screening
                if(bdcd.egger.NoS.U==1 & bdcd.egger.NoS.A==1 & bdcd.egger.NoS.R==1){matR["slr1_res1_mlr1", "cdEgger.NoS"] <- matR["slr1_res1_mlr1", "cdEgger.NoS"]+1}
                if(bdcd.egger.NoS.U==1 & bdcd.egger.NoS.A==1 & bdcd.egger.NoS.R==2){matR["slr1_res2_mlr1", "cdEgger.NoS"] <- matR["slr1_res2_mlr1", "cdEgger.NoS"]+1}
                if(bdcd.egger.NoS.U==1 & bdcd.egger.NoS.A==1 & bdcd.egger.NoS.R==3){matR["slr1_res3_mlr1", "cdEgger.NoS"] <- matR["slr1_res3_mlr1", "cdEgger.NoS"]+1}
                
                if(bdcd.egger.NoS.U==1 & bdcd.egger.NoS.A==2 & bdcd.egger.NoS.R==1){matR["slr1_res1_mlr2", "cdEgger.NoS"] <- matR["slr1_res1_mlr2", "cdEgger.NoS"]+1}
                if(bdcd.egger.NoS.U==1 & bdcd.egger.NoS.A==2 & bdcd.egger.NoS.R==2){matR["slr1_res2_mlr2", "cdEgger.NoS"] <- matR["slr1_res2_mlr2", "cdEgger.NoS"]+1}
                if(bdcd.egger.NoS.U==1 & bdcd.egger.NoS.A==2 & bdcd.egger.NoS.R==3){matR["slr1_res3_mlr2", "cdEgger.NoS"] <- matR["slr1_res3_mlr2", "cdEgger.NoS"]+1}
                
                if(bdcd.egger.NoS.U==1 & bdcd.egger.NoS.A==3 & bdcd.egger.NoS.R==1){matR["slr1_res1_mlr3", "cdEgger.NoS"] <- matR["slr1_res1_mlr3", "cdEgger.NoS"]+1}
                if(bdcd.egger.NoS.U==1 & bdcd.egger.NoS.A==3 & bdcd.egger.NoS.R==2){matR["slr1_res2_mlr3", "cdEgger.NoS"] <- matR["slr1_res2_mlr3", "cdEgger.NoS"]+1}
                if(bdcd.egger.NoS.U==1 & bdcd.egger.NoS.A==3 & bdcd.egger.NoS.R==3){matR["slr1_res3_mlr3", "cdEgger.NoS"] <- matR["slr1_res3_mlr3", "cdEgger.NoS"]+1}
                
                if(bdcd.egger.NoS.U==2 & bdcd.egger.NoS.A==1 & bdcd.egger.NoS.R==1){matR["slr2_res1_mlr1", "cdEgger.NoS"] <- matR["slr2_res1_mlr1", "cdEgger.NoS"]+1}
                if(bdcd.egger.NoS.U==2 & bdcd.egger.NoS.A==1 & bdcd.egger.NoS.R==2){matR["slr2_res2_mlr1", "cdEgger.NoS"] <- matR["slr2_res2_mlr1", "cdEgger.NoS"]+1}
                if(bdcd.egger.NoS.U==2 & bdcd.egger.NoS.A==1 & bdcd.egger.NoS.R==3){matR["slr2_res3_mlr1", "cdEgger.NoS"] <- matR["slr2_res3_mlr1", "cdEgger.NoS"]+1}
                
                if(bdcd.egger.NoS.U==2 & bdcd.egger.NoS.A==2 & bdcd.egger.NoS.R==1){matR["slr2_res1_mlr2", "cdEgger.NoS"] <- matR["slr2_res1_mlr2", "cdEgger.NoS"]+1}
                if(bdcd.egger.NoS.U==2 & bdcd.egger.NoS.A==2 & bdcd.egger.NoS.R==2){matR["slr2_res2_mlr2", "cdEgger.NoS"] <- matR["slr2_res2_mlr2", "cdEgger.NoS"]+1}
                if(bdcd.egger.NoS.U==2 & bdcd.egger.NoS.A==2 & bdcd.egger.NoS.R==3){matR["slr2_res3_mlr2", "cdEgger.NoS"] <- matR["slr2_res3_mlr2", "cdEgger.NoS"]+1}
                
                if(bdcd.egger.NoS.U==2 & bdcd.egger.NoS.A==3 & bdcd.egger.NoS.R==1){matR["slr2_res1_mlr3", "cdEgger.NoS"] <- matR["slr2_res1_mlr3", "cdEgger.NoS"]+1}
                if(bdcd.egger.NoS.U==2 & bdcd.egger.NoS.A==3 & bdcd.egger.NoS.R==2){matR["slr2_res2_mlr3", "cdEgger.NoS"] <- matR["slr2_res2_mlr3", "cdEgger.NoS"]+1}
                if(bdcd.egger.NoS.U==2 & bdcd.egger.NoS.A==3 & bdcd.egger.NoS.R==3){matR["slr2_res3_mlr3", "cdEgger.NoS"] <- matR["slr2_res3_mlr3", "cdEgger.NoS"]+1}
                
                if(bdcd.egger.NoS.U==3 & bdcd.egger.NoS.A==1 & bdcd.egger.NoS.R==1){matR["slr3_res1_mlr1", "cdEgger.NoS"] <- matR["slr3_res1_mlr1", "cdEgger.NoS"]+1}
                if(bdcd.egger.NoS.U==3 & bdcd.egger.NoS.A==1 & bdcd.egger.NoS.R==2){matR["slr3_res2_mlr1", "cdEgger.NoS"] <- matR["slr3_res2_mlr1", "cdEgger.NoS"]+1}
                if(bdcd.egger.NoS.U==3 & bdcd.egger.NoS.A==1 & bdcd.egger.NoS.R==3){matR["slr3_res3_mlr1", "cdEgger.NoS"] <- matR["slr3_res3_mlr1", "cdEgger.NoS"]+1}
                
                if(bdcd.egger.NoS.U==3 & bdcd.egger.NoS.A==2 & bdcd.egger.NoS.R==1){matR["slr3_res1_mlr2", "cdEgger.NoS"] <- matR["slr3_res1_mlr2", "cdEgger.NoS"]+1}
                if(bdcd.egger.NoS.U==3 & bdcd.egger.NoS.A==2 & bdcd.egger.NoS.R==2){matR["slr3_res2_mlr2", "cdEgger.NoS"] <- matR["slr3_res2_mlr2", "cdEgger.NoS"]+1}
                if(bdcd.egger.NoS.U==3 & bdcd.egger.NoS.A==2 & bdcd.egger.NoS.R==3){matR["slr3_res3_mlr2", "cdEgger.NoS"] <- matR["slr3_res3_mlr2", "cdEgger.NoS"]+1}
                
                if(bdcd.egger.NoS.U==3 & bdcd.egger.NoS.A==3 & bdcd.egger.NoS.R==1){matR["slr3_res1_mlr3", "cdEgger.NoS"] <- matR["slr3_res1_mlr3", "cdEgger.NoS"]+1}
                if(bdcd.egger.NoS.U==3 & bdcd.egger.NoS.A==3 & bdcd.egger.NoS.R==2){matR["slr3_res2_mlr3", "cdEgger.NoS"] <- matR["slr3_res2_mlr3", "cdEgger.NoS"]+1}
                if(bdcd.egger.NoS.U==3 & bdcd.egger.NoS.A==3 & bdcd.egger.NoS.R==3){matR["slr3_res3_mlr3", "cdEgger.NoS"] <- matR["slr3_res3_mlr3", "cdEgger.NoS"]+1}
            }
        }
      }

      
      
      if("CDcML" %in% colnames.mat1){
        ###################################
        # CDcML
        res.BDCDcML.u <- BiDirCDcML(b_X = bxU, b_Y = byU, se_X = sxU, se_Y = syU, n_X = n, n_Y = n, sig.cutoff = pFilterBDC, num_pert = 100)
        res.BDCDcML.a <- BiDirCDcML(b_X = bxA, b_Y = byA, se_X = sxA, se_Y = syA, n_X = n, n_Y = n, sig.cutoff = pFilterBDC, num_pert = 100)
        res.BDCDcML.r <- BiDirCDcML(b_X = bxR, b_Y = byR, se_X = sxR, se_Y = syR, n_X = n, n_Y = n, sig.cutoff = pFilterBDC, num_pert = 100)
        
        ############################################
        # Unadjusted - CDcCML
        ############################################
        u.S.DP<-NA
        
        # CIs - .S.DP
        u.ci.low.XY.S.DP <- res.BDCDcML.u$XtoY.est.S.DP - (qnorm(0.975)*res.BDCDcML.u$XtoY.se.S.DP)
        u.ci.upp.XY.S.DP <- res.BDCDcML.u$XtoY.est.S.DP + (qnorm(0.975)*res.BDCDcML.u$XtoY.se.S.DP)
        
        u.ci.low.YX.S.DP <- res.BDCDcML.u$YtoX.est.S.DP - (qnorm(0.975)*res.BDCDcML.u$YtoX.se.S.DP)
        u.ci.upp.YX.S.DP <- res.BDCDcML.u$YtoX.est.S.DP + (qnorm(0.975)*res.BDCDcML.u$YtoX.se.S.DP)
        
        # decision rules with CIs  --- .S.DP
        if((u.ci.low.XY.S.DP>0 & u.ci.upp.XY.S.DP<1)| (u.ci.low.XY.S.DP> (-1) & u.ci.upp.XY.S.DP<0)){
          if((u.ci.low.YX.S.DP>0 & u.ci.upp.YX.S.DP<1) | (u.ci.low.YX.S.DP> (-1) & u.ci.upp.YX.S.DP<0)){
            u.S.DP <- 3
          }else{
            u.S.DP <- 1
          }
        }else if((u.ci.low.YX.S.DP>0 & u.ci.upp.YX.S.DP<1) | (u.ci.low.YX.S.DP> (-1) & u.ci.upp.YX.S.DP<0)){
          if((u.ci.low.XY.S.DP>0 & u.ci.upp.XY.S.DP<1)| (u.ci.low.XY.S.DP> (-1) & u.ci.upp.XY.S.DP<0)){
            u.S.DP <- 3
          }else{
            u.S.DP <- 2
          }
        }else{
          u.S.DP <- 3
        }
        
        ############################################
        # Adjusted - CDcML
        ############################################
        a.S.DP<-NA
        
        # CIs - .S.DP
        a.ci.low.XY.S.DP <- res.BDCDcML.a$XtoY.est.S.DP - (qnorm(0.975)*res.BDCDcML.a$XtoY.se.S.DP)
        a.ci.upp.XY.S.DP <- res.BDCDcML.a$XtoY.est.S.DP + (qnorm(0.975)*res.BDCDcML.a$XtoY.se.S.DP)
        
        a.ci.low.YX.S.DP <- res.BDCDcML.a$YtoX.est.S.DP - (qnorm(0.975)*res.BDCDcML.a$YtoX.se.S.DP)
        a.ci.upp.YX.S.DP <- res.BDCDcML.a$YtoX.est.S.DP + (qnorm(0.975)*res.BDCDcML.a$YtoX.se.S.DP)
        
        # decision rules with CIs  --- .S.DP
        if((a.ci.low.XY.S.DP>0 & a.ci.upp.XY.S.DP<1)| (a.ci.low.XY.S.DP>-1 & a.ci.upp.XY.S.DP<0)){
          if((a.ci.low.YX.S.DP>0 & a.ci.upp.YX.S.DP<1) | (a.ci.low.YX.S.DP>-1 & a.ci.upp.YX.S.DP<0)){
            a.S.DP<- 3
          }else{
            a.S.DP<- 1
          }
        }else if((a.ci.low.YX.S.DP>0 & a.ci.upp.YX.S.DP<1) | (a.ci.low.YX.S.DP>-1 & a.ci.upp.YX.S.DP<0)){
          if((a.ci.low.XY.S.DP>0 & a.ci.upp.XY.S.DP<1)| (a.ci.low.XY.S.DP>-1 & a.ci.upp.XY.S.DP<0)){
            a.S.DP<- 3
          }else{
            a.S.DP<- 2
          }
        }else{
          a.S.DP<- 3
        }
        
        ############################################
        # Residual - CDcCML
        ############################################
        r.S.DP<-NA
        
        # CIs - .S.DP
        r.ci.low.XY.S.DP <- res.BDCDcML.r$XtoY.est.S.DP - (qnorm(0.975)*res.BDCDcML.r$XtoY.se.S.DP)
        r.ci.upp.XY.S.DP <- res.BDCDcML.r$XtoY.est.S.DP + (qnorm(0.975)*res.BDCDcML.r$XtoY.se.S.DP)
        
        r.ci.low.YX.S.DP <- res.BDCDcML.r$YtoX.est.S.DP - (qnorm(0.975)*res.BDCDcML.r$YtoX.se.S.DP)
        r.ci.upp.YX.S.DP <- res.BDCDcML.r$YtoX.est.S.DP + (qnorm(0.975)*res.BDCDcML.r$YtoX.se.S.DP)
        
        # decision rules with CIs  --- .S.DP
        if((r.ci.low.XY.S.DP>0 & r.ci.upp.XY.S.DP<1)| (r.ci.low.XY.S.DP> (-1) & r.ci.upp.XY.S.DP<0)){
          if((r.ci.low.YX.S.DP>0 & r.ci.upp.YX.S.DP<1) | (r.ci.low.YX.S.DP> (-1) & r.ci.upp.YX.S.DP<0)){
            r.S.DP <- 3
          }else{
            r.S.DP <- 1
          }
        }else if((r.ci.low.YX.S.DP>0 & r.ci.upp.YX.S.DP<1) | (r.ci.low.YX.S.DP> (-1) & r.ci.upp.YX.S.DP<0)){
          if((r.ci.low.XY.S.DP>0 & r.ci.upp.XY.S.DP<1)| (r.ci.low.XY.S.DP> (-1) & r.ci.upp.XY.S.DP<0)){
            r.S.DP <- 3
          }else{
            r.S.DP <- 2
          }
        }else{
          r.S.DP <- 3
        }
        
        
        if(u.S.DP==1 & a.S.DP==1 & r.S.DP==1){matR["slr1_res1_mlr1", "CDcML"] <- matR["slr1_res1_mlr1", "CDcML"]+1}
        if(u.S.DP==1 & a.S.DP==1 & r.S.DP==2){matR["slr1_res2_mlr1", "CDcML"] <- matR["slr1_res2_mlr1", "CDcML"]+1}
        if(u.S.DP==1 & a.S.DP==1 & r.S.DP==3){matR["slr1_res3_mlr1", "CDcML"] <- matR["slr1_res3_mlr1", "CDcML"]+1}
        
        if(u.S.DP==1 & a.S.DP==2 & r.S.DP==1){matR["slr1_res1_mlr2", "CDcML"] <- matR["slr1_res1_mlr2", "CDcML"]+1}
        if(u.S.DP==1 & a.S.DP==2 & r.S.DP==2){matR["slr1_res2_mlr2", "CDcML"] <- matR["slr1_res2_mlr2", "CDcML"]+1}
        if(u.S.DP==1 & a.S.DP==2 & r.S.DP==3){matR["slr1_res3_mlr2", "CDcML"] <- matR["slr1_res3_mlr2", "CDcML"]+1}
        
        if(u.S.DP==1 & a.S.DP==3 & r.S.DP==1){matR["slr1_res1_mlr3", "CDcML"] <- matR["slr1_res1_mlr3", "CDcML"]+1}
        if(u.S.DP==1 & a.S.DP==3 & r.S.DP==2){matR["slr1_res2_mlr3", "CDcML"] <- matR["slr1_res2_mlr3", "CDcML"]+1}
        if(u.S.DP==1 & a.S.DP==3 & r.S.DP==3){matR["slr1_res3_mlr3", "CDcML"] <- matR["slr1_res3_mlr3", "CDcML"]+1}
        
        if(u.S.DP==2 & a.S.DP==1 & r.S.DP==1){matR["slr2_res1_mlr1", "CDcML"] <- matR["slr2_res1_mlr1", "CDcML"]+1}
        if(u.S.DP==2 & a.S.DP==1 & r.S.DP==2){matR["slr2_res2_mlr1", "CDcML"] <- matR["slr2_res2_mlr1", "CDcML"]+1}
        if(u.S.DP==2 & a.S.DP==1 & r.S.DP==3){matR["slr2_res3_mlr1", "CDcML"] <- matR["slr2_res3_mlr1", "CDcML"]+1}
        
        if(u.S.DP==2 & a.S.DP==2 & r.S.DP==1){matR["slr2_res1_mlr2", "CDcML"] <- matR["slr2_res1_mlr2", "CDcML"]+1}
        if(u.S.DP==2 & a.S.DP==2 & r.S.DP==2){matR["slr2_res2_mlr2", "CDcML"] <- matR["slr2_res2_mlr2", "CDcML"]+1}
        if(u.S.DP==2 & a.S.DP==2 & r.S.DP==3){matR["slr2_res3_mlr2", "CDcML"] <- matR["slr2_res3_mlr2", "CDcML"]+1}
        
        if(u.S.DP==2 & a.S.DP==3 & r.S.DP==1){matR["slr2_res1_mlr3", "CDcML"] <- matR["slr2_res1_mlr3", "CDcML"]+1}
        if(u.S.DP==2 & a.S.DP==3 & r.S.DP==2){matR["slr2_res2_mlr3", "CDcML"] <- matR["slr2_res2_mlr3", "CDcML"]+1}
        if(u.S.DP==2 & a.S.DP==3 & r.S.DP==3){matR["slr2_res3_mlr3", "CDcML"] <- matR["slr2_res3_mlr3", "CDcML"]+1}
        
        if(u.S.DP==3 & a.S.DP==1 & r.S.DP==1){matR["slr3_res1_mlr1", "CDcML"] <- matR["slr3_res1_mlr1", "CDcML"]+1}
        if(u.S.DP==3 & a.S.DP==1 & r.S.DP==2){matR["slr3_res2_mlr1", "CDcML"] <- matR["slr3_res2_mlr1", "CDcML"]+1}
        if(u.S.DP==3 & a.S.DP==1 & r.S.DP==3){matR["slr3_res3_mlr1", "CDcML"] <- matR["slr3_res3_mlr1", "CDcML"]+1}
        
        if(u.S.DP==3 & a.S.DP==2 & r.S.DP==1){matR["slr3_res1_mlr2", "CDcML"] <- matR["slr3_res1_mlr2", "CDcML"]+1}
        if(u.S.DP==3 & a.S.DP==2 & r.S.DP==2){matR["slr3_res2_mlr2", "CDcML"] <- matR["slr3_res2_mlr2", "CDcML"]+1}
        if(u.S.DP==3 & a.S.DP==2 & r.S.DP==3){matR["slr3_res3_mlr2", "CDcML"] <- matR["slr3_res3_mlr2", "CDcML"]+1}
        
        if(u.S.DP==3 & a.S.DP==3 & r.S.DP==1){matR["slr3_res1_mlr3", "CDcML"] <- matR["slr3_res1_mlr3", "CDcML"]+1}
        if(u.S.DP==3 & a.S.DP==3 & r.S.DP==2){matR["slr3_res2_mlr3", "CDcML"] <- matR["slr3_res2_mlr3", "CDcML"]+1}
        if(u.S.DP==3 & a.S.DP==3 & r.S.DP==3){matR["slr3_res3_mlr3", "CDcML"] <- matR["slr3_res3_mlr3", "CDcML"]+1}
        
      }
      
      if("MRcML" %in% colnames.mat1){
        ##########################################
        # MRcML
        ##########################################
        res.BDMRcML.u <-BiDirMRcML(b_X = bxU,
                                   b_Y = byU,
                                   se_X = sxU,
                                   se_Y = syU,
                                   n_X = n,
                                   n_Y = n,
                                   sig.cutoff = pFilterBDC,
                                   num_pert = 100)
        res.BDMRcML.a <-BiDirMRcML(b_X = bxA,
                                   b_Y = byA,
                                   se_X = sxA,
                                   se_Y = syA,
                                   n_X = n,
                                   n_Y = n,
                                   sig.cutoff = pFilterBDC,
                                   num_pert = 100)
        res.BDMRcML.r <-BiDirMRcML(b_X = bxR,
                                   b_Y = byR,
                                   se_X = sxR,
                                   se_Y = syR,
                                   n_X = n,
                                   n_Y = n,
                                   sig.cutoff = pFilterBDC,
                                   num_pert = 100)
        
        ############################################
        # Unadjusted - MRcCML
        ############################################
        u.S.DP<-NA
        
        u.p.XY.BDMRcML.S.DP<- pnorm(-abs(res.BDMRcML.u$XtoY.est.S.DP / res.BDMRcML.u$XtoY.se.S.DP))*2
        u.p.YX.BDMRcML.S.DP<- pnorm(-abs(res.BDMRcML.u$YtoX.est.S.DP / res.BDMRcML.u$YtoX.se.S.DP))*2
        # decision rules
        if(u.p.XY.BDMRcML.S.DP <= 0.05 & u.p.YX.BDMRcML.S.DP > 0.05){
          u.S.DP<-1
        }#return case 1
        if(u.p.XY.BDMRcML.S.DP > 0.05 & u.p.YX.BDMRcML.S.DP <= 0.05){
          u.S.DP<-2
        }# return case 2
        if((u.p.XY.BDMRcML.S.DP >0.05 & u.p.YX.BDMRcML.S.DP > 0.05) |(u.p.XY.BDMRcML.S.DP <=0.05 & u.p.YX.BDMRcML.S.DP <= 0.05)){
          u.S.DP<-3
        }# return case 3
        
        
        
        ############################################
        # Adjusted - MRcCML
        ############################################
        a.S.DP<-NA
        
        a.p.XY.BDMRcML.S.DP<- pnorm(-abs(res.BDMRcML.a$XtoY.est.S.DP / res.BDMRcML.a$XtoY.se.S.DP))*2
        a.p.YX.BDMRcML.S.DP<- pnorm(-abs(res.BDMRcML.a$YtoX.est.S.DP / res.BDMRcML.a$YtoX.se.S.DP))*2
        # decision rules
        if(a.p.XY.BDMRcML.S.DP <= 0.05 & a.p.YX.BDMRcML.S.DP > 0.05){
          a.S.DP<-1
        }#return case 1
        if(a.p.XY.BDMRcML.S.DP > 0.05 & a.p.YX.BDMRcML.S.DP <= 0.05){
          a.S.DP<-2
        }# return case 2
        if((a.p.XY.BDMRcML.S.DP >0.05 & a.p.YX.BDMRcML.S.DP > 0.05) | (a.p.XY.BDMRcML.S.DP <=0.05 & a.p.YX.BDMRcML.S.DP <= 0.05)){
          a.S.DP<-3
        }# return case 3
        
        ############################################
        # residual - MRcCML
        ############################################
        r.S.DP<-NA
        
        r.p.XY.BDMRcML.S.DP<- pnorm(-abs(res.BDMRcML.r$XtoY.est.S.DP / res.BDMRcML.r$XtoY.se.S.DP))*2
        r.p.YX.BDMRcML.S.DP<- pnorm(-abs(res.BDMRcML.r$YtoX.est.S.DP / res.BDMRcML.r$YtoX.se.S.DP))*2
        # decision rules
        if(r.p.XY.BDMRcML.S.DP <= 0.05 & r.p.YX.BDMRcML.S.DP > 0.05){
          r.S.DP<-1
        }#return case 1
        if(r.p.XY.BDMRcML.S.DP > 0.05 & r.p.YX.BDMRcML.S.DP <= 0.05){
          r.S.DP<-2
        }# return case 2
        if((r.p.XY.BDMRcML.S.DP >0.05 & r.p.YX.BDMRcML.S.DP > 0.05) |(r.p.XY.BDMRcML.S.DP <=0.05 & r.p.YX.BDMRcML.S.DP <= 0.05)){
          r.S.DP<-3
        }# return case 3
        
        if(u.S.DP==1 & a.S.DP==1 & r.S.DP==1){matR["slr1_res1_mlr1", "MRcML"] <- matR["slr1_res1_mlr1", "MRcML"]+1}
        if(u.S.DP==1 & a.S.DP==1 & r.S.DP==2){matR["slr1_res2_mlr1", "MRcML"] <- matR["slr1_res2_mlr1", "MRcML"]+1}
        if(u.S.DP==1 & a.S.DP==1 & r.S.DP==3){matR["slr1_res3_mlr1", "MRcML"] <- matR["slr1_res3_mlr1", "MRcML"]+1}
        
        if(u.S.DP==1 & a.S.DP==2 & r.S.DP==1){matR["slr1_res1_mlr2", "MRcML"] <- matR["slr1_res1_mlr2", "MRcML"]+1}
        if(u.S.DP==1 & a.S.DP==2 & r.S.DP==2){matR["slr1_res2_mlr2", "MRcML"] <- matR["slr1_res2_mlr2", "MRcML"]+1}
        if(u.S.DP==1 & a.S.DP==2 & r.S.DP==3){matR["slr1_res3_mlr2", "MRcML"] <- matR["slr1_res3_mlr2", "MRcML"]+1}
        
        if(u.S.DP==1 & a.S.DP==3 & r.S.DP==1){matR["slr1_res1_mlr3", "MRcML"] <- matR["slr1_res1_mlr3", "MRcML"]+1}
        if(u.S.DP==1 & a.S.DP==3 & r.S.DP==2){matR["slr1_res2_mlr3", "MRcML"] <- matR["slr1_res2_mlr3", "MRcML"]+1}
        if(u.S.DP==1 & a.S.DP==3 & r.S.DP==3){matR["slr1_res3_mlr3", "MRcML"] <- matR["slr1_res3_mlr3", "MRcML"]+1}
        
        if(u.S.DP==2 & a.S.DP==1 & r.S.DP==1){matR["slr2_res1_mlr1", "MRcML"] <- matR["slr2_res1_mlr1", "MRcML"]+1}
        if(u.S.DP==2 & a.S.DP==1 & r.S.DP==2){matR["slr2_res2_mlr1", "MRcML"] <- matR["slr2_res2_mlr1", "MRcML"]+1}
        if(u.S.DP==2 & a.S.DP==1 & r.S.DP==3){matR["slr2_res3_mlr1", "MRcML"] <- matR["slr2_res3_mlr1", "MRcML"]+1}
        
        if(u.S.DP==2 & a.S.DP==2 & r.S.DP==1){matR["slr2_res1_mlr2", "MRcML"] <- matR["slr2_res1_mlr2", "MRcML"]+1}
        if(u.S.DP==2 & a.S.DP==2 & r.S.DP==2){matR["slr2_res2_mlr2", "MRcML"] <- matR["slr2_res2_mlr2", "MRcML"]+1}
        if(u.S.DP==2 & a.S.DP==2 & r.S.DP==3){matR["slr2_res3_mlr2", "MRcML"] <- matR["slr2_res3_mlr2", "MRcML"]+1}
        
        if(u.S.DP==2 & a.S.DP==3 & r.S.DP==1){matR["slr2_res1_mlr3", "MRcML"] <- matR["slr2_res1_mlr3", "MRcML"]+1}
        if(u.S.DP==2 & a.S.DP==3 & r.S.DP==2){matR["slr2_res2_mlr3", "MRcML"] <- matR["slr2_res2_mlr3", "MRcML"]+1}
        if(u.S.DP==2 & a.S.DP==3 & r.S.DP==3){matR["slr2_res3_mlr3", "MRcML"] <- matR["slr2_res3_mlr3", "MRcML"]+1}
        
        if(u.S.DP==3 & a.S.DP==1 & r.S.DP==1){matR["slr3_res1_mlr1", "MRcML"] <- matR["slr3_res1_mlr1", "MRcML"]+1}
        if(u.S.DP==3 & a.S.DP==1 & r.S.DP==2){matR["slr3_res2_mlr1", "MRcML"] <- matR["slr3_res2_mlr1", "MRcML"]+1}
        if(u.S.DP==3 & a.S.DP==1 & r.S.DP==3){matR["slr3_res3_mlr1", "MRcML"] <- matR["slr3_res3_mlr1", "MRcML"]+1}
        
        if(u.S.DP==3 & a.S.DP==2 & r.S.DP==1){matR["slr3_res1_mlr2", "MRcML"] <- matR["slr3_res1_mlr2", "MRcML"]+1}
        if(u.S.DP==3 & a.S.DP==2 & r.S.DP==2){matR["slr3_res2_mlr2", "MRcML"] <- matR["slr3_res2_mlr2", "MRcML"]+1}
        if(u.S.DP==3 & a.S.DP==2 & r.S.DP==3){matR["slr3_res3_mlr2", "MRcML"] <- matR["slr3_res3_mlr2", "MRcML"]+1}
        
        if(u.S.DP==3 & a.S.DP==3 & r.S.DP==1){matR["slr3_res1_mlr3", "MRcML"] <- matR["slr3_res1_mlr3", "MRcML"]+1}
        if(u.S.DP==3 & a.S.DP==3 & r.S.DP==2){matR["slr3_res2_mlr3", "MRcML"] <- matR["slr3_res2_mlr3", "MRcML"]+1}
        if(u.S.DP==3 & a.S.DP==3 & r.S.DP==3){matR["slr3_res3_mlr3", "MRcML"] <- matR["slr3_res3_mlr3", "MRcML"]+1}
      }
      
      ######################################
      # Start MR Steiger
      ######################################
      if(("MRS.IVW" %in% colnames.mat1) | ("MRS.wMedian" %in% colnames.mat1) | "MRS.Egger" %in% colnames.mat1){
        u.ivw <- NA
        a.ivw <- NA
        r.ivw <- NA
        
        u.wMed <- NA
        a.wMed <- NA
        r.wMed <- NA
        
        u.Egger <- NA
        a.Egger <- NA
        r.Egger <- NA
        
        #Sample sizes for p_exp, p_out
        n_exp<-n
        n_out<-n
        
        r_exp_SLR <- NULL 
        r_exp_MLR <- NULL
        r_exp_Res <- NULL
        
        r_out_SLR <- NULL 
        r_out_MLR <- NULL
        r_out_Res <- NULL
        # get vector of absolute correlations for snp-exp (r_exp) and snp-out (r_out)
        for(i in 1:ncol(matGX)){
          r_exp_SLR <- c(r_exp_SLR, abs( u.betaGX_X[i]/(sqrt( u.betaGX_X[i]^2+(n-2)*u.seGX_X[i]^2 )) ) )
          r_exp_MLR <- c(r_exp_MLR, abs( a.betaGX_X[i]/(sqrt( a.betaGX_X[i]^2+(n-2)*a.seGX_X[i]^2 )) ) )
          r_exp_Res <- c(r_exp_Res, abs( r.betaGX_X[i]/(sqrt( r.betaGX_X[i]^2+(n-2)*r.seGX_X[i]^2 )) ) )
          
          r_out_SLR <- c(r_out_SLR, abs( u.betaGX_Y[i]/(sqrt( u.betaGX_Y[i]^2+(n-2)*u.seGX_Y[i]^2 )) ) ) 
          r_out_MLR <- c(r_out_MLR, abs( a.betaGX_Y[i]/(sqrt( a.betaGX_Y[i]^2+(n-2)*a.seGX_Y[i]^2 )) ) )
          r_out_Res <- c(r_out_Res, abs( r.betaGX_Y[i]/(sqrt( r.betaGX_Y[i]^2+(n-2)*r.seGX_Y[i]^2 )) ) )
        }
        
        mrs.U <- mr_steiger(pxU, pyU, n_exp, n_out, r_exp_SLR, r_out_SLR, r_xxo = 1, r_yyo = 1)
        mrs.A <- mr_steiger(pxA, pyA, n_exp, n_out, r_exp_MLR, r_out_MLR, r_xxo = 1, r_yyo = 1)
        mrs.R <- mr_steiger(pxR, pyR, n_exp, n_out, r_exp_Res, r_out_Res, r_xxo = 1, r_yyo = 1)
          
        if("MRS.IVW" %in% colnames.mat1){
          u.tsmr.IVW.GX <- TwoSampleMR::mr_ivw(b_exp=u.betaGX_X, se_exp = u.seGX_X, b_out=u.betaGX_Y, se_out = u.seGX_Y, parameters = default_parameters())
          a.tsmr.IVW.GX <- TwoSampleMR::mr_ivw(b_exp=a.betaGX_X, se_exp = a.seGX_X, b_out=a.betaGX_Y, se_out = a.seGX_Y, parameters = default_parameters())
          r.tsmr.IVW.GX <- TwoSampleMR::mr_ivw(b_exp=r.betaGX_X, se_exp = r.seGX_X, b_out=r.betaGX_Y, se_out = r.seGX_Y, parameters = default_parameters())
          
          if(mrs.U$correct_causal_direction==TRUE & mrs.U$steiger_test <= sig.level & u.tsmr.IVW.GX$pval <= sig.level){
            u.ivw <- 1
            }else if(mrs.U$correct_causal_direction==FALSE & mrs.U$steiger_test <= sig.level & u.tsmr.IVW.GX$pval <= sig.level){
              u.ivw <- 2
            }else{u.ivw <- 3}
       
          if(mrs.A$correct_causal_direction==TRUE & mrs.A$steiger_test <= sig.level & a.tsmr.IVW.GX$pval <= sig.level){
            a.ivw <- 1
            }else if(mrs.A$correct_causal_direction==FALSE & mrs.A$steiger_test <= sig.level & a.tsmr.IVW.GX$pval <= sig.level){
            a.ivw <- 2
            }else{a.ivw <- 3}   
          
          if(mrs.R$correct_causal_direction==TRUE & mrs.R$steiger_test <= sig.level & r.tsmr.IVW.GX$pval <= sig.level){
            r.ivw <- 1
          }else if(mrs.R$correct_causal_direction==FALSE & mrs.R$steiger_test <= sig.level & r.tsmr.IVW.GX$pval <= sig.level){
            r.ivw <- 2
          }else{r.ivw <- 3}   
          
          if(u.ivw==1 & a.ivw==1 & r.ivw==1){matR["slr1_res1_mlr1", "MRS.IVW"] <- matR["slr1_res1_mlr1", "MRS.IVW"]+1}
          if(u.ivw==1 & a.ivw==1 & r.ivw==2){matR["slr1_res2_mlr1", "MRS.IVW"] <- matR["slr1_res2_mlr1", "MRS.IVW"]+1}
          if(u.ivw==1 & a.ivw==1 & r.ivw==3){matR["slr1_res3_mlr1", "MRS.IVW"] <- matR["slr1_res3_mlr1", "MRS.IVW"]+1}
          
          if(u.ivw==1 & a.ivw==2 & r.ivw==1){matR["slr1_res1_mlr2", "MRS.IVW"] <- matR["slr1_res1_mlr2", "MRS.IVW"]+1}
          if(u.ivw==1 & a.ivw==2 & r.ivw==2){matR["slr1_res2_mlr2", "MRS.IVW"] <- matR["slr1_res2_mlr2", "MRS.IVW"]+1}
          if(u.ivw==1 & a.ivw==2 & r.ivw==3){matR["slr1_res3_mlr2", "MRS.IVW"] <- matR["slr1_res3_mlr2", "MRS.IVW"]+1}
          
          if(u.ivw==1 & a.ivw==3 & r.ivw==1){matR["slr1_res1_mlr3", "MRS.IVW"] <- matR["slr1_res1_mlr3", "MRS.IVW"]+1}
          if(u.ivw==1 & a.ivw==3 & r.ivw==2){matR["slr1_res2_mlr3", "MRS.IVW"] <- matR["slr1_res2_mlr3", "MRS.IVW"]+1}
          if(u.ivw==1 & a.ivw==3 & r.ivw==3){matR["slr1_res3_mlr3", "MRS.IVW"] <- matR["slr1_res3_mlr3", "MRS.IVW"]+1}
          
          if(u.ivw==2 & a.ivw==1 & r.ivw==1){matR["slr2_res1_mlr1", "MRS.IVW"] <- matR["slr2_res1_mlr1", "MRS.IVW"]+1}
          if(u.ivw==2 & a.ivw==1 & r.ivw==2){matR["slr2_res2_mlr1", "MRS.IVW"] <- matR["slr2_res2_mlr1", "MRS.IVW"]+1}
          if(u.ivw==2 & a.ivw==1 & r.ivw==3){matR["slr2_res3_mlr1", "MRS.IVW"] <- matR["slr2_res3_mlr1", "MRS.IVW"]+1}
          
          if(u.ivw==2 & a.ivw==2 & r.ivw==1){matR["slr2_res1_mlr2", "MRS.IVW"] <- matR["slr2_res1_mlr2", "MRS.IVW"]+1}
          if(u.ivw==2 & a.ivw==2 & r.ivw==2){matR["slr2_res2_mlr2", "MRS.IVW"] <- matR["slr2_res2_mlr2", "MRS.IVW"]+1}
          if(u.ivw==2 & a.ivw==2 & r.ivw==3){matR["slr2_res3_mlr2", "MRS.IVW"] <- matR["slr2_res3_mlr2", "MRS.IVW"]+1}
          
          if(u.ivw==2 & a.ivw==3 & r.ivw==1){matR["slr2_res1_mlr3", "MRS.IVW"] <- matR["slr2_res1_mlr3", "MRS.IVW"]+1}
          if(u.ivw==2 & a.ivw==3 & r.ivw==2){matR["slr2_res2_mlr3", "MRS.IVW"] <- matR["slr2_res2_mlr3", "MRS.IVW"]+1}
          if(u.ivw==2 & a.ivw==3 & r.ivw==3){matR["slr2_res3_mlr3", "MRS.IVW"] <- matR["slr2_res3_mlr3", "MRS.IVW"]+1}
          
          if(u.ivw==3 & a.ivw==1 & r.ivw==1){matR["slr3_res1_mlr1", "MRS.IVW"] <- matR["slr3_res1_mlr1", "MRS.IVW"]+1}
          if(u.ivw==3 & a.ivw==1 & r.ivw==2){matR["slr3_res2_mlr1", "MRS.IVW"] <- matR["slr3_res2_mlr1", "MRS.IVW"]+1}
          if(u.ivw==3 & a.ivw==1 & r.ivw==3){matR["slr3_res3_mlr1", "MRS.IVW"] <- matR["slr3_res3_mlr1", "MRS.IVW"]+1}
          
          if(u.ivw==3 & a.ivw==2 & r.ivw==1){matR["slr3_res1_mlr2", "MRS.IVW"] <- matR["slr3_res1_mlr2", "MRS.IVW"]+1}
          if(u.ivw==3 & a.ivw==2 & r.ivw==2){matR["slr3_res2_mlr2", "MRS.IVW"] <- matR["slr3_res2_mlr2", "MRS.IVW"]+1}
          if(u.ivw==3 & a.ivw==2 & r.ivw==3){matR["slr3_res3_mlr2", "MRS.IVW"] <- matR["slr3_res3_mlr2", "MRS.IVW"]+1}
          
          if(u.ivw==3 & a.ivw==3 & r.ivw==1){matR["slr3_res1_mlr3", "MRS.IVW"] <- matR["slr3_res1_mlr3", "MRS.IVW"]+1}
          if(u.ivw==3 & a.ivw==3 & r.ivw==2){matR["slr3_res2_mlr3", "MRS.IVW"] <- matR["slr3_res2_mlr3", "MRS.IVW"]+1}
          if(u.ivw==3 & a.ivw==3 & r.ivw==3){matR["slr3_res3_mlr3", "MRS.IVW"] <- matR["slr3_res3_mlr3", "MRS.IVW"]+1}
        } 
        
        if("MRS.wMedian" %in% colnames.mat1){
          u.tsmr.wMed.GX <- TwoSampleMR::mr_weighted_median(b_exp=u.betaGX_X, se_exp = u.seGX_X, b_out=u.betaGX_Y, se_out = u.seGX_Y, parameters = default_parameters())
          a.tsmr.wMed.GX <- TwoSampleMR::mr_weighted_median(b_exp=a.betaGX_X, se_exp = a.seGX_X, b_out=a.betaGX_Y, se_out = a.seGX_Y, parameters = default_parameters())
          r.tsmr.wMed.GX <- TwoSampleMR::mr_weighted_median(b_exp=r.betaGX_X, se_exp = r.seGX_X, b_out=r.betaGX_Y, se_out = r.seGX_Y, parameters = default_parameters())
          
          if(mrs.U$correct_causal_direction==TRUE & mrs.U$steiger_test <= sig.level & u.tsmr.wMed.GX$pval <= sig.level){
            u.wMed <- 1
          }else if(mrs.U$correct_causal_direction==FALSE & mrs.U$steiger_test <= sig.level & u.tsmr.wMed.GX$pval <= sig.level){
            u.wMed <- 2
          }else{u.wMed <- 3}
          
          if(mrs.A$correct_causal_direction==TRUE & mrs.A$steiger_test <= sig.level & a.tsmr.wMed.GX$pval <= sig.level){
            a.wMed <- 1
          }else if(mrs.A$correct_causal_direction==FALSE & mrs.A$steiger_test <= sig.level & a.tsmr.wMed.GX$pval <= sig.level){
            a.wMed <- 2
          }else{a.wMed <- 3}   
          
          if(mrs.R$correct_causal_direction==TRUE & mrs.R$steiger_test <= sig.level & r.tsmr.wMed.GX$pval <= sig.level){
            r.wMed <- 1
          }else if(mrs.R$correct_causal_direction==FALSE & mrs.R$steiger_test <= sig.level & r.tsmr.wMed.GX$pval <= sig.level){
            r.wMed <- 2
          }else{r.wMed <- 3} 
          
          if(u.wMed==1 & a.wMed==1 & r.wMed==1){matR["slr1_res1_mlr1", "MRS.wMedian"] <- matR["slr1_res1_mlr1", "MRS.wMedian"]+1}
          if(u.wMed==1 & a.wMed==1 & r.wMed==2){matR["slr1_res2_mlr1", "MRS.wMedian"] <- matR["slr1_res2_mlr1", "MRS.wMedian"]+1}
          if(u.wMed==1 & a.wMed==1 & r.wMed==3){matR["slr1_res3_mlr1", "MRS.wMedian"] <- matR["slr1_res3_mlr1", "MRS.wMedian"]+1}
          
          if(u.wMed==1 & a.wMed==2 & r.wMed==1){matR["slr1_res1_mlr2", "MRS.wMedian"] <- matR["slr1_res1_mlr2", "MRS.wMedian"]+1}
          if(u.wMed==1 & a.wMed==2 & r.wMed==2){matR["slr1_res2_mlr2", "MRS.wMedian"] <- matR["slr1_res2_mlr2", "MRS.wMedian"]+1}
          if(u.wMed==1 & a.wMed==2 & r.wMed==3){matR["slr1_res3_mlr2", "MRS.wMedian"] <- matR["slr1_res3_mlr2", "MRS.wMedian"]+1}
          
          if(u.wMed==1 & a.wMed==3 & r.wMed==1){matR["slr1_res1_mlr3", "MRS.wMedian"] <- matR["slr1_res1_mlr3", "MRS.wMedian"]+1}
          if(u.wMed==1 & a.wMed==3 & r.wMed==2){matR["slr1_res2_mlr3", "MRS.wMedian"] <- matR["slr1_res2_mlr3", "MRS.wMedian"]+1}
          if(u.wMed==1 & a.wMed==3 & r.wMed==3){matR["slr1_res3_mlr3", "MRS.wMedian"] <- matR["slr1_res3_mlr3", "MRS.wMedian"]+1}
          
          if(u.wMed==2 & a.wMed==1 & r.wMed==1){matR["slr2_res1_mlr1", "MRS.wMedian"] <- matR["slr2_res1_mlr1", "MRS.wMedian"]+1}
          if(u.wMed==2 & a.wMed==1 & r.wMed==2){matR["slr2_res2_mlr1", "MRS.wMedian"] <- matR["slr2_res2_mlr1", "MRS.wMedian"]+1}
          if(u.wMed==2 & a.wMed==1 & r.wMed==3){matR["slr2_res3_mlr1", "MRS.wMedian"] <- matR["slr2_res3_mlr1", "MRS.wMedian"]+1}
          
          if(u.wMed==2 & a.wMed==2 & r.wMed==1){matR["slr2_res1_mlr2", "MRS.wMedian"] <- matR["slr2_res1_mlr2", "MRS.wMedian"]+1}
          if(u.wMed==2 & a.wMed==2 & r.wMed==2){matR["slr2_res2_mlr2", "MRS.wMedian"] <- matR["slr2_res2_mlr2", "MRS.wMedian"]+1}
          if(u.wMed==2 & a.wMed==2 & r.wMed==3){matR["slr2_res3_mlr2", "MRS.wMedian"] <- matR["slr2_res3_mlr2", "MRS.wMedian"]+1}
          
          if(u.wMed==2 & a.wMed==3 & r.wMed==1){matR["slr2_res1_mlr3", "MRS.wMedian"] <- matR["slr2_res1_mlr3", "MRS.wMedian"]+1}
          if(u.wMed==2 & a.wMed==3 & r.wMed==2){matR["slr2_res2_mlr3", "MRS.wMedian"] <- matR["slr2_res2_mlr3", "MRS.wMedian"]+1}
          if(u.wMed==2 & a.wMed==3 & r.wMed==3){matR["slr2_res3_mlr3", "MRS.wMedian"] <- matR["slr2_res3_mlr3", "MRS.wMedian"]+1}
          
          if(u.wMed==3 & a.wMed==1 & r.wMed==1){matR["slr3_res1_mlr1", "MRS.wMedian"] <- matR["slr3_res1_mlr1", "MRS.wMedian"]+1}
          if(u.wMed==3 & a.wMed==1 & r.wMed==2){matR["slr3_res2_mlr1", "MRS.wMedian"] <- matR["slr3_res2_mlr1", "MRS.wMedian"]+1}
          if(u.wMed==3 & a.wMed==1 & r.wMed==3){matR["slr3_res3_mlr1", "MRS.wMedian"] <- matR["slr3_res3_mlr1", "MRS.wMedian"]+1}
          
          if(u.wMed==3 & a.wMed==2 & r.wMed==1){matR["slr3_res1_mlr2", "MRS.wMedian"] <- matR["slr3_res1_mlr2", "MRS.wMedian"]+1}
          if(u.wMed==3 & a.wMed==2 & r.wMed==2){matR["slr3_res2_mlr2", "MRS.wMedian"] <- matR["slr3_res2_mlr2", "MRS.wMedian"]+1}
          if(u.wMed==3 & a.wMed==2 & r.wMed==3){matR["slr3_res3_mlr2", "MRS.wMedian"] <- matR["slr3_res3_mlr2", "MRS.wMedian"]+1}
          
          if(u.wMed==3 & a.wMed==3 & r.wMed==1){matR["slr3_res1_mlr3", "MRS.wMedian"] <- matR["slr3_res1_mlr3", "MRS.wMedian"]+1}
          if(u.wMed==3 & a.wMed==3 & r.wMed==2){matR["slr3_res2_mlr3", "MRS.wMedian"] <- matR["slr3_res2_mlr3", "MRS.wMedian"]+1}
          if(u.wMed==3 & a.wMed==3 & r.wMed==3){matR["slr3_res3_mlr3", "MRS.wMedian"] <- matR["slr3_res3_mlr3", "MRS.wMedian"]+1}
          
        }
        
        if("MRS.Egger" %in% colnames.mat1){
          u.tsmr.Egger.GX <- TwoSampleMR::mr_egger_regression(b_exp=u.betaGX_X, se_exp = u.seGX_X, b_out=u.betaGX_Y, se_out = u.seGX_Y, parameters = default_parameters())
          a.tsmr.Egger.GX <- TwoSampleMR::mr_egger_regression(b_exp=a.betaGX_X, se_exp = a.seGX_X, b_out=a.betaGX_Y, se_out = a.seGX_Y, parameters = default_parameters())
          r.tsmr.Egger.GX <- TwoSampleMR::mr_egger_regression(b_exp=r.betaGX_X, se_exp = r.seGX_X, b_out=r.betaGX_Y, se_out = r.seGX_Y, parameters = default_parameters())
          
          if(mrs.U$correct_causal_direction==TRUE & mrs.U$steiger_test <= sig.level & u.tsmr.Egger.GX$pval <= sig.level){
            u.Egger <- 1
          }else if(mrs.U$correct_causal_direction==FALSE & mrs.U$steiger_test <= sig.level & u.tsmr.Egger.GX$pval <= sig.level){
            u.Egger <- 2
          }else{u.Egger <- 3}
          
          if(mrs.A$correct_causal_direction==TRUE & mrs.A$steiger_test <= sig.level & a.tsmr.Egger.GX$pval <= sig.level){
            a.Egger <- 1
          }else if(mrs.A$correct_causal_direction==FALSE & mrs.A$steiger_test <= sig.level & a.tsmr.Egger.GX$pval <= sig.level){
            a.Egger <- 2
          }else{a.Egger <- 3}   
          
          if(mrs.R$correct_causal_direction==TRUE & mrs.R$steiger_test <= sig.level & r.tsmr.Egger.GX$pval <= sig.level){
            r.Egger <- 1
          }else if(mrs.R$correct_causal_direction==FALSE & mrs.R$steiger_test <= sig.level & r.tsmr.Egger.GX$pval <= sig.level){
            r.Egger <- 2
          }else{r.Egger <- 3} 
          
          if(u.Egger==1 & a.Egger==1 & r.Egger==1){matR["slr1_res1_mlr1", "MRS.Egger"] <- matR["slr1_res1_mlr1", "MRS.Egger"]+1}
          if(u.Egger==1 & a.Egger==1 & r.Egger==2){matR["slr1_res2_mlr1", "MRS.Egger"] <- matR["slr1_res2_mlr1", "MRS.Egger"]+1}
          if(u.Egger==1 & a.Egger==1 & r.Egger==3){matR["slr1_res3_mlr1", "MRS.Egger"] <- matR["slr1_res3_mlr1", "MRS.Egger"]+1}
          
          if(u.Egger==1 & a.Egger==2 & r.Egger==1){matR["slr1_res1_mlr2", "MRS.Egger"] <- matR["slr1_res1_mlr2", "MRS.Egger"]+1}
          if(u.Egger==1 & a.Egger==2 & r.Egger==2){matR["slr1_res2_mlr2", "MRS.Egger"] <- matR["slr1_res2_mlr2", "MRS.Egger"]+1}
          if(u.Egger==1 & a.Egger==2 & r.Egger==3){matR["slr1_res3_mlr2", "MRS.Egger"] <- matR["slr1_res3_mlr2", "MRS.Egger"]+1}
          
          if(u.Egger==1 & a.Egger==3 & r.Egger==1){matR["slr1_res1_mlr3", "MRS.Egger"] <- matR["slr1_res1_mlr3", "MRS.Egger"]+1}
          if(u.Egger==1 & a.Egger==3 & r.Egger==2){matR["slr1_res2_mlr3", "MRS.Egger"] <- matR["slr1_res2_mlr3", "MRS.Egger"]+1}
          if(u.Egger==1 & a.Egger==3 & r.Egger==3){matR["slr1_res3_mlr3", "MRS.Egger"] <- matR["slr1_res3_mlr3", "MRS.Egger"]+1}
          
          if(u.Egger==2 & a.Egger==1 & r.Egger==1){matR["slr2_res1_mlr1", "MRS.Egger"] <- matR["slr2_res1_mlr1", "MRS.Egger"]+1}
          if(u.Egger==2 & a.Egger==1 & r.Egger==2){matR["slr2_res2_mlr1", "MRS.Egger"] <- matR["slr2_res2_mlr1", "MRS.Egger"]+1}
          if(u.Egger==2 & a.Egger==1 & r.Egger==3){matR["slr2_res3_mlr1", "MRS.Egger"] <- matR["slr2_res3_mlr1", "MRS.Egger"]+1}
          
          if(u.Egger==2 & a.Egger==2 & r.Egger==1){matR["slr2_res1_mlr2", "MRS.Egger"] <- matR["slr2_res1_mlr2", "MRS.Egger"]+1}
          if(u.Egger==2 & a.Egger==2 & r.Egger==2){matR["slr2_res2_mlr2", "MRS.Egger"] <- matR["slr2_res2_mlr2", "MRS.Egger"]+1}
          if(u.Egger==2 & a.Egger==2 & r.Egger==3){matR["slr2_res3_mlr2", "MRS.Egger"] <- matR["slr2_res3_mlr2", "MRS.Egger"]+1}
          
          if(u.Egger==2 & a.Egger==3 & r.Egger==1){matR["slr2_res1_mlr3", "MRS.Egger"] <- matR["slr2_res1_mlr3", "MRS.Egger"]+1}
          if(u.Egger==2 & a.Egger==3 & r.Egger==2){matR["slr2_res2_mlr3", "MRS.Egger"] <- matR["slr2_res2_mlr3", "MRS.Egger"]+1}
          if(u.Egger==2 & a.Egger==3 & r.Egger==3){matR["slr2_res3_mlr3", "MRS.Egger"] <- matR["slr2_res3_mlr3", "MRS.Egger"]+1}
          
          if(u.Egger==3 & a.Egger==1 & r.Egger==1){matR["slr3_res1_mlr1", "MRS.Egger"] <- matR["slr3_res1_mlr1", "MRS.Egger"]+1}
          if(u.Egger==3 & a.Egger==1 & r.Egger==2){matR["slr3_res2_mlr1", "MRS.Egger"] <- matR["slr3_res2_mlr1", "MRS.Egger"]+1}
          if(u.Egger==3 & a.Egger==1 & r.Egger==3){matR["slr3_res3_mlr1", "MRS.Egger"] <- matR["slr3_res3_mlr1", "MRS.Egger"]+1}
          
          if(u.Egger==3 & a.Egger==2 & r.Egger==1){matR["slr3_res1_mlr2", "MRS.Egger"] <- matR["slr3_res1_mlr2", "MRS.Egger"]+1}
          if(u.Egger==3 & a.Egger==2 & r.Egger==2){matR["slr3_res2_mlr2", "MRS.Egger"] <- matR["slr3_res2_mlr2", "MRS.Egger"]+1}
          if(u.Egger==3 & a.Egger==2 & r.Egger==3){matR["slr3_res3_mlr2", "MRS.Egger"] <- matR["slr3_res3_mlr2", "MRS.Egger"]+1}
          
          if(u.Egger==3 & a.Egger==3 & r.Egger==1){matR["slr3_res1_mlr3", "MRS.Egger"] <- matR["slr3_res1_mlr3", "MRS.Egger"]+1}
          if(u.Egger==3 & a.Egger==3 & r.Egger==2){matR["slr3_res2_mlr3", "MRS.Egger"] <- matR["slr3_res2_mlr3", "MRS.Egger"]+1}
          if(u.Egger==3 & a.Egger==3 & r.Egger==3){matR["slr3_res3_mlr3", "MRS.Egger"] <- matR["slr3_res3_mlr3", "MRS.Egger"]+1}
          
          
        }
        
        
      }
      
      
    }# end sims loop
    
    
      # write tables to working directory
      write.table(matR, file=paste0(table.name,"_seed",SEED,"_matR.txt"), quote=F, row.names=T)
      write.table(corGXmr, file=paste0(table.name,"_seed",SEED,"_corGXmr.txt"), quote=F, row.names=F, col.names = F)
      write.table(corGYmr, file=paste0(table.name,"_seed",SEED,"_corGYmr.txt"), quote=F, row.names=F, col.names = F)
      write.table(corGXsr, file=paste0(table.name,"_seed",SEED,"_corGXsr.txt"), quote=F, row.names=F, col.names = F)
      write.table(corGYsr, file=paste0(table.name,"_seed",SEED,"_corGYsr.txt"), quote=F, row.names=F, col.names = F)
      write.table(corGXms, file=paste0(table.name,"_seed",SEED,"_corGXms.txt"), quote=F, row.names=F, col.names = F)
      write.table(corGYms, file=paste0(table.name,"_seed",SEED,"_corGYms.txt"), quote=F, row.names=F, col.names = F)
      
      
      if(plot.pdf){
          # create plot and save to working directory
          cols <- polychrome(n = 29)
          cols <- cols[c(1:4,7,5:6,8:19,21:27,20,28:29)]
          
          rownames1 <- rownames(matR)
          
          pdf(paste(table.name,"_seed",SEED,"_matR.pdf", sep = ""))
          par(mar=c(7, 3, 3, 7)+0.2, xpd=TRUE)
          barplot(as.matrix(matR),
                  col = cols,
                  las=2)
          legend("right", inset = c(-0.27,0), legend = rownames1,
                 fill = cols, box.lty = 0, cex = 0.7,xpd = T)
          dev.off()
      }
    
    
    return(list("matR" = matR, "corGXmr" = corGXmr, "corGYmr"=corGYmr, "corGXsr"=corGXsr,"corGYsr"=corGYsr, "corGXms"=corGXms, "corGYms"=corGYms))
    
  }
