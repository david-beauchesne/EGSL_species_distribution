# Run init.r before other scripts
rm(list=ls())
 # for use in R console.
 # set own relevant directory if working in R console, otherwise ignore if in terminal
setwd("/Users/davidbeauchesne/Dropbox/PhD/PhD_obj2/Structure_Comm_EGSL/EGSL_species_distribution/")
# -----------------------------------------------------------------------------
# PROJECT:
#    Evaluating the structure of the communities of the estuary
#       and gulf of St.Lawrence
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# DETAILS:
#   The goal of this script is to run JSDMs from the HMSC package on data from
#       the annual northern Gulf DFO trawl survey
# -----------------------------------------------------------------------------
library(reshape2)
library(tidyr)
library(dplyr)
library(magrittr)
# Dependencies
    library(Rcpp)
    library(RcppArmadillo)
    library(coda)

# Installing the package
    # library(devtools)
    # install_github("guiblanchet/HMSC")
    library(HMSC)

# Other libraries
    # install.packages('beanplot', dependencies = T)
    # install.packages('corrplot', dependencies = T)
    # install.packages('circlize', dependencies = T)
    library(beanplot)
    library(corrplot)
    library(circlize)

# Functions
    source('./Script/crossValid.r')

# Importing data required data from RawData
    northPluri <- readRDS('./RData/northPluriCor.rds')
    fig <- "../../../Wiki/docs/img/"

# Basic parameters needed for the analysis
    # Name of environmental covariables
        # , # add at some point
        envCov <- c('Prof','SSAL_MEAN','SalMoyMoy','TempMoyMoy','STEMMEAN','BTEMMEAN','O2_Sat_Mea','LaDeTow','LoDeTow')

    # Groups of environmental covariables, for variance partitioning
        # length(envGroup) == length(envCov) + 1; because of Intercept
        envGroup <- c('Intercept','Bathymetry','Salinity','Salinity','Temperature','Temperature','Temperature','Oxygen', 'Spatial', 'Spatial')

# -------------------------------------------------------------
# Remove NAs from environmental variables selected for analysis
# -------------------------------------------------------------

    NAs <- northPluri[, envCov] %>%
            lapply(X = ., FUN = function(x) which(is.na(x))) %>%
            unlist(.) %>%
            unique(.)

    northPluri <- northPluri[-NAs, ]

# ========================================
# 1. Formatting the data for HMSC analysis
# ========================================

    # Species list
        sp <- northPluri[, c('EspGen','N_EspSci')] %>%
                .[!duplicated(.), ] %>%
                .[order(.[, 'EspGen']), ]

    # ---------------------------------
    # Y: sample units by species matrix
    # ---------------------------------

        # The data has to be formatted so that lines are trawl sessions and columns are species
        # Creating wide version of dataset to have stations as rows and species captured as columns
        # Only for presence absence
        # Other fields could also be used for count or weight
        northPluri_wide <- dcast(northPluri,
                                 formula = No_Rel + No_Stn + DatDeTow + DatFiTow + HreDeb + HreFin + LaDeTow + LoDeTow + LaFiTow + LoFiTow + Prof + SSAL_MEAN + STEMMEAN + BTEMMEAN + HAB_C_E + O2_Sat_Mea + SalMoyMoy + TempMoyMoy + Megahabita + MHVar_3x3 ~ EspGen,
                                 value.var = 'EspGen',
                                 fun.aggregate = length)

        # Extracting only presence/absence data
        # This corresponds to the Y matrix for the HMSC package
        # A unique ID for each station is created by combining survey number (i.e. year) and station number
            Station <- paste(northPluri_wide[, 'No_Rel'], '_', northPluri_wide[, 'No_Stn'], sep = '')
            Y <- northPluri_wide[, as.character(sp[,'EspGen'])]
            rownames(Y) <- Station

    # -----------------------------------------
    # Pi: sample units by random effects matrix
    # -----------------------------------------

        # Create a dataframe for random effects, columns have to be factors
        # Using survey nomber, which correspond to years, as a random effect in the analysis
        # Ultimately, there is likely a correlation between stations done during a single year in a single strata
        # It would be a good thing to analyze spatial dependence between stations
            Pi <- data.frame(sampling_unit = as.factor(Station),
                             survey_number = as.factor(northPluri_wide[, 'No_Rel']))

    # ----------------------------------------------------
    # X: sampling units by environmental covariates matrix
    # ----------------------------------------------------

        # Create a matrix for the values of environmental covariates at each sampling unit location
        # The values have to be numeric
        # Bathymetry used, extracted from northPluri rather than epipelagic and benthic habitat data

            X <- matrix(ncol = length(envCov),
                        nrow = nrow(northPluri_wide))
            X <- northPluri_wide[, envCov]
            rownames(X) <- Station

    # ----------------------------
    # as.HMSCdata for HMSC package
    # ----------------------------
        # Creating HMSC dataset for analyses
            northPluri_HMSC <- as.HMSCdata(Y = Y, X = X, Random = Pi, interceptX = T, scaleX = T)
            saveRDS(northPluri_HMSC, file = './RData/northPluri_HMSC.rds')

# ===============================
# 2. Performing the MCMC sampling
# ===============================

    # # Sampling of posterior distribution
    # # Simpler version when uninformative priors are sufficient and a priori parameters do not need to be set
    #     model <- hmsc(northPluri_HMSC,
    #                   family = "probit",
    #                   niter = 10000,
    #                   nburn = 1000,
    #                   thin = 10)

        model <- hmsc(northPluri_HMSC,
                      family = "probit",
                      niter = 100000,
                      nburn = 1000,
                      thin = 100)

    # save model
        saveRDS(model, file = './RData/modelHSMC.rds')

# =========================================
# 3. Producing MCMC trace and density plots
# =========================================

    # Mixing objects
        mixingParamX <- as.mcmc(model, parameters = "paramX")
        mixingMeansParamX <- as.mcmc(model, parameters = "meansParamX")
        mixingMeansVarX <- as.mcmc(model, parameters = "varX")
        mixingParamLatent <- as.mcmc(model, parameters = "paramLatent")

    # Save meanParamX for trace and density plots
    saveRDS(mixingMeansParamX, file = './RData/mixingMeansParamX.rds')

    # Trace and density plots to visually diagnose mcmc chains
    # Another way to check for convergence is to use diagnostic tests such as Geweke's convergence diagnostic (geweke.diag function in coda) and the Gelman and Rubin's convergence diagnostic (gelman.diag function in coda).
        paramModel <- colnames(mixingMeansParamX)
        nParam <- length(paramModel)

        jpeg(paste(fig,'MCMCTracePlot.jpeg',sep=''), width = 6, height = (1.5*nParam), res = 150, units = 'in')
        par(mfrow = c(nParam, 2), mar = rep(2, 4))
        for(i in 1:ncol(mixingMeansParamX)) {
          traceplot(mixingMeansParamX[,i], col = "blue", main = paste('Trace of ', paramModel[i]))
          densplot(mixingMeansParamX[,i], col = "orange", main = paste('Density of ', paramModel[i]))
        }
        dev.off()

# ================================
# 4. Producing posterior summaries
# ================================
    # # Violin plot
    #     par(mar=c(6,4,1,1))
    #     mixingParamXDF <- as.data.frame(mixingParamX)
    #     beanplot(mixingParamXDF, las = 2)
    #     points(1:30, as.vector(param$param$paramX), pch=19, col="blue", cex=2)
    #
    # # Box plot
    #     par(mar=c(6,4,1,1))
    #     boxplot(mixingParamXDF, las = 2)
    #     points(1:30, as.vector(param$param$paramX), pch=19, col="blue", cex=2)
    #
    # Average
        average <- apply(model$results$estimation$paramX, 1:2, mean)
    # 95% confidence intervals
        CI.025 <- apply(model$results$estimation$paramX, 1:2, quantile, probs = 0.025)
        CI.975 <- apply(model$results$estimation$paramX, 1:2, quantile, probs = 0.975)

    # Summary table
        paramXCITable <- cbind(unlist(as.data.frame(average)),
                             unlist(as.data.frame(CI.025)),
                             unlist(as.data.frame(CI.975)))
        colnames(paramXCITable) <- c("average", "lowerCI", "upperCI")
        rownames(paramXCITable) <- paste(rep(colnames(average), each = nrow(average)), "_", rep(rownames(average), ncol(average)), sep="")

      # Save summary table
        saveRDS(paramXCITable, file = './RData/modelPostSumm.rds')

    # Credible intervals
        paramXCITable_Full <- paramXCITable

        beg <- seq(1,nrow(paramXCITable_Full), by = 124)
        end <- seq(124,nrow(paramXCITable_Full), by = 124)
        sign <- abs(as.numeric(paramXCITable_Full[, 'lowerCI'] <= 0 & paramXCITable_Full[, 'upperCI'] >=0) - 3)
        sp <- northPluri[, c('EspGen','N_EspSci')] %>%
                .[!duplicated(.), ] %>%
                .[order(.[, 'EspGen']), ]

    # Export figure
        jpeg(paste(fig,'credibleInterval.jpeg',sep=''), width = 6, height = (1.5*nParam), res = 150, units = 'in')
        par(mfrow = c((length(beg)+1),1))
        for(i in 1:length(beg)) {
            paramXCITable <- paramXCITable_Full[beg[i]:end[i], ]
            cols <- sign[beg[i]:end[i]]
            par(mar=c(1,2,1,1))
            plot(0, 0, xlim = c(1, nrow(paramXCITable)), ylim = round(range(paramXCITable)), type = "n", xlab = "", ylab = "", main=paste(colnames(mixingMeansParamX)[i]), xaxt="n", bty = 'n')
            axis(1,1:124,las=2, labels = rep('',124))
            abline(h = 0,col = 'grey')
            arrows(x0 = 1:nrow(paramXCITable), x1 = 1:nrow(paramXCITable), y0 = paramXCITable[, 2], y1 = paramXCITable[, 3], code = 3, angle = 90, length = 0.05, col = cols)
            points(1:nrow(paramXCITable), paramXCITable[,1], pch = 15, cex = 1, col = cols)
        }
        mtext(text = sp[,'N_EspSci'], side = 1, line = 1, outer = FALSE, at = 1:124, col = 1, las = 2, cex = 0.4)
        dev.off()

# =========================
# 5.1 Variance partitioning
# =========================

    # for parameter names: colnames(mixingMeansParamX)
    nGroup <- length(unique(envGroup)) + 2
    variationPart <- variPart(model, envGroup)
    saveRDS(variationPart, file = './RData/variPart.rds')

    Colour <- rainbow(n = nGroup, s = 1, v = 1, start = 0, end = max(1, nGroup - 1)/nGroup, alpha = 1)

    jpeg(paste(fig,'variancePartitioning.jpeg',sep=''), width = 6, height = 4, res = 150, units = 'in')
        par(mar = c(6,3,1,1))
        barplot(t(variationPart), col=Colour, names.arg = sp[, 'N_EspSci'], las = 2, cex.names = 0.4, cex.axis = 0.6)

        # Create legend elements
            legendVector <- character(nGroup)
            variPartLabel <- c(unique(envGroup), 'Random site', 'Random plot')

            for(i in 1:nGroup) {
                legendVector[i] <- paste(variPartLabel[i], ' (mean = ', round(mean(variationPart[, i]), 4)*100, "%)", sep="")
            }

        legend('bottomleft', legend = legendVector, fill = Colour, bg = 'white', cex = 0.5)
    dev.off()

# ========================================================================
# 5.2 Variance partitioning for individual parameters (species diagnotics)
# ========================================================================

    # Extract variance partitioning per parameters for individual taxa diagnostics
        nGroup <- length(envCov) + 2
        variationPart <- variPart(model, c('Intercept',envCov))
        saveRDS(variationPart, file = './RData/variPartInd.rds')

# =======================
# 6. Association networks
# =======================
    # Extract all estimated associatin matrix
        assoMat <- corRandomEff(model)

    # Average
        siteMean <- apply(assoMat[, , , 1], 1:2, mean)
        plotMean <- apply(assoMat[, , , 2], 1:2, mean)
        colnames(plotMean) <- rownames(plotMean) <- sp[,'N_EspSci']

    #--------------------
    ### Site level effect
    #--------------------
    # Build matrix of colours for chordDiagram
        siteDrawCol <- matrix(NA, nrow = nrow(siteMean), ncol = ncol(siteMean))
        siteDrawCol[which(siteMean > 0.4, arr.ind=TRUE)]<-"red"
        siteDrawCol[which(siteMean < -0.4, arr.ind=TRUE)]<-"blue"
    # Build matrix of "significance" for corrplot
        siteDraw <- siteDrawCol
        siteDraw[which(!is.na(siteDraw), arr.ind = TRUE)] <- 0
        siteDraw[which(is.na(siteDraw), arr.ind = TRUE)] <- 1
        siteDraw <- matrix(as.numeric(siteDraw), nrow = nrow(siteMean), ncol = ncol(siteMean))
    #--------------------
    ### Plot level effect
    #--------------------
    # Build matrix of colours for chordDiagram
        plotDrawCol <- matrix(NA, nrow = nrow(plotMean), ncol = ncol(plotMean))
        plotDrawCol[which(plotMean > 0.4, arr.ind=TRUE)]<-"red"
        plotDrawCol[which(plotMean < -0.4, arr.ind=TRUE)]<-"blue"
    # Build matrix of "significance" for corrplot
        plotDraw <- plotDrawCol
        plotDraw[which(!is.na(plotDraw), arr.ind = TRUE)] <- 0
        plotDraw[which(is.na(plotDraw), arr.ind = TRUE)] <- 1
        plotDraw <- matrix(as.numeric(plotDraw), nrow = nrow(plotMean), ncol = ncol(plotMean))

    # # plotDraw plots
        # par(mfrow=c(1,2))
        # # Matrix plot
        # Colour <- colorRampPalette(c("blue", "white", "red"))(200)
        # corrplot::corrplot(siteMean, method = "color", col = Colour, type = "lower", diag = FALSE, p.mat = siteDraw, tl.srt = 45)
        # # Chord diagram
        # circlize::chordDiagram(siteMean, symmetric = TRUE, annotationTrack = c("name", "grid"), grid.col = "grey", col = siteDrawCol)

    # siteDraw plots
        jpeg(paste(fig,'correlationPlot.jpeg',sep=''), width = 7, height = 7, res = 150, units = 'in')
        # Matrix plot
        Colour <- colorRampPalette(c("blue", "white", "red"))(200)
        corrplot::corrplot(plotMean, method = "color", col = Colour, type = "lower", diag = FALSE, p.mat = plotDraw, tl.srt = 65, tl.cex = 0.4, pch.cex = 0.3, tl.col = 'black')
        dev.off()

        jpeg(paste(fig,'chordDiagram.jpeg',sep=''), width = 7, height = 7, res = 150, units = 'in')
        # Chord diagram
        circlize::chordDiagram(plotMean, symmetric = TRUE, annotationTrack = 'grid', grid.col = "grey", col = plotDrawCol, preAllocateTracks = 1)
        circlize::circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
          xlim = circlize::get.cell.meta.data("xlim")
          ylim = circlize::get.cell.meta.data("ylim")
          sector.name = circlize::get.cell.meta.data("sector.index")
          circlize::circos.text(mean(xlim), ylim[1] + .1, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = 0.4)
        #   circlize::circos.axis(h = "top", labels.cex = 0.5, major.tick.percentage = 0.2, sector.index = sector.name, track.index = 2)
        }, bg.border = NA)
        dev.off()

# ===============================================
# 7. Computing the explanatory power of the model
# ===============================================
    # Prevalence
        prevSp <- colSums(northPluri_HMSC$Y)
    # Coefficient of multiple determination
        R2 <- Rsquared(model, averageSp = FALSE)
        R2comm <- Rsquared(model, averageSp = TRUE)

    # Save R^2 calculation for individual summaries
        saveRDS(R2, file = './RData/modelR2.rds')
        saveRDS(R2comm, file = './RData/modelR2comm.rds')

    # Draw figure
        jpeg(paste(fig,'r2summaries.jpeg',sep=''), width = 6, height = 5, res = 150, units = 'in')
        plot(prevSp, R2, xlab = "Prevalence", ylab = expression(R^2), cex = 0.8, pch=19, las=1, cex.lab = 1, main = 'Explanatory power of the model')
        abline(h = R2comm, col = "blue", lwd = 2)
        dev.off()

    # Extract all MCMC of paramX
        # model$results$estimation$paramX

    ### Full joint probability distribution
        # fullPost <- jposterior(model)

# =============================================
# 8. Generating predictions for validation data
# =============================================

    modelAUC <- crossValidation(data = northPluri_HMSC,  nCV = 20, validPct = 0.2)
    saveRDS(modelAUC, file = './RData/modelAUC.rds')

    meanAUC <- colMeans(modelAUC)
    sdAUC <- apply(modelAUC, MARGIN = 2, FUN = sd)

    jpeg(paste(fig,'crossValidation.jpeg',sep=''), width = 6, height = 5, res = 150, units = 'in')
        plot(prevSp, meanAUC, ylim = c(0,1), xlab = "Prevalence", ylab = 'AUC', cex = 0.8, pch=19, las=1, cex.lab = 1, main = 'Monte Carlo cross-validation with AUC of ROC curves')
        abline(h = 0.5, col = 'grey')
        arrows(x0 = prevSp, x1 = prevSp, y0 = (meanAUC - sdAUC), y1 = (meanAUC + sdAUC), code = 3, angle = 90, length = 0.05)
        points(prevSp, meanAUC, pch = 19, cex = 0.8)
    dev.off()

# # =============================================
# # 9. Generating predictions from complete model
# # =============================================
#
#     modelPredictions <- predict(model)
#     saveRDS(modelPredictions, file = './RData/modelPredictions.rds')

# =========================================
# 10. Generating predictions study area grid
# =========================================

    # Load HMSC model
        model <- readRDS('./RData/modelHSMC.rds')

    # Importing data for study grid that I wish to use for species distribution predictions
        egsl <- readRDS('./RData/egsl_grid.rds')

    # Use 'Bathy_Mean' instead of 'Prof', as 'Prof' comes from northPluri data
        colnames(egsl)[which(colnames(egsl) == 'Bathy_Mean')] <- 'Prof'

    # 'Prof' is positive, while 'Bathy_Mean' is negative, uniformize for model
        egsl[, 'Prof'] <- -egsl[, 'Prof']

    # 'x' and 'y' in egsl grid data, change to match environmental variables used
        colN <- colnames(egsl)
        colnames(egsl)[colN == 'x'] <- 'LoDeTow'
        colnames(egsl)[colN == 'y'] <- 'LaDeTow'

    # Remove NAs
        envCov <- c('Prof','SSAL_MEAN','SalMoyMoy','TempMoyMoy','STEMMEAN','BTEMMEAN','O2_Sat_Mea','LaDeTow','LoDeTow')
        NAs <- egsl[, envCov] %>%
                lapply(X = ., FUN = function(x) which(is.na(x))) %>%
                unlist(.) %>%
                unique(.)

        egsl <- egsl[-NAs, ]

    # New environmental covariables matrix (X matrix)
        Xnew <- egsl[, envCov]

    # New site- and plot-level random effect (Pi matrix)
        PiNew <- data.frame(sampling_unit = as.factor(egsl[, 'ID']),
                         survey_number = as.factor(rep(1, nrow(egsl))))

    # Organize the data into an HMSCdata object
        dataEGSL <- as.HMSCdata(X = Xnew, Random = PiNew, scaleX = T, interceptX = T)

    # Generate predictions
        predEGSL <- predict(model, newdata = dataEGSL)

    # Replace NAs in grid
        data <- matrix(nrow = (length(NAs) + nrow(predEGSL)), ncol = ncol(predEGSL), data = 0)
        dimnames(data) = list(paste('ID',seq(1:nrow(data)), sep = ''), colnames(predEGSL))
        data[NAs, ] <- NA
        data[seq(1:nrow(data))[-NAs], ] <- predEGSL
        predEGSL <- data

    # Save predictions
        saveRDS(predEGSL, file = './RData/predEGSL.rds')
