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


# Continued from './Script/3_northPluri_JSDM.r'

# =========================================
# 9. Generating predictions study area grid
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
        dataVal <- as.HMSCdata(X = Xnew, Random = PiNew, scaleX = T, interceptX = T)

    # Generate predictions
        predVal <- predict(model, newdata = dataVal)

# ==================================================
# 10. Generating diagnostics information per species
# ==================================================

    # Plot distributions, create individual pdf per taxa
        library(sp)
        library(rgdal)

    # Load grid
        egsl_grid <- readOGR(dsn = "../../../PhD_obj0/Study_Area/RData/", layer = "egsl_grid")

    # Species list
        northPluri <- readRDS('./RData/northPluriCor.rds')
        sp <- northPluri[, c('EspGen','N_EspSci')] %>%
                .[!duplicated(.), ] %>%
                .[order(.[, 'EspGen']), ]

    # Graphs
        rbPal <- colorRampPalette(c('#2f6eb9','#2aadba','#b45f5f'))
        # Have to make sure that they are correctly matched by ID in grid.
            taxa <- '4'
            data <- numeric(length(NAs) + nrow(predVal))
            data[NAs] <- NA
            data[seq(1:length(data))[-NAs]] <- predVal[, taxa]

        # Plots
            cols <- rbPal(50)[as.numeric(cut(data, breaks = 50))]
            # plot(egsl_grid, col = cols, border = 'grey')
            plot(egsl_grid, col = cols, border = '#dbdbdb', main = paste(sp[sp == taxa, 'N_EspSci']))
