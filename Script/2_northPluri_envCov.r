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
library(sp)
library(rgdal)
library(rgeos)
library(reshape2)
# source('../../../PhD_RawData/Script/Function/coordTOsppt.r')
source('../../MEGA/MEGAsync/PhDData/Unclassified/Scripts/script/Function/coordTOsppt.r')


# Importing data required data from RawData
    northPluri <- readRDS('./RData/northPluriCor.rds')

# Extracting environmental data for northPluri stations
    # At some point I will need to decide whether the extractions will be at the beginning, the end, the middle, or all along the trawl path
    # For now, I will take the middle of the path
    # For the code for the line, beginning or end, see script to format rawdata in './PhD/PhD_RawData/'
    matCoordMid <- northPluri[,c('LoDeTow','LoFiTow','LaDeTow','LaFiTow','ID')] # middle trawl session
    proj <- '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'
    spptmid <- coordTOsppt(matCoord = matCoordMid, data = northPluri, proj = proj, type = "2")

# Import environmental data
    # epipelagic <- readOGR(dsn = "../../../PhD_RawData/data/Epipelagic_habitats_DFO", layer = "St_Lawrence_coastal_and_epipelagic_habitats")
    # benthic <- readOGR(dsn = "../../../PhD_RawData/data/Benthic_habitats_DFO", layer = "Quadrillage_10X10km_Marin")
    epipelagic <- readOGR(dsn = "../../MEGA/MEGAsync/PhDData/Unclassified/Epipelagic_habitats_DFO", layer = "St_Lawrence_coastal_and_epipelagic_habitats")
    benthic <- readOGR(dsn = "../../MEGA/MEGAsync/PhDData/Unclassified/Benthic_habitats_DFO", layer = "Quadrillage_10X10km_Marin")

# Transform point projection to environmental data projection
    spptmid <- spTransform(spptmid, CRSobj = CRS(proj4string(epipelagic)))

# Intersect points with environmental data
    epiStations <- over(spptmid, epipelagic)
    bentStations <- over(spptmid, benthic)

# Bind environmental covariables with northPluri data
    epiVar <- c('SSTWK30','SSAL_MEAN','FWRINFSR','STEMMEAN','BTEMMEAN','SANDBEACH','MUDFLAT','MARSH','SEDIMENT_C','HAB_C_E', 'TURBIDMEAN', 'SSAL_MIN', 'SSAL_MAX', 'BSAL_MEAN', 'BSAL_MIN', 'BSAL_MAX', 'STEMMIN', 'STEMMAX', 'BTEMMIN', 'BTEMMAX')
    bentVar <- c('Bathy_Mean','Geomorph_1','O2_Sat_Mea','SalMoyMoy','TempMoyMoy','SS_Code','SS_Desc_An','Megahabita','MHVar_3x3','SalMinMoy', 'SalMaxMoy', 'TempMinMoy', 'TempMaxMoy')
    Prof <- rowMeans(northPluri[, c('Prof_1','Prof_2')])
    northPluri <- cbind(northPluri, Prof, epiStations[, epiVar], bentStations[, bentVar])

# Check for correlation between environmental covariables
    northPluri_wide <- dcast(northPluri,
                             formula = No_Rel + No_Stn + DatDeTow + DatFiTow + HreDeb +
                                       HreFin + LaDeTow + LoDeTow + LaFiTow + LoFiTow +
                                       Prof_1 + Prof_2 + Prof + SSTWK30 + SSAL_MEAN +
                                       FWRINFSR + STEMMEAN + BTEMMEAN + SANDBEACH + MUDFLAT +
                                       MARSH + SEDIMENT_C + HAB_C_E + Bathy_Mean +
                                       Geomorph_1 + O2_Sat_Mea + SalMoyMoy + TempMoyMoy +
                                       SS_Code + Megahabita + MHVar_3x3 + SalMinMoy + SalMaxMoy +
                                       TempMinMoy + TempMaxMoy + TURBIDMEAN + SSAL_MIN + SSAL_MAX +
                                       BSAL_MEAN + BSAL_MIN + BSAL_MAX + STEMMIN + STEMMAX + BTEMMIN +
                                       BTEMMAX ~ EspGen,
                             value.var = 'EspGen',
                             fun.aggregate = length)

    corVar <- c('SSTWK30','SSAL_MEAN','FWRINFSR','STEMMEAN','BTEMMEAN','SANDBEACH','MUDFLAT','MARSH','Bathy_Mean','O2_Sat_Mea','SalMoyMoy','TempMoyMoy','Prof','SalMinMoy', 'SalMaxMoy', 'TempMinMoy', 'TempMaxMoy','TURBIDMEAN', 'SSAL_MIN', 'SSAL_MAX', 'BSAL_MEAN', 'BSAL_MIN', 'BSAL_MAX', 'STEMMIN', 'STEMMAX', 'BTEMMIN', 'BTEMMAX')

    round(cor(northPluri_wide[, corVar], use = 'na.or.complete'), 2)

# After correlation Check
    # Bind environmental covariables with northPluri data
        northPluri <- readRDS('./RData/northPluriCor.rds')
        epiVar <- c('SSAL_MEAN','STEMMEAN','BTEMMEAN','HAB_C_E','STEMMIN','BTEMMIN')
        bentVar <- c('Bathy_Mean','O2_Sat_Mea','SalMoyMoy','TempMoyMoy','Megahabita','MHVar_3x3')
        Prof <- rowMeans(northPluri[, c('Prof_1','Prof_2')])
        northPluri <- cbind(northPluri, Prof, epiStations[, epiVar], bentStations[, bentVar])

# Save analysis dataset (overwrite pre-existing data, no need to )
    saveRDS(northPluri, file = './RData/northPluriCor.rds')

# Visual of variables to expert plus 'Bathy_Mean'
    rbPal <- colorRampPalette(c('#bfceda','#0b416c'))
    fig <- "../../../Wiki/docs/img/"

    # Benthic habitats
    jpeg(paste(fig,'benthic.jpeg',sep=''), width = 8, height = 5.33, res = 100, units = 'in')
    par(mfrow = c(2,3), mar = c(0,0,0,0),pin = c(4,2),pty = "m",xaxs = "i",xaxt = 'n',xpd = FALSE,yaxs = "i",yaxt = 'n')
    for(i in 1:length(bentVar)) {
        data <- benthic@data[, bentVar[i]]
        if(class(data) == 'factor') {
            nCol <- length(unique(data))
            cols <- rbPal(nCol)[data]
        } else {
            cols <- rbPal(50)[as.numeric(cut(data, breaks = 50))]
        }
        plot(benthic, col = cols, border = cols)
        text(x = 350000, y = 925000, labels = paste(bentVar[i]))
    }
    dev.off()

    # Epipelagic habitats
    jpeg(paste(fig,'epipelagic.jpeg',sep=''), width = 8, height = 5.33, res = 100, units = 'in')
    par(mfrow = c(2,3), mar = c(0,0,0,0),pin = c(4,2),pty = "m",xaxs = "i",xaxt = 'n',xpd = FALSE,yaxs = "i",yaxt = 'n')
    for(i in 1:length(epiVar)) {
        data <- epipelagic@data[, epiVar[i]]
        if(class(data) == 'factor') {
            nCol <- length(unique(data))
            cols <- rbPal(nCol)[data]
        } else {
            cols <- rbPal(50)[as.numeric(cut(data, breaks = 50))]
        }
        plot(epipelagic, col = cols, border = cols)
        text(x = 350000, y = 925000, labels = paste(epiVar[i]))
    }
    dev.off()

# Intersect study grid with environmental data
    egsl_grid <- readOGR(dsn = "../../../PhD_obj0/Study_Area/RData/", layer = "egsl_grid")

    benthic@data <- benthic@data[, bentVar]
    epipelagic@data <- epipelagic@data[, epiVar]

    egsl_gridBent <- aggregate(x = benthic, by = egsl_grid, FUN = mean, areaWeighted = T)
    egsl_gridEpi <- aggregate(x = epipelagic, by = egsl_grid, FUN = mean, areaWeighted = T)

# Get cell centroid
    egslCentroid <- gCentroid(egsl_grid, byid = T)
    proj <- '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'
    egslCentroid <- spTransform(egslCentroid, CRSobj = CRS(proj))

# Save data only, to be used to produce predictions of taxa distribution in the St. Lawrence
    egsl_grid@data <- cbind(egsl_grid@data, egslCentroid@coords, egsl_gridBent@data, egsl_gridEpi@data)
    saveRDS(egsl_grid@data, file = './RData/egsl_grid.rds')

# # Visualise grid data
#     data <- egsl_grid@data[, 'O2_Sat_Mea']
#     cols <- rbPal(50)[as.numeric(cut(data, breaks = 50))]
#     plot(egsl_grid, col = cols, border = cols)
