# -----------------------------------------------------------------------------
# PROJECT:
#    Evaluating the structure of the communities of the estuary
#    and gulf of St.Lawrence
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# REPOSITORY
#   EGSL species distribution
# -----------------------------------------------------------------------------

# /TODO: function to take data per species and make it per trawl survey (North)
# /TODO: function to take data per trawl and make it per species (South)
# /REVIEW: look at tidyr for this
# /REVIEW: review the problematic coordinates for fisheries data. I contacted S. Hurtubise on this subject.

# -----------------------------------------------------------------------------
# SCRIPT
#   Importing raw data in R
# -----------------------------------------------------------------------------

# Libraries
    library(sp)
    library(rgdal)
    library(rgeos)
    library(raster)
    library(magrittr)
    library(dplyr)

# Functions
    source('script/dmsTOdd.r')
    source('script/ptTOsplines.r')
    source('script/coordTOsppt.r')
    source('script/plotEGSL.r')

# Figure folder
    fig <- "../../../Wiki/docs/img/"

# RawData folders
    RawData <- "/Users/davidbeauchesne/Dropbox/PhD/PhD_RawData/"
    dataFolder1 <- paste(RawData, 'NordGolfe_RelevePluriSp_MPO/', sep = '') #Pluri-specific DFO surveys 2011-2015, North
    dataFolder2 <- paste(RawData, 'SudGolfe_RelevePluriSp_MPO/', sep = '') #Pluri-specific DFO surveys 2011-2015, South
    dataFolder3 <- paste(RawData, 'ObsMer_1999-2015_MPO/', sep = '') #Sea observers 1999-2015
    dataFolder4 <- paste(RawData, 'ZIF_Fisheries_2010-2015-MPO/', sep = '') #Fisheries logbooks 2010-2015
    dataFolder5 <- paste(RawData, 'Peche_sentinelle_2011-2015/', sep = '') #Peche sentinelles 2011-2015

# -----------------------------------------------------------------------------
# DATA: Pluri-specific DFO surveys, North
#
# INFORMATION:
    # Source ";"    /* Source de données, Alfred Needler=6, Lady Hammond=8, Teleost=16*/
    # No_Rel ";"    /* Numéro du relevé */
    # Nbpc ";"      /* Numéro de bateau de pêche commercial */
    # DatDeTow ";"  /* Date de début de trait Format aaaa-mm-jj */
    # DatFiTow ";"  /* Date de fin de trait Format aaaa-mm-jj */
    # No_Stn ";"    /* Numéro de la station */
    # EngGen ";"    /* Code d'engin de pêche utilisé */
    # Resul ";"     /* Code résultat opération */
    # EspGen ";"    /* Code d'espèce */
    # N_EspSci ";"  /* Nom scientifique de l'espèce*/
    # N_EspF ";"    /* Nom français de l'espèce */
    # HreDeb ";"    /* Heure début du trait, Format HH:MM:SS */
    # HreFin ";"    /* Heure fin de trait, Format HH:MM:SS */
    # LaDeTow ";"   /* Latitude positionnement début du trait, unité=ddmm.%% */
    # LoDeTow ";"   /* Longitude positionnement début du trait, unité=ddmm.%% */
    # LaFiTow ";"   /* Latitude positionnement fin du trait, unité=ddmm.%% */
    # LoFiTow ";"   /* Latitude début du trait, unité=ddmm.%% */
    # Prof_1 ";"    /* Profondeur début du trait, unité=m */
    # Prof_2 ";"    /* Profondeur fin de tait, unité=m */
    # WCapOri ";"   /* Poids de la capture, unité=kg */
    # SNb_Capt;     /* Nombre total individus capturés */
# -----------------------------------------------------------------------------

northPluri <- read.table(file = paste(dataFolder1, 'Te8-12.txt', sep = ''), header = TRUE, sep = ";", dec = '.')

# Transforming ddmm.%% to degree decimals
# See dmsTOdd() for more information
    northPluri[, "LaDeTow"] <- unlist(lapply(X = northPluri[, "LaDeTow"], FUN = dmsTOdd))
    northPluri[, "LoDeTow"] <- unlist(lapply(X = northPluri[, "LoDeTow"], FUN = dmsTOdd, type = 'long'))
    northPluri[, "LaFiTow"] <- unlist(lapply(X = northPluri[, "LaFiTow"], FUN = dmsTOdd))
    northPluri[, "LoFiTow"] <- unlist(lapply(X = northPluri[, "LoFiTow"], FUN = dmsTOdd, type = 'long'))
    ID <- seq(1,nrow(northPluri))
    northPluri <- cbind(ID, northPluri)

# Transforming data columns in correct type
    northPluri[,5] <- as.Date(northPluri[,5])
    northPluri[,6] <- as.Date(northPluri[,6])
    northPluri[,19] <- as.numeric(paste(northPluri[,19])) # coerces NAs in place of '. ', which is what I want
    northPluri[,20] <- as.numeric(paste(northPluri[,20])) # coerces NAs in place of '. ', which is what I want
    northPluri[,22] <- as.numeric(paste(northPluri[,22])) # coerces NAs in place of '.', which is what I want

# north plurispecific survey points and lines
    matCoordBeg <- northPluri[,c('LoDeTow','LaDeTow','ID')] # beginning trawl session
    matCoordEnd <- northPluri[,c('LoFiTow','LaFiTow','ID')] # end trawl session
    matCoordMid <- northPluri[,c('LoDeTow','LoFiTow','LaDeTow','LaFiTow','ID')] # middle trawl session
    matCoordLine <- northPluri[,c('LoDeTow','LoFiTow','LaDeTow','LaFiTow','ID')]
    proj <- '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'
    data <- northPluri

    spptmid <- coordTOsppt(matCoord = matCoordMid, data = northPluri, proj = proj, type = "2")
    # spptbeg <- coordTOsppt(matCoord = matCoordBeg, data = northPluri, proj = proj, type = "1")
    # spptend <- coordTOsppt(matCoord = matCoordEnd, data = northPluri, proj = proj, type = "1")
    # trawlLines <- ptTOsplines(matCoord = matCoordLine, data = northPluri, proj = proj)

# Export figures
    plotEGSL(data = spptmid, outputName = 'northpluri', outputFolder = fig, outputType = 'jpeg', exportFig = TRUE, size = 'lg')
    plotEGSL(data = spptmid, outputName = 'northpluri', outputFolder = fig, outputType = 'jpeg', exportFig = TRUE, size = 'sm', grid = FALSE)

# -----------------------------------------------------------------------------
# DATA: Pluri-specific DFO surveys, South
#
# INFORMATION:
#   All relevant information pertaining to species caught and their unique ID
#       is in the file paste(dataFolder2,'Gulf_Maritimes species codes.xls',sep='')
#
#   Data is already formatted in a way that allows me to create shapefile. However,
#       I will need to make usable at a species level, which it is not at the
#       moment.
#
#   This format allows me to see which species were NOT captured, which is not
#       possible with the format from the northern survey. I should make a function
#       to go back and forth between the two formats.
# -----------------------------------------------------------------------------

southPluri <- read.csv(file = paste(dataFolder2, 'DBeauchesne_MAR2016.csv', sep = ''), header = TRUE, sep = ",", dec = '.')
ID <- seq(1,nrow(southPluri))
southPluri <- cbind(ID, southPluri)

# Create shapefiles
    # south plurispecific survey points and lines
    matCoordBeg <- southPluri[,c('longitude_st','latitude_st','ID')] # beginning trawl session
    matCoordEnd <- southPluri[,c('longitude_end','latitude_end','ID')] # end trawl session
    matCoordMid <- southPluri[,c('longitude_st','longitude_end','latitude_st','latitude_end','ID')] # middle trawl session
    matCoordLine <- southPluri[,c('longitude_st','longitude_end','latitude_st','latitude_end','ID')]
    proj <- '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'
    data <- southPluri

    spptmid <- coordTOsppt(matCoord = matCoordMid, data = southPluri, proj = proj, type = "2")
    # spptbeg <- coordTOsppt(matCoord = matCoordBeg, data = southPluri, proj = proj, type = "1")
    # spptend <- coordTOsppt(matCoord = matCoordEnd, data = southPluri, proj = proj, type = "1")
    # trawlLines <- ptTOsplines(matCoord = matCoordLine, data = southPluri, proj = proj)

# Export figures
    plotEGSL(data = spptmid, outputName = 'southpluri', outputFolder = fig, outputType = 'jpeg', exportFig = TRUE, size = 'lg')
    plotEGSL(data = spptmid, outputName = 'southpluri', outputFolder = fig, outputType = 'jpeg', exportFig = TRUE, size = 'sm', grid = FALSE)

# -----------------------------------------------------------------------------
# DATA: Sea observer, data from 1999 to 2015
#
# INFORMATION:
#
# -----------------------------------------------------------------------------
fileNames <- dir(dataFolder3, pattern =".csv")
OBS <- read.csv(file = paste(dataFolder3, fileNames[1], sep = ''), header = TRUE, sep = ",", dec = '.')
OBS_engins <- read.csv(file = paste(dataFolder3, fileNames[2], sep = ''), header = TRUE, sep = ",", dec = '.')
OBS_species <- read.csv(file = paste(dataFolder3, fileNames[3], sep = ''), header = TRUE, sep = ",", dec = '.')

# Transform coordinates in degrees dedimal
    OBS[, "Latitude_deb"] <- unlist(lapply(X = OBS[, "Latitude_deb"], FUN = dmsTOdd2))
    OBS[, "Longitude_deb"] <- unlist(lapply(X = OBS[, "Longitude_deb"], FUN = dmsTOdd2, type = 'long'))
    OBS[, "Latitude_fin"] <- unlist(lapply(X = OBS[, "Latitude_fin"], FUN = dmsTOdd2))
    OBS[, "Longitude_fin"] <- unlist(lapply(X = OBS[, "Longitude_fin"], FUN = dmsTOdd2, type = 'long'))

# Add numeric ID
    ID <- seq(1,nrow(OBS))
    OBS <- cbind(ID, OBS)

# Create shapefile
    proj <- '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'
    matCoordBeg <- OBS[,c('Longitude_deb','Latitude_deb','ID')] # beginning trawl session
    spptbeg <- coordTOsppt(matCoord = matCoordBeg, data = OBS, proj = proj, type = "1")

# Export figures
    plotEGSL(data = spptbeg, outputName = 'seaobs', outputFolder = fig, outputType = 'jpeg', exportFig = TRUE, size = 'lg')
    plotEGSL(data = spptbeg, outputName = 'seaobs', outputFolder = fig, outputType = 'jpeg', exportFig = TRUE, size = 'sm', grid = FALSE)

# -----------------------------------------------------------------------------
# DATA: Fisheries logbook data (ZIF) 2010-2015
#
# INFORMATION:
# -----------------------------------------------------------------------------
fileNames <- dir(dataFolder4, pattern =".csv")
fisheries2010 <- read.csv(file = paste(dataFolder4, fileNames[1], sep = ''), header = TRUE, sep = ",", dec = '.')
fisheries2011 <- read.csv(file = paste(dataFolder4, fileNames[2], sep = ''), header = TRUE, sep = ",", dec = '.')
fisheries2012 <- read.csv(file = paste(dataFolder4, fileNames[3], sep = ''), header = TRUE, sep = ",", dec = '.')
fisheries2013 <- read.csv(file = paste(dataFolder4, fileNames[4], sep = ''), header = TRUE, sep = ",", dec = '.')
fisheries2014 <- read.csv(file = paste(dataFolder4, fileNames[5], sep = ''), header = TRUE, sep = ",", dec = '.')
fisheries2015 <- read.csv(file = paste(dataFolder4, fileNames[6], sep = ''), header = TRUE, sep = ",", dec = '.')
fisheries <- rbind(fisheries2010, fisheries2011, fisheries2012, fisheries2013, fisheries2014, fisheries2015)

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# Some data have coordinates that have an unknown format
# They only reprensent app. 0.01% of the dataset
# Until that is resolved, I will simply use the ones that are ok
problem <- unique(c(which(nchar(fisheries[,'Latitude']) != 6), which(nchar(fisheries[,'Longitude']) != 6)))
fisheries <- fisheries[-problem, ]
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# Transform coordinates
    fisheries[, "Latitude"] <- unlist(lapply(X = fisheries[, "Latitude"], FUN = dmsTOdd2))
    fisheries[, "Longitude"] <- unlist(lapply(X = fisheries[, "Longitude"], FUN = dmsTOdd2, type = 'long'))

# Individual numeric ID
    ID <- seq(1,nrow(fisheries))
    fisheries <- cbind(ID, fisheries)

# Creating shapefile
    proj <- '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'
    matCoordBeg <- fisheries[,c('Longitude','Latitude','ID')] # beginning trawl session
    spptbeg <- coordTOsppt(matCoord = matCoordBeg, data = fisheries, proj = proj, type = "1")

# Export figures
    plotEGSL(data = spptbeg, outputName = 'fish', outputFolder = fig, outputType = 'jpeg', exportFig = TRUE, size = 'lg')
    plotEGSL(data = spptbeg, outputName = 'fish', outputFolder = fig, outputType = 'jpeg', exportFig = TRUE, size = 'sm', grid = FALSE)

# -----------------------------------------------------------------------------
# DATA: Sentinel fisheries for fixed and mobile gear 2011-2015
#
# INFORMATION:
#
# /*************  Structure du fichier de données pour les pêches sentinelles avec engins fixes ********************************************************/
#
#     Source ";"     /* Source de données, Pêche sentinelle Qc=1, Pêche sentinelle TN=2*/
#     No_Rel ";"     /* Numéro du relevé */
#     NoStn ";"      /* Numéro de la station */
#     EngGen ";"     /* code de l'engin de pêche utilisé */
#     DatLevFx ";"   /* Date de l'activité de pêche */
#     EspGen ";"     /* Code d'espèce */
#     N_EspSci ";"   /* Nom scientifique de l'espèce*/
#     N_EspF ";"     /* Nom français de l'espèce */
#     LaLevFx ";"    /* Positionnement de l'engin de pêche, Latitude unité=ddmm.%% */
#     LoLevFX ";"    /* Positionnement de l'engin de pêche, Longitude unité=ddmm.%%  */
#     DpLevFx ";"    /* Profondeur, unité=m */
#     WCapOri ";"    /* Poids de la capture, unité=kg */
#     SNb_Capt;      /* Nombre individus capturés */
#
#
# /*************  Structure du fichier de données pour les pêches sentinelles avec engins mobiles ********************************************************/
#     Source ";"    /* Source de données, Alfred Needler=6, Lady Hammond=8, Teleost=16*/
#     No_Rel ";"    /* Numéro du relevé */
#     DatDeTow ";"  /* Date de début de trait Format aaaa-mm-jj */
#     DatFiTow ";"  /* Date de fin de trait Format aaaa-mm-jj */
#     NoStn ";"     /* Numéro de la station */
#     EngGen ";"    /* Code d'engin de pêche utilisé */
#     Resul ";"     /* Code résultat opération */
#     EspGen ";"    /* Code d'espèce */
#     N_EspSci ";"  /* Nom scientifique de l'espèce*/
#     N_EspF ";"    /* Nom français de l'espèce */
#     HreDeb ";"    /* Heure début du trait, Format HH:MM:SS */
#     HreFin ";"    /* Heure fin de trait, Format HH:MM:SS */
#     LaDeTow ";"   /* Latitude positionnement début du trait, unité=ddmm.%% */
#     LoDeTow ";"   /* Longitude positionnement début du trait, unité=ddmm.%% */
#     LaFiTow ";"   /* Latitude positionnement fin du trait, unité=ddmm.%% */
#     LoFiTow ";"   /* Latitude début du trait, unité=dd.%% */
#     Prof_1 ";"    /* Profondeur début du trait, unité=m */
#     Prof_2 ";"    /* Profondeur fin de tait, unité=m */
#     WCapOri ";"   /* Poids de la capture, unité=kg */
#     SNb_Capt;     /* Nombre individus capturés */
# # -----------------------------------------------------------------------------
sentFix <- read.table(file = paste(dataFolder5, 'PSF_QC-TN_Rel18-22.txt', sep = ''), header = TRUE, sep = ";", dec = '.')
sentMobQC <- read.table(file = paste(dataFolder5, 'PSM_QC_42-50.txt', sep = ''), header = TRUE, sep = ";", dec = '.')
sentMobTN <- read.table(file = paste(dataFolder5, 'PSM_TN_41-50.txt', sep = ''), header = TRUE, sep = ";", dec = '.')

# # Transforming data columns in correct type
    sentFix[, 13] <- as.numeric(paste(sentFix[,13]))
    sentMobTN[, 8] <- as.numeric(paste(sentMobTN[,8]))
    sentMobTN[, 15] <- as.numeric(paste(sentMobTN[,15]))
    sentMobTN[, 16] <- as.numeric(paste(sentMobTN[,16]))
    sentMobTN[, 19] <- as.numeric(paste(sentMobTN[,19]))
    sentMobTN[, 20] <- as.numeric(paste(sentMobTN[,20]))
    sentMobQC[, 5] <- as.numeric(paste(sentMobQC[,5]))
    sentMobQC[, 8] <- as.numeric(paste(sentMobQC[,8]))
    sentMobQC[, 19] <- as.numeric(paste(sentMobQC[,19]))
    sentMobQC[, 20] <- as.numeric(paste(sentMobQC[,20]))

# Transforming ddmm.%% to degree decimals
# See dmsTOdd() for more information
    sentFix[, "LaLevFx"] <- unlist(lapply(X = sentFix[, "LaLevFx"], FUN = dmsTOdd))
    sentFix[, "LoLevFX"] <- unlist(lapply(X = sentFix[, "LoLevFX"], FUN = dmsTOdd, type = 'long'))
    ID <- seq(1,nrow(sentFix))
    sentFix <- cbind(ID, sentFix)

    sentMobQC[, "LaDeTow"] <- unlist(lapply(X = sentMobQC[, "LaDeTow"], FUN = dmsTOdd))
    sentMobQC[, "LoDeTow"] <- unlist(lapply(X = sentMobQC[, "LoDeTow"], FUN = dmsTOdd, type = 'long'))
    sentMobQC[, "LaFiTow"] <- unlist(lapply(X = sentMobQC[, "LaFiTow"], FUN = dmsTOdd))
    sentMobQC[, "LoFiTow"] <- unlist(lapply(X = sentMobQC[, "LoFiTow"], FUN = dmsTOdd, type = 'long'))
    Region <- rep("QC",nrow(sentMobQC))
    sentMobQC <- cbind(Region, sentMobQC)

    sentMobTN[, "LaDeTow"] <- unlist(lapply(X = sentMobTN[, "LaDeTow"], FUN = dmsTOdd))
    sentMobTN[, "LoDeTow"] <- unlist(lapply(X = sentMobTN[, "LoDeTow"], FUN = dmsTOdd, type = 'long'))
    sentMobTN[, "LaFiTow"] <- unlist(lapply(X = sentMobTN[, "LaFiTow"], FUN = dmsTOdd))
    sentMobTN[, "LoFiTow"] <- unlist(lapply(X = sentMobTN[, "LoFiTow"], FUN = dmsTOdd, type = 'long'))
    Region <- rep("TN",nrow(sentMobTN))
    sentMobTN <- cbind(Region, sentMobTN)

    sentMob <- rbind(sentMobQC, sentMobTN)
    ID <- seq(1,nrow(sentMob))
    sentMob <- cbind(ID, sentMob)

# Create shapefiles
    proj <- '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'
    #Fixed
    matCoordBeg <- sentFix[,c('LoLevFX','LaLevFx','ID')]
    spptFix <- coordTOsppt(matCoord = matCoordBeg, data = sentFix, proj = proj)
    #Mobile
    matCoordBeg <- sentMob[,c('LoDeTow','LaDeTow','ID')]
    spptbegMob <- coordTOsppt(matCoord = matCoordBeg, data = sentMob, proj = proj)

# Export figures
    # Fixed
    plotEGSL(data = spptFix, outputName = 'sentfix', outputFolder = fig, outputType = 'jpeg', exportFig = TRUE, size = 'lg')
    plotEGSL(data = spptFix, outputName = 'sentfix', outputFolder = fig, outputType = 'jpeg', exportFig = TRUE, size = 'sm', grid = FALSE)

    # Mobile
    plotEGSL(data = spptbegMob, outputName = 'sentmob', outputFolder = fig, outputType = 'jpeg', exportFig = TRUE, size = 'lg')
    plotEGSL(data = spptbegMob, outputName = 'sentmob', outputFolder = fig, outputType = 'jpeg', exportFig = TRUE, size = 'sm', grid = FALSE)
