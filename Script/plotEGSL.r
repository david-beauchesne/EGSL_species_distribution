#' @name
#'
#' @title
#'
#' @description
#'
#' @param
#'
#' @return
#'
#' @author
#' David Beauchesne
#'
#' @references
#'
#' @importFrom
#'
#' @example
#'
#' @rdname
#'
#' @export

# Function to create points for trawl survey data

# /TODO: Add polygons
# /TODO: Add constant extent for the graphs
# /TODO: Add multiple data on graph

plotEGSL <- function(data,
                     outputName = '',
                     outputFolder = '',
                     outputType = 'jpeg',
                     exportFig = FALSE,
                     colFeature = '#707141',
                     outline = TRUE,
                     grid = TRUE,
                     size = 'lg',
                     type = 'points',
                     gSimp = 1000) {

    # function to plot EGSL with data overlaid on theme map
    #
    # Parameters:
    #   data            = feature(s) to add to theme map
    #   outputName      = name of file to generate
    #   outputFolder    = folder to export figure
    #   outputType      = type of figure output c('jpeg','pdf')
    #   exportFig       = TRUE or FALSE, export figure
    #   colFeature      = color of plotted features
    #   outline         = TRUE or FALSE, plot EGSL outline
    #   grid            = TRUE or FALSE, plot EGSL grid
    #   size            = size of figure output c('sm','lg')
    #   type            = type of data to overlay c('points','lines','polygons')
    #   gSimp           = simplify value for gSimplify()

    # Function:
    if(size == 'lg') {
        wd <- 10
        ht <- 8
        resol <- 300
    } else if(size == 'sm') {
        wd <- 3
        ht <- 2
        resol <- 75
    }

    if(exportFig == TRUE) {
        if(outputType == 'jpeg') {
            jpeg(paste(outputFolder,outputName,'_',size,'.',outputType,sep = ''),
                 width=wd,
                 height=ht,
                 res = resol,
                 units = 'in')
        } else if(outputType == 'pdf') {
            pdf(paste(outputFolder,outputName,'_',size,'.',outputType,sep = ''),
                 width=wd,
                 height=ht)
        }
    }

    par(mar = c(0,0,0,0),
        pin = c(4,2),
        pty = "m",
        xaxs = "i",
        xaxt = 'n',
        xpd = FALSE,
        yaxs = "i",
        yaxt = 'n')

    # PLOTS
        # plot EGSL grid
        if(grid == TRUE) {
            egsl_grid <- readOGR(dsn = "../../../PhD_obj0/Study_area/RData", layer = "egsl_grid")
            plot(egsl_grid, border = '#dddddd', col = '#A9C2E4')
        }

        # plot EGSL contour
        # add to plot if grid == TRUE
        # create new plot if grid == FALSE
        if(outline == TRUE & grid == TRUE) {
            egsl <- readOGR(dsn = "../../../PhD_RawData/EGSL", layer = "egsl")
            egsl <- gSimplify(egsl, gSimp, topologyPreserve=TRUE)
            plot(egsl, border = 'black', col = 'transparent', add = TRUE)
        } else if(outline == TRUE & grid == FALSE) {
            egsl <- readOGR(dsn = "../../../PhD_RawData/EGSL", layer = "egsl")
            egsl <- gSimplify(egsl, gSimp, topologyPreserve=TRUE)
            plot(egsl, border = 'black', col = 'transparent')
        }

        # Transform data to the right projection
        if(grid == TRUE) {
            data <- spTransform(data, CRSobj = CRS(proj4string(egsl_grid)))
        } else {
            data <- spTransform(data, CRSobj = CRS(proj4string(egsl)))
        }

        # Overlay data on EGSL theme map
        if(type == 'points') {
            points(data, pch = 21, col = 'transparent', bg = colFeature, cex = 0.5)
        } else if(type == 'lines') {
            lines(data, lwd = 2, col = colFeature)
        }

        if(exportFig == TRUE) {
            dev.off()
        }
} # plotEGSL
