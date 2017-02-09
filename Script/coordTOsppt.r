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

coordTOsppt <- function(matCoord, data, proj, type = '1') {
    # data = the data frame of data corresponding to the points
    # type = type of points to create
    #   '1' = beginning or end (user specified)
    #   '2' = middle of trawl session
    #   '3' = points every X distance (to do at some point)

    if(type == '1') {
        # Matrix of coordinates with columns 1:2 = (lon, lat) and  column 3 = ID
        sppt <- SpatialPointsDataFrame(coords = matCoord[,1:2], data = data, proj4string = CRS(proj), match.ID = TRUE)
    } else if(type == '2') {
        # Matrix of coordinates with columns 1:4 = (lonBeg, lonEnd, latBeg, latEnd) and column 5 = ID

        lon <- -(abs(matCoord[,1]) + abs(matCoord[,2]))/2
        lat <- (matCoord[,3] + matCoord[,4])/2
        sppt <- SpatialPointsDataFrame(coords = data.frame(lon, lat), data = data, proj4string = CRS(proj), match.ID = TRUE)
    }

    return(sppt)
} # coordTOsppt
