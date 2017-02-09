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

# Function to create lines for trawl survey data

ptTOlines <- function(coordinates) {
    # Provide a vector with elements 1:4 = (lonBeg, lonEnd, latBeg, latEnd) and element 5 = ID
    coords <- matrix(nrow = 2, ncol = 2, data = as.numeric(coordinates[1:4]))
    line <- Line(coords)
    lines <- Lines(list(line), ID = coordinates[5])
    return(lines)
}

ptTOsplines <- function(matCoord, data, proj) {
    # Matrix of coordinates with columns 1:4 = (lonBeg, lonEnd, latBeg, latEnd) and column 5 = ID
    # data = the data frame of data corresponding to the lines
    lines <- apply(matCoord, 1, ptTOlines)
    lineSHP <- SpatialLines(lines, proj4string = CRS(proj))
    lineSHP <- SpatialLinesDataFrame(lineSHP, data, match.ID = TRUE)
    return(lineSHP)
}
