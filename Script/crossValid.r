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

crossValidation <- function(data, nCV, validPct, as.strata = TRUE, family = "probit", niter = 10000, nburn = 1000, thin = 10) {

    # ======================================
    # Monte Carlo cross-validation procedure
    # ======================================

        # Matrix to store cross-validation results
        modelAUC <- matrix(ncol = ncol(data$Y), nrow = nCV, data = 0, dimnames = list(paste('iter_', seq(1,nCV), sep = ''), colnames(data$Y)))

        for(cv in 1:nCV){
            # ---------------------------------------
            # Random sampling of dataset observations
            # ---------------------------------------
                # Index of validation observations randomly selected
                if(as.strata == TRUE) {
                    strata <- unique(data$Random[, 'random2'])
                    nStrata <- length(strata)
                    perStrata <- aggregate(data$Random, by = list(c(data$Random[, 'random2'])), FUN = length)[,2]
                    nValid <- round(validPct * perStrata)

                    randomSample <- character()
                    for(i in 1:nStrata) {
                        randomSample <- c(randomSample, as.character(sample(data$Random[data$Random[, 'random2'] == strata[i],'random1'], size = nValid[i], replace = FALSE)))
                    }
                } else {
                    randomSample <- as.character(sample(data$Random[, 'random1'], size = round(validPct * nrow(data$Random)), replace = FALSE))
                }

            # ----------------------------------------------------
            # Divide dataset into training and validation datasets
            # ----------------------------------------------------
                # Identify which rows should be set aside for cross validation
                    removeSample <- which(rownames(X) %in% randomSample)

                # Validation dataset
                    dataValid <- data
                    dataValid[[1]] <- dataValid[[1]][removeSample, ]
                    dataValid[[2]] <- dataValid[[2]][removeSample, ]
                    dataValid[[3]] <- dataValid[[3]][removeSample, ]
                    dataValid[[3]][, 1] <- as.factor(as.character(dataValid[[3]][, 1])) # Need the proper number of levels, build new index

                # Training dataset
                    # HMSC dataset minus sampling units set aside for cross validation
                    dataTrain <- data
                    dataTrain[[1]] <- dataTrain[[1]][-removeSample, ]
                    dataTrain[[2]] <- dataTrain[[2]][-removeSample, ]
                    dataTrain[[3]] <- dataTrain[[3]][-removeSample, ]
                    dataTrain[[3]][, 1] <- as.factor(as.character(dataTrain[[3]][, 1])) # Need the proper number of levels, build new index

            # ------------------------------------
            # Evaluate model with training dataset
            # ------------------------------------
                model <- HMSC::hmsc(dataTrain, family = family, niter = niter, nburn = nburn, thin = thin)

            # -----------------------------------
            # Predictions for validation datasets
            # -----------------------------------
                predVal <- predict(model, newdata = dataValid)

            #------------------------------------
            # Cross-validation: AUC of ROC curves
            #------------------------------------
                for(i in 1:ncol(predVal)) {
                    modelAUC[cv, i] <- ModelMetrics::auc(dataValid$Y[, i], predVal[, i])
                }
        }

        return(modelAUC)
}
