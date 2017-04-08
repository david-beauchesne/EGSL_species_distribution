# Run init.r before other scripts
rm(list=ls())
 # for use in R console.
 # set own relevant directory if working in R console, otherwise ignore if in terminal
setwd("/Users/davidbeauchesne/Dropbox/PhD/PhD_obj2/Structure_Comm_EGSL/EGSL_species_distribution")
# -----------------------------------------------------------------------------
# PROJECT:
#    Evaluating the structure of the communities of the estuary
#    and gulf of St.Lawrence
# -----------------------------------------------------------------------------

# Studying the HMSC package, contact is G. Blanchet for inquiries

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

# First part of the workflow:
    # 1. Setting up an HMSC object by dening the model structure (e.g. whether traits or random effects are to be included) and organizing the data (typically imported from les with standard R commands).
    # 2. Dening the priors required for Bayesian inference.
    # 3. Initiating the model parameters
    # 4. Running the actual estimation scheme. This involves setting the Markov chain Monte Carlo (MCMC) sampler (e.g. by dening the number of iterations and how they will be thinned). As running the estimation scheme may take some time, we recommend the user to save the HMSC object to a le outside of R (e.g. called \model.RData"). This object includes the model structure, the data, and the full posterior distribution.

# Second part of the workflow:
    # No set order, but the user might typically wish to do the following:
        # - Assess the convergence of the MCMC chains,
        # - Produce posterior summaries (e.g. posterior means and quantiles) such as tables or plots,
        # - Measure the explanatory power of the model,
        # - Perform a variance partitioning among the fixed and random effects,
        # - Make predictions.

# Structure of the HMSC model
    # For ecologically meaningful analyses, the minimal set of data needed are either of the following:
        # 1. The occurrence matrix for one species as well as the environmental covariate matrix, in which case the HMSC model corresponds to a traditional single-species model.
        # 2. The occurrence matrix for several species, in which case the HMSC model corresponds to a model-based ordination.

# Formatting data
    # species: numeric values
    # environmental covariates: numeric values
    # random effects: factor if single random effect, data frame if multiple random effects
    # traits: numeric values
    # phylogeny: square symmetric correlation or covariance matrix
    # autocorrelated random effect: either a data frame with the first column a factor while the other columns are spatial or temporal coordinates or a list of data frame following the structure presented previously if there are multiple autocorrelated random effects.

    # Typically highly recommended to scale (center and divide by the standard deviation) the environmental covariates and the traits so that their mean is zero and their variance one.
        # this removes the potential eects units can have on the parameter estimation
        # the default priors are compatible with the scaled covariates

# Example 1
    # 1. Reading the data
        # Community matrix
        spComm <- read.csv("./documentation/HMSC-data/simulated/Y.csv")
        # Environmental covariates
        env <- read.csv("./documentation/HMSC-data/simulated/X.csv")
        # Random effects
        sitePlot <- read.csv("./documentation/HMSC-data/simulated/Pi.csv")

    # 2. Creating HMSC data
        # Convert all columns of Pi to a factor
        for(i in 1:ncol(sitePlot)) {
            sitePlot[,i] <- as.factor(sitePlot[,i])
        }

        simulEx1 <- as.HMSCdata(Y = spComm, X = env, Random = sitePlot, interceptX = FALSE, scaleX = FALSE)

        # Alternatively, load the data directly. Available for this example, but I will have to do it myself for my own data
        data("simulEx1")

    # 3. Defining prior distribution
        # We are using uninformative priors for this part. However, the user can modify this if desired.
        simulEx1prior <- as.HMSCprior(simulEx1)

    # 4. Setting initial model parameters
        simulEx1param <- as.HMSCparam(simulEx1, simulEx1prior)

    # 5. Setting the values of the true parameters"
        # This is only possible with simulated data for which we know the actual parameter values. Otherwise, it is impossible.
        data("simulParamEx1")

    # 6. Performing the MCMC sampling
        model <- hmsc(simulEx1,
                      param = simulEx1param,
                      priors = simulEx1prior,
                      family = "probit",
                      niter = 10000,
                      nburn = 1000,
                      thin = 10)

        # Simpler version when uninformative priors are sufficient and a priori parameters do not need to be set
        model <- hmsc(simulEx1,
                      family = "probit",
                      niter = 10000,
                      nburn = 1000,
                      thin = 10)

    # 7. Producing MCMC trace and density plots
        # Mixing objects
        mixingParamX <- as.mcmc(model, parameters = "paramX")
        mixingMeansParamX <- as.mcmc(model, parameters = "meansParamX")
        mixingMeansVarX <- as.mcmc(model, parameters = "varX")
        mixingParamLatent <- as.mcmc(model, parameters = "paramLatent")

        plot(mixingMeansParamX, col = "blue")

    # 8. Producing posterior summaries
        # Violin plot
            par(mar=c(6,4,1,1))
            mixingParamXDF <- as.data.frame(mixingParamX)
            beanplot(mixingParamXDF, las = 2)
            points(1:30, as.vector(simulParamEx1$param$paramX), pch=19, col="blue", cex=2)

        # Box plot
            par(mar=c(6,4,1,1))
            boxplot(mixingParamXDF, las = 2)
            points(1:30, as.vector(simulParamEx1$param$paramX), pch=19, col="blue", cex=2)

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

            # Print summary table
            paramXCITable

        # Credible intervals
            par(mar=c(7,4,1,1))
            plot(0, 0, xlim = c(1, nrow(paramXCITable)), ylim = range(paramXCITable), type = "n", xlab = "", ylab = "", main="paramX", xaxt="n")
            axis(1,1:30,rownames(paramXCITable),las=2)
            abline(h = 0,col = "grey")
            arrows(x0 = 1:nrow(paramXCITable), x1 = 1:nrow(paramXCITable), y0 = paramXCITable[, 2], y1 = paramXCITable[, 3], code = 3, angle = 90, length = 0.05)
            points(1:nrow(paramXCITable), paramXCITable[,1], pch = 15, cex = 2)
            points(1:nrow(paramXCITable), as.vector(simulParamEx1$param$paramX),col = "blue", pch = 19, cex = 2)

    # 9. Variance partitioning
        variationPart <- variPart(model, c(rep("climate", 2), "habitat"))

        Colour <- c("orange", "blue", "darkgreen", "purple")
        barplot(t(variationPart), col=Colour, las=1)
        legend("bottomleft",
                legend = c(paste("Fixed climate (mean = ", round(mean(variationPart[, 1]), 4)*100, "%)", sep=""),
                           paste("Fixed habitat (mean = ", round(mean(variationPart[, 2]), 4)*100, "%)", sep=""),
                           paste("Random site (mean = ", round(mean(variationPart[, 3]), 4)*100, "%)", sep=""),
                           paste("Random plot (mean = ", round(mean(variationPart[, 4]), 4)*100, "%)", sep="")),
               fill = Colour,
               bg="white")

    # 10. Association networks
        # Extract all estimated associatin matrix
            assoMat <- corRandomEff(model)
        # Average
            siteMean <- apply(assoMat[, , , 1], 1:2, mean)
            plotMean <- apply(assoMat[, , , 2], 1:2, mean)

        #=======================
        ### Associations to draw
        #=======================
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

        # plotDraw plots
            par(mfrow=c(1,2))
            # Matrix plot
            Colour <- colorRampPalette(c("blue", "white", "red"))(200)
            corrplot::corrplot(siteMean, method = "color", col = Colour, type = "lower", diag = FALSE, p.mat = siteDraw, tl.srt = 45)
            # Chord diagram
            circlize::chordDiagram(siteMean, symmetric = TRUE, annotationTrack = c("name", "grid"), grid.col = "grey", col = siteDrawCol)

        # siteDraw plots
            par(mfrow=c(1,2))
            # Matrix plot
            Colour <- colorRampPalette(c("blue", "white", "red"))(200)
            corrplot::corrplot(plotMean, method = "color", col = Colour, type = "lower", diag = FALSE, p.mat = plotDraw, tl.srt = 45)
            # Chord diagram
            circlize::chordDiagram(plotMean, symmetric = TRUE, annotationTrack = c("name", "grid"), grid.col = "grey", col = plotDrawCol)

    # 11. Computing the explanatory power of the model
        # Prevalence
            prevSp <- colSums(simulEx1$Y)
        # Coefficient of multiple determination
            R2 <- Rsquared(model, averageSp = FALSE)
            R2comm <- Rsquared(model, averageSp = TRUE)
        # Draw figure
            par(mar=c(5,6,0.5,0.5))
            plot(prevSp, R2, xlab = "Prevalence", ylab = expression(R^2), pch=19, las=1,cex.lab = 2)
            abline(h = R2comm, col = "blue", lwd = 2)

        # Extract all MCMC of paramX
            model$results$estimation$paramX

        ### Full joint probability distribution
            fullPost <- jposterior(model)

    # 12. Generating predictions for training data
        # predictions that are not conditional on the occurrences of other species
            predTrain <- predict(model)

    #13. Generating predictions for new data
        # New environmental covariates
            newPred <- 100
            nEnv <- ncol(simulEx1$X)
            Xnew <- matrix(nrow = newPred, ncol = nEnv)
            colnames(Xnew) <- colnames(simulEx1$X)
            Xnew[, 1] <- 1
            Xnew[, 2] <- mean(simulEx1$X[, 2])
            Xnew[, 3] <- seq(min(simulEx1$X[, 3]), max(simulEx1$X[, 3]), length = newPred)

        # New site- and plot-level random effect
            RandomSel <- sample(200, 100)
            RandomNew <- simulEx1$Random[RandomSel, ]
            for(i in 1:ncol(RandomNew)) {
                RandomNew[, i] <- as.factor(as.character(RandomNew[, i]))
            }
            colnames(RandomNew) <- colnames(simulEx1$Random)

        # Organize the data into an HMSCdata object
            dataVal <- as.HMSCdata(X = Xnew, Random = RandomNew, scaleX = FALSE, interceptX = FALSE)

        # Generate predictions
            predVal <- predict(model, newdata = dataVal)

        # Plot predictions
            plot(0, 0, type="n", xlim = range(dataVal$X[, 3]), ylim = c(0, 1),
                 xlab = "Environmental covariate 3",
                 ylab = "Probability of occurrence")
            Colours <- rainbow(ncol(predVal))
            for(i in 1:ncol(predVal)) {
                lines(dataVal$X[, 3], predVal[, i], col = Colours[i], lwd=3)
            }
