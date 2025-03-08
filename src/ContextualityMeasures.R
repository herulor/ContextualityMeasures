################################################################################
# # Measures of contextuality and noncontextuality
# #
# # Functions for computing the measures from
# # arXiv:1903.07170v5
# #
# # Víctor H Cervantes
# # 2019
################################################################################



# # Functions


#' Verify values taken by variables with same name in different contexts
#'
#' checkConnections(x) returns TRUE if the variables with the same name
#' have the same values in all prop.tables of x
#' @param x A list of prop.table
#' @return Logical.
#' @examples
checkConnections <- function (x) {

    namesProperties <- unique(unlist(lapply(x,
                                            function (y) names(dimnames(y)))))

    test <- TRUE

    for (iiName in namesProperties) {

        contextsIiName  <- sapply(x,
                                  function (y) iiName %in% names(dimnames(y)))

        dimensionIiName <- sapply(x[contextsIiName],
                                  function (y) which(names(dimnames(y)) %in% iiName))

        valuesIiName    <- sapply(lapply(x[contextsIiName],
                                         function (y) dimnames(y)),
                                  function(z) z[iiName])

        for (iiContext in seq(length(valuesIiName))) {
            if (iiContext < length(valuesIiName)) {
                for (jjContext in seq(from = iiContext + 1, to = length(valuesIiName))) {

                    test <- identical(sort(valuesIiName[[jjContext]]),
                                      sort(valuesIiName[[iiContext]]))

                    if (!test) {
                        message('Variable ', iiName, ' does not take the same values in all contexts.')
                        return(test)
                    }
                }
            }
        }

    }

    return(test)
}

#' Check if the object is of class RVSystem
#'
#' isRVSystem(x) returns TRUE if object has the appropriate structure
#'
#' @param x A list of prop.table
#' @return Logical
#' @examples
isRVSystem <- function (x) {

    test <- TRUE

    if (!is.list(x)) {
        message('x is not a list. An RVSystem must be a list of prop.tables with named dimensions.')
        test <- FALSE
    }

    xNames <- unlist(lapply(x, function(y) names(dimnames(y))))

    if (any(sapply(xNames, is.null))) {
        message('Some dimension of x has no name. An RVSystem must be a list of prop.tables with named dimensions')
        test <- FALSE
    }

    if (any(sapply(xNames, function (y) y == ''))) {
        message('Some dimension of x has no name. An RVSystem must be a list of prop.tables with named dimensions')
        test <- FALSE
    }

    if (!all(sapply(x, is.table))) {
        message('Some element of x is not a table. An RVSystem must be a list of prop.tables with named dimensions')
        test <- FALSE
    }

    if (!all(sapply(sapply(sapply(x, sum), function (y) all.equal(y, 1)), isTRUE))) {
        message('Some table of x is not a prop.table. An RVSystem must be a list of prop.tables with named dimensions')
        test <- FALSE
    }

    if (class(x) != 'RVSystem') {
        class(x) <- 'RVSystem'
    }

    test <- test & checkConnections(x)

    return(test)
}




#' Check if the object is of class DRVSystem
#'
#' isDRVSystem(x) returns TRUE if object has the appropriate structure
#'
#' @param x A list of prop.table with all tables of size 2x...x2
#' @return Logical
#' @examples
isDRVSystem <- function (x) {

    test <- isRVSystem(x)

    if (test) {
        if (!all(unlist(lapply(sapply(x, dim), function (y) y <= 2)))) {
            test <- FALSE
            message('x contains some prop.table for a variable with more than two values. A dichotomous RVSystem contains only dichotomous (and maybe deterministic) variables.')
        }
    }

    return(test)

}


#' A constructor function for the class RVSystem
#'
#' RVSystem(...) returns an object of class RVSystem if a set of prop.tables that
#' passes the check for RVSystem is given as the arguments
#'
#' @param ... one or more prop.tables to construct the system
#' @return An RVSystem if the set passes the test. NULL otherwise.
#' @examples
RVSystem <- function (x, ...) {

    if (!is.list(x)) {
        out <- list(x, ...)
    } else {
        out <- x
    }

    class(out) <- 'RVSystem'
    if (isRVSystem(out)) {
        if (is.null(names(out)) || any(names(out) == "")) {
            names(out) <- paste0('Ctx', seq(length(out)))
        }
        return(out)
    }
}




#' Construct the reduced Boolean Matrix for the linear programming
#' for computing the contextuality measures (CNT1, CNT2, and CNT3)
#' and the noncontextuality measure NCNT2
#'
#' cntMatrices(x) returns a list containing:
#'   cntMatrixL, cntMatrixB, cntMatrixC: boolean (sparse) matrices with the appropriate incidences
#'     of elements of the coupling vector (identified by the columns) to the vector of probabilities
#'     describing the system for the low-order marginals, bunch, and connection probabilities, respectively.
#'   lowMargins, bunchProbabilities, connectionProbabilities: the vectors of probabilities that describe the system.
#'
#'
#' @param x an RVSystem
#' @return A cntMatrices object
#' @examples
cntMatrices <- function (x, withMatrices = TRUE) {

    require('slam')
    require('abind')

    if (suppressMessages(isDRVSystem(x))) {

        if (is.null(names(x)) || any(names(x) == "")) {
            names(x) <- paste0('Ctx', seq(length(x)))
        }

        nContexts     <- length(x)
        nVariables    <- sum(sapply(x, function (y) length(dim(y))))

        nColumns <- 2 ^ nVariables

        namesProperties <- unique(unlist(lapply(x,
                                                function (y) names(dimnames(y)))))

        refValues        <- character(length = length(namesProperties))
        names(refValues) <- namesProperties

        for (iiName in namesProperties) {

            for (jjContext in names(x)) {

                if (iiName %in% names(dimnames(x[[jjContext]]))) {

                    refValues[iiName] <-  dimnames(x[[jjContext]])[[iiName]][2]

                    break
                }
            }
        }

        # Low-order marginals

        if (withMatrices) {
            cntMatrixL <- slam::simple_triplet_zero_matrix(nrow = nVariables, ncol = nColumns,
                                                           mode = 'integer')

            cntMatrixL           <- rbind(t(rep(1, nColumns)), cntMatrixL)
        } else {
            cntMatrixL <- matrix(0, nrow = nVariables + 1)
        }
        rownames(cntMatrixL) <- rep('', nVariables + 1)

        lowMargins        <- c(1, numeric(nVariables))
        names(lowMargins) <- rep('', nVariables + 1)

        iiRow <- 2

        for (jjContext in names(x)) {

            for (kkVariable in seq(dimnames(x[[jjContext]]))) {

                kkVariableName <- names(dimnames(x[[jjContext]]))[kkVariable]

                if (withMatrices) {
                    cntMatrixL[iiRow, ]         <- rep(c(1, 0), each = 2 ** (nVariables - iiRow + 1))
                }
                rownames(cntMatrixL)[iiRow] <- paste0(jjContext, '.', kkVariableName)

                lowMargins[iiRow]        <- apply(X = x[[jjContext]],
                                                  MARGIN = kkVariable,
                                                  FUN = sum)[refValues[kkVariableName]]

                names(lowMargins)[iiRow] <- paste0(jjContext, '.', kkVariableName)

                iiRow <- iiRow + 1

            }
        }

        # Bunch probabilities

        bunchesVariables    <- lapply(x, function (y) names(dimnames(y)))
        nBunchesVariables   <- sapply(x, function (y) length(dim(y)))
        nBunchProbabilities <- sum((2 ^ nBunchesVariables) - nBunchesVariables - 1)

        if (withMatrices) {
            cntMatrixB <- slam::simple_triplet_zero_matrix(nrow = nBunchProbabilities,
                                                           ncol = nColumns,
                                                           mode = 'integer')
        } else {
            cntMatrixB <- matrix(0, nrow = nBunchProbabilities)
        }

        rownames(cntMatrixB)      <- rep('', nBunchProbabilities)

        bunchProbabilities        <- integer(nBunchProbabilities)
        names(bunchProbabilities) <- rep('', nBunchProbabilities)

        iiRow <- 1
        for (jjContext in names(x)) {

            if (nBunchesVariables[jjContext] > 1) {
                for (kkMargins in 2:nBunchesVariables[jjContext]) {

                    mmMarginNames <- combn(x = seq_along(bunchesVariables[[jjContext]]),
                                           m = kkMargins)

                    for (kkkMargin in seq(ncol(mmMarginNames))) {

                        indexVariables  <- mmMarginNames[, kkkMargin]
                        namesVariables  <- bunchesVariables[[jjContext]][indexVariables]
                        valuesVariables <- refValues[namesVariables]

                        marginsTable <- apply(x[[jjContext]], indexVariables, sum)

                        bunchProbabilities[iiRow] <- abind::asub(marginsTable,
                                                                 as.list(valuesVariables))

                        if (withMatrices) {
                            cntMatrixB[iiRow, ] <- apply(cntMatrixL[paste0(jjContext, '.', namesVariables), ], 2, prod)
                        }

                        rownames(cntMatrixB)[iiRow]          <-
                            names(bunchProbabilities)[iiRow] <-
                            paste(c(jjContext, namesVariables), collapse = '.')

                        iiRow <- iiRow + 1
                    }
                }
            }

        }

        # Connection probabilities

        contextsNames <- names(x)

        connectionBunches   <- lapply(namesProperties,
                                      function (z)
                                          contextsNames[sapply(
                                              lapply(x,
                                                     function (y)
                                                         names(dimnames(y)) %in% z), any)]
        )
        names(connectionBunches) <- namesProperties

        nConnectionVariables     <- sapply(connectionBunches, length)
        nConnectionProbabilities <- sum(choose(nConnectionVariables, 2))

        if (withMatrices) {
            cntMatrixC <- slam::simple_triplet_zero_matrix(nrow = nConnectionProbabilities,
                                                           ncol = nColumns,
                                                           mode = 'integer')
        } else {
            cntMatrixC <- matrix(0, nrow = nConnectionProbabilities)
        }

        rownames(cntMatrixC)           <- rep('', nConnectionProbabilities)

        connectionProbabilities        <- integer(nConnectionProbabilities)
        names(connectionProbabilities) <- rep('', nConnectionProbabilities)

        iiRow <- 1
        for (jjProperty in namesProperties) {

            if (nConnectionVariables[jjProperty] > 1) {
                mmContextPairs <- combn(x = seq_along(connectionBunches[[jjProperty]]),
                                        m = 2)

                for (kkkPair in seq(ncol(mmContextPairs))) {

                    indexContexts   <- mmContextPairs[, kkkPair]
                    namesContexts   <- connectionBunches[[jjProperty]][indexContexts]
                    valuesVariables <- rep(refValues[jjProperty], length(indexContexts))

                    connectionProbabilities[iiRow] <- min(lowMargins[paste0(connectionBunches[[jjProperty]][mmContextPairs[, kkkPair]],
                                                                            '.', jjProperty)])

                    if (withMatrices) {
                        cntMatrixC[iiRow, ] <- apply(cntMatrixL[paste0(connectionBunches[[jjProperty]][mmContextPairs[, kkkPair]],
                                                                       '.', jjProperty), ], 2, prod)
                    }

                    rownames(cntMatrixC)[iiRow]               <-
                        names(connectionProbabilities)[iiRow] <-
                        paste(c(jjProperty, connectionBunches[[jjProperty]][mmContextPairs[, kkkPair]]), collapse = '.')

                    iiRow <- iiRow + 1
                }
            }
        }

        out <- list(cntMatrixL = cntMatrixL, cntMatrixB = cntMatrixB, cntMatrixC = cntMatrixC,
                    lowMargins = lowMargins, bunchProbabilities = bunchProbabilities,
                    connectionProbabilities = connectionProbabilities)

        class(out) <- 'cntMatrices'

        return(out)
    } else {
        stop('x is not a dichotomous RVSystem')
    }
}



#' Check if the object is of class cntMatrices
#'
#' isCntMatrices(x) returns TRUE if object has the appropriate structure
#'
#' @param x The object to be tested
#' @return Logical
#' @examples
isCntMatrices <- function (x) {

    test <- TRUE

    cntNames <- c('cntMatrixL', 'cntMatrixB', 'cntMatrixC',
                  'lowMargins', 'bunchProbabilities', 'connectionProbabilities')

    nameTest <- cntNames %in% names(x)
    if (!all(nameTest)) {
        notFound <- ifelse(sum(!nameTest) == 1, cntNames[!nameTest], paste(cntNames[!nameTest], collapse = ', '))
        message(paste(cntNames[!(cntNames %in% names(x))], collapse = ', '), ' not found as an element of x.')
        test <- FALSE
        return(test)
    }

    if (any(class(x[['cntMatrixL']]) == 'simple_triplet_matrix')) {
        require('slam')
    }

    if (any(class(x[['cntMatrixL']]) == 'Matrix')) {
        require('Matrix')
    }

    if (sum(diff(sapply(x[c('cntMatrixL', 'cntMatrixB', 'cntMatrixC')], ncol))) != 0) {
        message('The matrices cntMatrixL, cntMatrixB, and cntMatrixC must all have the same number of columns')
        test <- FALSE
    }

    if (sum(x[['cntMatrixL']] == 0, x[['cntMatrixL']] == 1) != prod(dim(x[['cntMatrixL']]))) {
        message('cntMatrixL must be a Boolean matrix')
        test <- FALSE
    }

    if (sum(x[['cntMatrixB']] == 0, x[['cntMatrixB']] == 1) != prod(dim(x[['cntMatrixB']]))) {
        message('cntMatrixB must be a Boolean matrix')
        test <- FALSE
    }

    if (sum(x[['cntMatrixC']] == 0, x[['cntMatrixC']] == 1) != prod(dim(x[['cntMatrixC']]))) {
        message('cntMatrixC must be a Boolean matrix')
        test <- FALSE
    }

    if (any(abs(x[['lowMargins']] - 0.5) > 0.5)) {
        message('lowMargins must contain probabilities (i.e. real numbers between 0 and 1)')
        test <- FALSE
    }

    if (any(abs(x[['bunchesProbabilities']] - 0.5) > 0.5)) {
        message('bunchesProbabilities must contain probabilities (i.e. real numbers between 0 and 1)')
        test <- FALSE
    }

    if (any(abs(x[['connectionProbabilities']] - 0.5) > 0.5)) {
        message('connectionProbabilities must contain probabilities (i.e. real numbers between 0 and 1)')
        test <- FALSE
    }

    if (nrow(x[['cntMatrixL']]) != length(x[['lowMargins']])) {
        message('The number of elements of lowMargins must be equal to the number of rows of cntMatrixL')
        test <- FALSE
    }

    if (nrow(x[['cntMatrixB']]) != length(x[['bunchProbabilities']])) {
        message('The number of elements of bunchProbabilities must be equal to the number of rows of cntMatrixB')
        test <- FALSE
    }

    if (nrow(x[['cntMatrixC']]) != length(x[['connectionProbabilities']])) {
        message('The number of elements of connectionProbabilities must be equal to the number of rows of cntMatrixC')
        test <- FALSE
    }

    if (any(rownames(x[['cntMatrixL']]) != names(x[['lowMargins']]))) {
        message('The elements of lowMargins and the rows of cntMatrixL must have the same names')
        test <- FALSE
    }

    if (any(rownames(x[['cntMatrixB']]) != names(x[['bunchProbabilities']]))) {
        message('The elements of bunchProbabilities and the rows of cntMatrixB must have the same names')
        test <- FALSE
    }

    if (any(rownames(x[['cntMatrixC']]) != names(x[['connectionProbabilities']]))) {
        message('The elements of connectionProbabilities and the rows of cntMatrixC must have the same names')
        test <- FALSE
    }

    return(test)
}



#' Compute the contextuality measure CNT1
#'
#' CNT1(x) returns a list containing:
#'     CNT1: The value of the CNT1 contextuality measure
#'     solution: The solution found for optimal coupling
#'     status: Whether the solution was successfully found in the linear programming task.
#'         It will return 0 for the optimal solution being found, and non-zero otherwise.
#'
#'
#' @param x a cntMatrices object
#' @return A Cnt.Meas object with the CNT1 measure
CNT1 <- function (x) {
    if (suppressMessages(isCntMatrices(x))) {
        require('Rglpk')

        if (class(x[['cntMatrixL']]) == 'simple_triplet_matrix') {
            obj <- slam::col_sums(x[['cntMatrixC']])
        } else {
            obj <- colSums(x[['cntMatrixC']])
        }

        mat <- rbind(x[['cntMatrixL']], x[['cntMatrixB']])
        rhs <- c(x[['lowMargins']], x[['bunchProbabilities']])
        dir <- rep('==', length(rhs))


        optimalCoupling <- Rglpk_solve_LP(obj = obj, mat = mat, dir = dir, rhs = rhs,
                                          max = TRUE)


        if (class(x[['cntMatrixL']]) == 'simple_triplet_matrix') {
            fittedLow <- t(slam::matprod_simple_triplet_matrix(x = x[['cntMatrixL']],
                                                             y = optimalCoupling[['solution']]))
            fittedCon <- t(slam::matprod_simple_triplet_matrix(x = x[['cntMatrixC']],
                                                             y = optimalCoupling[['solution']]))
            fittedBun <- t(slam::matprod_simple_triplet_matrix(x = x[['cntMatrixB']],
                                                             y = optimalCoupling[['solution']]))
        } else {
            fittedLow <- t(x[['cntMatrixL']] %*% optimalCoupling[['solution']])
            fittedCon <- t(x[['cntMatrixC']] %*% optimalCoupling[['solution']])
            fittedBun <- t(x[['cntMatrixB']] %*% optimalCoupling[['solution']])
        }

        CNT1 <- sum(abs(x[['connectionProbabilities']] - fittedCon))

        if (CNT1 < (200 * .Machine$double.eps)) {
            CNT1 <- 0
        }

        output <- list(CNT1 = CNT1, solution = optimalCoupling[['solution']],
                       status = optimalCoupling[['status']],
                       fittedLow = fittedLow, fittedCon = fittedCon, fittedBun = fittedBun)
        class(output) <- 'Cnt.Meas'

        if (output[['status']] != 0) {
            warning('Warning!. Linear programming failed to find an optimal solution')
        }

        return(output)
    } else {
        stop('x is not a cntMatrices object')
    }
}



#' Compute the contextuality measure CNT2
#'
#' CNT2(x) returns a list containing:
#'     CNT2: The value of the CNT2 contextuality measure
#'     solution: The solution found for optimal coupling
#'     status: Whether the solution was successfully found in the linear programming task.
#'         It will return 0 for the optimal solution being found, and non-zero otherwise.
#'
#'
#' @param x a cntMatrices object
#' @return A Cnt.Meas object with the CNT2 measure
CNT2 <- function (x) {
    if (suppressMessages(isCntMatrices(x))) {
        require('Rglpk')


        minusLimit <- cbind(x[['cntMatrixB']],
                            -diag(nrow = nrow(x[['cntMatrixB']])))

        plusLimit  <- cbind(x[['cntMatrixB']],
                            diag(nrow = nrow(x[['cntMatrixB']])))

        x[['cntMatrixL']] <- cbind(x[['cntMatrixL']],
                                   matrix(rep(0, nrow(plusLimit) * nrow(x[['cntMatrixL']])),
                                          nrow = nrow(x[['cntMatrixL']])))

        x[['cntMatrixC']] <- cbind(x[['cntMatrixC']],
                                   matrix(rep(0, nrow(plusLimit) * nrow(x[['cntMatrixC']])),
                                          nrow = nrow(x[['cntMatrixC']])))

        obj <- c(rep(0, ncol(x[['cntMatrixB']])), rep(1, nrow(plusLimit)))
        mat <- rbind(x[['cntMatrixL']], x[['cntMatrixC']],
                     minusLimit, plusLimit)
        dir <- c(rep('==', nrow(x[['cntMatrixL']]) + nrow(x[['cntMatrixC']])),
                 rep('<=', nrow(minusLimit)),
                 rep('>=', nrow(plusLimit)))
        rhs <- c(x[['lowMargins']],
                 x[['connectionProbabilities']],
                 x[['bunchProbabilities']],
                 x[['bunchProbabilities']])

        optimalCoupling <- Rglpk_solve_LP(obj = obj, mat = mat, dir = dir, rhs = rhs,
                                          max = FALSE)


        if (class(x[['cntMatrixB']]) == 'simple_triplet_matrix') {
            fittedLow <- t(slam::matprod_simple_triplet_matrix(x = x[['cntMatrixL']][, 1:ncol(x[['cntMatrixB']])],
                                                             y = optimalCoupling[['solution']][1:ncol(x[['cntMatrixB']])]))
            fittedCon <- t(slam::matprod_simple_triplet_matrix(x = x[['cntMatrixC']][, 1:ncol(x[['cntMatrixB']])],
                                                             y = optimalCoupling[['solution']][1:ncol(x[['cntMatrixB']])]))
            fittedBun <- t(slam::matprod_simple_triplet_matrix(x = x[['cntMatrixB']],
                                                             y = optimalCoupling[['solution']][1:ncol(x[['cntMatrixB']])]))
        } else {
            fittedLow <- t(x[['cntMatrixL']][, 1:ncol(x[['cntMatrixB']])] %*% optimalCoupling[['solution']][1:ncol(x[['cntMatrixB']])])
            fittedCon <- t(x[['cntMatrixC']][, 1:ncol(x[['cntMatrixB']])] %*% optimalCoupling[['solution']][1:ncol(x[['cntMatrixB']])])
            fittedBun <- t(x[['cntMatrixB']][, 1:ncol(x[['cntMatrixB']])] %*% optimalCoupling[['solution']][1:ncol(x[['cntMatrixB']])])
        }

        CNT2 <- sum(abs(x[['bunchProbabilities']] - fittedBun))

        if (abs(CNT2) < (200 * .Machine$double.eps)) {
            CNT2 <- 0
        }

        output <- list(CNT2 = CNT2, solution = optimalCoupling[['solution']],
                       status = optimalCoupling[['status']],
                       fittedLow = fittedLow, fittedCon = fittedCon, fittedBun = fittedBun)
        class(output) <- 'Cnt.Meas'

        if (output[['status']] != 0) {
            warning('Warning!. Linear programming failed to find an optimal solution')
        }

        return(output)
    } else {
        stop('x is not a cntMatrices object')
    }
}


#' Compute the noncontextuality measure NCNT2
#'
#' NCNT2(x) returns a list containing:
#'     NCNT2: The value of the NCNT2 noncontextuality measure
#'     solution: The solution found for optimal coupling
#'     status: Whether the solution was successfully found in the linear programming task.
#'         It will return 0 for the optimal solution being found, and non-zero otherwise.
#'
#' @param x a cntMatrices object
#' @param coupling optional linear programming solution of a coupling for the system (Rglpk_solve_LP or FindCoupling)
#' @return A Cnt.Meas object with the NCNT2 measure
NCNT2 <- function (x, coupling = NULL) {

    if (suppressMessages(isCntMatrices(x))) {
        require('Rglpk')

        if (is.null(coupling)) {
            coupling <- FindCoupling(x)
        }

        optimalCoupling <- coupling

        if ((coupling[['status']] == 0) & (nrow(x[['cntMatrixB']]) > 0)) {
            dMatrix <- diag(nrow = nrow(x[['cntMatrixB']]))

            x[['cntMatrixL']] <- cbind(x[['cntMatrixL']],
                                       matrix(rep(0, nrow(x[['cntMatrixL']])),
                                              nrow = nrow(x[['cntMatrixL']])))

            x[['cntMatrixC']] <- cbind(x[['cntMatrixC']],
                                       matrix(rep(0, nrow(x[['cntMatrixC']])),
                                              nrow = nrow(x[['cntMatrixC']])))

            obj <- c(rep(0, ncol(x[['cntMatrixB']])), 1)
            dir <- c(rep('==', nrow(x[['cntMatrixL']]) + nrow(x[['cntMatrixB']]) + nrow(x[['cntMatrixC']])))

            rhs <- c(x[['lowMargins']],
                     x[['bunchProbabilities']],
                     x[['connectionProbabilities']])

            optimalCouplingPlus <- optimalCouplingMinus <- list()
            optimalDs <- numeric(nrow(dMatrix))

            if (nrow(dMatrix) > 0) {
                for (iiD in seq(nrow(dMatrix))) {

                    matPlus  <- rbind(x[['cntMatrixL']],
                                      cbind(x[['cntMatrixB']], dMatrix[, iiD]),
                                      x[['cntMatrixC']])

                    matMinus <- rbind(x[['cntMatrixL']],
                                      cbind(x[['cntMatrixB']], -dMatrix[, iiD]),
                                      x[['cntMatrixC']])


                    optimalCouplingPlus[[iiD]]  <- Rglpk_solve_LP(obj = obj, mat = matPlus, dir = dir, rhs = rhs,
                                                                  max = TRUE)

                    optimalCouplingMinus[[iiD]] <- Rglpk_solve_LP(obj = obj, mat = matMinus, dir = dir, rhs = rhs,
                                                                  max = TRUE)


                    optimalDs[iiD] <- min(optimalCouplingPlus[[iiD]][['solution']][length(obj)],
                                          optimalCouplingMinus[[iiD]][['solution']][length(obj)])
                }

                x[['cntMatrixL']] <- x[['cntMatrixL']][, -ncol(x[['cntMatrixL']])]
                x[['cntMatrixC']] <- x[['cntMatrixC']][, -ncol(x[['cntMatrixC']])]


                NCNT2 <- min(optimalDs)
            } else {
                NCNT2 <- 0
            }

        } else {
            NCNT2 <- 0
        }

        if (abs(NCNT2) < (200 * .Machine$double.eps)) {
            NCNT2 <- 0
        }

        if (class(x[['cntMatrixB']]) == 'simple_triplet_matrix') {
            fittedLow <- t(slam::matprod_simple_triplet_matrix(x = x[['cntMatrixL']],
                                                             y = optimalCoupling[['solution']]))
            fittedCon <- t(slam::matprod_simple_triplet_matrix(x = x[['cntMatrixC']],
                                                             y = optimalCoupling[['solution']]))
            fittedBun <- t(slam::matprod_simple_triplet_matrix(x = x[['cntMatrixB']],
                                                             y = optimalCoupling[['solution']]))
        } else {
            fittedLow <- t(x[['cntMatrixL']] %*% optimalCoupling[['solution']])
            fittedCon <- t(x[['cntMatrixC']] %*% optimalCoupling[['solution']])
            fittedBun <- t(x[['cntMatrixB']] %*% optimalCoupling[['solution']])
        }

        output <- list(NCNT2 = NCNT2, solution = coupling[['solution']],
                       status = coupling[['status']],
                       fittedLow = fittedLow, fittedCon = fittedCon, fittedBun = fittedBun)
        class(output) <- 'Cnt.Meas'

        return(output)
    } else {
        stop('x is not a cntMatrices object')
    }
}


#' Compute the hierarchical contextuality measure CNT2h
#'
#' CNT2h(x) returns a list containing:
#'     CNT2h: The value of the CNT2h contextuality measure
#'     solution: The solution found for optimal coupling
#'     status: Whether the solution was successfully found in the linear programming task.
#'         It will return 0 for the optimal solution being found, and non-zero otherwise.
#'
#'
#' @param x a cntMatrices object
#' @return A Cnt.Meas object with the CNT2 measure
CNT2h <- function (x) {
    if (suppressMessages(isCntMatrices(x))) {
        require('Rglpk')

        if (length(x[['bunchProbabilities']]) >= 1) {

            namesProperties <- unique(unlist(lapply(strsplit(names(x[['lowMargins']])[-1], '\\.'),
                                                    function (y) y[2])
            ))

            bunchLevel <- unlist(lapply(lapply(strsplit(rownames(x[['cntMatrixB']]), '\\.'),
                                               function (y)
                                                   namesProperties %in% y), sum))

            opt <- 0
            for (llLevel in unique(bunchLevel)) {

                minusLimit <- cbind(x[['cntMatrixB']][bunchLevel == llLevel, ],
                                    -diag(nrow = sum(bunchLevel == llLevel)))

                plusLimit  <- cbind(x[['cntMatrixB']][bunchLevel == llLevel, ],
                                    diag(nrow = sum(bunchLevel == llLevel)))

                tempLowMat <- cbind(x[['cntMatrixL']],
                                    matrix(rep(0, nrow(plusLimit) * nrow(x[['cntMatrixL']])),
                                           nrow = nrow(x[['cntMatrixL']])))
                if (sum(bunchLevel < llLevel) == 0) {
                    if (class(x[['cntMatrixB']]) == 'simple_triplet_matrix') {
                        tempBunMat <- x[['cntMatrixB']][bunchLevel < llLevel, ]
                        tempBunMat$ncol <- ncol(plusLimit)
                    } else {
                        tempBunMat <- matrix(numeric(), ncol = ncol(plusLimit))
                    }
                } else {
                    tempBunMat <- cbind(x[['cntMatrixB']][bunchLevel < llLevel, ],
                                        matrix(rep(0, max(nrow(plusLimit) * nrow(x[['cntMatrixB']][bunchLevel < llLevel, ]),
                                                          nrow(plusLimit))),
                                               nrow = nrow(x[['cntMatrixB']][bunchLevel < llLevel, ])))
                }

                tempConMat <- cbind(x[['cntMatrixC']],
                                    matrix(rep(0, nrow(plusLimit) * nrow(x[['cntMatrixC']])),
                                           nrow = nrow(x[['cntMatrixC']])))

                obj <- c(rep(0, ncol(x[['cntMatrixB']])), rep(1, nrow(plusLimit)))
                mat <- rbind(tempLowMat, tempConMat, tempBunMat,
                             minusLimit, plusLimit)
                dir <- c(rep('==', nrow(tempLowMat) + nrow(tempConMat) + nrow(tempBunMat)),
                         rep('<=', nrow(minusLimit)),
                         rep('>=', nrow(plusLimit)))
                rhs <- c(x[['lowMargins']],
                         x[['connectionProbabilities']],
                         x[['bunchProbabilities']][bunchLevel <  llLevel],
                         x[['bunchProbabilities']][bunchLevel == llLevel],
                         x[['bunchProbabilities']][bunchLevel == llLevel])

                optimalCoupling <- Rglpk_solve_LP(obj = obj, mat = mat, dir = dir, rhs = rhs,
                                                  max = FALSE)

                opt <- optimalCoupling$optimum

                if (opt > 200 * .Machine$double.eps) {
                    break
                }

            }

            CNT2h <- sum(opt)

            if (class(x[['cntMatrixB']]) == 'simple_triplet_matrix') {
                fittedLow <- t(slam::matprod_simple_triplet_matrix(x = x[['cntMatrixL']][, 1:ncol(x[['cntMatrixB']])],
                                                                 y = optimalCoupling[['solution']][1:ncol(x[['cntMatrixB']])]))
                fittedCon <- t(slam::matprod_simple_triplet_matrix(x = x[['cntMatrixC']][, 1:ncol(x[['cntMatrixB']])],
                                                                 y = optimalCoupling[['solution']][1:ncol(x[['cntMatrixB']])]))
                fittedBun <- t(slam::matprod_simple_triplet_matrix(x = x[['cntMatrixB']],
                                                                 y = optimalCoupling[['solution']][1:ncol(x[['cntMatrixB']])]))
            } else {
                fittedLow <- t(x[['cntMatrixL']][, 1:ncol(x[['cntMatrixB']])] %*% optimalCoupling[['solution']][1:ncol(x[['cntMatrixB']])])
                fittedCon <- t(x[['cntMatrixC']][, 1:ncol(x[['cntMatrixB']])] %*% optimalCoupling[['solution']][1:ncol(x[['cntMatrixB']])])
                fittedBun <- t(x[['cntMatrixB']][, 1:ncol(x[['cntMatrixB']])] %*% optimalCoupling[['solution']][1:ncol(x[['cntMatrixB']])])
            }


            if (abs(CNT2h) < (200 * .Machine$double.eps)) {
                CNT2h <- 0
            }



            upperCNT2h <- upperLevels <- numeric(length = max(bunchLevel) - llLevel)
            upperFits <- list()

            if (llLevel < max(bunchLevel)) {

                iiUpper <- 1
                for (mmLevel in seq(from = llLevel + 1, to = max(bunchLevel))) {

                    minusLimit <- cbind(x[['cntMatrixB']][bunchLevel >= llLevel & bunchLevel <= mmLevel, ],
                                        -diag(nrow = sum(bunchLevel >= llLevel & bunchLevel <= mmLevel)))

                    plusLimit  <- cbind(x[['cntMatrixB']][bunchLevel >= llLevel & bunchLevel <= mmLevel, ],
                                        diag(nrow = sum(bunchLevel >= llLevel & bunchLevel <= mmLevel)))

                    tempLowMat <- cbind(x[['cntMatrixL']],
                                        matrix(rep(0, nrow(plusLimit) * nrow(x[['cntMatrixL']])),
                                               nrow = nrow(x[['cntMatrixL']])))
                    if (sum(bunchLevel < llLevel) == 0) {
                        if (class(x[['cntMatrixB']]) == 'simple_triplet_matrix') {
                            tempBunMat <- x[['cntMatrixB']][bunchLevel < llLevel, ]
                            tempBunMat$ncol <- ncol(plusLimit)
                        } else {
                            tempBunMat <- matrix(numeric(), ncol = ncol(plusLimit))
                        }
                    } else {
                        tempBunMat <- cbind(x[['cntMatrixB']][bunchLevel < llLevel, ],
                                            matrix(rep(0, max(nrow(plusLimit) * nrow(x[['cntMatrixB']][bunchLevel < llLevel, ]),
                                                              nrow(plusLimit))),
                                                   nrow = nrow(x[['cntMatrixB']][bunchLevel < llLevel, ])))
                    }

                    tempConMat <- cbind(x[['cntMatrixC']],
                                        matrix(rep(0, nrow(plusLimit) * nrow(x[['cntMatrixC']])),
                                               nrow = nrow(x[['cntMatrixC']])))

                    obj <- c(rep(0, ncol(x[['cntMatrixB']])), rep(1, nrow(plusLimit)))
                    mat <- rbind(tempLowMat, tempConMat, tempBunMat,
                                 minusLimit, plusLimit)
                    dir <- c(rep('==', nrow(tempLowMat) + nrow(tempConMat) + nrow(tempBunMat)),
                             rep('<=', nrow(minusLimit)),
                             rep('>=', nrow(plusLimit)))
                    rhs <- c(x[['lowMargins']],
                             x[['connectionProbabilities']],
                             x[['bunchProbabilities']][bunchLevel <  llLevel],
                             x[['bunchProbabilities']][bunchLevel >= llLevel & bunchLevel <= mmLevel],
                             x[['bunchProbabilities']][bunchLevel >= llLevel & bunchLevel <= mmLevel])

                    optimalCoupling <- Rglpk_solve_LP(obj = obj, mat = mat, dir = dir, rhs = rhs,
                                                      max = FALSE)

                    opt <- optimalCoupling$optimum

                    upperCNT2h[iiUpper]  <- sum(opt)
                    upperLevels[iiUpper] <- mmLevel
                    upperFits[[iiUpper]] <- list()

                    if (class(x[['cntMatrixB']]) == 'simple_triplet_matrix') {
                        upperFits[[iiUpper]][['fittedLow']] <- t(slam::matprod_simple_triplet_matrix(x = x[['cntMatrixL']][, 1:ncol(x[['cntMatrixB']])],
                                                                         y = optimalCoupling[['solution']][1:ncol(x[['cntMatrixB']])]))
                        upperFits[[iiUpper]][['fittedCon']] <- t(slam::matprod_simple_triplet_matrix(x = x[['cntMatrixC']][, 1:ncol(x[['cntMatrixB']])],
                                                                         y = optimalCoupling[['solution']][1:ncol(x[['cntMatrixB']])]))
                        upperFits[[iiUpper]][['fittedBun']] <- t(slam::matprod_simple_triplet_matrix(x = x[['cntMatrixB']],
                                                                         y = optimalCoupling[['solution']][1:ncol(x[['cntMatrixB']])]))
                    } else {
                        upperFits[[iiUpper]][['fittedLow']] <- t(x[['cntMatrixL']][, 1:ncol(x[['cntMatrixB']])] %*% optimalCoupling[['solution']][1:ncol(x[['cntMatrixB']])])
                        upperFits[[iiUpper]][['fittedCon']] <- t(x[['cntMatrixC']][, 1:ncol(x[['cntMatrixB']])] %*% optimalCoupling[['solution']][1:ncol(x[['cntMatrixB']])])
                        upperFits[[iiUpper]][['fittedBun']] <- t(x[['cntMatrixB']][, 1:ncol(x[['cntMatrixB']])] %*% optimalCoupling[['solution']][1:ncol(x[['cntMatrixB']])])
                    }

                    iiUpper <- iiUpper + 1
                }
            }

        } else {
            CNT2h <- 0
            optimalCoupling <- FindCoupling(x)
            fittedLow <- x[['lowMargins']]
            fittedCon <- x[['connectionProbabilities']]
            fittedBun <- x[['bunchProbabilities']]
        }

        output <- list(CNT2h = CNT2h, level = llLevel, solution = optimalCoupling[['solution']],
                       status = optimalCoupling[['status']],
                       fittedLow = fittedLow, fittedCon = fittedCon, fittedBun = fittedBun)


        if (exists('upperCNT2h', inherits = FALSE) && length(upperCNT2h) > 0) {
            output[['upperCNT2h']]  <- upperCNT2h
            output[['upperLevels']] <- upperLevels
            output[['upperFits']]   <- upperFits
        }


        class(output) <- 'Cnt.Meas'

        return(output)
    } else {
        stop('x is not a cntMatrices object')
    }
}

#' Compute the hierarchical noncontextuality measure NCNT2h
#'
#' NCNT2h(x) returns a list containing:
#'     NCNT2h: The value of the NCNT2h contextuality measure
#'     solution: The solution found for optimal coupling
#'     status: Whether the solution was successfully found in the linear programming task.
#'         It will return 0 for the optimal solution being found, and non-zero otherwise.
#'
#'
#' @param x a cntMatrices object
#' @return A Cnt.Meas object with the NCNT2 measure
NCNT2h <- function (x) {
    if (suppressMessages(isCntMatrices(x))) {
        require('Rglpk')

        NCNT2h <- 0
        optimalCoupling <- FindCoupling(x)

        if (length(x[['bunchProbabilities']]) >= 1) {

            namesProperties <- unique(unlist(lapply(strsplit(names(x[['lowMargins']])[-1], '\\.'),
                                                    function (y) y[2])
            ))

            bunchLevel <- unlist(lapply(lapply(strsplit(rownames(x[['cntMatrixB']]), '\\.'),
                                               function (y)
                                                   namesProperties %in% y), sum))

            tempLowMat <- cbind(x[['cntMatrixL']],
                                matrix(rep(0, nrow(x[['cntMatrixL']])),
                                       nrow = nrow(x[['cntMatrixL']])))

            tempConMat <- cbind(x[['cntMatrixC']],
                                matrix(rep(0, nrow(x[['cntMatrixC']])),
                                       nrow = nrow(x[['cntMatrixC']])))

            obj <- c(rep(0, ncol(x[['cntMatrixB']])), 1)

            lowerNCNT2h <- numeric(length(unique(bunchLevel)))
            lowerLevels <- numeric(length(unique(bunchLevel)))
            lowerFits   <- list()

            opt <- 0
            for (llLevel in unique(bunchLevel)) {

                y <- x
                y[['cntMatrixB']] <- x[['cntMatrixB']][bunchLevel <= llLevel, ]
                y[['bunchProbabilities']] <- x[['bunchProbabilities']][bunchLevel <= llLevel]

                optimalCouplingCheck  <- FindCoupling(y)

                status <- optimalCouplingCheck$status

                if (status > 0) {
                    llLevel <- llLevel - 1
                    break
                }

                dMatrix  <- diag(nrow = sum(bunchLevel == llLevel))

                optimalCouplingPlus <- optimalCouplingMinus <- list()
                optimalDs <- numeric(nrow(dMatrix))

                if (nrow(dMatrix) > 0) {
                    for (iiD in seq(nrow(dMatrix))) {

                        matPlus  <- cbind(x[['cntMatrixB']][bunchLevel == llLevel, ],  dMatrix[, iiD])
                        matMinus <- cbind(x[['cntMatrixB']][bunchLevel == llLevel, ], -dMatrix[, iiD])

                        if (sum(bunchLevel < llLevel) == 0) {
                            if (class(x[['cntMatrixB']]) == 'simple_triplet_matrix') {
                                tempBunMat <- x[['cntMatrixB']][bunchLevel < llLevel, ]
                                tempBunMat$ncol <- ncol(matPlus)
                            } else {
                                tempBunMat <- matrix(numeric(), ncol = ncol(matPlus))
                            }
                        } else {
                            tempBunMat <- cbind(x[['cntMatrixB']][bunchLevel < llLevel, ],
                                                matrix(rep(0, nrow(x[['cntMatrixB']][bunchLevel < llLevel, ])),
                                                       nrow = nrow(x[['cntMatrixB']][bunchLevel < llLevel, ])))
                        }

                        matP <- rbind(tempLowMat, tempConMat, tempBunMat,
                                      matPlus)

                        matM <- rbind(tempLowMat, tempConMat, tempBunMat,
                                      matMinus)

                        dir <- c(rep('==', nrow(tempLowMat) + nrow(tempConMat) + nrow(tempBunMat) + nrow(matPlus)))
                        rhs <- c(x[['lowMargins']],
                                 x[['connectionProbabilities']],
                                 x[['bunchProbabilities']][bunchLevel <  llLevel],
                                 x[['bunchProbabilities']][bunchLevel == llLevel])


                        optimalCouplingPlus[[iiD]]  <- Rglpk_solve_LP(obj = obj, mat = matP, dir = dir, rhs = rhs,
                                                                      max = TRUE)

                        optimalCouplingMinus[[iiD]] <- Rglpk_solve_LP(obj = obj, mat = matM, dir = dir, rhs = rhs,
                                                                      max = TRUE)


                        optimalDs[iiD] <- min(optimalCouplingPlus[[iiD]][['solution']][length(obj)],
                                              optimalCouplingMinus[[iiD]][['solution']][length(obj)])
                    }

                    NCNT2h <- min(optimalDs)

                    if (abs(NCNT2h) < (200 * .Machine$double.eps)) {
                        NCNT2h <- 0
                    }

                    lowerNCNT2h[llLevel - 1] <- NCNT2h
                    lowerLevels[llLevel - 1] <- llLevel
                    lowerFits[[llLevel - 1]]   <- list()

                    optimalCoupling <- optimalCouplingCheck

                    if (class(x[['cntMatrixB']]) == 'simple_triplet_matrix') {
                        lowerFits[[llLevel - 1]][['fittedLow']] <- t(slam::matprod_simple_triplet_matrix(x = x[['cntMatrixL']],
                                                                                                       y = optimalCoupling[['solution']]))
                        lowerFits[[llLevel - 1]][['fittedCon']] <- t(slam::matprod_simple_triplet_matrix(x = x[['cntMatrixC']],
                                                                                                       y = optimalCoupling[['solution']]))
                        lowerFits[[llLevel - 1]][['fittedBun']] <- t(slam::matprod_simple_triplet_matrix(x = x[['cntMatrixB']],
                                                                                                       y = optimalCoupling[['solution']]))
                    } else {
                        lowerFits[[llLevel - 1]][['fittedLow']] <- t(x[['cntMatrixL']] %*% optimalCoupling[['solution']])
                        lowerFits[[llLevel - 1]][['fittedCon']] <- t(x[['cntMatrixC']] %*% optimalCoupling[['solution']])
                        lowerFits[[llLevel - 1]][['fittedBun']] <- t(x[['cntMatrixB']] %*% optimalCoupling[['solution']])
                    }

                }

            }

            if (any(lowerLevels == 0)) {
                keep   <- which(lowerLevels > 0)
                lowerNCNT2h <- lowerNCNT2h[keep]
                lowerLevels <- lowerLevels[keep]
            }

            if (class(x[['cntMatrixB']]) == 'simple_triplet_matrix') {
                fittedLow <- t(slam::matprod_simple_triplet_matrix(x = x[['cntMatrixL']],
                                                                 y = optimalCoupling[['solution']]))
                fittedCon <- t(slam::matprod_simple_triplet_matrix(x = x[['cntMatrixC']],
                                                                 y = optimalCoupling[['solution']]))
                fittedBun <- t(slam::matprod_simple_triplet_matrix(x = x[['cntMatrixB']],
                                                                 y = optimalCoupling[['solution']]))
            } else {
                fittedLow <- t(x[['cntMatrixL']] %*% optimalCoupling[['solution']])
                fittedCon <- t(x[['cntMatrixC']] %*% optimalCoupling[['solution']])
                fittedBun <- t(x[['cntMatrixB']] %*% optimalCoupling[['solution']])
            }


        } else {
            llLevel <- 1
        }

        output <- list(NCNT2h = NCNT2h, level = max(get0('lowerLevels'), 1), solution = optimalCoupling[['solution']],
                       status = optimalCoupling[['status']],
                       fittedLow = fittedLow, fittedCon = fittedCon, fittedBun = fittedBun)

        if (exists('lowerNCNT2h', inherits = FALSE) && length(lowerNCNT2h) > 0) {
            output[['lowerNCNT2h']]  <- lowerNCNT2h
            output[['lowerLevels']] <- lowerLevels
            output[['lowerFits']]   <- lowerFits
        }

        class(output) <- 'Cnt.Meas'

        return(output)
    } else {
        stop('x is not a cntMatrices object')
    }
}




#' Compute the contextuality measure CNT3
#'
#' CNT3(x) returns a list containing:
#'     CNT3: The value of the CNT3 contextuality measure
#'     solution: The solution found for optimal coupling
#'     status: Whether the solution was successfully found in the linear programming task.
#'         It will return 0 for the optimal solution being found, and non-zero otherwise.
#'
#'
#' @param x a cntMatrices object
#' @return A Cnt.Meas object with the CNT3 measure
CNT3 <- function (x) {
    if (suppressMessages(isCntMatrices(x))) {
        require('Rglpk')


        obj <- c(rep(0, ncol(x[['cntMatrixB']])), rep(1, ncol(x[['cntMatrixB']])))
        mat <- rbind(x[['cntMatrixL']], x[['cntMatrixB']], x[['cntMatrixC']])
        mat <- cbind(mat, -mat)
        dir <- c(rep('==', nrow(x[['cntMatrixL']]) + nrow(x[['cntMatrixB']]) + nrow(x[['cntMatrixC']])))
        rhs <- c(x[['lowMargins']],
                 x[['bunchProbabilities']],
                 x[['connectionProbabilities']])

        optimalCoupling <- Rglpk_solve_LP(obj = obj, mat = mat, dir = dir, rhs = rhs,
                                          max = FALSE)


        CNT3 <- sum(abs(optimalCoupling[['solution']][1:ncol(x[['cntMatrixL']])] -
                            optimalCoupling[['solution']][(ncol(x[['cntMatrixL']]) + 1):length(optimalCoupling[['solution']])]
        )) - 1

        if (abs(CNT3) < (200 * .Machine$double.eps)) {
            CNT3 <- 0
        }

        if (class(x[['cntMatrixB']]) == 'simple_triplet_matrix') {
            fittedLow <- t(slam::matprod_simple_triplet_matrix(x = cbind(x[['cntMatrixL']], -x[['cntMatrixL']]),
                                                             y = optimalCoupling[['solution']]))
            fittedCon <- t(slam::matprod_simple_triplet_matrix(x = cbind(x[['cntMatrixC']], -x[['cntMatrixC']]),
                                                             y = optimalCoupling[['solution']]))
            fittedBun <- t(slam::matprod_simple_triplet_matrix(x = cbind(x[['cntMatrixB']], -x[['cntMatrixB']]),
                                                             y = optimalCoupling[['solution']]))
        } else {
            fittedLow <- t(x[['cntMatrixL']] %*% optimalCoupling[['solution']])
            fittedCon <- t(x[['cntMatrixC']] %*% optimalCoupling[['solution']])
            fittedBun <- t(x[['cntMatrixB']] %*% optimalCoupling[['solution']])
        }

        output <- list(CNT3 = CNT3, solution = optimalCoupling[['solution']],
                       status = optimalCoupling[['status']],
                       fittedLow = fittedLow, fittedCon = fittedCon, fittedBun = fittedBun)
        class(output) <- 'Cnt.Meas'

        return(output)
    } else {
        stop('x is not a cntMatrices object')
    }
}



#' Compute the contextuality measure CNFR - Contextual fraction
#' under the consistification of systems
#'
#' CNFR(x) returns a list containing:
#'     CNFR: The value of the CNFR contextuality measure
#'     solution: The solution found for optimal coupling
#'     status: Whether the solution was successfully found in the linear programming task.
#'         It will return 0 for the optimal solution being found, and non-zero otherwise.
#'
#'
#' @param x a cntMatrices object
#' @return A Cnt.Meas object with the CNFR measure
CNFR <- function (x) {
    if (any(class(x) == 'cntFracMatrices') & suppressMessages(isCntMatrices(x))) {
        require('Rglpk')

        obj <- rep(1, ncol(x[['cntMatrixL']]))
        mat <- rbind(x[['cntMatrixB']], x[['cntMatrixC']])
        dir <- c(rep('<=', nrow(x[['cntMatrixB']]) + nrow(x[['cntMatrixC']])))
        rhs <- c(x[['bunchProbabilities']],
                 x[['connectionProbabilities']])

        optimalCoupling <- Rglpk_solve_LP(obj = obj, mat = mat, dir = dir, rhs = rhs,
                                          max = TRUE)


        CNFR <- 1 - optimalCoupling[['optimum']]

        if (abs(CNFR) < (200 * .Machine$double.eps)) {
            CNFR <- 0
        }

        output <- list(CNFR = CNFR, solution = optimalCoupling[['solution']],
                       status = optimalCoupling[['status']])
        class(output) <- 'Cnt.Meas'

        return(output)
    } else {
        stop('x is not a cntFracMatrices object')
    }
}






#' Find a coupling for a cntMatrices representation of a RVSystem
#'
#' FindCoupling(x) returns a list containing:
#'     solution: The solution found for optimal coupling
#'     status: Whether the solution was successfully found in the linear programming task.
#'         It will return 0 for the optimal solution being found, and non-zero otherwise.
#'
#' @param x a cntMatrices object
#' @return List with solution and status from Rglpk_solve_LP optimization
FindCoupling <- function (x) {

    if (suppressMessages(isCntMatrices(x))) {
        obj <- rep(1, ncol(x[['cntMatrixL']]))
        mat <- rbind(x[['cntMatrixL']], x[['cntMatrixB']], x[['cntMatrixC']])
        dir <- c(rep('==', nrow(x[['cntMatrixL']]) + nrow(x[['cntMatrixB']]) + nrow(x[['cntMatrixC']])))
        rhs <- c(x[['lowMargins']],
                 x[['bunchProbabilities']],
                 x[['connectionProbabilities']])

        coupling <- Rglpk_solve_LP(obj = obj, mat = mat, dir = dir, rhs = rhs,
                                   max = FALSE)

        return(list(solution = coupling[['solution']], status = coupling[['status']], optimum = coupling[['optimum']]))
    } else {
        stop('x is not a cntMatrices object')
    }
}





#' Find the maximum value of the sum of signed elements of a vector
#' with an odd number of sign changes.
#'
#' S1Function(x) returns the maximal sum with odd number of sign changes
#'
#' @param x a numeric vector
#' @return numeric
S1Function <- function (x)  {

    if (any(is.na(x))) {
        stop('NA values not allowed')
    }

    nNegative <- sum(x <= 0)

    if (nNegative == 0) {
        x <- sort(x)
        x[1] <- -x[1]
        return(sum(x))
    }

    if ((nNegative %% 2) == 1) {
        return(sum(abs(x)))
    }

    posValues <- sort(x[x > 0])
    negValues <- sort(x[x <= 0], decreasing = TRUE)

    if ((nNegative == length(x))) {
        x <- sort(x, decreasing = TRUE)
        x[-1] <- abs(x[-1])
        return(sum(x))
    }

    if ((posValues[1] + negValues[1]) > 0) {
        negValues[-1] <- abs(negValues[-1])
    } else {
        posValues[1] <- -posValues[1]
        negValues    <- abs(negValues)
    }
    return(sum(c(negValues, posValues)))

}




#' Simple Marginal Selectivity differences sum for binary cyclic systems
#'
#' SimpleMSBinaryCyclic(x) returns the sum of the absolute differences in the expected values
#' for cyclic DRVSystem objects
#'
#' @param x a DRVSystem object or a cntMatrices object
#' @param values if NULL (Default), random variables are assumed to take values 1 (for the reference category)
#' and -1 for the other category
#' @return numeric
SimpleMSBinaryCyclic <- function(x, values = NULL) {

    if (!is.null(values)) {
        if (!is.numeric(values) || (length(values) != 2)) {
            warning('values should be a numeric vector with two values. Replaced with NULL')
            values <- NULL
        }
    }

    if (suppressMessages(isDRVSystem(x))) {
        x <- cntMatrices(x)
    }

    if (suppressMessages(isCntMatrices(x))) {

        properties <- unique(sapply(strsplit(names(x[['lowMargins']])[-1],
                                             split = '\\.'),
                                    function (y) y[2]))

        expectationsList <- list()
        if (is.null(values)) {
            values <- c(1, -1)
        }

        for (iiProperty in properties) {
            expectationsList[[iiProperty]] <- x[['lowMargins']][grep(
                pattern = paste0('.', iiProperty, '(\\.|$)'),
                perl = TRUE,
                x = names(x[["lowMargins"]])
            )]
            if (length(expectationsList[[iiProperty]]) != 2) {
                stop(iiProperty, ' appears in more (or less) than two contexts.\nSystem cannot be cyclic.')
            }

            expectationsList[[iiProperty]] <- (values[1] * expectationsList[[iiProperty]]) +
                (values[2] * (1 - expectationsList[[iiProperty]]))
        }


        output <- sum(abs(sapply(expectationsList, diff)))

        return(output)
    } else {
        stop('x is not a DRVSystem or a cntMatrices object')
    }
}




#' Rectangle limits for binary cyclic systems
#'
#' RectangleProduct(x) returns the limits for the n-dim rectangle of admissible products
#' for cyclic DRVSystem objects
#'
#' @param x a DRVSystem object or a cntMatrices object
#' @param values if NULL (Default), random variables are assumed to take values 1 (for the reference category)
#' and -1 for the other category
#' @return List with limits for n-dimensional rectangle
RectangleProduct <- function(x, values = NULL) {

    if (!is.null(values)) {
        if (!is.numeric(values) || (length(values) != 2)) {
            warning('values should be a numeric vector with two values. Replaced with NULL')
            values <- NULL
        }
    }

    if (suppressMessages(isDRVSystem(x))) {
        x <- cntMatrices(x)
    }

    if (suppressMessages(isCntMatrices(x))) {

        expectationsList <- list()
        limitsList <- list()
        if (is.null(values)) {
            values <- c(1, -1)
        }

        contexts <- unique(sapply(strsplit(names(x[['bunchProbabilities']]),
                                           split = '\\.'),
                                  function (y) y[1]))

        for (iiContext in seq_along(contexts)) {
            jjMargins <- grep(pattern = contexts[iiContext],
                              x = names(x[['lowMargins']]))
            margins <- x[['lowMargins']][jjMargins]

            expectationsList[[iiContext]] <- c(max(0, sum(margins) - 1),
                                               min(margins))

            limitsList[[iiContext]] <- (expectationsList[[iiContext]] * (diff(values)^2)) +
                (prod(values) * sum(margins)) +
                ((values[2]^2) * (1 - sum(margins)))

        }

        names(limitsList) <- contexts

        return(limitsList)
    } else {
        stop('x is not a DRVSystem or a cntMatrices object')
    }
}



#' Distance of binary cyclic system to its rectangle limits
#'
#' RectangleProductDistance(x) returns the L1 distance to the rectangle limits
#' for cyclic DRVSystem objects
#'
#' @param x a DRVSystem object or a cntMatrices object
#' @param values if NULL (Default), random variables are assumed to take values 1 (for the reference category)
#' and -1 for the other category
#' @return numeric
RectangleProductDistance <- function(x, values = NULL) {

    if (!is.null(values)) {
        if (!is.numeric(values) || (length(values) != 2)) {
            warning('values should be a numeric vector with two values. Replaced with NULL')
            values <- NULL
        }
    }

    if (suppressMessages(isDRVSystem(x))) {
        x <- cntMatrices(x)
    }

    if (suppressMessages(isCntMatrices(x))) {

        if (is.null(values)) {
            values <- c(1, -1)
        }

        limitsList <- RectangleProduct(x, values)

        contexts <- unique(sapply(strsplit(names(x[['bunchProbabilities']]),
                                           split = '\\.'),
                                  function (y) y[1]))
        distances <- productExpectations <- numeric(length = length(contexts))

        for (iiContext in seq_along(contexts)) {
            jjMargins <- grep(pattern = paste0('^', contexts[iiContext], '\\.'),
                              x = names(x[['lowMargins']]))

            cell11 <- x[['bunchProbabilities']][
                grep(pattern = paste0('^', contexts[iiContext], '\\.'),
                     x = names(x[['bunchProbabilities']]))
                ]

            margins <- x[['lowMargins']][jjMargins]

            productExpectations[iiContext] <-(cell11 * (diff(values)^2)) +
                (prod(values) * sum(margins)) +
                ((values[2]^2) * (1 - sum(margins)))

            distances[iiContext] <- min(abs(productExpectations[iiContext] -
                                                limitsList[[iiContext]]))
        }

        output <- min(distances)

        return(output)
    } else {
        stop('x is not a DRVSystem or a cntMatrices object')
    }
}




#' Compute the inequality for cyclic DRVSystems proved in Contextuality-by-Default
#' for determining whether a system is contextual or not
#'
#' @param x a DRVSystem object or a cntMatrices object
#' @param values if NULL (Default), random variables are assumed to take values 1 (for the reference category)
#' and -1 for the other category
#' @return numeric
DeltaCFunction <- function(x, values = NULL) {

    if (suppressMessages(isDRVSystem(x))) {
        x <- cntMatrices(x, withMatrices = FALSE)
    }

    if (is.null(values)) {
        values <- c(1, -1)
    }
    if (suppressMessages(isCntMatrices(x))) {
        sMS <- SimpleMSBinaryCyclic(x, values)


        contexts <- unique(sapply(strsplit(names(x[['bunchProbabilities']]),
                                           split = '\\.'),
                                  function (y) y[1]))
        diagonalExpectations <- numeric(length = length(contexts))
        for (iiContext in seq_along(contexts)) {
            jjMargins <- grep(pattern = paste0('^', contexts[iiContext], '\\.'),
                              x = names(x[['lowMargins']]))
            margins <- x[['lowMargins']][jjMargins]

            if (length(margins) != 2) {
                stop(contexts[iiContext], ' contains more (or less) than two variables.\nSystem cannot be cyclic.')
            }

            cell11 <- x[['bunchProbabilities']][
                grep(pattern = paste0('^', contexts[iiContext], '\\.'),
                     x = names(x[['bunchProbabilities']]))
                ]
            diagonalExpectations[iiContext] <- ((values[1]^2) * cell11) +
                ((values[2]^2) * (1 - sum(margins) + cell11)) +
                (prod(values) * (sum(margins) - (2 * cell11)))
        }

        lambdaC <- S1Function(diagonalExpectations)

        Delta  <- min(sMS + length(x[['bunchProbabilities']]) - 2,
                      length(x[['bunchProbabilities']]))

        deltaC <- lambdaC - Delta

        mA     <- RectangleProductDistance(x, values)

        output <- list()
        output[['lambdaC']] <- lambdaC
        output[['sMSdiff']] <- sMS
        output[['deltaC']]  <- deltaC
        output[['mA']]      <- mA

        class(output) <- 'Cnt.Meas'
        return(output)
    } else {
        stop('x is not a DRVSystem or a cntMatrices object')
    }
}



#' Construct the Boolean Matrix for the linear programming
#' for computing the contextual fraction measure
#'
#' cntFracMatrices(x) returns a list containing:
#'   cntMatrixL, cntMatrixB, cntMatrixC: boolean (sparse) matrices with the appropriate incidences
#'     of elements of the coupling vector (identified by the columns) to the vector of probabilities
#'     describing the system for the low-order marginals, bunch, and connection probabilities, respectively.
#'   lowMargins, bunchProbabilities, connectionProbabilities: the vectors of probabilities that describe the system.
#'
#'
#' @param x an RVSystem
#' @return a cntMatrices object
#' @examples
cntFracMatrices <- function (x) {

    require('slam')
    require('abind')

    if (suppressMessages(isDRVSystem(x))) {

        if (is.null(names(x)) || any(names(x) == "")) {
            names(x) <- paste0('Ctx', seq(length(x)))
        }

        nContexts     <- length(x)
        nVariables    <- sum(sapply(x, function (y) length(dim(y))))

        nColumns <- 2 ^ nVariables

        namesProperties <- unique(unlist(lapply(x,
                                                function (y) names(dimnames(y)))))

        refValues        <- character(length = length(namesProperties))
        ferValues        <- character(length = length(namesProperties))
        names(ferValues) <- names(refValues) <- namesProperties

        for (iiName in namesProperties) {
            for (jjContext in names(x)) {
                if (iiName %in% names(dimnames(x[[jjContext]]))) {
                    ferValues[iiName] <-  dimnames(x[[jjContext]])[[iiName]][1]
                    refValues[iiName] <-  dimnames(x[[jjContext]])[[iiName]][2]
                }
            }
        }

        # Low-order marginals

        cntMatrixL <- slam::simple_triplet_zero_matrix(nrow = 2 * nVariables, ncol = nColumns,
                                                       mode = 'integer')

        cntMatrixL           <- rbind(t(rep(1, nColumns)), cntMatrixL)
        rownames(cntMatrixL) <- rep('', (2 * nVariables) + 1)

        lowMargins        <- c(1, numeric(2 * nVariables))
        names(lowMargins) <- rep('', (2 * nVariables) + 1)

        iiRow <- 2

        for (jjContext in names(x)) {

            for (kkVariable in seq(dimnames(x[[jjContext]]))) {

                kkVariableName <- names(dimnames(x[[jjContext]]))[kkVariable]

                cntMatrixL[iiRow, ]             <- rep(c(1, 0), each = 2 ** (nVariables - (iiRow / 2)))
                cntMatrixL[iiRow + 1, ]         <- rep(c(0, 1), each = 2 ** (nVariables - (iiRow / 2)))

                rownames(cntMatrixL)[iiRow + 1] <- paste0(jjContext, '.',
                                                          kkVariableName, '.',
                                                          ferValues[kkVariableName])
                rownames(cntMatrixL)[iiRow]     <- paste0(jjContext, '.',
                                                          kkVariableName, '.',
                                                          refValues[kkVariableName])

                lowMargins[(iiRow + 1):iiRow]  <- apply(X = x[[jjContext]],
                                                        MARGIN = kkVariable,
                                                        FUN = sum)

                names(lowMargins)[iiRow + 1]   <- paste0(jjContext, '.',
                                                         kkVariableName, '.',
                                                         ferValues[kkVariableName])
                names(lowMargins)[iiRow]       <- paste0(jjContext, '.',
                                                         kkVariableName, '.',
                                                         refValues[kkVariableName])

                iiRow <- iiRow + 2

            }
        }

        # Bunch probabilities

        bunchesVariables    <- lapply(x, function (y) names(dimnames(y)))
        nBunchesVariables   <- sapply(x, function (y) length(dim(y)))
        nBunchProbabilities <- sum(2 ^ nBunchesVariables)

        cntMatrixB <- slam::simple_triplet_zero_matrix(nrow = nBunchProbabilities,
                                                       ncol = nColumns,
                                                       mode = 'integer')
        rownames(cntMatrixB) <- rep('', nBunchProbabilities)

        bunchProbabilities <- integer(nBunchProbabilities)
        names(bunchProbabilities) <- rep('', nBunchProbabilities)

        if (nBunchProbabilities > 1) {
            iiRow <- 1
            for (jjContext in names(x)) {

                contextProbabilities <- as.data.frame(x[[jjContext]])
                contextVariables     <- names(contextProbabilities)[
                    -grep('Freq', names(contextProbabilities))
                    ]

                for (kkVariable in contextVariables) {
                    contextProbabilities[, kkVariable] <- paste0(kkVariable, '.',
                                                                 contextProbabilities[, kkVariable])
                }

                for (mmCell in seq(nrow(contextProbabilities))) {
                    bunchProbabilities[iiRow]        <- contextProbabilities[mmCell, 'Freq']

                    cntMatrixB[iiRow, ] <- apply(cntMatrixL[paste0(jjContext, '.',
                                                                   contextProbabilities[mmCell,
                                                                                        -grep('Freq', names(contextProbabilities))])
                                                            , ], 2, prod)

                    rownames(cntMatrixB)[iiRow] <- names(bunchProbabilities)[iiRow] <- paste0(jjContext, '.',
                                                                                              paste(contextProbabilities[mmCell,
                                                                                                                         -grep('Freq', names(contextProbabilities))
                                                                                                                         ], collapse = '.'))

                    iiRow <- iiRow + 1
                }
            }
        }

        # Connection probabilities

        contextsNames <- names(x)

        connectionBunches   <- lapply(namesProperties,
                                      function (z)
                                          contextsNames[sapply(
                                              lapply(x,
                                                     function (y)
                                                         names(dimnames(y)) %in% z), any)]
        )
        names(connectionBunches) <- namesProperties

        nConnectionVariables     <- sapply(connectionBunches, length)
        nConnectionProbabilities <- 4 * sum(choose(nConnectionVariables, 2))

        cntMatrixC <- slam::simple_triplet_zero_matrix(nrow = nConnectionProbabilities,
                                                       ncol = nColumns,
                                                       mode = 'integer')
        rownames(cntMatrixC) <- rep('', nConnectionProbabilities)

        connectionProbabilities <- integer(nConnectionProbabilities)
        names(connectionProbabilities) <- rep('', nConnectionProbabilities)

        if (nConnectionProbabilities > 1) {
            iiRow <- 1
            for (jjProperty in namesProperties) {

                mmContextPairs <- combn(x = seq_along(connectionBunches[[jjProperty]]),
                                        m = 2)

                for (kkkPair in seq(ncol(mmContextPairs))) {

                    indexContexts  <- mmContextPairs[, kkkPair]
                    namesContexts  <- connectionBunches[[jjProperty]][indexContexts]

                    valuesVar <- c(refValues[jjProperty], ferValues[jjProperty])

                    for (iiValues in seq_along(valuesVar)) {
                        for (jjValues in seq_along(valuesVar)) {

                            if (iiValues == jjValues) {
                                connectionProbabilities[iiRow] <- min(lowMargins[
                                    paste0(connectionBunches[[jjProperty]][indexContexts],
                                           '.', jjProperty,
                                           '.', valuesVar[c(iiValues, jjValues)])
                                    ])
                            } else {
                                connectionProbabilities[iiRow] <- lowMargins[
                                    paste0(connectionBunches[[jjProperty]][indexContexts][1],
                                           '.', jjProperty, '.', valuesVar[iiValues])
                                    ] -
                                    min(lowMargins[
                                        paste0(connectionBunches[[jjProperty]][indexContexts],
                                               '.', jjProperty, '.', valuesVar[iiValues])
                                        ])
                            }

                            cntMatrixC[iiRow, ] <- apply(cntMatrixL[paste0(connectionBunches[[jjProperty]][indexContexts],
                                                                           '.', jjProperty,
                                                                           '.', valuesVar[c(iiValues, jjValues)]), ], 2, prod)

                            rownames(cntMatrixC)[iiRow] <-
                                names(connectionProbabilities)[iiRow] <-
                                paste(c(jjProperty, connectionBunches[[jjProperty]][indexContexts], valuesVar[c(iiValues, jjValues)]), collapse = '.')

                            iiRow <- iiRow + 1
                        }
                    }
                }
            }
        }

        out <- list(cntMatrixL = cntMatrixL, cntMatrixB = cntMatrixB, cntMatrixC = cntMatrixC,
                    lowMargins = lowMargins, bunchProbabilities = bunchProbabilities,
                    connectionProbabilities = connectionProbabilities)

        class(out) <- c('cntMatrices', 'cntFracMatrices')

        return(out)
    } else {
        stop('x is not a dichotomous RVSystem')
    }
}





#' Print the value of the (non)contextuality measure from a Cnt.Meas object
#'
#' print.Cnt.Meas(x) prints the value of a (non)contextuality measure
#'
#'
#' @param x a Cnt.Meas object
print.Cnt.Meas <- function (x) {
    if (class(x) == 'Cnt.Meas') {
        measure <- names(x)[grep('^CNT|^NCNT|^HCNT|^HNCNT|CNFR|delta', names(x))]
        print(paste0(measure, ': ', x[[measure]]))
    }
}



#' Compute the CbD (non)contextuality measures for a DRVSystem or a cntMatrices object
#'
#' CbDMeasures(x) returns a list containnig the values of the (non)contextuality measures for x
#'
#' @param x A DRVSystem or cntMatrices object
#' @param hierLevels Logical. If TRUE, upper levels for CNT2h and lower levels for NCNT2h when available. Defaults to FALSE
#' @param as.df Logical. If TRUE, return as data.frame. Defaults to FALSE
#' @return A list containing the values of the measures CNT1, CNT2, CNT2h, CNT3, CNFR, NCNT2, NCNT2h.
#'  If x is cyclic systems, also includes the Delta (CHSH-type) inequality and
#'  the distance from the system to the n-dimensional rectangle.
#'  If x is a DRVSystem object (rather than a ctMatrices object), the Contextual Fraction
#'  measure is also computed.
CbDMeasures <- function (x, hierLevels = FALSE, as.df = FALSE) {
    if (suppressMessages(isDRVSystem(x))) {
        y <- cntFracMatrices(x)
        x <- cntMatrices(x)
    }

    if (suppressMessages(isCntMatrices(x))) {

        deltaX <- suppressWarnings(try(DeltaCFunction(x), silent = TRUE))

        CNT1   <- CNT1(x)
        CNT2   <- CNT2(x)
        CNT2h  <- CNT2h(x)
        CNT3   <- CNT3(x)
        NCNT2  <- NCNT2(x)
        NCNT2h <- NCNT2h(x)

        if (exists('y', inherits = FALSE)) {
            CNFR <- CNFR(y)
        } else {
        }

        output <- list(
            CNT1   = CNT1[['CNT1']],
            CNT2   = CNT2[['CNT2']],
            CNT3   = CNT3[['CNT3']],
            NCNT2  = NCNT2[['NCNT2']]
        )

        if (exists('CNFR', inherits = FALSE)) {
            output[['CNFR']] <- CNFR[['CNFR']]
        }

        if (length(grep('Error', deltaX)) < 1) {
            output[['Delta']] <- deltaX[['deltaC']]
            output[['mA']] <- deltaX[['mA']]
        }

        if (hierLevels) {
            if (hasName(CNT2h, 'upperCNT2h')) {
                CNT2hs <- c(CNT2h[['CNT2h']], CNT2h[['upperCNT2h']])
                names(CNT2hs) <- seq(from = CNT2h[['level']], to = max(CNT2h[['upperLevels']]))
                for (iiLevels in names(CNT2hs)) {
                    output[[paste0('CNT2h.', iiLevels)]] <- CNT2hs[iiLevels]
                }
            } else {
                output[[paste0('CNT2h.', CNT2h[['level']])]] <- CNT2h[['CNT2h']]
            }
            if (hasName(NCNT2h, 'lowerNCNT2h')) {
                NCNT2hs <- NCNT2h[['lowerNCNT2h']]
                names(NCNT2hs) <- NCNT2h[['lowerLevels']]
                for (iiLevels in names(NCNT2hs)) {
                    output[[paste0('NCNT2h', '.', iiLevels)]] <- NCNT2hs[iiLevels]
                }
            }
        } else {
            output[['CNT2h']]  <- CNT2h[['CNT2h']]
            output[['levelC']] <- CNT2h[['level']]
            output[['NCNT2h']] <- NCNT2h[['NCNT2h']]
            output[['levelN']] <- NCNT2h[['level']]
        }

        if (as.df) {
            output <- data.frame(output)
        }

        return(output)

    } else {
        stop('x is not a DRVSystem or a cntMatrices object')
    }
}



#' Generate cyclic systems with the underlying structure of a maximally contextual
#' cyclic system of the given rank
#'
#' cyclicRVSystem(n) returns an RVSystem object that is dichotomous and cyclic
#'
#' @param n integer greater or equal to 2
#' @param refTable
#' @param negTable
#' @return a cyclic DRVSystem of rank n
cyclicRVSystem <- function (n, refTable = NULL, negTable = NULL) {

    if (is.null(refTable)) {
        refTable <- as.table(diag(2)) / 2
    }

    if (is.null(negTable)) {
        negTable <- refTable[2:1, ]
    }

    varNames <- LETTERS[c(1:n, 1)]

    rvSystem <- do.call(RVSystem, lapply(as.list(
        c(rep('refTable', n - 1),
          'negTable')),
        as.name))

    for (iiContext in seq_along(rvSystem)) {
        dimnames(rvSystem[[iiContext]])[[1]] <-
            dimnames(rvSystem[[iiContext]])[[2]] <-
            c('-1', '1')
        names(dimnames(rvSystem[[iiContext]])) <- varNames[c(iiContext, iiContext + 1)]
    }

    return(rvSystem)

}



#' Generates two 2x2 prop.tables with given margins
#' (for a reference category) and a given probability of
#' the conjunction (of both reference categories)
#'
#' generateRefFerTables(product, margin1, margin2)
#'
#' @param product numeric between 0 and 1
#' @param margin1 numeric between 0 and 1
#' @param margin2 numeric between 0 and 1
#'
#' @return A list with prop.tables refTable and ferTable
generateRefFerTables <- function(product = 0.5, margin1 = 0.5, margin2 = 0.5) {
    if (round(product, 10) > round(margin1, 10) || round(product, 10) > round(margin2, 10)) {
        stop("product cell needs to have a smaller probability than each margin")
    }
    if (round(product, 10) > 1 || round(product, 10) < 0 ||
        round(margin1, 10) > 1 || round(margin1, 10) < 0 ||
        round(margin2, 10) > 1 || round(margin2, 10) < 0) {
        stop("product and margins need to be probabilities")
    }

    if (round(product, 10) < round(margin1 + margin2 - 1, 10)) {
        stop("product cell needs to have a larger probability than the sum of the margins minus 1")
    }

    maxProd <- min(margin1, margin2)
    minProd <- max(0, margin1 + margin2 - 1)

    epsilon <- maxProd - product

    ferTable <- as.table(matrix(c(1 - margin1 - margin2 + maxProd - epsilon,
                                  margin1 - maxProd + epsilon,
                                  margin2 - maxProd + epsilon,
                                  maxProd - epsilon),
                                ncol = 2))
    refTable <-  as.table(matrix(c(1 - margin1 - margin2 + minProd + epsilon,
                                   margin1 - minProd - epsilon,
                                   margin2 - minProd - epsilon,
                                   minProd + epsilon),
                                 ncol = 2))

    output <- list(refTable = refTable, ferTable = ferTable)

    return(output)
}



#' checks whether x gives the structure of a PR box
#'
#' isPRStructure(x)
#'
#' @param x a binary table
#' @return logical
isPRStructure <- function (x) {
    test <- TRUE
    if (!is.table(x) & !is.matrix(x)) {
        message('x should be a table or a matrix')
        test <- FALSE
    }

    if (nrow(x) != ncol(x)) {
        message('x should be a square table')
        test <- FALSE
    }

    if (!all(x %in% c(0, 1))) {
        message('x should contain only 0s or 1s')
        test <- FALSE
    }

    if (sum(diag(x[-1, ])) != 0) {
        message('x should contain only 0s in the main subdiagonal')
        test <- FALSE
    }

    n <- ncol(x)

    if (sum(x) != (n^2 - n + 1)) {
        message('x should contain ', n^2 - n + 1, ' 1s')
        test <- FALSE
    }

    if (!is.null(rownames(x))) {
        namesDim <- unique(c(rownames(x), colnames(x)))
        if (length(namesDim) != (nrow(x) + ncol(x))) {
            message('Each row and column of x should have a unique name')
            test <- FALSE
        }
    }

    return(test)
}


#' Generates a nxn table that gives the structure of an n-PR box
#' arXiv:1903.07170v5 # for several non-cyclic systems based on PR boxes
#' by Brunner, Scarani & Gisin (2006) doi:10.1063/1.2352857
#'
#' generatePRBoxNStructure(n)
#'
#' @param n an integer equal or greater than 2 (e.g. 2L, 4L)
#' @return nxn table of zeros and ones with the structure for a n-PR box
generatePRBoxNStructure <- function (n) {

    if (!is.integer(n) || (n < 2)) {
        stop('n should be at least 2')
    }


    x <- matrix(data = 1, nrow = n, ncol = n)
    y <- matrix(data = seq(n^2), nrow = n, ncol = n)

    x[diag(y[-1, , drop = FALSE])] <- 0

    rownames(x) <- paste0('A', seq(nrow(x)))
    colnames(x) <- paste0('B', seq(nrow(x)))

    return(x)
}



#' Generates a RVSystem with the structure of an n-PR box
#' populated by joint distributions given by refTable, negTable and chgTable
#' or perfectly (anti)correlated joint distributions (i.e, the PR box system) if all are NULL
#'
#' @param PRBoxStructure a table giving a PRBoxStruvture (see isPRStructure)
#' @param posTable a prop.table that will be used for all entries of PRBoxStructure equal to 1
#' @param negTable a prop.table that will be used for all entries of PRBoxStructure equal to 0 (except when chgTable is not NULL)
#' @param chgTable a prop.table that will be used for the first entry of PRBoxStructure that equals 0
generatePRBoxLikeSystem <- function(PRBoxStructure, posTable = NULL, negTable = NULL, chgTable = NULL) {

    if(suppressMessages(isPRStructure(PRBoxStructure))) {

        if (is.null(posTable)) {
            posTable <- as.table(diag(2)) / 2
        }

        if (is.null(negTable)) {
            negTable <- posTable[2:1, ]
        }

        if (is.null(chgTable)) {
            chgTable <- negTable
        }

        listJoints <- list()
        rows <- row(PRBoxStructure)
        cols <- col(PRBoxStructure)


        jjOddCtx <- 1
        for (iiContext in seq_along(PRBoxStructure)) {
            if (PRBoxStructure[iiContext] == 1) {
                listJoints[[iiContext]] <- posTable
            } else {
                if (jjOddCtx == 1) {
                    listJoints[[iiContext]] <- chgTable
                    jjOddCtx <- 0
                } else {
                    listJoints[[iiContext]] <- negTable
                }
            }

            names(dimnames(listJoints[[iiContext]])) <-
                c(rownames(PRBoxStructure)[rows[iiContext]],
                  colnames(PRBoxStructure)[cols[iiContext]])

            dimnames(listJoints[[iiContext]])[[1]] <-
                dimnames(listJoints[[iiContext]])[[2]] <-
                c('-1', '1')
        }

        out <- do.call(RVSystem, listJoints)

        return(out)
    }
}

#' Create a prop.table from a set of patterns
#'
#' MakeJointFromPatterns(patterns) returns a prop.table with named dimensions
#'
#' @param patterns a matrix or data.frame of deterministic patterns. A deterministic pattern is given by an integer vector.
#' @param namesVec optional list of names for the columns of patterns.
MakeJointFromPatterns <- function (patterns, namesVec = NULL) {

    if (!(is.matrix(patterns) | is.data.frame(patterns))) {
        stop('patterns must be a matrix.')
    }

    if (!is.null(namesVec)) {
        if (length(namesVec) != ncol(patterns)) {
            warning('The number of names in namesVec should be the same as the number of columns of patterns.')
        } else {
            colnames(patterns) <- namesVec
        }
    } else {
        if (is.null(colnames(patterns))) {
            colnames(patterns) <- paste0('V.', seq(ncol(patterns)))
        }
    }

    jointTable <- prop.table(table(as.data.frame(patterns)))

    return(jointTable)

}
