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
RVSystem <- function (...) {
    out <- list(...)
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
#' @return
#' @examples
cntMatrices <- function (x) {

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

        cntMatrixL <- slam::simple_triplet_zero_matrix(nrow = nVariables, ncol = nColumns,
                                                       mode = 'integer')

        cntMatrixL           <- rbind(t(rep(1, nColumns)), cntMatrixL)
        rownames(cntMatrixL) <- rep('', nVariables + 1)

        lowMargins        <- c(1, numeric(nVariables))
        names(lowMargins) <- rep('', nVariables + 1)

        iiRow <- 2

        for (jjContext in names(x)) {

            for (kkVariable in seq(dimnames(x[[jjContext]]))) {

                kkVariableName <- names(dimnames(x[[jjContext]]))[kkVariable]

                cntMatrixL[iiRow, ]         <- rep(c(1, 0), each = 2 ** (nVariables - iiRow + 1))
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

        cntMatrixB <- slam::simple_triplet_zero_matrix(nrow = nBunchProbabilities,
                                                       ncol = nColumns,
                                                       mode = 'integer')

        rownames(cntMatrixB)      <- rep('', nBunchProbabilities)

        bunchProbabilities        <- integer(nBunchProbabilities)
        names(bunchProbabilities) <- rep('', nBunchProbabilities)

        iiRow <- 1
        for (jjContext in names(x)) {

            for (kkMargins in 2:nBunchesVariables[jjContext]) {

                mmMarginNames <- combn(x = seq_along(bunchesVariables[[jjContext]]),
                                       m = kkMargins)

                for (kkkMargin in seq(ncol(mmMarginNames))) {

                    indexVariables  <- mmMarginNames[, kkkMargin]
                    namesVariables  <- bunchesVariables[[jjContext]][indexVariables]
                    valuesVariables <- refValues[namesVariables]

                    marginsTable <- apply(x[[jjContext]], indexVariables, sum)

                    bunchProbabilities[iiRow] <- abind::asub(marginsTable,
                                                             as.list(valuesVariables),
                                                             indexVariables)

                    cntMatrixB[iiRow, ] <- apply(cntMatrixL[paste0(jjContext, '.', namesVariables), ], 2, prod)

                    rownames(cntMatrixB)[iiRow]          <-
                        names(bunchProbabilities)[iiRow] <-
                        paste(c(jjContext, namesVariables), collapse = '.')

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
        nConnectionProbabilities <- sum(choose(nConnectionVariables, 2))

        cntMatrixC <- slam::simple_triplet_zero_matrix(nrow = nConnectionProbabilities,
                                                       ncol = nColumns,
                                                       mode = 'integer')

        rownames(cntMatrixC)           <- rep('', nConnectionProbabilities)

        connectionProbabilities        <- integer(nConnectionProbabilities)
        names(connectionProbabilities) <- rep('', nConnectionProbabilities)

        iiRow <- 1
        for (jjProperty in namesProperties) {

            mmContextPairs <- combn(x = seq_along(connectionBunches[[jjProperty]]),
                                    m = 2)

            for (kkkPair in seq(ncol(mmContextPairs))) {

                indexContexts   <- mmContextPairs[, kkkPair]
                namesContexts   <- connectionBunches[[jjProperty]][indexContexts]
                valuesVariables <- rep(refValues[jjProperty], length(indexContexts))

                connectionProbabilities[iiRow] <- min(lowMargins[paste0(connectionBunches[[jjProperty]],
                                                                        '.', jjProperty)])

                cntMatrixC[iiRow, ] <- apply(cntMatrixL[paste0(connectionBunches[[jjProperty]],
                                                               '.', jjProperty), ], 2, prod)

                rownames(cntMatrixC)[iiRow]               <-
                    names(connectionProbabilities)[iiRow] <-
                    paste(c(jjProperty, connectionBunches[[jjProperty]]), collapse = '.')

                iiRow <- iiRow + 1
            }
        }

        out <- list(cntMatrixL = cntMatrixL, cntMatrixB = cntMatrixB, cntMatrixC = cntMatrixC,
                    lowMargins = lowMargins, bunchProbabilities = bunchProbabilities,
                    connectionProbabilities = connectionProbabilities)

        class(out) <- 'cntMatrices'

        return(out)
    } else {
        stop('x is not a divhotomous RVSystem')
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

    if (class(x[['cntMatrixL']]) == 'simple_triplet_matrix') {
        require('slam')
    }

    if (class(x[['cntMatrixL']]) == 'Matrix') {
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
            CNT1 <- sum(x[['connectionProbabilities']] -
                            slam::matprod_simple_triplet_matrix(x = x[['cntMatrixC']],
                                                                y = optimalCoupling[['solution']]))
        } else {
            CNT1 <- sum(x[['connectionProbabilities']] - (x[['cntMatrixC']] %*% optimalCoupling[['solution']]))
        }

        if (abs(CNT1) < (200 * .Machine$double.eps)) {
            CNT1 <- 0
        }

        output <- list(CNT1 = CNT1, solution = optimalCoupling[['solution']],
                       status = optimalCoupling[['status']])
        class(output) <- 'Cnt.Meas'

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


        if (class(x[['cntMatrixL']]) == 'simple_triplet_matrix') {
            CNT2 <- sum(abs(x[['bunchProbabilities']] -
                                slam::matprod_simple_triplet_matrix(x = x[['cntMatrixB']],
                                                                    y = optimalCoupling[['solution']][1:ncol(x[['cntMatrixB']])])
            ))
        } else {
            CNT2 <- sum(abs(x[['connectionProbabilities']] -
                                (x[['cntMatrixB']] %*%
                                     optimalCoupling[['solution']][1:ncol(x[['cntMatrixB']])])))
        }

        if (abs(CNT2) < (200 * .Machine$double.eps)) {
            CNT2 <- 0
        }

        output <- list(CNT2 = CNT2, solution = optimalCoupling[['solution']],
                       status = optimalCoupling[['status']])
        class(output) <- 'Cnt.Meas'

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
#' @param y optional linear programming solution of a coupling for the system
NCNT2 <- function (x, coupling = NULL) {

    if (suppressMessages(isCntMatrices(x))) {
        require('Rglpk')

        if (is.null(coupling)) {
            coupling <- FindCoupling(x)
        }

        if (coupling[['status']] == 0) {
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

            NCNT2 <- min(optimalDs)
        } else {
            NCNT2 <- 0
        }

        if (abs(NCNT2) < (200 * .Machine$double.eps)) {
            NCNT2 <- 0
        }


        output <- list(NCNT2 = NCNT2, solution = coupling[['solution']],
                       status = coupling[['status']])
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

        output <- list(CNT3 = CNT3, solution = optimalCoupling[['solution']],
                       status = optimalCoupling[['status']])
        class(output) <- 'Cnt.Meas'

        return(output)
    } else {
        stop('x is not a cntMatrices object')
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

        return(list(solution = coupling[['solution']], status = coupling[['status']]))
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
                pattern = paste0('.', iiProperty),
                x = names(x[["lowMargins"]])
            )]
            if (length(expectationsList[[iiProperty]]) != 2) {
                stop(iiProperty, ' appears in more (or less) than two contexts.\\System cannot be cyclic.')
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
            jjMargins <- grep(pattern = contexts[iiContext],
                              x = names(x[['lowMargins']]))

            cell11 <- x[['bunchProbabilities']][
                grep(pattern = contexts[iiContext],
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
DeltaCFunction <- function(x, values = NULL) {

    if (suppressMessages(isDRVSystem(x))) {
        x <- cntMatrices(x)
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
            jjMargins <- grep(pattern = contexts[iiContext],
                              x = names(x[['lowMargins']]))
            margins <- x[['lowMargins']][jjMargins]

            cell11 <- x[['bunchProbabilities']][
                grep(pattern = contexts[iiContext],
                     x = names(x[['bunchProbabilities']]))
                ]
            diagonalExpectations[iiContext] <- ((values[1]^2) * cell11) +
                ((values[2]^2) * (1 - sum(margins) + cell11)) +
                (prod(values) * (sum(margins) - (2 * cell11)))
        }

        lambdaC <- S1Function(diagonalExpectations)

        deltaC <- lambdaC - sMS - length(x[['bunchProbabilities']]) + 2

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



#' Construct the reduced Boolean Matrix for the linear programming
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
#' @return
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

                    break
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
        nBunchProbabilities <- sum(2 * (2 ^ nBunchesVariables))

        cntMatrixB <- slam::simple_triplet_zero_matrix(nrow = nBunchProbabilities,
                                                       ncol = nColumns,
                                                       mode = 'integer')
        rownames(cntMatrixB) <- rep('', nBunchProbabilities)

        bunchProbabilities <- integer(nBunchProbabilities)
        names(bunchProbabilities) <- rep('', nBunchProbabilities)

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

        cntMatrixC <- slam::simple_triplet_zero_matrix(nrow = nConnectionProbabilities,
                                                       ncol = nColumns,
                                                       mode = 'integer')
        rownames(cntMatrixC) <- rep('', nConnectionProbabilities)

        connectionProbabilities <- integer(nConnectionProbabilities)
        names(connectionProbabilities) <- rep('', nConnectionProbabilities)

        iiRow <- 1
        for (jjProperty in namesProperties) {

            mmContextPairs <- combn(x = seq_along(connectionBunches[[jjProperty]]),
                                    m = 2)

            for (kkkPair in seq(ncol(mmContextPairs))) {

                indexContexts  <- mmContextPairs[, kkkPair]
                namesContexts  <- connectionBunches[[jjProperty]][indexContexts]
                valuesVariables <- rep(refValues[jjProperty], length(indexContexts))

                connectionProbabilities[iiRow] <- min(lowMargins[paste0(connectionBunches[[jjProperty]],
                                                                        '.', jjProperty)])

                cntMatrixC[iiRow, ] <- apply(cntMatrixL[paste0(connectionBunches[[jjProperty]],
                                                               '.', jjProperty), ], 2, prod)

                rownames(cntMatrixC)[iiRow] <-
                    names(connectionProbabilities)[iiRow] <-
                    paste(c(jjProperty, connectionBunches[[jjProperty]]), collapse = '.')

                iiRow <- iiRow + 1
            }
        }

        out <- list(cntMatrixL = cntMatrixL, cntMatrixB = cntMatrixB, cntMatrixC = cntMatrixC,
                    lowMargins = lowMargins, bunchProbabilities = bunchProbabilities,
                    connectionProbabilities = connectionProbabilities)

        class(out) <- 'cntMatrices'

        return(out)
    } else {
        stop('x is not a divhotomous RVSystem')
    }
}




print.Cnt.Meas <- function (x) {
    if (class(x) == 'Cnt.Meas') {
        measure <- names(x)[grep('CNT|delta', names(x))]
        print(paste0(measure, ': ', x[[measure]]))
    }
}



CbDMeasures <- function (x) {
    if (suppressMessages(isDRVSystem(x))) {
        x <- cntMatrices(x)
    }

    if (suppressMessages(isCntMatrices(x))) {

        deltaX <- suppressWarnings(try(DeltaCFunction(x), silent = TRUE))

        output <- list(
            CNT1  = CNT1(x)[['CNT1']],
            CNT2  = CNT2(x)[['CNT2']],
            CNT3  = CNT3(x)[['CNT3']],
            NCNT2 = NCNT2(x)[['NCNT2']]
        )

        if (!grep('Error', deltaX)) {
            output[['Delta']] <- deltaX[['deltaC']]
            output[['mA']] <- deltaX[['mA']]
        }


        return(output)

    } else {
        stop('x is not a DRVSystem or a cntMatrices object')
    }
}

