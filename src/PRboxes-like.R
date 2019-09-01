###############################################################################
# # Measures of contextuality and noncontextuality
# # Computations of the measures from
# # arXiv:1903.07170v5 # for several non-cyclic systems based on PR boxes
# # by Brunner, Scarani & Gisin (2006) doi:10.1063/1.2352857
# #
# # Víctor H Cervantes
# # 2019
###############################################################################


library("ggplot2")
library("dplyr")

source("./ContextualityMeasures.R")


cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
                "#D55E00", "#CC79A7")


measureLabels <- c('CNT1', 'CNT2', 'CNT3', 'NCNT2',
                   'CNT1/CNT2', 'CNT1/CNT3', '2*CNT2/CNT3')



generateRefFerTables <- function(product = 0.5, margin1 = 0.5, margin2 = 0.5) {
  if (product > margin1 || product > margin2) {
    stop("product cell needs to have a smaller probability than each margin")
  }
  if (product > 1 || product < 0 || margin1 > 1 || margin1 < 0 || margin2 > 1 ||
      margin2 < 0) {
    stop("product and margins need to be probabilities")
  }

  refTable <- as.table(matrix(c(product, margin1 - product, margin2 - product,
                                1 - margin1 - margin2 + product), ncol = 2))
  ferTable <- refTable[2:1, ]

  output <- list(refTable = refTable, ferTable = ferTable)

  return(output)
}



isPRStructure <- function (x) {
  test <- TRUE
  if (!is.table(x) & !is.matrix(x)) {
    message('x should be a table or a matrix')
    test <- FALSE
  }

  if (!all(x %in% c(0, 1))) {
    message('x should contain only 0s or 1s')
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



################################################################################
# #  Constant |cor| varying cor from 1 in refTable (-1 in ferTable) to -1 (1)
################################################################################

nCors <- 41

listReFerTables <- list()
products <- seq(from = 0.5, to = 0, length.out = nCors)

for (iiProduct in seq_along(products)) {
    listReFerTables[[iiProduct]] <- generateRefFerTables(product = products[iiProduct])
    names(listReFerTables)[iiProduct] <- sprintf('%.3f', ((4 * products[iiProduct]) - 1))

}


measuresVaryingCor <- list()


    for (jjCor in seq(nCors)) {

      print(jjCor)

        refTable <- listReFerTables[[jjCor]][['refTable']]
        negTable <- listReFerTables[[jjCor]][['ferTable']]

        measuresVaryingCor[[jjCor]] <- sapply(CbDMeasures(
            generatePRBoxLikeSystem(PRBoxStructure = generatePRBoxNStructure(n = 3L),
                           posTable = refTable,
                           negTable = negTable,
                           chgTable = NULL)),
            function (y) y[1])
    }

    names(measuresVaryingCor) <- paste0('Value.', seq(nCors))

    measuresVaryingCor <- as.data.frame(measuresVaryingCor)

    measuresVaryingCor['CNT1/CNT2', ] <-
        measuresVaryingCor['CNT1', ] /
        measuresVaryingCor['CNT2', ]

    measuresVaryingCor['CNT1/CNT3', ] <-
        measuresVaryingCor['CNT1', ] /
        measuresVaryingCor['CNT3', ]

    measuresVaryingCor['2*CNT2/CNT3', ] <-
        2 * measuresVaryingCor['CNT2', ] /
        measuresVaryingCor['CNT3', ]

    measuresVaryingCor[, 'Measure'] <- rownames(measuresVaryingCor)
    measuresVaryingCor[, 'Measure'] <- ordered(measuresVaryingCor[, 'Measure'],
                                                         levels = measureLabels)

    measuresVaryingCor <- reshape(measuresVaryingCor, direction = 'long',
                                            varying = grep('Value', names(measuresVaryingCor)),
                                            timevar = 'Expectation')

    measuresVaryingCor[, 'Expectation'] <- names(listReFerTables)[measuresVaryingCor[, 'Expectation']]

    measuresVaryingCor[abs(measuresVaryingCor[, 'Value']) == 'Inf', 'Value'] <- NaN
    measuresVaryingCor[(!is.nan(measuresVaryingCor[, 'Value'])) &
                                     (abs(measuresVaryingCor[, 'Value']) < (10 * .Machine$double.eps)), 'Value'] <- 0


figureData <- filter(measuresVaryingCor,
                     grepl('CNT./', measuresVaryingCor[, 'Measure']))

varCorFig <- ggplot(figureData, aes(x = as.numeric(Expectation), y = Value)) +
    geom_point(aes(colour = Measure), size = 0.7) +
    facet_grid(Measure ~ ., scales = 'free_y') +
    xlab('Common expectation of the product') +
    scale_colour_brewer(palette="Dark2")

ggsave(filename = '../output/PR3constantAbsCorRatios.png', plot = varCorFig,
       width = 7, height = 5)

figureData <- filter(measuresVaryingCor,
                     grepl('^N?CNT.$', measuresVaryingCor[, 'Measure']))

measCorFig <- ggplot(figureData, aes(x = as.numeric(Expectation), y = Value)) +
    geom_point(aes(colour = Measure), size = 0.7) +
    facet_grid(Measure ~ ., scales = 'free_y') +
    xlab('Common expectation of the product') +
    scale_colour_brewer(palette="Dark2")

ggsave(filename = '../output/PR3constantAbsCorMeasures.png', plot = measCorFig,
       width = 7, height = 5)





################################################################################
# #  Maximal but 1 varying 1 cor from 1 in refTable (-1 in ferTable) to -1 (1)
################################################################################

nCors <- 41

listReFerTables <- list()
products <- seq(from = 0.5, to = 0, length.out = nCors)

for (iiProduct in seq_along(products)) {
    listReFerTables[[iiProduct]] <- generateRefFerTables(product = products[iiProduct])
    names(listReFerTables)[iiProduct] <- sprintf('%.3f', ((4 * products[iiProduct]) - 1))

}


measuresVaryingCor <- list()


    for (jjCor in seq(nCors)) {

      print(jjCor)

        refTable <- listReFerTables[[1]][['refTable']]
        negTable <- listReFerTables[[1]][['ferTable']]
        chgTable <- listReFerTables[[jjCor]][['ferTable']]

        measuresVaryingCor[[jjCor]] <- sapply(CbDMeasures(
            generatePRBoxLikeSystem(PRBoxStructure = generatePRBoxNStructure(n = 3L),
                           posTable = refTable,
                           negTable = negTable,
                           chgTable = chgTable)),
            function (y) y[1])
    }

    names(measuresVaryingCor) <- paste0('Value.', seq(nCors))

    measuresVaryingCor <- as.data.frame(measuresVaryingCor)

    measuresVaryingCor['CNT1/CNT2', ] <-
        measuresVaryingCor['CNT1', ] /
        measuresVaryingCor['CNT2', ]

    measuresVaryingCor['CNT1/CNT3', ] <-
        measuresVaryingCor['CNT1', ] /
        measuresVaryingCor['CNT3', ]

    measuresVaryingCor['2*CNT2/CNT3', ] <-
        2 * measuresVaryingCor['CNT2', ] /
        measuresVaryingCor['CNT3', ]

    measuresVaryingCor[, 'Measure'] <- rownames(measuresVaryingCor)
    measuresVaryingCor[, 'Measure'] <- ordered(measuresVaryingCor[, 'Measure'],
                                                         levels = measureLabels)


    ### SCATTERPLOT


    measuresVaryingCor <- reshape(measuresVaryingCor, direction = 'long',
                                            varying = grep('Value', names(measuresVaryingCor)),
                                            timevar = 'Expectation')

    measuresVaryingCor[, 'Expectation'] <- names(listReFerTables)[measuresVaryingCor[, 'Expectation']]

    measuresVaryingCor[abs(measuresVaryingCor[, 'Value']) == 'Inf', 'Value'] <- NaN
    measuresVaryingCor[(!is.nan(measuresVaryingCor[, 'Value'])) &
                                     (abs(measuresVaryingCor[, 'Value']) < (10 * .Machine$double.eps)), 'Value'] <- 0


figureData <- filter(measuresVaryingCor,
                     grepl('CNT./', measuresVaryingCor[, 'Measure']))

varCorFig <- ggplot(figureData, aes(x = as.numeric(Expectation), y = round(Value, 10))) +
    geom_point(aes(colour = Measure), size = 0.7) +
    facet_grid(Measure ~ ., scales = 'free_y') +
    xlab('Expectation of the product in varying context') +
    ylab('Value') +
    scale_colour_brewer(palette="Dark2")

ggsave(filename = '../output/PR3Varying1Ratios.png', plot = varCorFig,
       width = 7, height = 5)

figureData <- filter(measuresVaryingCor,
                     grepl('^N?CNT.$', measuresVaryingCor[, 'Measure']))

measCorFig <- ggplot(figureData, aes(x = as.numeric(Expectation), y = round(Value, 10))) +
    geom_point(aes(colour = Measure), size = 0.7) +
    facet_grid(Measure ~ ., scales = 'free_y') +
    xlab('Expectation of the product in varying context') +
    ylab('Value') +
    scale_colour_brewer(palette="Dark2")

ggsave(filename = '../output/PR3Varying1Measures.png', plot = measCorFig,
       width = 7, height = 5)

