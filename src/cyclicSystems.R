################################################################################
# # Measures of contextuality and noncontextuality
# #
# # Computations of the measures from
# # arXiv:1903.07170v5
# # for several cyclic systems
# #
# # Víctor H Cervantes
# # 2019
################################################################################

library('ggplot2')
library('dplyr')

source('./ContextualityMeasures.R')


cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

measureLabels <- c('Delta', 'mA', 'CNT1', 'CNT2', 'CNT3', 'NCNT2',
                   'Delta/CNT2', '-Delta/NCNT2', 'mADelta/NCNT2', 'minDelta/NCNT2', 'CNT1/CNT2', '2*CNT2/CNT3')

################################################################################
# # Generate cyclic system of rank n
################################################################################

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



generateRefFerTables <- function (product = 0.5, margin1 = 0.5, margin2 = 0.5) {
    if (product > margin1 || product > margin2) {
        stop('product cell needs to have a smaller probability than each margin')
    }
    if (product > 1 || product < 0 || margin1 > 1 || margin1 < 0 || margin2 > 1 || margin2 < 0) {
        stop('product and margins need to be probabilities')
    }

    refTable <- as.table(matrix(c(product, margin1 - product, margin2 - product, 1 - margin1 - margin2 + product), ncol = 2))
    ferTable <- refTable[2:1, ]

    output <- list(refTable = refTable, ferTable = ferTable)

    return(output)
}



################################################################################
# #  Maximally contextual systems CbD measures
################################################################################

maximallyContextualSystems <- list()

maxRank <- 6 # actually maxRank - 1

for (iiRank in seq(maxRank)) {
    maximallyContextualSystems[[iiRank]] <- sapply(CbDMeasures(cyclicRVSystem(n = iiRank + 1)),
                                                   function (y) y[1])
}

names(maximallyContextualSystems) <- paste0('Value.', seq_along(maximallyContextualSystems) + 1)

maximallyContextualSystems <- as.data.frame(maximallyContextualSystems)

maximallyContextualSystems['CNT1/CNT2', ] <-
    maximallyContextualSystems['CNT1', ] /
    maximallyContextualSystems['CNT2', ]

maximallyContextualSystems['2*CNT2/CNT3', ] <-
    2 * maximallyContextualSystems['CNT2', ] /
    maximallyContextualSystems['CNT3', ]

maximallyContextualSystems['Delta/CNT2', ] <-
    maximallyContextualSystems['Delta', ] /
    maximallyContextualSystems['CNT2', ]

maximallyContextualSystems['Delta/NCNT2', ] <-
    maximallyContextualSystems['Delta', ] /
    maximallyContextualSystems['NCNT2', ]


maximallyContextualSystems[, 'Measure'] <- rownames(maximallyContextualSystems)
maximallyContextualSystems[, 'Measure'] <- ordered(maximallyContextualSystems[, 'Measure'],
                                                   levels = measureLabels)


maximallyContextualSystemsLong <- reshape(maximallyContextualSystems, direction = 'long',
                                          varying = grep('Value', names(maximallyContextualSystems)),
                                          timevar = 'Rank')



figureData <- filter(maximallyContextualSystemsLong,
                     grepl('/C|Delta$', maximallyContextualSystemsLong[, 'Measure']))

maximalFig <- ggplot(figureData, aes(x = Rank, y = Value)) +
    geom_line(aes(colour = Measure, linetype = Measure)) +
    scale_x_continuous(breaks = seq(2, maxRank + 1)) +
    scale_y_continuous(breaks = seq(1, maxRank), limits = c(0, maxRank)) +
    scale_colour_brewer(palette="Dark2")


ggsave(filename = '../output/maximallyContextual.png',
       plot = maximalFig,
       width = 6, height = 4)


################################################################################
# #  Constant |cor| varying cor from 1 in refTable (-1 in ferTable) to -1 (1)
################################################################################

nCors <- 35

listReFerTables <- list()
products <- seq(from = 0.5, to = 0, length.out = nCors)

for (iiProduct in seq_along(products)) {
    listReFerTables[[iiProduct]] <- generateRefFerTables(product = products[iiProduct])
    names(listReFerTables)[iiProduct] <- sprintf('%.3f', ((4 * products[iiProduct]) - 1))

}


measuresVaryingCor <- list()

for (iiRank in seq(maxRank)) {

    measuresVaryingCor[[iiRank]] <- list()
    for (jjCor in seq(nCors)) {

        refTable <- listReFerTables[[jjCor]][['refTable']]
        if ((iiRank %% 2) == 0) {
            negTable <- listReFerTables[[jjCor]][['refTable']]
        } else {
            negTable <- listReFerTables[[jjCor]][['ferTable']]
        }

        measuresVaryingCor[[iiRank]][[jjCor]] <- sapply(CbDMeasures(
            cyclicRVSystem(n = iiRank + 1,
                           refTable = refTable,
                           negTable = negTable)),
            function (y) y[1])
    }

    names(measuresVaryingCor[[iiRank]]) <- paste0('Value.', seq(nCors))

    measuresVaryingCor[[iiRank]] <- as.data.frame(measuresVaryingCor[[iiRank]])

    measuresVaryingCor[[iiRank]]['CNT1/CNT2', ] <-
        measuresVaryingCor[[iiRank]]['CNT1', ] /
        measuresVaryingCor[[iiRank]]['CNT2', ]

    measuresVaryingCor[[iiRank]]['2*CNT2/CNT3', ] <-
        2 * measuresVaryingCor[[iiRank]]['CNT2', ] /
        measuresVaryingCor[[iiRank]]['CNT3', ]

    measuresVaryingCor[[iiRank]]['Delta/CNT2', ] <-
        measuresVaryingCor[[iiRank]]['Delta', ] /
        measuresVaryingCor[[iiRank]]['CNT2', ]

    measuresVaryingCor[[iiRank]]['-Delta/NCNT2', ] <-
        -measuresVaryingCor[[iiRank]]['Delta', ] /
        measuresVaryingCor[[iiRank]]['NCNT2', ]

    measuresVaryingCor[[iiRank]]['mADelta/NCNT2', ] <-
        measuresVaryingCor[[iiRank]]['mA', ] /
        measuresVaryingCor[[iiRank]]['NCNT2', ]

    measuresVaryingCor[[iiRank]]['minDelta/NCNT2', ] <-
        apply(measuresVaryingCor[[iiRank]][c('mADelta/NCNT2', '-Delta/NCNT2'), ],
              2, min)


    measuresVaryingCor[[iiRank]][, 'Measure'] <- rownames(measuresVaryingCor[[iiRank]])
    measuresVaryingCor[[iiRank]][, 'Measure'] <- ordered(measuresVaryingCor[[iiRank]][, 'Measure'],
                                                         levels = measureLabels)


    measuresVaryingCor[[iiRank]] <- reshape(measuresVaryingCor[[iiRank]], direction = 'long',
                                            varying = grep('Value', names(measuresVaryingCor[[iiRank]])),
                                            timevar = 'Expectation')

    measuresVaryingCor[[iiRank]][, 'Expectation'] <- names(listReFerTables)[measuresVaryingCor[[iiRank]][, 'Expectation']]

    measuresVaryingCor[[iiRank]][abs(measuresVaryingCor[[iiRank]][, 'Value']) == 'Inf', 'Value'] <- NaN
    measuresVaryingCor[[iiRank]][(!is.nan(measuresVaryingCor[[iiRank]][, 'Value'])) &
                                     (abs(measuresVaryingCor[[iiRank]][, 'Value']) < (10 * .Machine$double.eps)), 'Value'] <- 0

    measuresVaryingCor[[iiRank]][, 'Rank'] <- iiRank + 1

}

measuresVaryingCorLong <- do.call('rbind', measuresVaryingCor)

figureData <- filter(measuresVaryingCorLong,
                     grepl('CNT./|^Delta|min', measuresVaryingCorLong[, 'Measure']))

varCorFig <- ggplot(figureData, aes(x = as.numeric(Expectation), y = Value)) +
    geom_point(aes(colour = Measure), size = 0.7) +
    facet_grid(Measure ~ Rank, scales = 'free_y') +
    xlab('Common expectation of the product') +
    scale_colour_brewer(palette="Dark2")

ggsave(filename = '../output/constantAbsCorRatios.png', plot = varCorFig,
       width = 11, height = 6.5)

figureData <- filter(measuresVaryingCorLong,
                     grepl('Delta$|^.?CNT.$', measuresVaryingCorLong[, 'Measure']))

measCorFig <- ggplot(figureData, aes(x = as.numeric(Expectation), y = Value)) +
    geom_point(aes(colour = Measure), size = 0.7) +
    facet_grid(Measure ~ Rank, scales = 'free_y') +
    xlab('Common expectation of the product') +
    scale_colour_brewer(palette="Dark2")

ggsave(filename = '../output/constantAbsCorMeasures.png', plot = measCorFig,
       width = 11, height = 6.5)


################################################################################
# #  Consistently connected. Varying one context from perfect contextual system
################################################################################

nCors <- 35

listReFerTables <- list()
products <- seq(from = 0.5, to = 0, length.out = nCors)

for (iiProduct in seq_along(products)) {
    listReFerTables[[iiProduct]] <- generateRefFerTables(product = products[iiProduct])
    names(listReFerTables)[iiProduct] <- sprintf('%.3f', ((4 * products[iiProduct]) - 1))

}


measuresVaryingCor <- list()

for (iiRank in seq(maxRank)) {
    measuresVaryingCor[[iiRank]] <- list()
    for (jjCor in seq(nCors)) {

        refTable <- listReFerTables[[jjCor]][[1]]
        if (iiRank %% 2) {
            negTable <- listReFerTables[[jjCor]][[1]]
        } else {
            negTable <- listReFerTables[[jjCor]][['ferTable']]
        }

        measuresVaryingCor[[iiRank]][[jjCor]] <- sapply(CbDMeasures(
            cyclicRVSystem(n = iiRank + 1,
                           refTable = refTable,
                           negTable = negTable)),
            function (y) y[1])
    }

    names(measuresVaryingCor[[iiRank]]) <- paste0('Value.', seq(nCors))

    measuresVaryingCor[[iiRank]] <- as.data.frame(measuresVaryingCor[[iiRank]])

    measuresVaryingCor[[iiRank]]['CNT1/CNT2', ] <-
        measuresVaryingCor[[iiRank]]['CNT1', ] /
        measuresVaryingCor[[iiRank]]['CNT2', ]

    measuresVaryingCor[[iiRank]]['2*CNT2/CNT3', ] <-
        2 * measuresVaryingCor[[iiRank]]['CNT2', ] /
        measuresVaryingCor[[iiRank]]['CNT3', ]

    measuresVaryingCor[[iiRank]]['Delta/CNT2', ] <-
        measuresVaryingCor[[iiRank]]['Delta', ] /
        measuresVaryingCor[[iiRank]]['CNT2', ]

    measuresVaryingCor[[iiRank]]['-Delta/NCNT2', ] <-
        -measuresVaryingCor[[iiRank]]['Delta', ] /
        measuresVaryingCor[[iiRank]]['NCNT2', ]

    measuresVaryingCor[[iiRank]]['mADelta/NCNT2', ] <-
        measuresVaryingCor[[iiRank]]['mA', ] /
        measuresVaryingCor[[iiRank]]['NCNT2', ]

    measuresVaryingCor[[iiRank]]['minDelta/NCNT2', ] <-
        apply(measuresVaryingCor[[iiRank]][c('mADelta/NCNT2', '-Delta/NCNT2'), ],
              2, min)


    measuresVaryingCor[[iiRank]][, 'Measure'] <- rownames(measuresVaryingCor[[iiRank]])
    measuresVaryingCor[[iiRank]][, 'Measure'] <- ordered(measuresVaryingCor[[iiRank]][, 'Measure'],
                                                         levels = measureLabels)


    measuresVaryingCor[[iiRank]] <- reshape(measuresVaryingCor[[iiRank]], direction = 'long',
                                            varying = grep('Value', names(measuresVaryingCor[[iiRank]])),
                                            timevar = 'Expectation')

    measuresVaryingCor[[iiRank]][, 'Expectation'] <- names(listReFerTables)[measuresVaryingCor[[iiRank]][, 'Expectation']]

    measuresVaryingCor[[iiRank]][abs(measuresVaryingCor[[iiRank]][, 'Value']) == 'Inf', 'Value'] <- NaN
    measuresVaryingCor[[iiRank]][(!is.nan(measuresVaryingCor[[iiRank]][, 'Value'])) &
                                     (abs(measuresVaryingCor[[iiRank]][, 'Value']) < (10 * .Machine$double.eps)), 'Value'] <- 0

    measuresVaryingCor[[iiRank]][, 'Rank'] <- iiRank + 1

}

measuresVaryingCorLong <- do.call('rbind', measuresVaryingCor)

figureData <- filter(measuresVaryingCorLong,
                     grepl('CNT./|^Delta|min', measuresVaryingCorLong[, 'Measure']))

varCorFig <- ggplot(figureData, aes(x = as.numeric(Expectation), y = Value)) +
    geom_point(aes(colour = Measure), size = 0.7) +
    facet_grid(Measure ~ Rank, scales = 'free_y') +
    xlab('Expectation of the product in varying context') +
    scale_colour_brewer(palette="Dark2")

ggsave(filename = '../output/consistentlyConnectedVarying1Ratios.png', plot = varCorFig,
       width = 11, height = 6.5)

figureData <- filter(measuresVaryingCorLong,
                     grepl('Delta$|^.?CNT.$', measuresVaryingCorLong[, 'Measure']))

measCorFig <- ggplot(figureData, aes(x = as.numeric(Expectation), y = Value)) +
    geom_point(aes(colour = Measure), size = 0.7) +
    facet_grid(Measure ~ Rank, scales = 'free_y') +
    xlab('Expectation of the product in varying context') +
    scale_colour_brewer(palette="Dark2")

ggsave(filename = '../output/consistentlyConnectedVarying1Measures.png', plot = measCorFig,
       width = 11, height = 6.5)



################################################################################
# #  Inconsistently connected systems - Margins .3 and .7 a
################################################################################

nCors <- 35

listReFerTables <- list()
margin1 <- 0.3
margin2 <- 0.7
products <- seq(from = 0.3, to = 0, length.out = nCors)

for (iiProduct in seq_along(products)) {
    listReFerTables[[iiProduct]] <- generateRefFerTables(product = products[iiProduct],
                                                         margin1 = margin1,
                                                         margin2 = margin2)
    expectation <- (4 * products[iiProduct]) - (2 * (margin1 + margin2)) + 1
    names(listReFerTables)[iiProduct] <- sprintf('%.3f', expectation)

}


measuresVaryingCor <- list()

for (iiRank in seq(maxRank)) {
    measuresVaryingCor[[iiRank]] <- list()
    for (jjCor in seq(nCors)) {

        refTable <- listReFerTables[[jjCor]][['refTable']]
        if (iiRank %% 2) {
            negTable <- listReFerTables[[jjCor]][['refTable']]
        } else {
            negTable <- listReFerTables[[jjCor]][['ferTable']]
        }

        measuresVaryingCor[[iiRank]][[jjCor]] <- sapply(CbDMeasures(
            cyclicRVSystem(n = iiRank + 1,
                           refTable = refTable,
                           negTable = negTable)),
            function (y) y[1])
    }

    names(measuresVaryingCor[[iiRank]]) <- paste0('Value.', seq(nCors))

    measuresVaryingCor[[iiRank]] <- as.data.frame(measuresVaryingCor[[iiRank]])

    measuresVaryingCor[[iiRank]]['CNT1/CNT2', ] <-
        measuresVaryingCor[[iiRank]]['CNT1', ] /
        measuresVaryingCor[[iiRank]]['CNT2', ]

    measuresVaryingCor[[iiRank]]['2*CNT2/CNT3', ] <-
        2 * measuresVaryingCor[[iiRank]]['CNT2', ] /
        measuresVaryingCor[[iiRank]]['CNT3', ]

    measuresVaryingCor[[iiRank]]['Delta/CNT2', ] <-
        measuresVaryingCor[[iiRank]]['Delta', ] /
        measuresVaryingCor[[iiRank]]['CNT2', ]

    measuresVaryingCor[[iiRank]]['-Delta/NCNT2', ] <-
        -measuresVaryingCor[[iiRank]]['Delta', ] /
        measuresVaryingCor[[iiRank]]['NCNT2', ]

    measuresVaryingCor[[iiRank]]['mADelta/NCNT2', ] <-
        measuresVaryingCor[[iiRank]]['mA', ] /
        measuresVaryingCor[[iiRank]]['NCNT2', ]

    measuresVaryingCor[[iiRank]]['minDelta/NCNT2', ] <-
        apply(measuresVaryingCor[[iiRank]][c('mADelta/NCNT2', '-Delta/NCNT2'), ],
              2, min)


    measuresVaryingCor[[iiRank]][, 'Measure'] <- rownames(measuresVaryingCor[[iiRank]])
    measuresVaryingCor[[iiRank]][, 'Measure'] <- ordered(measuresVaryingCor[[iiRank]][, 'Measure'],
                                                         levels = measureLabels)


    measuresVaryingCor[[iiRank]] <- reshape(measuresVaryingCor[[iiRank]], direction = 'long',
                                            varying = grep('Value', names(measuresVaryingCor[[iiRank]])),
                                            timevar = 'Expectation')

    measuresVaryingCor[[iiRank]][, 'Expectation'] <- names(listReFerTables)[measuresVaryingCor[[iiRank]][, 'Expectation']]

    measuresVaryingCor[[iiRank]][abs(measuresVaryingCor[[iiRank]][, 'Value']) == 'Inf', 'Value'] <- NaN
    measuresVaryingCor[[iiRank]][(!is.nan(measuresVaryingCor[[iiRank]][, 'Value'])) &
                                     (abs(measuresVaryingCor[[iiRank]][, 'Value']) < (10 * .Machine$double.eps)), 'Value'] <- 0

    measuresVaryingCor[[iiRank]][, 'Rank'] <- iiRank + 1

}

measuresVaryingCorLong <- do.call('rbind', measuresVaryingCor)

figureData <- filter(measuresVaryingCorLong,
                     grepl('CNT./|^Delta|min', measuresVaryingCorLong[, 'Measure']))

varCorFig <- ggplot(figureData, aes(x = as.numeric(Expectation), y = Value)) +
    geom_point(aes(colour = Measure), size = 0.7) +
    facet_grid(Measure ~ Rank, scales = 'free_y') +
    xlab('Common expectation of the product') +
    scale_colour_brewer(palette="Dark2")

ggsave(filename = '../output/inconsistent3v7cyclicRatios.png', plot = varCorFig,
       width = 11, height = 6.5)

figureData <- filter(measuresVaryingCorLong,
                     grepl('Delta$|^.?CNT.$', measuresVaryingCorLong[, 'Measure']))

measCorFig <- ggplot(figureData, aes(x = as.numeric(Expectation), y = Value)) +
    geom_point(aes(colour = Measure), size = 0.7) +
    facet_grid(Measure ~ Rank, scales = 'free_y') +
    xlab('Common expectation of the product') +
    scale_colour_brewer(palette="Dark2")

ggsave(filename = '../output/inconsistent3v7cyclicMeasures.png', plot = measCorFig,
       width = 11, height = 6.5)



################################################################################
# #  Inconsistently connected systems - Margins .3 and .7 b
################################################################################

nCors <- 35

listReFerTables <- list()
margin1 <- 0.3
margin2 <- 0.3
products <- seq(from = 0.3, to = 0, length.out = nCors)

for (iiProduct in seq_along(products)) {
    listReFerTables[[iiProduct]] <- generateRefFerTables(product = products[iiProduct],
                                                         margin1 = margin1,
                                                         margin2 = margin2)
    expectation <- (4 * products[iiProduct]) - (2 * (margin1 + margin2)) + 1
    names(listReFerTables)[iiProduct] <- sprintf('%.3f', expectation)

}


measuresVaryingCor <- list()

for (iiRank in seq(maxRank)) {
    measuresVaryingCor[[iiRank]] <- list()
    for (jjCor in seq(nCors)) {

        refTable <- listReFerTables[[jjCor]][['refTable']]
        if (iiRank %% 2) {
            negTable <- listReFerTables[[jjCor]][['refTable']]
        } else {
            negTable <- listReFerTables[[jjCor]][['ferTable']]
        }

        measuresVaryingCor[[iiRank]][[jjCor]] <- sapply(CbDMeasures(
            cyclicRVSystem(n = iiRank + 1,
                           refTable = refTable,
                           negTable = negTable)),
            function (y) y[1])
    }

    names(measuresVaryingCor[[iiRank]]) <- paste0('Value.', seq(nCors))

    measuresVaryingCor[[iiRank]] <- as.data.frame(measuresVaryingCor[[iiRank]])

    measuresVaryingCor[[iiRank]]['CNT1/CNT2', ] <-
        measuresVaryingCor[[iiRank]]['CNT1', ] /
        measuresVaryingCor[[iiRank]]['CNT2', ]

    measuresVaryingCor[[iiRank]]['2*CNT2/CNT3', ] <-
        2 * measuresVaryingCor[[iiRank]]['CNT2', ] /
        measuresVaryingCor[[iiRank]]['CNT3', ]

    measuresVaryingCor[[iiRank]]['Delta/CNT2', ] <-
        measuresVaryingCor[[iiRank]]['Delta', ] /
        measuresVaryingCor[[iiRank]]['CNT2', ]

    measuresVaryingCor[[iiRank]]['-Delta/NCNT2', ] <-
        -measuresVaryingCor[[iiRank]]['Delta', ] /
        measuresVaryingCor[[iiRank]]['NCNT2', ]

    measuresVaryingCor[[iiRank]]['mADelta/NCNT2', ] <-
        measuresVaryingCor[[iiRank]]['mA', ] /
        measuresVaryingCor[[iiRank]]['NCNT2', ]

    measuresVaryingCor[[iiRank]]['minDelta/NCNT2', ] <-
        apply(measuresVaryingCor[[iiRank]][c('mADelta/NCNT2', '-Delta/NCNT2'), ],
              2, min)


    measuresVaryingCor[[iiRank]][, 'Measure'] <- rownames(measuresVaryingCor[[iiRank]])
    measuresVaryingCor[[iiRank]][, 'Measure'] <- ordered(measuresVaryingCor[[iiRank]][, 'Measure'],
                                                         levels = measureLabels)


    measuresVaryingCor[[iiRank]] <- reshape(measuresVaryingCor[[iiRank]], direction = 'long',
                                            varying = grep('Value', names(measuresVaryingCor[[iiRank]])),
                                            timevar = 'Expectation')

    measuresVaryingCor[[iiRank]][, 'Expectation'] <- names(listReFerTables)[measuresVaryingCor[[iiRank]][, 'Expectation']]

    measuresVaryingCor[[iiRank]][abs(measuresVaryingCor[[iiRank]][, 'Value']) == 'Inf', 'Value'] <- NaN
    measuresVaryingCor[[iiRank]][(!is.nan(measuresVaryingCor[[iiRank]][, 'Value'])) &
                                     (abs(measuresVaryingCor[[iiRank]][, 'Value']) < (10 * .Machine$double.eps)), 'Value'] <- 0

    measuresVaryingCor[[iiRank]][, 'Rank'] <- iiRank + 1

}

measuresVaryingCorLong <- do.call('rbind', measuresVaryingCor)

figureData <- filter(measuresVaryingCorLong,
                     grepl('CNT./|^Delta|min', measuresVaryingCorLong[, 'Measure']))

varCorFig <- ggplot(figureData, aes(x = as.numeric(Expectation), y = Value)) +
    geom_point(aes(colour = Measure), size = 0.7) +
    facet_grid(Measure ~ Rank, scales = 'free_y') +
    xlab('Common expectation of the product') +
    scale_colour_brewer(palette="Dark2")

ggsave(filename = '../output/inconsistent3v7cyclic3Ratios.png', plot = varCorFig,
       width = 11, height = 6.5)

figureData <- filter(measuresVaryingCorLong,
                     grepl('Delta$|^.?CNT.$', measuresVaryingCorLong[, 'Measure']))

measCorFig <- ggplot(figureData, aes(x = as.numeric(Expectation), y = Value)) +
    geom_point(aes(colour = Measure), size = 0.7) +
    facet_grid(Measure ~ Rank, scales = 'free_y') +
    xlab('Common expectation of the product') +
    scale_colour_brewer(palette="Dark2")

ggsave(filename = '../output/inconsistent3v7cyclic3Measures.png', plot = measCorFig,
       width = 11, height = 6.5)





################################################################################
# #  Inconsistently connected systems - Margins .53 and .47 a
################################################################################

nCors <- 35

listReFerTables <- list()
margin1 <- 0.53
margin2 <- 0.47
products <- seq(from = 0.3, to = 0, length.out = nCors)

for (iiProduct in seq_along(products)) {
    listReFerTables[[iiProduct]] <- generateRefFerTables(product = products[iiProduct],
                                                         margin1 = margin1,
                                                         margin2 = margin2)
    expectation <- (4 * products[iiProduct]) - (2 * (margin1 + margin2)) + 1
    names(listReFerTables)[iiProduct] <- sprintf('%.3f', expectation)

}


measuresVaryingCor <- list()

for (iiRank in seq(maxRank)) {
    measuresVaryingCor[[iiRank]] <- list()
    for (jjCor in seq(nCors)) {

        refTable <- listReFerTables[[jjCor]][['refTable']]
        if (iiRank %% 2) {
            negTable <- listReFerTables[[jjCor]][['refTable']]
        } else {
            negTable <- listReFerTables[[jjCor]][['ferTable']]
        }

        measuresVaryingCor[[iiRank]][[jjCor]] <- sapply(CbDMeasures(
            cyclicRVSystem(n = iiRank + 1,
                           refTable = refTable,
                           negTable = negTable)),
            function (y) y[1])
    }


    names(measuresVaryingCor[[iiRank]]) <- paste0('Value.', seq(nCors))

    measuresVaryingCor[[iiRank]] <- as.data.frame(measuresVaryingCor[[iiRank]])

    measuresVaryingCor[[iiRank]]['CNT1/CNT2', ] <-
        measuresVaryingCor[[iiRank]]['CNT1', ] /
        measuresVaryingCor[[iiRank]]['CNT2', ]

    measuresVaryingCor[[iiRank]]['2*CNT2/CNT3', ] <-
        2 * measuresVaryingCor[[iiRank]]['CNT2', ] /
        measuresVaryingCor[[iiRank]]['CNT3', ]

    measuresVaryingCor[[iiRank]]['Delta/CNT2', ] <-
        measuresVaryingCor[[iiRank]]['Delta', ] /
        measuresVaryingCor[[iiRank]]['CNT2', ]

    measuresVaryingCor[[iiRank]]['-Delta/NCNT2', ] <-
        -measuresVaryingCor[[iiRank]]['Delta', ] /
        measuresVaryingCor[[iiRank]]['NCNT2', ]

    measuresVaryingCor[[iiRank]]['mADelta/NCNT2', ] <-
        measuresVaryingCor[[iiRank]]['mA', ] /
        measuresVaryingCor[[iiRank]]['NCNT2', ]

    measuresVaryingCor[[iiRank]]['minDelta/NCNT2', ] <-
        apply(measuresVaryingCor[[iiRank]][c('mADelta/NCNT2', '-Delta/NCNT2'), ],
              2, min)


    measuresVaryingCor[[iiRank]][, 'Measure'] <- rownames(measuresVaryingCor[[iiRank]])
    measuresVaryingCor[[iiRank]][, 'Measure'] <- ordered(measuresVaryingCor[[iiRank]][, 'Measure'],
                                                         levels = measureLabels)


    measuresVaryingCor[[iiRank]] <- reshape(measuresVaryingCor[[iiRank]], direction = 'long',
                                            varying = grep('Value', names(measuresVaryingCor[[iiRank]])),
                                            timevar = 'Expectation')

    measuresVaryingCor[[iiRank]][, 'Expectation'] <- names(listReFerTables)[measuresVaryingCor[[iiRank]][, 'Expectation']]

    measuresVaryingCor[[iiRank]][abs(measuresVaryingCor[[iiRank]][, 'Value']) == 'Inf', 'Value'] <- NaN
    measuresVaryingCor[[iiRank]][(!is.nan(measuresVaryingCor[[iiRank]][, 'Value'])) &
                                     (abs(measuresVaryingCor[[iiRank]][, 'Value']) < (10 * .Machine$double.eps)), 'Value'] <- 0

    measuresVaryingCor[[iiRank]][, 'Rank'] <- iiRank + 1

}

measuresVaryingCorLong <- do.call('rbind', measuresVaryingCor)

figureData <- filter(measuresVaryingCorLong,
                     grepl('CNT./|^Delta|min', measuresVaryingCorLong[, 'Measure']))

varCorFig <- ggplot(figureData, aes(x = as.numeric(Expectation), y = Value)) +
    geom_point(aes(colour = Measure), size = 0.7) +
    facet_grid(Measure ~ Rank, scales = 'free_y') +
    xlab('Common expectation of the product') +
    scale_colour_brewer(palette="Dark2")

ggsave(filename = '../output/inconsistent53v47cyclicRatios.png', plot = varCorFig,
       width = 11, height = 6.5)

figureData <- filter(measuresVaryingCorLong,
                     grepl('Delta$|^.?CNT.$', measuresVaryingCorLong[, 'Measure']))

measCorFig <- ggplot(figureData, aes(x = as.numeric(Expectation), y = Value)) +
    geom_point(aes(colour = Measure), size = 0.7) +
    facet_grid(Measure ~ Rank, scales = 'free_y') +
    xlab('Common expectation of the product') +
    scale_colour_brewer(palette="Dark2")

ggsave(filename = '../output/inconsistent53v47cyclicMeasures.png', plot = measCorFig,
       width = 11, height = 6.5)





################################################################################
# #  Inconsistently connected systems - Margins .53 and .47 b
################################################################################

nCors <- 35

listReFerTables <- list()
margin1 <- 0.47
margin2 <- 0.47
products <- seq(from = 0.3, to = 0, length.out = nCors)

for (iiProduct in seq_along(products)) {
    listReFerTables[[iiProduct]] <- generateRefFerTables(product = products[iiProduct],
                                                         margin1 = margin1,
                                                         margin2 = margin2)
    expectation <- (4 * products[iiProduct]) - (2 * (margin1 + margin2)) + 1
    names(listReFerTables)[iiProduct] <- sprintf('%.3f', expectation)

}


measuresVaryingCor <- list()

for (iiRank in seq(maxRank)) {
    measuresVaryingCor[[iiRank]] <- list()
    for (jjCor in seq(nCors)) {

        refTable <- listReFerTables[[jjCor]][['refTable']]
        if (iiRank %% 2) {
            negTable <- listReFerTables[[jjCor]][['refTable']]
        } else {
            negTable <- listReFerTables[[jjCor]][['ferTable']]
        }

        measuresVaryingCor[[iiRank]][[jjCor]] <- sapply(CbDMeasures(
            cyclicRVSystem(n = iiRank + 1,
                           refTable = refTable,
                           negTable = negTable)),
            function (y) y[1])
    }

    names(measuresVaryingCor[[iiRank]]) <- paste0('Value.', seq(nCors))

    measuresVaryingCor[[iiRank]] <- as.data.frame(measuresVaryingCor[[iiRank]])

    measuresVaryingCor[[iiRank]]['CNT1/CNT2', ] <-
        measuresVaryingCor[[iiRank]]['CNT1', ] /
        measuresVaryingCor[[iiRank]]['CNT2', ]

    measuresVaryingCor[[iiRank]]['2*CNT2/CNT3', ] <-
        2 * measuresVaryingCor[[iiRank]]['CNT2', ] /
        measuresVaryingCor[[iiRank]]['CNT3', ]

    measuresVaryingCor[[iiRank]]['Delta/CNT2', ] <-
        measuresVaryingCor[[iiRank]]['Delta', ] /
        measuresVaryingCor[[iiRank]]['CNT2', ]

    measuresVaryingCor[[iiRank]]['-Delta/NCNT2', ] <-
        -measuresVaryingCor[[iiRank]]['Delta', ] /
        measuresVaryingCor[[iiRank]]['NCNT2', ]

    measuresVaryingCor[[iiRank]]['mADelta/NCNT2', ] <-
        measuresVaryingCor[[iiRank]]['mA', ] /
        measuresVaryingCor[[iiRank]]['NCNT2', ]

    measuresVaryingCor[[iiRank]]['minDelta/NCNT2', ] <-
        apply(measuresVaryingCor[[iiRank]][c('mADelta/NCNT2', '-Delta/NCNT2'), ],
              2, min)


    measuresVaryingCor[[iiRank]][, 'Measure'] <- rownames(measuresVaryingCor[[iiRank]])
    measuresVaryingCor[[iiRank]][, 'Measure'] <- ordered(measuresVaryingCor[[iiRank]][, 'Measure'],
                                                         levels = measureLabels)


    measuresVaryingCor[[iiRank]] <- reshape(measuresVaryingCor[[iiRank]], direction = 'long',
                                            varying = grep('Value', names(measuresVaryingCor[[iiRank]])),
                                            timevar = 'Expectation')

    measuresVaryingCor[[iiRank]][, 'Expectation'] <- names(listReFerTables)[measuresVaryingCor[[iiRank]][, 'Expectation']]

    measuresVaryingCor[[iiRank]][abs(measuresVaryingCor[[iiRank]][, 'Value']) == 'Inf', 'Value'] <- NaN
    measuresVaryingCor[[iiRank]][(!is.nan(measuresVaryingCor[[iiRank]][, 'Value'])) &
                                     (abs(measuresVaryingCor[[iiRank]][, 'Value']) < (10 * .Machine$double.eps)), 'Value'] <- 0

    measuresVaryingCor[[iiRank]][, 'Rank'] <- iiRank + 1

}

measuresVaryingCorLong <- do.call('rbind', measuresVaryingCor)

figureData <- filter(measuresVaryingCorLong,
                     grepl('CNT./|^Delta|min', measuresVaryingCorLong[, 'Measure']))

varCorFig <- ggplot(figureData, aes(x = as.numeric(Expectation), y = Value)) +
    geom_point(aes(colour = Measure), size = 0.7) +
    facet_grid(Measure ~ Rank, scales = 'free_y') +
    xlab('Common expectation of the product') +
    scale_colour_brewer(palette="Dark2")

ggsave(filename = '../output/inconsistent53v47cyclic3Ratios.png', plot = varCorFig,
       width = 11, height = 6.5)

figureData <- filter(measuresVaryingCorLong,
                     grepl('Delta$|^.?CNT.$', measuresVaryingCorLong[, 'Measure']))

measCorFig <- ggplot(figureData, aes(x = as.numeric(Expectation), y = Value)) +
    geom_point(aes(colour = Measure), size = 0.7) +
    facet_grid(Measure ~ Rank, scales = 'free_y') +
    xlab('Common expectation of the product') +
    scale_colour_brewer(palette="Dark2")

ggsave(filename = '../output/inconsistent53v47cyclic3Measures.png', plot = measCorFig,
       width = 11, height = 6.5)






################################################################################
# #  Maximal but one inconsistently connected systems - Margins .3 and .7 a
################################################################################

nCors <- 35

listReFerTables <- list()
margin1 <- 0.3
margin2 <- 0.7
products <- seq(from = 0.3, to = 0, length.out = nCors)

for (iiProduct in seq_along(products)) {
    listReFerTables[[iiProduct]] <- generateRefFerTables(product = products[iiProduct],
                                                         margin1 = margin1,
                                                         margin2 = margin2)
    expectation <- (4 * products[iiProduct]) - (2 * (margin1 + margin2)) + 1
    names(listReFerTables)[iiProduct] <- sprintf('%.3f', expectation)

}


measuresVaryingCor <- list()

for (iiRank in seq(maxRank)) {
    measuresVaryingCor[[iiRank]] <- list()
    for (jjCor in seq(nCors)) {

        refTable <- listReFerTables[[1]][['refTable']]
        if (iiRank %% 2) {
            negTable <- listReFerTables[[jjCor]][['refTable']]
        } else {
            negTable <- listReFerTables[[jjCor]][['ferTable']]
        }

        measuresVaryingCor[[iiRank]][[jjCor]] <- sapply(CbDMeasures(
            cyclicRVSystem(n = iiRank + 1,
                           refTable = refTable,
                           negTable = negTable)),
            function (y) y[1])
    }

    names(measuresVaryingCor[[iiRank]]) <- paste0('Value.', seq(nCors))

    measuresVaryingCor[[iiRank]] <- as.data.frame(measuresVaryingCor[[iiRank]])

    measuresVaryingCor[[iiRank]]['CNT1/CNT2', ] <-
        measuresVaryingCor[[iiRank]]['CNT1', ] /
        measuresVaryingCor[[iiRank]]['CNT2', ]

    measuresVaryingCor[[iiRank]]['2*CNT2/CNT3', ] <-
        2 * measuresVaryingCor[[iiRank]]['CNT2', ] /
        measuresVaryingCor[[iiRank]]['CNT3', ]

    measuresVaryingCor[[iiRank]]['Delta/CNT2', ] <-
        measuresVaryingCor[[iiRank]]['Delta', ] /
        measuresVaryingCor[[iiRank]]['CNT2', ]

    measuresVaryingCor[[iiRank]]['-Delta/NCNT2', ] <-
        -measuresVaryingCor[[iiRank]]['Delta', ] /
        measuresVaryingCor[[iiRank]]['NCNT2', ]

    measuresVaryingCor[[iiRank]]['mADelta/NCNT2', ] <-
        measuresVaryingCor[[iiRank]]['mA', ] /
        measuresVaryingCor[[iiRank]]['NCNT2', ]

    measuresVaryingCor[[iiRank]]['minDelta/NCNT2', ] <-
        apply(measuresVaryingCor[[iiRank]][c('mADelta/NCNT2', '-Delta/NCNT2'), ],
              2, min)


    measuresVaryingCor[[iiRank]][, 'Measure'] <- rownames(measuresVaryingCor[[iiRank]])
    measuresVaryingCor[[iiRank]][, 'Measure'] <- ordered(measuresVaryingCor[[iiRank]][, 'Measure'],
                                                         levels = measureLabels)


    measuresVaryingCor[[iiRank]] <- reshape(measuresVaryingCor[[iiRank]], direction = 'long',
                                            varying = grep('Value', names(measuresVaryingCor[[iiRank]])),
                                            timevar = 'Expectation')

    measuresVaryingCor[[iiRank]][, 'Expectation'] <- names(listReFerTables)[measuresVaryingCor[[iiRank]][, 'Expectation']]

    measuresVaryingCor[[iiRank]][abs(measuresVaryingCor[[iiRank]][, 'Value']) == 'Inf', 'Value'] <- NaN
    measuresVaryingCor[[iiRank]][(!is.nan(measuresVaryingCor[[iiRank]][, 'Value'])) &
                                     (abs(measuresVaryingCor[[iiRank]][, 'Value']) < (10 * .Machine$double.eps)), 'Value'] <- 0

    measuresVaryingCor[[iiRank]][, 'Rank'] <- iiRank + 1

}

measuresVaryingCorLong <- do.call('rbind', measuresVaryingCor)

figureData <- filter(measuresVaryingCorLong,
                     grepl('CNT./|^Delta|min', measuresVaryingCorLong[, 'Measure']))

varCorFig <- ggplot(figureData, aes(x = as.numeric(Expectation), y = Value)) +
    geom_point(aes(colour = Measure), size = 0.7) +
    facet_grid(Measure ~ Rank, scales = 'free_y') +
    xlab('Expectation of the product in varying context') +
    scale_colour_brewer(palette="Dark2")

ggsave(filename = '../output/inconsistent3v7cyclicVarying1Ratios.png', plot = varCorFig,
       width = 11, height = 6.5)

figureData <- filter(measuresVaryingCorLong,
                     grepl('Delta$|^.?CNT.$', measuresVaryingCorLong[, 'Measure']))

measCorFig <- ggplot(figureData, aes(x = as.numeric(Expectation), y = Value)) +
    geom_point(aes(colour = Measure), size = 0.7) +
    facet_grid(Measure ~ Rank, scales = 'free_y') +
    xlab('Expectation of the product in varying context') +
    scale_colour_brewer(palette="Dark2")

ggsave(filename = '../output/inconsistent3v7cyclicVarying1Measures.png', plot = measCorFig,
       width = 11, height = 6.5)



################################################################################
# #  Maximal but one Inconsistently connected systems - Margins .3 and .7 b
################################################################################

nCors <- 35

listReFerTables <- list()
margin1 <- 0.3
margin2 <- 0.3
products <- seq(from = 0.3, to = 0, length.out = nCors)

for (iiProduct in seq_along(products)) {
    listReFerTables[[iiProduct]] <- generateRefFerTables(product = products[iiProduct],
                                                         margin1 = margin1,
                                                         margin2 = margin2)
    expectation <- (4 * products[iiProduct]) - (2 * (margin1 + margin2)) + 1
    names(listReFerTables)[iiProduct] <- sprintf('%.3f', expectation)

}


measuresVaryingCor <- list()

for (iiRank in seq(maxRank)) {
    measuresVaryingCor[[iiRank]] <- list()
    for (jjCor in seq(nCors)) {

        refTable <- listReFerTables[[1]][[1]]
        if (iiRank %% 2) {
            negTable <- listReFerTables[[jjCor]][[1]]
        } else {
            negTable <- listReFerTables[[jjCor]][['ferTable']]
        }

        measuresVaryingCor[[iiRank]][[jjCor]] <- sapply(CbDMeasures(
            cyclicRVSystem(n = iiRank + 1,
                           refTable = refTable,
                           negTable = negTable)),
            function (y) y[1])
    }

    names(measuresVaryingCor[[iiRank]]) <- paste0('Value.', seq(nCors))

    measuresVaryingCor[[iiRank]] <- as.data.frame(measuresVaryingCor[[iiRank]])

    measuresVaryingCor[[iiRank]]['CNT1/CNT2', ] <-
        measuresVaryingCor[[iiRank]]['CNT1', ] /
        measuresVaryingCor[[iiRank]]['CNT2', ]

    measuresVaryingCor[[iiRank]]['2*CNT2/CNT3', ] <-
        2 * measuresVaryingCor[[iiRank]]['CNT2', ] /
        measuresVaryingCor[[iiRank]]['CNT3', ]

    measuresVaryingCor[[iiRank]]['Delta/CNT2', ] <-
        measuresVaryingCor[[iiRank]]['Delta', ] /
        measuresVaryingCor[[iiRank]]['CNT2', ]

    measuresVaryingCor[[iiRank]]['-Delta/NCNT2', ] <-
        -measuresVaryingCor[[iiRank]]['Delta', ] /
        measuresVaryingCor[[iiRank]]['NCNT2', ]

    measuresVaryingCor[[iiRank]]['mADelta/NCNT2', ] <-
        measuresVaryingCor[[iiRank]]['mA', ] /
        measuresVaryingCor[[iiRank]]['NCNT2', ]

    measuresVaryingCor[[iiRank]]['minDelta/NCNT2', ] <-
        apply(measuresVaryingCor[[iiRank]][c('mADelta/NCNT2', '-Delta/NCNT2'), ],
              2, min)


    measuresVaryingCor[[iiRank]][, 'Measure'] <- rownames(measuresVaryingCor[[iiRank]])
    measuresVaryingCor[[iiRank]][, 'Measure'] <- ordered(measuresVaryingCor[[iiRank]][, 'Measure'],
                                                         levels = measureLabels)


    measuresVaryingCor[[iiRank]] <- reshape(measuresVaryingCor[[iiRank]], direction = 'long',
                                            varying = grep('Value', names(measuresVaryingCor[[iiRank]])),
                                            timevar = 'Expectation')

    measuresVaryingCor[[iiRank]][, 'Expectation'] <- names(listReFerTables)[measuresVaryingCor[[iiRank]][, 'Expectation']]

    measuresVaryingCor[[iiRank]][abs(measuresVaryingCor[[iiRank]][, 'Value']) == 'Inf', 'Value'] <- NaN
    measuresVaryingCor[[iiRank]][(!is.nan(measuresVaryingCor[[iiRank]][, 'Value'])) &
                                     (abs(measuresVaryingCor[[iiRank]][, 'Value']) < (10 * .Machine$double.eps)), 'Value'] <- 0

    measuresVaryingCor[[iiRank]][, 'Rank'] <- iiRank + 1

}

measuresVaryingCorLong <- do.call('rbind', measuresVaryingCor)

figureData <- filter(measuresVaryingCorLong,
                     grepl('CNT./|^Delta|min', measuresVaryingCorLong[, 'Measure']))

varCorFig <- ggplot(figureData, aes(x = as.numeric(Expectation), y = Value)) +
    geom_point(aes(colour = Measure), size = 0.7) +
    facet_grid(Measure ~ Rank, scales = 'free_y') +
    xlab('Expectation of the product in varying context') +
    scale_colour_brewer(palette="Dark2")

ggsave(filename = '../output/inconsistent3v7cyclic3Varying1Ratios.png', plot = varCorFig,
       width = 11, height = 6.5)

figureData <- filter(measuresVaryingCorLong,
                     grepl('Delta$|^.?CNT.$', measuresVaryingCorLong[, 'Measure']))

measCorFig <- ggplot(figureData, aes(x = as.numeric(Expectation), y = Value)) +
    geom_point(aes(colour = Measure), size = 0.7) +
    facet_grid(Measure ~ Rank, scales = 'free_y') +
    xlab('Expectation of the product in varying context') +
    scale_colour_brewer(palette="Dark2")

ggsave(filename = '../output/inconsistent3v7cyclic3Varying1Measures.png', plot = measCorFig,
       width = 11, height = 6.5)





################################################################################
# #  Maximal but one Inconsistently connected systems - Margins .53 and .47 a
################################################################################

nCors <- 35

listReFerTables <- list()
margin1 <- 0.53
margin2 <- 0.47
products <- seq(from = 0.3, to = 0, length.out = nCors)

for (iiProduct in seq_along(products)) {
    listReFerTables[[iiProduct]] <- generateRefFerTables(product = products[iiProduct],
                                                         margin1 = margin1,
                                                         margin2 = margin2)
    expectation <- (4 * products[iiProduct]) - (2 * (margin1 + margin2)) + 1
    names(listReFerTables)[iiProduct] <- sprintf('%.3f', expectation)

}


measuresVaryingCor <- list()

for (iiRank in seq(maxRank)) {
    measuresVaryingCor[[iiRank]] <- list()
    for (jjCor in seq(nCors)) {

        refTable <- listReFerTables[[1]][[1]]
        if (iiRank %% 2) {
            negTable <- listReFerTables[[jjCor]][[1]]
        } else {
            negTable <- listReFerTables[[jjCor]][['ferTable']]
        }

        measuresVaryingCor[[iiRank]][[jjCor]] <- sapply(CbDMeasures(
            cyclicRVSystem(n = iiRank + 1,
                           refTable = refTable,
                           negTable = negTable)),
            function (y) y[1])
    }

    names(measuresVaryingCor[[iiRank]]) <- paste0('Value.', seq(nCors))

    measuresVaryingCor[[iiRank]] <- as.data.frame(measuresVaryingCor[[iiRank]])

    measuresVaryingCor[[iiRank]]['CNT1/CNT2', ] <-
        measuresVaryingCor[[iiRank]]['CNT1', ] /
        measuresVaryingCor[[iiRank]]['CNT2', ]

    measuresVaryingCor[[iiRank]]['2*CNT2/CNT3', ] <-
        2 * measuresVaryingCor[[iiRank]]['CNT2', ] /
        measuresVaryingCor[[iiRank]]['CNT3', ]

    measuresVaryingCor[[iiRank]]['Delta/CNT2', ] <-
        measuresVaryingCor[[iiRank]]['Delta', ] /
        measuresVaryingCor[[iiRank]]['CNT2', ]

    measuresVaryingCor[[iiRank]]['-Delta/NCNT2', ] <-
        -measuresVaryingCor[[iiRank]]['Delta', ] /
        measuresVaryingCor[[iiRank]]['NCNT2', ]

    measuresVaryingCor[[iiRank]]['mADelta/NCNT2', ] <-
        measuresVaryingCor[[iiRank]]['mA', ] /
        measuresVaryingCor[[iiRank]]['NCNT2', ]

    measuresVaryingCor[[iiRank]]['minDelta/NCNT2', ] <-
        apply(measuresVaryingCor[[iiRank]][c('mADelta/NCNT2', '-Delta/NCNT2'), ],
              2, min)


    measuresVaryingCor[[iiRank]][, 'Measure'] <- rownames(measuresVaryingCor[[iiRank]])
    measuresVaryingCor[[iiRank]][, 'Measure'] <- ordered(measuresVaryingCor[[iiRank]][, 'Measure'],
                                                         levels = measureLabels)


    measuresVaryingCor[[iiRank]] <- reshape(measuresVaryingCor[[iiRank]], direction = 'long',
                                            varying = grep('Value', names(measuresVaryingCor[[iiRank]])),
                                            timevar = 'Expectation')

    measuresVaryingCor[[iiRank]][, 'Expectation'] <- names(listReFerTables)[measuresVaryingCor[[iiRank]][, 'Expectation']]

    measuresVaryingCor[[iiRank]][abs(measuresVaryingCor[[iiRank]][, 'Value']) == 'Inf', 'Value'] <- NaN
    measuresVaryingCor[[iiRank]][(!is.nan(measuresVaryingCor[[iiRank]][, 'Value'])) &
                                     (abs(measuresVaryingCor[[iiRank]][, 'Value']) < (10 * .Machine$double.eps)), 'Value'] <- 0

    measuresVaryingCor[[iiRank]][, 'Rank'] <- iiRank + 1

}

measuresVaryingCorLong <- do.call('rbind', measuresVaryingCor)

figureData <- filter(measuresVaryingCorLong,
                     grepl('CNT./|^Delta|min', measuresVaryingCorLong[, 'Measure']))

varCorFig <- ggplot(figureData, aes(x = as.numeric(Expectation), y = Value)) +
    geom_point(aes(colour = Measure), size = 0.7) +
    facet_grid(Measure ~ Rank, scales = 'free_y') +
    xlab('Expectation of the product in varying context') +
    scale_colour_brewer(palette="Dark2")

ggsave(filename = '../output/inconsistent53v47cyclicVarying1Ratios.png', plot = varCorFig,
       width = 11, height = 6.5)

figureData <- filter(measuresVaryingCorLong,
                     grepl('Delta$|^.?CNT.$', measuresVaryingCorLong[, 'Measure']))

measCorFig <- ggplot(figureData, aes(x = as.numeric(Expectation), y = Value)) +
    geom_point(aes(colour = Measure), size = 0.7) +
    facet_grid(Measure ~ Rank, scales = 'free_y') +
    xlab('Expectation of the product in varying context') +
    scale_colour_brewer(palette="Dark2")

ggsave(filename = '../output/inconsistent53v47cyclicVarying1Measures.png', plot = measCorFig,
       width = 11, height = 6.5)





################################################################################
# #  Maximal but one Inconsistently connected systems - Margins .53 and .47 b
################################################################################

nCors <- 35

listReFerTables <- list()
margin1 <- 0.47
margin2 <- 0.47
products <- seq(from = 0.3, to = 0, length.out = nCors)

for (iiProduct in seq_along(products)) {
    listReFerTables[[iiProduct]] <- generateRefFerTables(product = products[iiProduct],
                                                         margin1 = margin1,
                                                         margin2 = margin2)
    expectation <- (4 * products[iiProduct]) - (2 * (margin1 + margin2)) + 1
    names(listReFerTables)[iiProduct] <- sprintf('%.3f', expectation)

}


measuresVaryingCor <- list()

for (iiRank in seq(maxRank)) {
    measuresVaryingCor[[iiRank]] <- list()
    for (jjCor in seq(nCors)) {

        refTable <- listReFerTables[[1]][[1]]
        if (iiRank %% 2) {
            negTable <- listReFerTables[[jjCor]][[1]]
        } else {
            negTable <- listReFerTables[[jjCor]][['ferTable']]
        }

        measuresVaryingCor[[iiRank]][[jjCor]] <- sapply(CbDMeasures(
            cyclicRVSystem(n = iiRank + 1,
                           refTable = refTable,
                           negTable = negTable)),
            function (y) y[1])
    }

    names(measuresVaryingCor[[iiRank]]) <- paste0('Value.', seq(nCors))

    measuresVaryingCor[[iiRank]] <- as.data.frame(measuresVaryingCor[[iiRank]])

    measuresVaryingCor[[iiRank]]['CNT1/CNT2', ] <-
        measuresVaryingCor[[iiRank]]['CNT1', ] /
        measuresVaryingCor[[iiRank]]['CNT2', ]

    measuresVaryingCor[[iiRank]]['2*CNT2/CNT3', ] <-
        2 * measuresVaryingCor[[iiRank]]['CNT2', ] /
        measuresVaryingCor[[iiRank]]['CNT3', ]

    measuresVaryingCor[[iiRank]]['Delta/CNT2', ] <-
        measuresVaryingCor[[iiRank]]['Delta', ] /
        measuresVaryingCor[[iiRank]]['CNT2', ]

    measuresVaryingCor[[iiRank]]['-Delta/NCNT2', ] <-
        -measuresVaryingCor[[iiRank]]['Delta', ] /
        measuresVaryingCor[[iiRank]]['NCNT2', ]

    measuresVaryingCor[[iiRank]]['mADelta/NCNT2', ] <-
        measuresVaryingCor[[iiRank]]['mA', ] /
        measuresVaryingCor[[iiRank]]['NCNT2', ]

    measuresVaryingCor[[iiRank]]['minDelta/NCNT2', ] <-
        apply(measuresVaryingCor[[iiRank]][c('mADelta/NCNT2', '-Delta/NCNT2'), ],
              2, min)


    measuresVaryingCor[[iiRank]][, 'Measure'] <- rownames(measuresVaryingCor[[iiRank]])
    measuresVaryingCor[[iiRank]][, 'Measure'] <- ordered(measuresVaryingCor[[iiRank]][, 'Measure'],
                                                         levels = measureLabels)


    measuresVaryingCor[[iiRank]] <- reshape(measuresVaryingCor[[iiRank]], direction = 'long',
                                            varying = grep('Value', names(measuresVaryingCor[[iiRank]])),
                                            timevar = 'Expectation')

    measuresVaryingCor[[iiRank]][, 'Expectation'] <- names(listReFerTables)[measuresVaryingCor[[iiRank]][, 'Expectation']]

    measuresVaryingCor[[iiRank]][abs(measuresVaryingCor[[iiRank]][, 'Value']) == 'Inf', 'Value'] <- NaN
    measuresVaryingCor[[iiRank]][(!is.nan(measuresVaryingCor[[iiRank]][, 'Value'])) &
                                     (abs(measuresVaryingCor[[iiRank]][, 'Value']) < (10 * .Machine$double.eps)), 'Value'] <- 0

    measuresVaryingCor[[iiRank]][, 'Rank'] <- iiRank + 1

}

measuresVaryingCorLong <- do.call('rbind', measuresVaryingCor)

figureData <- filter(measuresVaryingCorLong,
                     grepl('CNT./|^Delta|min', measuresVaryingCorLong[, 'Measure']))

varCorFig <- ggplot(figureData, aes(x = as.numeric(Expectation), y = Value)) +
    geom_point(aes(colour = Measure), size = 0.7) +
    facet_grid(Measure ~ Rank, scales = 'free_y') +
    xlab('Expectation of the product in varying context') +
    scale_colour_brewer(palette="Dark2")

ggsave(filename = '../output/inconsistent53v47cyclic3Varying1Ratios.png', plot = varCorFig,
       width = 11, height = 6.5)

figureData <- filter(measuresVaryingCorLong,
                     grepl('Delta$|^.?CNT.$', measuresVaryingCorLong[, 'Measure']))

measCorFig <- ggplot(figureData, aes(x = as.numeric(Expectation), y = Value)) +
    geom_point(aes(colour = Measure), size = 0.7) +
    facet_grid(Measure ~ Rank, scales = 'free_y') +
    xlab('Expectation of the product in varying context') +
    scale_colour_brewer(palette="Dark2")

ggsave(filename = '../output/inconsistent53v47cyclic3Varying1Measures.png', plot = measCorFig,
       width = 11, height = 6.5)



