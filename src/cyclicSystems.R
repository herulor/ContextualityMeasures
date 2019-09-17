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

source('./src/ContextualityMeasures.R')


################################################################################
# #  Functions
################################################################################

generateSequenceCyclicSystemsCNT2 <- function (nCors, maxRank = 2,
                                               margin1 = NULL, margin2 = NULL,
                                               varying = 'all') {


    if (is.null(margin1)) {
        margin1 <- 0.5
    }
    if (is.null(margin2)) {
        margin2 <- 0.5
    }

    listReFerTables <- list()

    if (length(nCors) == 1) {
        products <- seq(to   = min(margin1, margin2),
                        from = max(0, margin1 + margin2 - 1),
                        length.out = nCors)
    } else {
        minCoup <- max(0, margin1 + margin2 - 1)
        maxCoup <- min(margin1, margin2)
        products <- (maxCoup - minCoup) * nCors + minCoup
    }

for (iiProduct in seq_along(products)) {
    listReFerTables[[iiProduct]] <- generateRefFerTables(product = products[iiProduct],
                                                         margin1 = margin1,
                                                         margin2 = margin2)
    traveledDistance <- (products[iiProduct] - min(products)) / (max(products) - min(products))
    names(listReFerTables)[iiProduct] <- sprintf('%.3f', traveledDistance)
}

measuresVaryingCor <- list()

for (iiRank in seq(maxRank)) {

    measuresVaryingCor[[iiRank]] <- list()
    for (jjCor in seq(nCors)) {

        if (varying == 'one') {
            jjRef <- 1
            jjFer <- 1
        } else {
            jjRef <- jjCor
            jjFer <- jjCor
        }

        if ((iiRank %% 2) == 0) {
            # Odd rank systems
            refTable <- listReFerTables[[jjFer]][['ferTable']]
            negTable <- listReFerTables[[jjCor]][['ferTable']]
        } else {
            # Even rank systems
            refTable <- listReFerTables[[jjRef]][['ferTable']]
            negTable <- listReFerTables[[jjCor]][['refTable']]
        }

        measuresVaryingCor[[iiRank]][[jjCor]] <- CNT2Cyclic(
            cyclicRVSystem(n = iiRank + 1,
                           refTable = refTable,
                           negTable = negTable)
        )
    }

    measuresVaryingCor[[iiRank]] <- do.call('rbind', measuresVaryingCor[[iiRank]])
    measuresVaryingCor[[iiRank]][, 'Distance'] <- names(listReFerTables)
    measuresVaryingCor[[iiRank]][, 'Rank'] <- iiRank + 1
}

measuresVaryingCorLong <- do.call('rbind', measuresVaryingCor)

return(measuresVaryingCorLong)
}


#' Compute the noncontextuality measure NCNT2
#'
#' CNT2Cyclic(x) returns a list containing:
#'     NCNT2: The value of the NCNT2 noncontextuality measure
#'     solution: The solution found for optimal coupling
#'     status: Whether the solution was successfully found in the linear programming task.
#'         It will return 0 for the optimal solution being found, and non-zero otherwise.
#'
#' @param x a cntMatrices object
#' @param coupling optional linear programming solution of a coupling for the system (Rglpk_solve_LP or FindCoupling)
#' @return A Cnt.Meas object with the NCNT2 measure
CNT2Cyclic <- function(x) {

    if (suppressMessages(isDRVSystem(x))) {
        x <- cntMatrices(x)
    }

    if (suppressMessages(isCntMatrices(x))) {
        deltaX <- suppressWarnings(try(DeltaCFunction(x), silent = TRUE))
    }

    deltaX <- lapply(deltaX[c('deltaC', 'mA')], round, 8)

    CNT2 <- round(ifelse(deltaX[['deltaC']] > 0, deltaX[['deltaC']],
                   max(deltaX[['deltaC']],
                       -deltaX[['mA']])), 8)

    if (CNT2 <= 0) {
        if (CNT2 == deltaX[['deltaC']]) {
            if (deltaX[['mA']] == deltaX[['deltaC']]) {
                systemType <- 'Noncontextual - Close to face and boundary'
            } else {
                systemType <- 'Noncontextual - Close to face'
            }
        } else {
            systemType <- 'Noncontextual - Close to boundary'
        }
    } else {
        systemType <- 'Contextual'
    }

    output <- data.frame('N_CNT2' = CNT2 / 4,
                         systemType = systemType)

    return(output)
}



plotTravel <- function (x, parity = 'odd', ncnt = 'together') {

    measuresVaryingCorLong <- x
    if (parity == 'odd') {
        parity <- 1
    } else {
        parity <- 0
    }

    figureData <- filter(measuresVaryingCorLong,
                         Rank %% 2 == parity)


    figureData <- mutate(figureData,
                         'Type_pos'      = c(as.character(figureData[-1, 'systemType']), ''),
                         'Type_prev'     = c('', as.character(figureData[-nrow(figureData), 'systemType'])),
                         Changed         = ifelse(systemType == Type_prev,
                                                  'same',
                                                  'different'),
                         Changes         = ifelse(systemType == Type_pos,
                                                  'same',
                                                  'different'),
                         N_CNT2          = round(N_CNT2, 4))

    group_by(figureData, Rank) %>% mutate(minimal = N_CNT2 == min(N_CNT2)) ->
        figureData

    figureData <- filter(as.data.frame(unclass(figureData)),
                         Changed == 'different' | Changes == 'different' |
                             minimal |
                             Distance == '0.000' | Distance == '1.000')

    figureData <- mutate(ungroup(figureData),
                         'N_CNT2.end'    = c(figureData[-1, 'N_CNT2'], NA),
                         'Distance.end'  = c(as.character(figureData[-1, 'Distance']), NA))

    figureData <- figureData[-nrow(figureData), ]
    figureData <- figureData[!(figureData[, 'Distance.end'] == '0.000'), ]

    figureData <- mutate(figureData,
                         Type = ifelse((N_CNT2 * N_CNT2.end) < 0,
                                       'No plot',
                                       ifelse(!((N_CNT2 <= 0) | (N_CNT2.end <= 0)),
                                       'Contextual',
                                       ifelse(N_CNT2.end < N_CNT2,
                                              as.character(Type_pos),
                                              as.character(systemType)))))

    figureData <- mutate(figureData,
                         Measure = ifelse((N_CNT2 <= 0) & (N_CNT2.end <= 0),
                                          'NCNT2', 'CNT2')) %>%
        select(Rank, Type, Measure, N_CNT2, N_CNT2.end, Distance, Distance.end)

    figureTypeCon   <- filter(figureData, Type == 'Contextual')
    figureTypeFace  <- filter(figureData, grepl('face', Type))
    figureTypeBound <- filter(figureData, grepl('boundary', Type) & !grepl('face', Type))


    if (ncnt == 'together') {
        cnt2Fig <- ggplot(figureData, aes(x = as.numeric(as.character(Distance)), y = N_CNT2)) +
            ylim(-0.25, 0.5) + xlim(0, 1)
        if (nrow(figureTypeCon) > 0) {
            cnt2Fig <- cnt2Fig + geom_segment(data = figureTypeCon,
                                              aes(x = as.numeric(as.character(Distance)),
                                                  xend = as.numeric(as.character(Distance.end)),
                                                  y = N_CNT2, yend = N_CNT2.end),
                                              linetype = 1)
        }
        if (nrow(figureTypeFace) > 0) {
            cnt2Fig <- cnt2Fig + geom_segment(data = figureTypeFace,
                                              aes(x = as.numeric(as.character(Distance)),
                                                  xend = as.numeric(as.character(Distance.end)),
                                                  y = N_CNT2, yend = N_CNT2.end),
                                              linetype = 3)
        }
        if (nrow(figureTypeBound) > 0) {
            cnt2Fig <- cnt2Fig + geom_segment(data = figureTypeBound,
                                              aes(x = as.numeric(as.character(Distance)),
                                                  xend = as.numeric(as.character(Distance.end)),
                                                  y = N_CNT2, yend = N_CNT2.end),
                                              linetype = 6)
        }

        cnt2Fig <- cnt2Fig + facet_wrap(. ~ Rank,
                                        labeller = labeller(Rank = label_both))
    } else {
        cnt2Fig <- ggplot(figureData, aes(x = as.numeric(as.character(Distance)), y = abs(N_CNT2))) +
            ylim(0, 0.5) + xlim(0, 1)
        if (nrow(figureTypeCon) > 0) {
            cnt2Fig <- cnt2Fig + geom_segment(data = figureTypeCon,
                                              aes(x = as.numeric(as.character(Distance)),
                                                  xend = as.numeric(as.character(Distance.end)),
                                                  y = abs(N_CNT2), yend = abs(N_CNT2.end)),
                                              linetype = 1)
        }
        if (nrow(figureTypeFace) > 0) {
            cnt2Fig <- cnt2Fig + geom_segment(data = figureTypeFace,
                                              aes(x = as.numeric(as.character(Distance)),
                                                  xend = as.numeric(as.character(Distance.end)),
                                                  y = abs(N_CNT2), yend = abs(N_CNT2.end)),
                                              linetype = 3)
        }
        if (nrow(figureTypeBound) > 0) {
            cnt2Fig <- cnt2Fig + geom_segment(data = figureTypeBound,
                                              aes(x = as.numeric(as.character(Distance)),
                                                  xend = as.numeric(as.character(Distance.end)),
                                                  y = abs(N_CNT2), yend = abs(N_CNT2.end)),
                                              linetype = 6)
        }

        cnt2Fig <- cnt2Fig + facet_grid(Measure ~ Rank, scales = 'free_y',
                                        labeller = labeller(Rank = label_both))
    }

    cnt2Fig <- cnt2Fig +
        xlab('Proportion of distance traversed') +
        ylab('Value') +
        scale_colour_brewer(palette="Dark2")

}


################################################################################
# #  Constant |cor| varying cor from 1 in refTable (-1 in ferTable) to -1 (1)
################################################################################

epsi <- c(-1.2e-3, 1.2e-3)
nCors <- sort(c(0, 1/2 + epsi, 1/3 + epsi, 1/4 + epsi, 1/5 + epsi, 1/6 + epsi, 1/7 + epsi, 2/3, 3/4 + epsi, 4/5, 5/6 + epsi, 1))
maxRank <- 6 #Actually maxRank - 1. So this goes 2 - 7

measuresVaryingCorLongCons <- generateSequenceCyclicSystemsCNT2(nCors = nCors, maxRank = maxRank,
                                                                margin1 = 0.5, margin2 = 0.5, varying = 'all')

togetherEvenFig <- plotTravel(x = measuresVaryingCorLongCons, parity = 'even', ncnt = 'together')
separateEvenFig <- plotTravel(x = measuresVaryingCorLongCons, parity = 'even', ncnt = 'separate')

togetherOddFig <- plotTravel(x = measuresVaryingCorLongCons, parity = 'odd', ncnt = 'together')
separateOddFig <- plotTravel(x = measuresVaryingCorLongCons, parity = 'odd', ncnt = 'separate')


ggsave(filename = './output/cyclicSystems_Consistent_Constant_correlation_Lines_together_Even_rank.png',
       plot = togetherEvenFig,
       width = 6.5, height = 4.0)

ggsave(filename = './output/cyclicSystems_Consistent_Constant_correlation_Lines_separate_Even_rank.png',
       plot = separateEvenFig,
       width = 7.5, height = 4.5)

ggsave(filename = './output/cyclicSystems_Consistent_Constant_correlation_Lines_together_Odd_rank.png',
       plot = togetherOddFig,
       width = 6.5, height = 4.0)

ggsave(filename = './output/cyclicSystems_Consistent_Constant_correlation_Lines_separate_Odd_rank.png',
       plot = separateOddFig,
       width = 7.5, height = 4.5)



################################################################################
# #  Inconsistently connected systems - Margins .55 and .4  - Constant correlations
################################################################################

nCors <- sort(c(0, 1/32 + epsi, 1/24, 5/16 + epsi, 11/16 + epsi, 1/2, 1/4, 1/6 + epsi, 1))

measuresVaryingCorLongInc <- generateSequenceCyclicSystemsCNT2(nCors = nCors, maxRank = maxRank,
                                                               margin1 = 0.55, margin2 = 0.4, varying = 'all')

togetherEvenFig <- plotTravel(x = measuresVaryingCorLongInc, parity = 'even', ncnt = 'together')
separateEvenFig <- plotTravel(x = measuresVaryingCorLongInc, parity = 'even', ncnt = 'separate')

togetherOddFig <- plotTravel(x = measuresVaryingCorLongInc, parity = 'odd', ncnt = 'together')
separateOddFig <- plotTravel(x = measuresVaryingCorLongInc, parity = 'odd', ncnt = 'separate')


ggsave(filename = './output/cyclicSystems_Inconsistent_0.55-0.4_Constant_correlation_Lines_together_Even_rank.png',
       plot = togetherEvenFig,
       width = 6.5, height = 4.0)

ggsave(filename = './output/cyclicSystems_Inconsistent_0.55-0.4_Constant_correlation_Lines_separate_Even_rank.png',
       plot = separateEvenFig,
       width = 7.5, height = 4.5)

ggsave(filename = './output/cyclicSystems_Inconsistent_0.55-0.4_Constant_correlation_Lines_together_Odd_rank.png',
       plot = togetherOddFig,
       width = 6.5, height = 4.0)

ggsave(filename = './output/cyclicSystems_Inconsistent_0.55-0.4_Constant_correlation_Lines_separate_Odd_rank.png',
       plot = separateOddFig,
       width = 7.5, height = 4.5)

