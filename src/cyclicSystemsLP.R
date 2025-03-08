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


cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

measureLabels <- c('Delta', 'mA', 'CNT1', 'CNT2', 'CNT2h', 'levelC', 'CNT3', 'CNFR', 'NCNT2', 'NCNT2h', 'levelN',
                   'Delta/CNT2', '-Delta/NCNT2', 'mADelta/NCNT2', 'minDelta/NCNT2', 'CNT1/CNT2', '2*CNT2/CNT3')




################################################################################
# #  Functions
################################################################################


generateSequenceCyclicSystemsMeasures <- function (nCors, maxRank = 2,
                                                   margin1 = NULL, margin2 = NULL,
                                                   varying = 'all') {


    if (is.null(margin1)) {
        margin1 <- 0.5
    }
    if (is.null(margin2)) {
        margin2 <- 0.5
    }

    listReFerTables <- list()

    products <- seq(to   = min(margin1, margin2),
                    from = max(0, margin1 + margin2 - 1),
                    length.out = nCors)

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

            measuresVaryingCor[[iiRank]][[jjCor]] <- #sapply(CbDMeasures(
                cyclicRVSystem(n = iiRank + 1,
                               refTable = refTable,
                               negTable = negTable)#),
#                function (y) y[1])
        }

       # names(measuresVaryingCor[[iiRank]]) <- paste0('Value.', seq(nCors))

       # measuresVaryingCor[[iiRank]] <- as.data.frame(measuresVaryingCor[[iiRank]])
#
#        measuresVaryingCor[[iiRank]]['CNT1/CNT2', ] <-
#            measuresVaryingCor[[iiRank]]['CNT1', ] /
#            measuresVaryingCor[[iiRank]]['CNT2', ]
#
#        measuresVaryingCor[[iiRank]]['2*CNT2/CNT3', ] <-
#            2 * measuresVaryingCor[[iiRank]]['CNT2', ] /
#            measuresVaryingCor[[iiRank]]['CNT3', ]
#
#        measuresVaryingCor[[iiRank]]['Delta/CNT2', ] <-
#            measuresVaryingCor[[iiRank]]['Delta', ] /
#            measuresVaryingCor[[iiRank]]['CNT2', ]
#
#        measuresVaryingCor[[iiRank]]['-Delta/NCNT2', ] <-
#            -measuresVaryingCor[[iiRank]]['Delta', ] /
#            measuresVaryingCor[[iiRank]]['NCNT2', ]
#
#        measuresVaryingCor[[iiRank]]['mADelta/NCNT2', ] <-
#            measuresVaryingCor[[iiRank]]['mA', ] /
#            measuresVaryingCor[[iiRank]]['NCNT2', ]
#
#        measuresVaryingCor[[iiRank]]['minDelta/NCNT2', ] <-
#            apply(measuresVaryingCor[[iiRank]][c('mADelta/NCNT2', '-Delta/NCNT2'), ],
#                  2, min)
#
#
#        measuresVaryingCor[[iiRank]][, 'Measure'] <- rownames(measuresVaryingCor[[iiRank]])
#        measuresVaryingCor[[iiRank]][, 'Measure'] <- ordered(measuresVaryingCor[[iiRank]][, 'Measure'],
#                                                             levels = measureLabels)
##
#
#        measuresVaryingCor[[iiRank]] <- reshape(measuresVaryingCor[[iiRank]], direction = 'long',
#                                                varying = grep('Value', names(measuresVaryingCor[[iiRank]])),
#                                                timevar = 'Expectation')
#
#        measuresVaryingCor[[iiRank]][, 'Expectation'] <- names(listReFerTables)[measuresVaryingCor[[iiRank]][, 'Expectation']]
#
#        measuresVaryingCor[[iiRank]][abs(measuresVaryingCor[[iiRank]][, 'Value']) == 'Inf', 'Value'] <- NaN
#        measuresVaryingCor[[iiRank]][(!is.nan(measuresVaryingCor[[iiRank]][, 'Value'])) &
#                                         (abs(measuresVaryingCor[[iiRank]][, 'Value']) < (10 * .Machine$double.eps)), 'Value'] <- 0
#
#        measuresVaryingCor[[iiRank]][, 'Rank'] <- iiRank + 1
    }

    measuresVaryingCorLong <- do.call('rbind', measuresVaryingCor)

    return(measuresVaryingCorLong)
}


################################################################################
# #  Maximally contextual systems CbD measures
################################################################################

maximallyContextualSystems <- list()

maxRank <- 6 # actually maxRank - 1

## Computations

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


## Figure

figureData <- filter(maximallyContextualSystemsLong,
                     grepl('/C|Delta$', maximallyContextualSystemsLong[, 'Measure']))

maximalFig <- ggplot(figureData, aes(x = Rank, y = Value)) +
    geom_line(aes(colour = Measure, linetype = Measure)) +
    scale_x_continuous(breaks = seq(2, maxRank + 1)) +
    scale_y_continuous(breaks = seq(1, maxRank), limits = c(0, maxRank)) +
    scale_colour_brewer(palette="Dark2")


ggsave(filename = './output/maximallyContextual.png',
       plot = maximalFig,
       width = 6, height = 4)





################################################################################
# #  Constant |cor| varying cor from 1 in refTable (-1 in ferTable) to -1 (1)
################################################################################

nCors <- 41

measuresVaryingCorLong <- generateSequenceCyclicSystemsMeasures(nCors = nCors, maxRank = maxRank,
                                                                margin1 = 0.5, margin2 = 0.5, varying = 'all')

## Save results

write.table(measuresVaryingCorLong, file = './output/cyclicSystems_Consistent_Constant_correlation.txt',
            sep = '\t', col.names = TRUE, row.names = FALSE)


## Figures

figureData <- filter(measuresVaryingCorLong,
                     grepl('CNT./|^Delta|min', measuresVaryingCorLong[, 'Measure'])
                     & Rank %% 2 == 0)

varCorFig <- ggplot(figureData, aes(x = as.numeric(Expectation), y = Value)) +
    geom_point(aes(colour = Measure), size = 0.7) +
    facet_grid(Measure ~ Rank, scales = 'free_y') +
    xlab('Proportion of distance traversed') +
    scale_colour_brewer(palette="Dark2")

ggsave(filename = './output/cyclicSystems_Consistent_Constant_correlation_Ratios_Even_rank.png', plot = varCorFig,
       width = 11, height = 6.5)



figureData <- filter(measuresVaryingCorLong,
                     grepl('CNT./|^Delta|min', measuresVaryingCorLong[, 'Measure'])
                     & Rank %% 2 == 1)

varCorFig <- ggplot(figureData, aes(x = as.numeric(Expectation), y = Value)) +
    geom_point(aes(colour = Measure), size = 0.7) +
    facet_grid(Measure ~ Rank, scales = 'free_y') +
    xlab('Proportion of distance traversed') +
    scale_colour_brewer(palette="Dark2")

ggsave(filename = './output/cyclicSystems_Consistent_Constant_correlation_Ratios_Odd_rank.png', plot = varCorFig,
       width = 11, height = 6.5)




measuresVaryingCorLong <- filter(measuresVaryingCorLong, Measure == 'CNT2' | Measure == 'NCNT2')


figureData <- filter(measuresVaryingCorLong,
                     grepl('Delta$|^.?CNT.$', measuresVaryingCorLong[, 'Measure'])
                     & Rank %% 2 == 0)

measCorFig <- ggplot(figureData, aes(x = as.numeric(Expectation), y = Value)) +
    geom_point(aes(colour = Measure), size = 0.7) +
    facet_grid(Measure ~ Rank, scales = 'free_y') +
    xlab('Proportion of distance traversed') +
    scale_colour_brewer(palette="Dark2")

ggsave(filename = './output/cyclicSystems_Consistent_Constant_correlation_Measures_Even_rank.png', plot = measCorFig,
       width = 7, height = 4.5)


figureData <- filter(measuresVaryingCorLong,
                     grepl('Delta$|^.?CNT.$', measuresVaryingCorLong[, 'Measure'])
                     & Rank %% 2 == 1)

measCorFig <- ggplot(figureData, aes(x = as.numeric(Expectation), y = Value)) +
    geom_point(aes(colour = Measure), size = 0.7) +
    facet_grid(Measure ~ Rank, scales = 'free_y') +
    xlab('Proportion of distance traversed') +
    scale_colour_brewer(palette="Dark2")

ggsave(filename = './output/cyclicSystems_Consistent_Constant_correlation_Measures_Odd_rank.png', plot = measCorFig,
       width = 7, height = 4.5)





################################################################################
# #  Consistently connected. Varying one context from perfect contextual system
################################################################################

measuresVaryingCorLong <- generateSequenceCyclicSystemsMeasures(nCors = nCors, maxRank = maxRank,
                                                                margin1 = 0.5, margin2 = 0.5, varying = 'one')


## Save results

write.table(measuresVaryingCorLong, file = './output/cyclicSystems_Consistent_Varying_one_context.txt',
            sep = '\t', col.names = TRUE, row.names = FALSE)


## Figures

figureData <- filter(measuresVaryingCorLong,
                     grepl('CNT./|^Delta|min', measuresVaryingCorLong[, 'Measure'])
                     & Rank %% 2 == 0)

varCorFig <- ggplot(figureData, aes(x = as.numeric(Expectation), y = Value)) +
    geom_point(aes(colour = Measure), size = 0.7) +
    facet_grid(Measure ~ Rank, scales = 'free_y') +
    xlab('Proportion of distance traversed') +
    scale_colour_brewer(palette="Dark2")

ggsave(filename = './output/cyclicSystems_Consistent_Varying_one_context_Ratios_Even_rank.png', plot = varCorFig,
       width = 11, height = 6.5)


figureData <- filter(measuresVaryingCorLong,
                     grepl('CNT./|^Delta|min', measuresVaryingCorLong[, 'Measure'])
                     & Rank %% 2 == 1)

varCorFig <- ggplot(figureData, aes(x = as.numeric(Expectation), y = Value)) +
    geom_point(aes(colour = Measure), size = 0.7) +
    facet_grid(Measure ~ Rank, scales = 'free_y') +
    xlab('Proportion of distance traversed') +
    scale_colour_brewer(palette="Dark2")

ggsave(filename = './output/cyclicSystems_Consistent_Varying_one_context_Ratios_Odd_rank.png', plot = varCorFig,
       width = 11, height = 6.5)




measuresVaryingCorLong <- filter(measuresVaryingCorLong, Measure == 'CNT2' | Measure == 'NCNT2')


figureData <- filter(measuresVaryingCorLong,
                     grepl('Delta$|^.?CNT.$', measuresVaryingCorLong[, 'Measure'])
                     & Rank %% 2 == 0)

measCorFig <- ggplot(figureData, aes(x = as.numeric(Expectation), y = Value)) +
    geom_point(aes(colour = Measure), size = 0.7) +
    facet_grid(Measure ~ Rank, scales = 'free_y') +
    xlab('Proportion of distance traversed') +
    scale_colour_brewer(palette="Dark2")

ggsave(filename = './output/cyclicSystems_Consistent_Varying_one_context_Measures_Even_rank.png', plot = measCorFig,
       width = 7, height = 4.5)




figureData <- filter(measuresVaryingCorLong,
                     grepl('Delta$|^.?CNT.$', measuresVaryingCorLong[, 'Measure'])
                     & Rank %% 2 == 1)

measCorFig <- ggplot(figureData, aes(x = as.numeric(Expectation), y = Value)) +
    geom_point(aes(colour = Measure), size = 0.7) +
    facet_grid(Measure ~ Rank, scales = 'free_y') +
    xlab('Proportion of distance traversed') +
    scale_colour_brewer(palette="Dark2")

ggsave(filename = './output/cyclicSystems_Consistent_Varying_one_context_Measures_Odd_rank.png', plot = measCorFig,
       width = 7, height = 4.5)



################################################################################
# #  Inconsistently connected systems - Margins .3 and .7 a - Constant correlations
################################################################################

measuresVaryingCorLong <- generateSequenceCyclicSystemsMeasures(nCors = nCors, maxRank = maxRank,
                                                                margin1 = 0.3, margin2 = 0.6, varying = 'all')


## Save results

write.table(measuresVaryingCorLong, file = './output/cyclicSystems_Inconsistent_0.3-0.6_Constant_correlation.txt',
            sep = '\t', col.names = TRUE, row.names = FALSE)


## Figures

figureData <- filter(measuresVaryingCorLong,
                     grepl('CNT./|^Delta|min', measuresVaryingCorLong[, 'Measure'])
                     & Rank %% 2 == 0)

varCorFig <- ggplot(figureData, aes(x = as.numeric(Expectation), y = Value)) +
    geom_point(aes(colour = Measure), size = 0.7) +
    facet_grid(Measure ~ Rank, scales = 'free_y') +
    xlab('Proportion of distance traversed') +
    scale_colour_brewer(palette="Dark2")

ggsave(filename = './output/cyclicSystems_Inconsistent_0.3-0.6_Constant_correlation_Ratios_Even_rank.png', plot = varCorFig,
       width = 11, height = 6.5)


figureData <- filter(measuresVaryingCorLong,
                     grepl('CNT./|^Delta|min', measuresVaryingCorLong[, 'Measure'])
                     & Rank %% 2 == 1)

varCorFig <- ggplot(figureData, aes(x = as.numeric(Expectation), y = Value)) +
    geom_point(aes(colour = Measure), size = 0.7) +
    facet_grid(Measure ~ Rank, scales = 'free_y') +
    xlab('Proportion of distance traversed') +
    scale_colour_brewer(palette="Dark2")

ggsave(filename = './output/cyclicSystems_Inconsistent_0.3-0.6_Constant_correlation_Ratios_Odd_rank.png', plot = varCorFig,
       width = 11, height = 6.5)



measuresVaryingCorLong <- filter(measuresVaryingCorLong, Measure == 'CNT2' | Measure == 'NCNT2')



figureData <- filter(measuresVaryingCorLong,
                     grepl('Delta$|^.?CNT.$', measuresVaryingCorLong[, 'Measure'])
                     & Rank %% 2 == 0)

measCorFig <- ggplot(figureData, aes(x = as.numeric(Expectation), y = Value)) +
    geom_point(aes(colour = Measure), size = 0.7) +
    facet_grid(Measure ~ Rank, scales = 'free_y') +
    xlab('Proportion of distance traversed') +
    scale_colour_brewer(palette="Dark2")


ggsave(filename = './output/cyclicSystems_Inconsistent_0.3-0.6_Constant_correlation_Measures_Even_rank.png', plot = measCorFig,
       width = 7, height = 4.5)


figureData <- filter(measuresVaryingCorLong,
                     grepl('Delta$|^.?CNT.$', measuresVaryingCorLong[, 'Measure'])
                     & Rank %% 2 == 1)

measCorFig <- ggplot(figureData, aes(x = as.numeric(Expectation), y = Value)) +
    geom_point(aes(colour = Measure), size = 0.7) +
    facet_grid(Measure ~ Rank, scales = 'free_y') +
    xlab('Proportion of distance traversed') +
    scale_colour_brewer(palette="Dark2")


ggsave(filename = './output/cyclicSystems_Inconsistent_0.3-0.6_Constant_correlation_Measures_Odd_rank.png', plot = measCorFig,
       width = 7, height = 4.5)




################################################################################
# #  Inconsistently connected systems - Margins .55 and .4  - Constant correlations
################################################################################

measuresVaryingCorLong <- generateSequenceCyclicSystemsMeasures(nCors = nCors, maxRank = maxRank,
                                                                margin1 = 0.55, margin2 = 0.4, varying = 'all')



## Save results

write.table(measuresVaryingCorLong, file = './output/cyclicSystems_Inconsistent_0.55-0.4_Constant_correlation.txt',
            sep = '\t', col.names = TRUE, row.names = FALSE)


## Figures

figureData <- filter(measuresVaryingCorLong,
                     grepl('CNT./|^Delta|min', measuresVaryingCorLong[, 'Measure'])
                     & Rank %% 2 == 0)

varCorFig <- ggplot(figureData, aes(x = as.numeric(Expectation), y = Value)) +
    geom_point(aes(colour = Measure), size = 0.7) +
    facet_grid(Measure ~ Rank, scales = 'free_y') +
    xlab('Proportion of distance traversed') +
    scale_colour_brewer(palette="Dark2")

ggsave(filename = './output/cyclicSystems_Inconsistent_0.55-0.4_Constant_correlation_Ratios_Even_rank.png', plot = varCorFig,
       width = 11, height = 6.5)


figureData <- filter(measuresVaryingCorLong,
                     grepl('CNT./|^Delta|min', measuresVaryingCorLong[, 'Measure'])
                     & Rank %% 2 == 1)

varCorFig <- ggplot(figureData, aes(x = as.numeric(Expectation), y = Value)) +
    geom_point(aes(colour = Measure), size = 0.7) +
    facet_grid(Measure ~ Rank, scales = 'free_y') +
    xlab('Proportion of distance traversed') +
    scale_colour_brewer(palette="Dark2")

ggsave(filename = './output/cyclicSystems_Inconsistent_0.55-0.4_Constant_correlation_Ratios_Odd_rank.png', plot = varCorFig,
       width = 11, height = 6.5)




measuresVaryingCorLong <- filter(measuresVaryingCorLong, Measure == 'CNT2' | Measure == 'NCNT2')


figureData <- filter(measuresVaryingCorLong,
                     grepl('Delta$|^.?CNT.$', measuresVaryingCorLong[, 'Measure'])
                     & Rank %% 2 == 0)

measCorFig <- ggplot(figureData, aes(x = as.numeric(Expectation), y = Value)) +
    geom_point(aes(colour = Measure), size = 0.7) +
    facet_grid(Measure ~ Rank, scales = 'free_y') +
    xlab('Proportion of distance traversed') +
    scale_colour_brewer(palette="Dark2")


ggsave(filename = './output/cyclicSystems_Inconsistent_0.55-0.4_Constant_correlation_Measures_Even_rank.png', plot = measCorFig,
       width = 7, height = 4.5)


figureData <- filter(measuresVaryingCorLong,
                     grepl('Delta$|^.?CNT.$', measuresVaryingCorLong[, 'Measure'])
                     & Rank %% 2 == 1)

measCorFig <- ggplot(figureData, aes(x = as.numeric(Expectation), y = Value)) +
    geom_point(aes(colour = Measure), size = 0.7) +
    facet_grid(Measure ~ Rank, scales = 'free_y') +
    xlab('Proportion of distance traversed') +
    scale_colour_brewer(palette="Dark2")


ggsave(filename = './output/cyclicSystems_Inconsistent_0.55-0.4_Constant_correlation_Measures_Odd_rank.png', plot = measCorFig,
       width = 7, height = 4.5)





################################################################################
# #  Maximal but one inconsistently connected systems - Margins .3 and .6
################################################################################

measuresVaryingCorLong <- generateSequenceCyclicSystemsMeasures(nCors = nCors, maxRank = maxRank,
                                                                margin1 = 0.3, margin2 = 0.6, varying = 'one')



## Save results

write.table(measuresVaryingCorLong, file = './output/cyclicSystems_Inconsistent_0.3-0.6_Varying_one.txt',
            sep = '\t', col.names = TRUE, row.names = FALSE)


## Figures

figureData <- filter(measuresVaryingCorLong,
                     grepl('CNT./|^Delta|min', measuresVaryingCorLong[, 'Measure'])
                     & Rank %% 2 == 0)

varCorFig <- ggplot(figureData, aes(x = as.numeric(Expectation), y = Value)) +
    geom_point(aes(colour = Measure), size = 0.7) +
    facet_grid(Measure ~ Rank, scales = 'free_y') +
    xlab('Proportion of distance traveled') +
    scale_colour_brewer(palette="Dark2")

ggsave(filename = './output/cyclicSystems_Inconsistent_0.3-0.6_Varying_one_context_Ratios_Even_rank.png', plot = varCorFig,
       width = 11, height = 6.5)


figureData <- filter(measuresVaryingCorLong,
                     grepl('CNT./|^Delta|min', measuresVaryingCorLong[, 'Measure'])
                     & Rank %% 2 == 1)

varCorFig <- ggplot(figureData, aes(x = as.numeric(Expectation), y = Value)) +
    geom_point(aes(colour = Measure), size = 0.7) +
    facet_grid(Measure ~ Rank, scales = 'free_y') +
    xlab('Proportion of distance traveled') +
    scale_colour_brewer(palette="Dark2")

ggsave(filename = './output/cyclicSystems_Inconsistent_0.3-0.6_Varying_one_context_Ratios_Odd_rank.png', plot = varCorFig,
       width = 11, height = 6.5)




measuresVaryingCorLong <- filter(measuresVaryingCorLong, Measure == 'CNT2' | Measure == 'NCNT2')


figureData <- filter(measuresVaryingCorLong,
                     grepl('Delta$|^.?CNT.$', measuresVaryingCorLong[, 'Measure'])
                     & Rank %% 2 == 0)

measCorFig <- ggplot(figureData, aes(x = as.numeric(Expectation), y = Value)) +
    geom_point(aes(colour = Measure), size = 0.7) +
    facet_grid(Measure ~ Rank, scales = 'free_y') +
    xlab('Proportion of distance traveled') +
    scale_colour_brewer(palette="Dark2")


ggsave(filename = './output/cyclicSystems_Inconsistent_0.3-0.6_Varying_one_context_Measures_Even_rank.png', plot = measCorFig,
       width = 7, height = 4.5)


figureData <- filter(measuresVaryingCorLong,
                     grepl('Delta$|^.?CNT.$', measuresVaryingCorLong[, 'Measure'])
                     & Rank %% 2 == 1)

measCorFig <- ggplot(figureData, aes(x = as.numeric(Expectation), y = Value)) +
    geom_point(aes(colour = Measure), size = 0.7) +
    facet_grid(Measure ~ Rank, scales = 'free_y') +
    xlab('Proportion of distance traveled') +
    scale_colour_brewer(palette="Dark2")


ggsave(filename = './output/cyclicSystems_Inconsistent_0.3-0.6_Varying_one_context_Measures_Odd_rank.png', plot = measCorFig,
       width = 7, height = 4.5)




################################################################################
# #  Maximal but one Inconsistently connected systems - Margins .53 and .47
################################################################################

measuresVaryingCorLong <- generateSequenceCyclicSystemsMeasures(nCors = nCors, maxRank = maxRank,
                                                                margin1 = 0.55, margin2 = 0.4, varying = 'one')


## Save results

write.table(measuresVaryingCorLong, file = './output/cyclicSystems_Inconsistent_0.55-0.4_Varying_one.txt',
            sep = '\t', col.names = TRUE, row.names = FALSE)


## Figures

figureData <- filter(measuresVaryingCorLong,
                     grepl('CNT./|^Delta|min', measuresVaryingCorLong[, 'Measure'])
                     & Rank %% 2 == 0)

varCorFig <- ggplot(figureData, aes(x = as.numeric(Expectation), y = Value)) +
    geom_point(aes(colour = Measure), size = 0.7) +
    facet_grid(Measure ~ Rank, scales = 'free_y') +
    xlab('Proportion of distance traveled') +
    scale_colour_brewer(palette="Dark2")

ggsave(filename = './output/cyclicSystems_Inconsistent_0.55-0.4_Varying_one_context_Ratios_Even_rank.png', plot = varCorFig,
       width = 11, height = 6.5)


figureData <- filter(measuresVaryingCorLong,
                     grepl('CNT./|^Delta|min', measuresVaryingCorLong[, 'Measure'])
                     & Rank %% 2 == 1)

varCorFig <- ggplot(figureData, aes(x = as.numeric(Expectation), y = Value)) +
    geom_point(aes(colour = Measure), size = 0.7) +
    facet_grid(Measure ~ Rank, scales = 'free_y') +
    xlab('Proportion of distance traveled') +
    scale_colour_brewer(palette="Dark2")

ggsave(filename = './output/cyclicSystems_Inconsistent_0.55-0.4_Varying_one_context_Ratios_Odd_rank.png', plot = varCorFig,
       width = 11, height = 6.5)




measuresVaryingCorLong <- filter(measuresVaryingCorLong, Measure == 'CNT2' | Measure == 'NCNT2')


figureData <- filter(measuresVaryingCorLong,
                     grepl('Delta$|^.?CNT.$', measuresVaryingCorLong[, 'Measure'])
                     & Rank %% 2 == 0)

measCorFig <- ggplot(figureData, aes(x = as.numeric(Expectation), y = Value)) +
    geom_point(aes(colour = Measure), size = 0.7) +
    facet_grid(Measure ~ Rank, scales = 'free_y') +
    xlab('Proportion of distance traveled') +
    scale_colour_brewer(palette="Dark2")


ggsave(filename = './output/cyclicSystems_Inconsistent_0.55-0.4_Varying_one_context_Measures_Even_rank.png', plot = measCorFig,
       width = 7, height = 4.5)


figureData <- filter(measuresVaryingCorLong,
                     grepl('Delta$|^.?CNT.$', measuresVaryingCorLong[, 'Measure'])
                     & Rank %% 2 == 1)

measCorFig <- ggplot(figureData, aes(x = as.numeric(Expectation), y = Value)) +
    geom_point(aes(colour = Measure), size = 0.7) +
    facet_grid(Measure ~ Rank, scales = 'free_y') +
    xlab('Proportion of distance traveled') +
    scale_colour_brewer(palette="Dark2")


ggsave(filename = './output/cyclicSystems_Inconsistent_0.55-0.4_Varying_one_context_Measures_Odd_rank.png', plot = measCorFig,
       width = 7, height = 4.5)


