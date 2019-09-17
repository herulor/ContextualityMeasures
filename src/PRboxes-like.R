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
library("GGally")
library("dplyr")

source("./src/ContextualityMeasures.R")


cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
                "#D55E00", "#CC79A7")


measureLabels <- c('CNT1', 'CNT2', 'CNT3', 'NCNT2',
                   'CNT1/CNT2', '2*CNT1/CNT3', '2*CNT2/CNT3')




################################################################################
# # Functions
################################################################################


generateSequencePRBoxlikeMeasures <- function (nCors, PRrank = 3L,
                                               margin1 = NULL, margin2 = NULL,
                                               varying = 'all') {


  if (is.null(margin1)) {
    margin1 <- 0.5
  }
  if (is.null(margin2)) {
    margin2 <- 0.5
  }

  listReFerTables <- list()

  products <- seq(from = min(margin1, margin2),
                  to   = max(0, margin1 + margin2 - 1),
                  length.out = nCors)


  for (iiProduct in seq(nCors)) {
    listReFerTables[[iiProduct]] <- generateRefFerTables(product = products[iiProduct])
    expectation <- 2 * sum(diag(listReFerTables[[iiProduct]][['refTable']])) - 1
    names(listReFerTables)[iiProduct] <- sprintf('%.3f', expectation)

  }

  measuresVaryingCor <- list()


  for (jjCor in seq(nCors)) {

    if (varying == 'one') {
      jjRef <- 1
      jjFer <- 1
    } else {
      jjRef <- jjCor
      jjFer <- jjCor
    }

    refTable <- listReFerTables[[jjRef]][['refTable']]
    negTable <- listReFerTables[[jjFer]][['refTable']][2:1, ]
    chgTable <- listReFerTables[[jjCor]][['refTable']][2:1, ]

    measuresVaryingCor[[jjCor]] <- sapply(CbDMeasures(
      generatePRBoxLikeSystem(PRBoxStructure = generatePRBoxNStructure(n = PRrank),
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

  measuresVaryingCor['2*CNT1/CNT3', ] <-
    2 * measuresVaryingCor['CNT1', ] /
    measuresVaryingCor['CNT3', ]

  measuresVaryingCor['2*CNT2/CNT3', ] <-
    2 * measuresVaryingCor['CNT2', ] /
    measuresVaryingCor['CNT3', ]

  measuresVaryingCor[, 'Measure'] <- rownames(measuresVaryingCor)
  measuresVaryingCor[, 'Measure'] <- ordered(measuresVaryingCor[, 'Measure'],
                                             levels = measureLabels)


  measuresVaryingCorLong <- reshape(measuresVaryingCor, direction = 'long',
                                    varying = grep('Value', names(measuresVaryingCor)),
                                    timevar = 'Expectation')

  measuresVaryingCorLong[, 'Expectation'] <- names(listReFerTables)[measuresVaryingCorLong[, 'Expectation']]

  measuresVaryingCorLong[abs(measuresVaryingCorLong[, 'Value']) == 'Inf', 'Value'] <- NaN
  measuresVaryingCorLong[(!is.nan(measuresVaryingCorLong[, 'Value'])) &
                           (abs(measuresVaryingCorLong[, 'Value']) < (10 * .Machine$double.eps)), 'Value'] <- 0


  return(measuresVaryingCorLong)
}




################################################################################
# #  Constant |cor| varying cor from 1 in refTable (-1 in ferTable) to -1 (1)
################################################################################

nCors <- 31

measuresVaryingCorLong <- generateSequencePRBoxlikeMeasures(nCors = nCors,
                                                            PRrank = 3L,
                                                            margin1 = 0.5,
                                                            margin2 = 0.5,
                                                            varying = 'all')

## Save results

#write.table(measuresVaryingCorLong, file = './output/PRlikeSystem_Consistent_Constant_correlation.txt',
#            sep = '\t', col.names = TRUE, row.names = FALSE)

measuresVaryingCorLong <- read.table('./output/PRlikeSystem_Consistent_Constant_correlation.txt',
                                      header = TRUE, sep = '\t')


# # Transform to wide table for scatterplots

dplyr::inner_join(dplyr::inner_join(
  dplyr::filter(measuresVaryingCorLong, Measure == 'CNT1') %>%
    mutate(CNT1 = Value) %>% select('Expectation', 'CNT1'),

  dplyr::filter(measuresVaryingCorLong, Measure == 'CNT2') %>%
    mutate(CNT2 = Value) %>% select('Expectation', 'CNT2')),

  dplyr::filter(measuresVaryingCorLong, Measure == 'CNT3') %>%
    mutate(CNT3 = Value) %>% select('Expectation', 'CNT3')) %>%

  mutate(Type = 'Consistent', Varying = 'All') ->

  measuresWideCA



## Figures

figureData <- filter(measuresVaryingCorLong,
                     grepl('CNT./', measuresVaryingCorLong[, 'Measure']))

varCorFig <- ggplot(figureData, aes(x = as.numeric(Expectation), y = Value)) +
  geom_point(aes(colour = Measure), size = 0.7) +
  facet_grid(Measure ~ ., scales = 'free_y') +
  xlab('Common expectation of the product') +
  scale_colour_brewer(palette="Dark2")

ggsave(filename = './output/PRlikeSystem_Consistent_Constant_correlation_Ratios.png', plot = varCorFig,
       width = 7, height = 5)

figureData <- filter(measuresVaryingCorLong,
                     grepl('^N?CNT.$', measuresVaryingCorLong[, 'Measure']))

measCorFig <- ggplot(figureData, aes(x = as.numeric(Expectation), y = Value)) +
  geom_point(aes(colour = Measure), size = 0.7) +
  facet_grid(Measure ~ ., scales = 'free_y') +
  xlab('Common expectation of the product') +
  scale_colour_brewer(palette="Dark2")

ggsave(filename = './output/PRlikeSystem_Consistent_Constant_correlation_Measures.png', plot = measCorFig,
       width = 7, height = 5)





################################################################################
# #  Maximal but 1 varying 1 cor from 1 in refTable (-1 in ferTable) to -1 (1)
################################################################################


measuresVaryingCorLong <- generateSequencePRBoxlikeMeasures(nCors = nCors,
                                                            PRrank = 3L,
                                                            margin1 = 0.5,
                                                            margin2 = 0.5,
                                                            varying = 'one')

## Save results

write.table(measuresVaryingCorLong, file = './output/PRlikeSystem_Consistent_Varying_one_context.txt',
            sep = '\t', col.names = TRUE, row.names = FALSE)

#measuresVaryingCorLong <- read.table(file = './output/PRlikeSystem_Consistent_Varying_one_context.txt',
#           header = TRUE, sep = '\t')

# # Transform to wide table for scatterplots

dplyr::inner_join(dplyr::inner_join(
  dplyr::filter(measuresVaryingCorLong, Measure == 'CNT1') %>%
    mutate(CNT1 = Value) %>% select('Expectation', 'CNT1'),

  dplyr::filter(measuresVaryingCorLong, Measure == 'CNT2') %>%
    mutate(CNT2 = Value) %>% select('Expectation', 'CNT2')),

  dplyr::filter(measuresVaryingCorLong, Measure == 'CNT3') %>%
    mutate(CNT3 = Value) %>% select('Expectation', 'CNT3')) %>%

  mutate(Type = 'Consistent', Varying = 'One') ->

  measuresWideCO



## Figures

figureData <- filter(measuresVaryingCorLong,
                     grepl('CNT./', measuresVaryingCorLong[, 'Measure']))

varCorFig <- ggplot(figureData, aes(x = as.numeric(Expectation), y = round(Value, 10))) +
  geom_point(aes(colour = Measure), size = 0.7) +
  facet_grid(Measure ~ ., scales = 'free_y') +
  xlab('Expectation of the product in varying context') +
  ylab('Value') +
  scale_colour_brewer(palette="Dark2")

ggsave(filename = './output/PRlikeSystem_Consistent_Varying_one_context_Ratios.png', plot = varCorFig,
       width = 7, height = 5)

figureData <- filter(measuresVaryingCorLong,
                     grepl('^N?CNT.$', measuresVaryingCorLong[, 'Measure']))

measCorFig <- ggplot(figureData, aes(x = as.numeric(Expectation), y = round(Value, 10))) +
  geom_point(aes(colour = Measure), size = 0.7) +
  facet_grid(Measure ~ ., scales = 'free_y') +
  xlab('Expectation of the product in varying context') +
  ylab('Value') +
  scale_colour_brewer(palette="Dark2")

ggsave(filename = './output/PRlikeSystem_Consistent_Varying_one_context_Measures.png', plot = measCorFig,
       width = 7, height = 5)


################################################################################
# # Inconsistent 0.55 - 0.40. - Constant correlation
################################################################################

measuresVaryingCorLong <- generateSequencePRBoxlikeMeasures(nCors = nCors,
                                                            PRrank = 3L,
                                                            margin1 = 0.55,
                                                            margin2 = 0.40,
                                                            varying = 'all')


## Save results

write.table(measuresVaryingCorLong, file = './output/PRlikeSystem_Inconsistent_0.55-0.40_Constant_correlation.txt',
            sep = '\t', col.names = TRUE, row.names = FALSE)


#measuresVaryingCorLong <- read.table('./output/PRlikeSystem_Inconsistent_0.55-0.40_Constant_correlation.txt',
#                                     header = TRUE, sep = '\t')

# # Transform to wide table for scatterplots

dplyr::inner_join(dplyr::inner_join(
  dplyr::filter(measuresVaryingCorLong, Measure == 'CNT1') %>%
    mutate(CNT1 = Value) %>% select('Expectation', 'CNT1'),

  dplyr::filter(measuresVaryingCorLong, Measure == 'CNT2') %>%
    mutate(CNT2 = Value) %>% select('Expectation', 'CNT2')),

  dplyr::filter(measuresVaryingCorLong, Measure == 'CNT3') %>%
    mutate(CNT3 = Value) %>% select('Expectation', 'CNT3')) %>%

  mutate(Type = 'Inconsistent', Varying = 'All') ->

  measuresWideIA


## Figures

figureData <- filter(measuresVaryingCorLong,
                     grepl('CNT./', measuresVaryingCorLong[, 'Measure']))

varCorFig <- ggplot(figureData, aes(x = as.numeric(Expectation), y = Value)) +
  geom_point(aes(colour = Measure), size = 0.7) +
  facet_grid(Measure ~ ., scales = 'free_y') +
  xlab('Common expectation of the product') +
  scale_colour_brewer(palette="Dark2")

ggsave(filename = './output/PRlikeSystem_Inconsistent_0.55-0.40_Constant_correlation_Ratios.png', plot = varCorFig,
       width = 7, height = 5)

figureData <- filter(measuresVaryingCorLong,
                     grepl('^N?CNT.$', measuresVaryingCorLong[, 'Measure']))

measCorFig <- ggplot(figureData, aes(x = as.numeric(Expectation), y = Value)) +
  geom_point(aes(colour = Measure), size = 0.7) +
  facet_grid(Measure ~ ., scales = 'free_y') +
  xlab('Common expectation of the product') +
  scale_colour_brewer(palette="Dark2")

ggsave(filename = './output/PRlikeSystem_Inconsistent_0.55-0.40_Constant_correlation_Measures.png', plot = measCorFig,
       width = 7, height = 5)





################################################################################
# # Inconsistent 0.55 - 0.40. Varying one context
################################################################################


measuresVaryingCorLong <- generateSequencePRBoxlikeMeasures(nCors = nCors,
                                                            PRrank = 3L,
                                                            margin1 = 0.55,
                                                            margin2 = 0.40,
                                                            varying = 'one')

## Save results

write.table(measuresVaryingCorLong, file = './output/PRlikeSystem_Inconsistent_0.55-0.40_Varying_one_context.txt',
            sep = '\t', col.names = TRUE, row.names = FALSE)

#measuresVaryingCorLong <- read.table('./output/PRlikeSystem_Inconsistent_0.55-0.40_Varying_one_context.txt',
#                                     header = TRUE, sep = '\t')

# # Transform to wide table for scatterplots

dplyr::inner_join(dplyr::inner_join(
  dplyr::filter(measuresVaryingCorLong, Measure == 'CNT1') %>%
    mutate(CNT1 = Value) %>% select('Expectation', 'CNT1'),

  dplyr::filter(measuresVaryingCorLong, Measure == 'CNT2') %>%
    mutate(CNT2 = Value) %>% select('Expectation', 'CNT2')),

  dplyr::filter(measuresVaryingCorLong, Measure == 'CNT3') %>%
    mutate(CNT3 = Value) %>% select('Expectation', 'CNT3')) %>%

  mutate(Type = 'Inconsistent', Varying = 'One') ->

  measuresWideIO


## Figures

figureData <- filter(measuresVaryingCorLong,
                     grepl('CNT./', measuresVaryingCorLong[, 'Measure']))

varCorFig <- ggplot(figureData, aes(x = as.numeric(Expectation), y = round(Value, 10))) +
  geom_point(aes(colour = Measure), size = 0.7) +
  facet_grid(Measure ~ ., scales = 'free_y') +
  xlab('Expectation of the product in varying context') +
  ylab('Value') +
  scale_colour_brewer(palette="Dark2")

ggsave(filename = './output/PRlikeSystem_Inconsistent_0.55-0.40_Varying_one_context_Ratios.png', plot = varCorFig,
       width = 7, height = 5)

figureData <- filter(measuresVaryingCorLong,
                     grepl('^N?CNT.$', measuresVaryingCorLong[, 'Measure']))

measCorFig <- ggplot(figureData, aes(x = as.numeric(Expectation), y = round(Value, 10))) +
  geom_point(aes(colour = Measure), size = 0.7) +
  facet_grid(Measure ~ ., scales = 'free_y') +
  xlab('Expectation of the product in varying context') +
  ylab('Value') +
  scale_colour_brewer(palette="Dark2")

ggsave(filename = './output/PRlikeSystem_Inconsistent_0.55-0.40_Varying_one_context_Measures.png', plot = measCorFig,
       width = 7, height = 5)




# # Scatterplot and alike


measuresWide <- filter(dplyr::bind_rows(measuresWideCA, measuresWideCO,
                                        measuresWideIA, measuresWideIO),
                       CNT1 > 0)



measuresWide <- mutate(measuresWide,
                       `CNT1/CNT2` = round(CNT1 / CNT2, 5),
                       `CNT1/CNT3` = round(CNT1 / CNT3, 5),
                       `CNT2/CNT1` = round(CNT2 / CNT1, 5),
                       `CNT2/CNT3` = round(CNT2 / CNT3, 5),
                       `CNT3/CNT1` = round(CNT3 / CNT2, 5),
                       sign = ifelse(Varying == 'All', as.numeric(Expectation) > 0, FALSE),
                       group = paste(sign, Type, Varying),
                       Cases = paste(Type, Varying))


scatterCnt <- ggpairs(measuresWide, mapping = aes(colour = group),
                      columns = c('CNT1', 'CNT2', 'CNT3'),
                      diag = NULL,
                      upper = list(continuous = 'points'))


constantCNT1 <- data.frame(CNT1.start = measuresWide[c(14, 15), 'CNT1'],
                           CNT1.end   = measuresWide[c(44, 38), 'CNT1'],
                           CNT2.start = measuresWide[c(14, 15), 'CNT2'],
                           CNT2.end   = measuresWide[c(44, 38), 'CNT2'])

constantCNT2 <- data.frame(CNT1.start = measuresWide[c(6, 16), 'CNT1'],
                           CNT1.end   = measuresWide[c(11, 17), 'CNT1'],
                           CNT2.start = measuresWide[c(6, 16), 'CNT2'],
                           CNT2.end   = measuresWide[c(11, 17), 'CNT2'])


measuresWidePlot <- filter(measuresWide, Type == 'Consistent')
measuresWidePlot <- mutate(measuresWidePlot,
                           Case = recode(Cases,
                                          'Consistent All' = 'b',
                                          'Consistent One' = 'a'))

linesCnt <- ggplot(measuresWidePlot, aes(y = CNT1, x = CNT2))  +
  geom_point(aes(shape = Case), size = 2, alpha = 0.7) +
  labs(shape = '') +
  scale_shape_manual(values = c(19, 5))
slopes <- unique(measuresWide[, 'CNT1/CNT2'])

linesCnt <- linesCnt +
  geom_segment(aes(x = 0, y = 0,
                   xend = max(CNT2),
                   yend = max(CNT2) * slopes[1]),
               linetype = 1,
               alpha = 0.2) +
  geom_segment(aes(x = 0, y = 0,
                   xend = max(CNT2),
                   yend = max(CNT2) * slopes[2]),
               linetype = 2,
               alpha = 0.2)

linesCnt <- linesCnt +
  geom_segment(data = constantCNT1,
               aes(x = CNT2.start, y = CNT1.start,
                   xend = CNT2.end, yend = CNT1.end),
               linetype = 3,
               alpha = 0.6) +
  geom_segment(data = constantCNT2,
               aes(x = CNT2.start, y = CNT1.start,
                   xend = CNT2.end, yend = CNT1.end),
               linetype = 4,
               alpha = 0.6)


ggsave(filename = './output/PRlikeSystem_CNT_Scatterplots.png', plot = scatterCnt,
       width = 7, height = 5)


ggsave(filename = './output/PRlikeSystem_CNT_Lineplots.png', plot = linesCnt,
       width = 6, height = 4.5)

