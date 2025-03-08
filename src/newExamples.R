################################################################################
# # Measures of contextuality and noncontextuality
# #
# # Additional examples of
# # C-C systems and contextuality analysis
# #
# # Víctor H Cervantes
# # 2020
################################################################################

library(ggplot2)
library(ggrepel)



source('./src/ContextualityMeasures.R')



# # Simplest
cornerPatterns <- matrix(c(0L, 1L, 0L, 1L), nrow = 2)
context1 <- context2 <- MakeJointFromPatterns(cornerPatterns)

names(dimnames(context1)) <- names(dimnames(context2)) <- paste0("V", c(1, 2))

jointPatternsSimplest <- expand.grid(V1 = gl(2, 1, labels = c(0L, 1L)),
                             V2 = gl(2, 1, labels = c(0L, 1L)))


nPatterns <- nrow(jointPatternsSimplest)
matSimplest <- list()
cbdSimplest <- data.frame(patterns = 0, CNT1 = 0, CNT2 = 0, CNT3 = 0, NCNT2 = 0)

row <- 1
for (iiSize in seq(nPatterns / 2)) {

    mmPatterns <- combn(nPatterns / 2, iiSize)

    for (jjPattern in seq(ncol(mmPatterns))) {

        jjPatternRows <- sort(c(mmPatterns[, jjPattern], 1 + nPatterns - mmPatterns[, jjPattern]))
        context3 <- MakeJointFromPatterns(jointPatternsSimplest[jjPatternRows, ])

        sysSimplest <- RVSystem(context1, context2, context3)
        matTemp     <- cntMatrices(sysSimplest)
        cbdTemp     <- CbDMeasures(matTemp)

        matSimplest[[row]]   <- matTemp
        cbdSimplest[row,  1] <- paste(rownames(jointPatternsSimplest)[jjPatternRows], collapse = '.')
        cbdSimplest[row, -1] <- as.vector(cbdTemp[c('CNT1', 'CNT2', 'CNT3', 'NCNT2')])
        row <- row + 1
    }

}





# # Uncycling cylcic Rank3
cornerPatterns <- matrix(c(0L, 1L, 0L, 1L), nrow = 2)
context1 <- context2 <- MakeJointFromPatterns(cornerPatterns)

names(dimnames(context1)) <- paste0("V", c(1, 2))
names(dimnames(context2)) <- paste0("V", c(2, 3))

jointPatternsUnCyc3 <- expand.grid(V1 = gl(2, 1, labels = c(0L, 1L)),
                             V2 = gl(2, 1, labels = c(0L, 1L)),
                             V3 = gl(2, 1, labels = c(0L, 1L)))

nPatterns <- nrow(jointPatternsUnCyc3)
matUnCyc3 <- list()
cbdUnCyc3 <- data.frame(patterns = 0, CNT1 = 0, CNT2 = 0, CNT3 = 0, NCNT2 = 0)

row <- 1
for (iiSize in seq(nPatterns / 2)) {

    mmPatterns <- combn(nPatterns / 2, iiSize)

    for (jjPattern in seq(ncol(mmPatterns))) {

        jjPatternRows <- sort(c(mmPatterns[, jjPattern], 1 + nPatterns - mmPatterns[, jjPattern]))
        context3 <- MakeJointFromPatterns(jointPatternsUnCyc3[jjPatternRows, ])

        sysUnCyc3 <- RVSystem(context1, context2, context3)
        matTemp     <- cntMatrices(sysUnCyc3)
        cbdTemp     <- CbDMeasures(matTemp)

        matUnCyc3[[row]]   <- matTemp
        cbdUnCyc3[row,  1] <- paste(rownames(jointPatternsUnCyc3)[jjPatternRows], collapse = '.')
        cbdUnCyc3[row, -1] <- as.vector(cbdTemp[c('CNT1', 'CNT2', 'CNT3', 'NCNT2')])
        row <- row + 1
    }

}



# # Ascending-Descending Staircase

cornerPatterns <- matrix(c(0L, 1L, 0L, 1L), nrow = 2)
context1 <- context2 <- context3 <- context4 <- MakeJointFromPatterns(cornerPatterns)

names(dimnames(context1)) <- paste0("V", c(1, 2))
names(dimnames(context2)) <- paste0("V", c(2, 3))
names(dimnames(context3)) <- paste0("V", c(3, 4))
names(dimnames(context4)) <- paste0("V", c(4, 1))

jointPatterns <- expand.grid(V1 = gl(2, 1, labels = c(0L, 1L)),
                             V2 = gl(2, 1, labels = c(0L, 1L)),
                             V3 = gl(2, 1, labels = c(0L, 1L)),
                             V4 = gl(2, 1, labels = c(0L, 1L)))

#jointPatterns <- jointPatterns[c(1, 2, 4, 6, 8, 9, 11, 13, 15, 16), ]
jointPatternsFS <- jointPatterns

nPatterns <- nrow(jointPatterns)

ADStaircaseFSSys <- list()
ADStaircaseFSCBD <- data.frame(patterns = 0, CNT1 = 0, CNT2 = 0, CNT2h = 0, level = 0, CNT3 = 0, NCNT2 = 0)

row <- 1
for (iiSize in seq(nPatterns / 2)) {

    mmPatterns <- combn(nPatterns / 2, iiSize)

    for (jjPattern in seq(ncol(mmPatterns))) {

        jjPatternRows <- sort(c(mmPatterns[, jjPattern], 1 + nPatterns - mmPatterns[, jjPattern]))
        context5 <- MakeJointFromPatterns(jointPatterns[jjPatternRows, ])

        ADStaircase <- RVSystem(context1, context2, context3, context4, context5)
        matADStaircase <- cntMatrices(ADStaircase)
        cbdADStaircase <- CbDMeasures(matADStaircase)

        ADStaircaseFSSys[[row]]   <- matADStaircase
        ADStaircaseFSCBD[row,  1] <- paste(rownames(jointPatterns)[jjPatternRows], collapse = '.')
        ADStaircaseFSCBD[row, -1] <- as.vector(cbdADStaircase[c('CNT1', 'CNT2', 'CNT2h', 'levelC', 'CNT3', 'NCNT2')])
        row <- row + 1
    }

}




# # Plots Ascending-Descending staircase
tabADStaircase <- ADStaircaseFSCBD
tabADStaircase[, 'patterns'] <- gsub('(9|10|11|12|13|14|15|16)', '', tabADStaircase[, 'patterns'])
tabADStaircase[, 'patterns'] <- gsub('\\.\\.+', '.', tabADStaircase[, 'patterns'])
tabADStaircase[, 'patterns'] <- gsub('\\.$', '', tabADStaircase[, 'patterns'])
tabADStaircase[, 'patterns'] <- gsub('\\.', ', ', tabADStaircase[, 'patterns'])


# # Full plot
slopes <- sort(unique(round(tabADStaircase$CNT1 / tabADStaircase$CNT2, 6))[-1])

uniquePairsCNTx <- as.data.frame(unique(cbind(
    CNT1 = round(tabADStaircase$CNT1, 5),
    CNT2 = round(tabADStaircase$CNT2, 5)
)))

plotADSFull <- ggplot(uniquePairsCNTx, aes(x = CNT2, y = CNT1)) +
    geom_point(size = 2, alpha = 0.9) +

    geom_segment(aes(x = 0.5, y = 1,
                     xend = 1, yend = 1),
                 linetype = 3,
                 alpha = 0.6) +

    geom_segment(aes(x = 1, y = 1,
                     xend = 1, yend = 2),
                 linetype = 3,
                 alpha = 0.6)


ggsave('./output/newComparisonLines.png', plotADSFull,
       width = 7, height = 5)


plotADSFull <- plotADSFull +
        geom_segment(aes(x = 0, y = 0,
                         xend = 2 / slopes[1],
                         yend = 2),
                     linetype = 1,
                     size = .3,
                     alpha = 0.05) +

        geom_segment(aes(x = 0, y = 0,
                         xend = 2 / slopes[7],
                         yend = 2),
                     linetype = 1,
                     size = .3,
                     alpha = 0.05) +

        geom_segment(aes(x = 0, y = 0,
                         xend = 2 / slopes[9],
                         yend = 2),
                     linetype = 1,
                     size = .3,
                     alpha = 0.05) +

        geom_segment(aes(x = 0, y = 0,
                         xend = 2 / slopes[11],
                         yend = 2),
                     linetype = 1,
                     size = .3,
                     alpha = 0.05) +

ggsave('./output/newComparisonLines.png', plotADSFull,
       width = 7, height = 5)

plotADSFull <- plotADSFull +
    geom_text_repel(data = tabADStaircase[c(3, 7), ], aes(label = patterns),
                    nudge_x = .01, nudge_y = -.01) +

    geom_text_repel(data = tabADStaircase[c(14), ], aes(label = patterns),
                    nudge_x = -.01, nudge_y = .01)

ggsave('./output/newComparisonLinesWithLabels.png', plotADSFull,
       width = 7, height = 5)


# # Full plot - shades
slopes <- sort(unique(round(tabADStaircase$CNT1 / tabADStaircase$CNT2, 6))[-1])
plotADSFull <- ggplot(uniquePairsCNTx, aes(x = CNT2, y = CNT1)) +
    geom_point(size = 2, alpha = 0.2) +
    geom_point(data = uniquePairsCNTx[c(2, 3, 6), ],
               size = 2, alpha = 0.9) +

        geom_segment(aes(x = 0, y = 0,
                         xend = 2 / slopes[1],
                         yend = 2),
                     linetype = 1,
                     size = .3,
                     alpha = 0.05) +

        geom_segment(aes(x = 0, y = 0,
                         xend = 2 / slopes[7],
                         yend = 2),
                     linetype = 1,
                     size = .3,
                     alpha = 0.05) +

        geom_segment(aes(x = 0, y = 0,
                         xend = 2 / slopes[9],
                         yend = 2),
                     linetype = 1,
                     size = .3,
                     alpha = 0.05) +

        geom_segment(aes(x = 0, y = 0,
                         xend = 2 / slopes[11],
                         yend = 2),
                     linetype = 1,
                     size = .3,
                     alpha = 0.05) +

    geom_segment(aes(x = 0.5, y = 1,
                     xend = 1, yend = 1),
                 linetype = 3,
                 alpha = 0.6) +

    geom_segment(aes(x = 1, y = 1,
                     xend = 1, yend = 2),
                 linetype = 3,
                 alpha = 0.6)

ggsave('./output/newComparisonLinesShades.png', plotADSFull,
       width = 7, height = 5)

plotADSFull <- plotADSFull +
    geom_text_repel(data = tabADStaircase[c(3, 7), ], aes(label = patterns),
                    nudge_x = .01, nudge_y = -.01) +

    geom_text_repel(data = tabADStaircase[c(14), ], aes(label = patterns),
                    nudge_x = -.01, nudge_y = .01)

ggsave('./output/newComparisonLinesWithLabelsShades.png', plotADSFull,
       width = 7, height = 5)



# # Plot slides subset
tabADSlides <- tabADStaircase[c(1, 6, 7, 86, 254, 255), ]
plotADSSlides <- ggplot(tabADSlides, aes(x = CNT2, y = CNT1)) +
    geom_point(size = 2, alpha = 0.7) +

        geom_segment(aes(x = 0, y = 0,
                         xend = 2 / slopes[1],
                         yend = 2),
                     linetype = 1,
                     size = .3,
                     alpha = 0.05) +

        geom_segment(aes(x = 0, y = 0,
                         xend = 2 / slopes[5],
                         yend = 2),
                     linetype = 1,
                     size = .3,
                     alpha = 0.05) +

        geom_segment(aes(x = 0, y = 0,
                         xend = 2 / slopes[9],
                         yend = 2),
                     linetype = 1,
                     size = .3,
                     alpha = 0.05) +

        geom_segment(aes(x = 0, y = 0,
                         xend = 2 / slopes[11],
                         yend = 2),
                     linetype = 1,
                     size = .3,
                     alpha = 0.05) +

    geom_segment(aes(x = 1, y = 1.25,
                     xend = 1, yend = 2),
                 linetype = 3,
                 alpha = 0.6) +

    geom_segment(aes(x = 1, y = 2,
                     xend = 2, yend = 2),
                 linetype = 3,
                 alpha = 0.6)

ggsave('./output/newComparisonLinesSubset.png', plotADSSlides,
       width = 7, height = 5)

plotADSSlides <- plotADSSlides +
    geom_text_repel(data = tabADSlides[5:6, ],
                    aes(label = patterns),
                    nudge_x = .05, nudge_y = 0) +

    geom_text_repel(data = tabADSlides[c(2, 4), ],
                    aes(label = patterns),
                    nudge_x = .02, nudge_y = -.02) +

    geom_text_repel(data = tabADSlides[3, ],
                    aes(label = patterns),
                    nudge_x = -.02, nudge_y = .02)


ggsave('./output/newComparisonLinesSubsetWithLabels.png', plotADSSlides,
       width = 7, height = 5)


# # Plot slides subset - shades
tabADSlides <- tabADStaircase[c(1, 6, 7, 86, 254, 255), ]
plotADSSlides <- ggplot(uniquePairsCNTx, aes(x = CNT2, y = CNT1)) +
    geom_point(size = 2, alpha = 0.2) +
    geom_point(data = tabADSlides,
               size = 2, alpha = 0.9) +

        geom_segment(aes(x = 0, y = 0,
                         xend = 2 / slopes[1],
                         yend = 2),
                     linetype = 1,
                     size = .3,
                     alpha = 0.05) +

        geom_segment(aes(x = 0, y = 0,
                         xend = 2 / slopes[5],
                         yend = 2),
                     linetype = 1,
                     size = .3,
                     alpha = 0.05) +

        geom_segment(aes(x = 0, y = 0,
                         xend = 2 / slopes[9],
                         yend = 2),
                     linetype = 1,
                     size = .3,
                     alpha = 0.05) +

        geom_segment(aes(x = 0, y = 0,
                         xend = 2 / slopes[11],
                         yend = 2),
                     linetype = 1,
                     size = .3,
                     alpha = 0.05) +

    geom_segment(aes(x = 1, y = 1.25,
                     xend = 1, yend = 2),
                 linetype = 3,
                 alpha = 0.6) +

    geom_segment(aes(x = 1, y = 2,
                     xend = 2, yend = 2),
                 linetype = 3,
                 alpha = 0.6)

ggsave('./output/newComparisonLinesSubsetShades.png', plotADSSlides,
       width = 7, height = 5)

plotADSSlides <- plotADSSlides +
    geom_text_repel(data = tabADSlides[5:6, ],
                    aes(label = patterns),
                    nudge_x = .05, nudge_y = 0) +

    geom_text_repel(data = tabADSlides[c(2, 4), ],
                    aes(label = patterns),
                    nudge_x = .02, nudge_y = -.02) +

    geom_text_repel(data = tabADSlides[3, ],
                    aes(label = patterns),
                    nudge_x = -.02, nudge_y = .02)


ggsave('./output/newComparisonLinesSubsetWithLabelsShades.png', plotADSSlides,
       width = 7, height = 5)







###########   Cyclic systems ----

ccEx2 <- sys2Ex

ccEx2$Ctx1[1, 2] <- .15
ccEx2$Ctx1[2, 2] <- .35

ccEx2$Ctx2[1, 1] <- .35
ccEx2$Ctx2[1, 2] <- .05
ccEx2$Ctx2[2, 1] <- .15
ccEx2$Ctx2[2, 2] <- .45

CbDMeasures(ccEx2, as.df = TRUE)
unclass(DeltaCFunction(ccEx2))



##### Other Examples non-contextual inconsistently connected rank 2 system

icEx2 <- sys2Ex

icEx2$Ctx1[1, 1] <- .25
icEx2$Ctx1[1, 2] <- .10
icEx2$Ctx1[2, 1] <- .30
icEx2$Ctx1[2, 2] <- .35

icEx2$Ctx2[1, 1] <- .05
icEx2$Ctx2[1, 2] <- .35
icEx2$Ctx2[2, 1] <- .15
icEx2$Ctx2[2, 2] <- .45

CbDMeasures(icEx2, as.df = TRUE)
unclass(DeltaCFunction(icEx2))

