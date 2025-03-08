################################################################################
# # Measures of contextuality and noncontextuality
# #
# # Additional examples of
# # C-C systems and contextuality analysis
# #
# # Víctor H Cervantes
# # 2020
################################################################################

set.seed(43)

# Three way

jointPatterns3W <- expand.grid(V1 = gl(2, 1, labels = c(0L, 1L)),
                               V2 = gl(2, 1, labels = c(0L, 1L)),
                               V3 = gl(2, 1, labels = c(0L, 1L)))

context1 <- MakeJointFromPatterns(jointPatterns3W)
context2 <- MakeJointFromPatterns(jointPatterns3W[c(2, 3, 5, 8), ])
context3 <- MakeJointFromPatterns(jointPatterns3W[c(1, 4, 6, 7), ])


sysEx83  <- RVSystem(context1, context2)
sysEx83B <- RVSystem(context2, context3)
sysEx83I <- RVSystem(context1, context1)


cbdEx83 <- bind_rows(lapply(list(sysEx83I, sysEx83, sysEx83B), CbDMeasures))
cbdEx83[, 'k'] <- 3
cbdEx83[, 'Structure'] <- c('I-I', 'I-P', 'P-P')

# Four way

jointPatterns4W <- expand.grid(V1 = gl(2, 1, labels = c(0L, 1L)),
                               V2 = gl(2, 1, labels = c(0L, 1L)),
                               V3 = gl(2, 1, labels = c(0L, 1L)),
                               V4 = gl(2, 1, labels = c(0L, 1L)))

context1 <- MakeJointFromPatterns(jointPatterns4W)
context2 <- MakeJointFromPatterns(jointPatterns4W[c(2, 3, 5, 8, 9, 12, 14, 15),  ]) # Pair-triple indep. not quadruple
context3 <- MakeJointFromPatterns(jointPatterns4W[-c(2, 3, 5, 8, 9, 12, 14, 15), ])
context4 <- MakeJointFromPatterns(jointPatterns4W[c(2, 3, 5, 8, 10, 11, 13, 16),  ]) # Only pair-indep
context5 <- MakeJointFromPatterns(jointPatterns4W[-c(2, 3, 5, 8, 10, 11, 13, 16),  ])


sysEx84I <- RVSystem(context1, context1)
sysEx84A <- RVSystem(context1, context2) # Indep v Triple-indep
sysEx84B <- RVSystem(context1, context4) # Indep v Pair-indep
sysEx84C <- RVSystem(context2, context3) # Triple v -Triple
sysEx84D <- RVSystem(context4, context5) # Pair v -Pair



cbdEx84 <- bind_rows(lapply(list(sysEx84I, sysEx84A, sysEx84C, sysEx84B, sysEx84D), CbDMeasures))

cbdEx84[, 'k'] <- 4
cbdEx84[, 'Structure'] <- c('I-I', 'I-T', 'T-T', 'I-P', 'P-P')



cbdExTab <- bind_rows(cbdEx83, cbdEx84)

initCols <- grep('(k|St)', names(cbdExTab))
cbdExTab <- bind_cols(cbdExTab[, initCols], cbdExTab[, -initCols])





# # Random systems

nSamples <- 1000

sysRandom3k <- list()
cbdRandom3k <- data.frame(CbDMeasures(sysEx83I))

sysRandom4k <- list()
cbdRandom4k <- data.frame(CbDMeasures(sysEx84I))

for (iiSample in seq(nSamples)) {
    context1 <- MakeJointFromPatterns(jointPatterns3W[sample(8, 100, replace = TRUE), ])
    context2 <- MakeJointFromPatterns(jointPatterns3W[sample(8, 100, replace = TRUE), ])
    sysRandom3k[[iiSample]] <- RVSystem(context1, context2)
    cbdRandom3k[iiSample, ] <- data.frame(CbDMeasures(sysRandom3k[[iiSample]]))

    context1 <- MakeJointFromPatterns(jointPatterns4W[sample(16, 100, replace = TRUE), ])
    context2 <- MakeJointFromPatterns(jointPatterns4W[sample(16, 100, replace = TRUE), ])
    sysRandom4k[[iiSample]] <- RVSystem(context1, context2)
    cbdRandom4k[iiSample, ] <- data.frame(CbDMeasures(sysRandom4k[[iiSample]]))

}

cbdRandom3k[, 'k'] <- 3
cbdRandom4k[, 'k'] <- 4
cbdRandomTab <- bind_rows(cbdRandom3k, cbdRandom4k)

cbdRandomTab %>%
    mutate('CNT1/CNT2' = round(CNT1 / CNT2, 5),
           'CNT3/CNT2' = round(CNT3 / CNT2, 5),
           'CNFR/CNT2' = round(CNFR / CNT2, 5),
           'CNFR/CNT3' = round(CNT3 / CNT2, 5)) %>%
    mutate('CNT1/CNT2h' = round(CNT1 / CNT2h, 5),
           'CNT2/CNT2h' = round(CNT2 / CNT2h, 5),
           'CNT3/CNT2h' = round(CNT3 / CNT2h, 5),
           'CNFR/CNT2h' = round(CNFR / CNT2h, 5),
           'CNFR/CNT3' = round(CNT3 / CNT2, 5)) %>%
   group_by(k, levelC, levelN) %>%
    select(k, levelC, levelN, matches('/')) %>%
    filter(!is.nan(.data[['CNT1/CNT2']])) %>%
    summarise_all(list(min = min, max = max, n = n_distinct)) %>%
    pivot_longer(matches('/'),
                 names_pattern = '(.*)_(.*)',
                 names_to = c('Stat', '.value')) ->
    cbdRandomContRatios

cbdRandomTab %>%
    mutate('NCNT2/NCNT2h' = round(NCNT2 / NCNT2h, 5),
           'NCNT2/NCNT2' = NCNT2/ NCNT2) %>%
   group_by(k, levelC, levelN) %>%
    select(k, levelC, levelN, matches('/')) %>%
    filter(!is.nan(.data[['NCNT2/NCNT2h']])) %>%
    summarise_all(list(min = min, max = max, n = n_distinct)) %>%
    pivot_longer(matches('/'),
                 names_pattern = '(.*)_(.*)',
                 names_to = c('Stat', '.value')) %>%
    filter(Stat != 'NCNT2/NCNT2') ->
    cbdRandomNonContRatios


cbdRandomRatios <- merge(cbdRandomContRatios, cbdRandomNonContRatios, all = TRUE)
