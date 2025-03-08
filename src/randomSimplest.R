################################################################################
# # Measures of contextuality and noncontextuality
# #
# # Additional examples of
# # C-C systems and contextuality analysis
# #
# # Víctor H Cervantes
# # 2020
################################################################################



Fill2by2Table <- function (myTable, myVector) {
    if (!is.table(myTable)) {
        stop('myTable must be a table')
    }
    if (!all(dim(myTable) == 2)) {
        stop('myTable must be a 2 by 2 table')
    }
    if (length(myVector) != 4) {
        stop('myVector must have 4 elements')
    }

    myTable[1, 1] <- as.numeric(myVector[1])
    myTable[1, 2] <- as.numeric(myVector[2])
    myTable[2, 1] <- as.numeric(myVector[3])
    myTable[2, 2] <- as.numeric(myVector[4])

    return(myTable)
}

# # Simplest
cornerPatterns <- matrix(c(0L, 1L, 0L, 1L), nrow = 2)
context1 <- MakeJointFromPatterns(cornerPatterns)
names(dimnames(context1)) <- paste0("V", c(1, 2))



probVectors <- expand.grid(A = seq(0, 1, by = 0.01),
                           B = seq(0, 1, by = 0.01),
                           C = seq(0, 1, by = 0.01))
probVectors[, 'D'] <- round(1 - rowSums(probVectors), 3)
probVectors <- probVectors[probVectors[, 'D'] >= 0, ]

set.seed(2048)


nDist <- nrow(probVectors)
nSamples <- 500

cbdRandom <- list()
sysRandom <- list()

for (iiRank in seq(from = 2, to = 6)) {
    if (iiRank == 2) {
        cbdRandom[[as.character(iiRank)]] <- data.frame(CbDMeasures(RVSystem(context1, context1)))
    } else {
        cbdRandom[[as.character(iiRank)]] <- data.frame(CbDMeasures(RVSystem(context1, context1, context1)))
    }
    iiRow <- 1
    for (iiRow in seq(nSamples)) {

        contextList <- list()

        for (jjContext in seq(iiRank)) {
            contextList[[jjContext]] <- Fill2by2Table(context1, probVectors[sample(nDist, 1), ])
        }

        sysRandom[[as.character(iiRank)]][[iiRow]] <- RVSystem(contextList)
        cbdRandom[[as.character(iiRank)]][iiRow, ] <- data.frame(CbDMeasures(sysRandom[[as.character(iiRank)]][[iiRow]]))
    }
    # # # Ratios

    cbdRandom[[as.character(iiRank)]][, 'CNT1/CNT2'] <-
        round(cbdRandom[[as.character(iiRank)]][, 'CNT1'] /
                  cbdRandom[[as.character(iiRank)]][, 'CNT2'], 4)

    cbdRandom[[as.character(iiRank)]][, 'CNT3/CNT2'] <-
        round(cbdRandom[[as.character(iiRank)]][, 'CNT3'] /
                  cbdRandom[[as.character(iiRank)]][, 'CNT2'], 4)

    cbdRandom[[as.character(iiRank)]][, 'CNFR/CNT2'] <-
        round(cbdRandom[[as.character(iiRank)]][, 'CNFR'] /
                  cbdRandom[[as.character(iiRank)]][, 'CNT2'], 4)

    cbdRandom[[as.character(iiRank)]][, 'CNFR/CNT3'] <-
        round(cbdRandom[[as.character(iiRank)]][, 'CNFR'] /
                  cbdRandom[[as.character(iiRank)]][, 'CNT3'], 4)

    cbdRandom[[as.character(iiRank)]][, 'Rank'] <- iiRank

    cbdRandom[[as.character(iiRank)]] %>%
        select(matches("Rank") | matches("/")) ->
        cbdRandom[[as.character(iiRank)]]

}

cbdRandom <- bind_rows(cbdRandom)

cbdRandom %>% group_by_all() %>%
    summarise(.groups = 'keep', n= n()) ->
    cbdRandomN

cbdRandomN %>% group_by(Rank) %>%
    transmute(n = n, contextual = !is.nan(.data[["CNT1/CNT2"]])) %>%
    group_by(Rank, contextual) %>%
    summarise(.groups = 'drop', n = sum(n)) %>%
    group_by(Rank) %>%
    summarise(.groups = 'keep',
              contextual = 100 * last(n) / sum(n)) ->
    cbdRandomP

cbdRandomN %>% group_by(Rank) %>%
    filter(!is.nan(.data[["CNT1/CNT2"]])) %>%
    select(!n) %>%
    summarise_all(list(min = min, max = max, n = n_distinct)) %>%
    select(Rank, matches("CNT1/"), matches("CNT3/"), matches("CNFR/CNT2"), matches("CNFR/CNT3")) ->
    cbdRandomS

cbdRandomT <- merge(cbdRandomP, cbdRandomS)
