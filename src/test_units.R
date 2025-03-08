################################################################################
# # Measures of contextuality and noncontextuality
# #
# # Example data and test units
# #
# # Víctor H Cervantes
# # 2019
################################################################################

source('./src/ContextualityMeasures.R')

# Structure data for tests

data("Titanic")
data("iris")

badData <- list(prop.table(Titanic),
                prop.table(table(Sepal = iris$Sepal.Length, Petal = iris$Petal.Length)),
                prop.table(Titanic),
                prop.table(table(Petal = iris$Petal.Width, Sepal = iris$Sepal.Width)))

class(badData) <- 'RVSystem'

goodData <- list(prop.table(Titanic),
                prop.table(table(Sepal = iris$Sepal.Length, Petal = iris$Petal.Length)),
                prop.table(Titanic),
                prop.table(table(Petal = iris$Petal.Length, Sepal = iris$Sepal.Length)))

class(goodData) <- 'RVSystem'



test_that('badData fails structure test',
          expect_false(isRVSystem(x = badData)))

test_that('goodData passes structure test',
          expect_true(isRVSystem(x = goodData)))




# # Dichotomous

maxAllowed <- as.table(matrix(c(cos(pi / 2), 1 - cos(pi / 2), 1 - cos(pi / 2), cos(pi / 2)) / 2, ncol = 2))
dimnames(maxAllowed)[[1]] <- dimnames(maxAllowed)[[2]] <- c('0', '1')

Ctx1 <- maxAllowed
Ctx2 <- maxAllowed
Ctx3 <- maxAllowed
Ctx4 <- maxAllowed

dimnames(Ctx4)[[1]] <- dimnames(Ctx4)[[2]][2:1] <- c('0', '1')

names(dimnames(Ctx1)) <- c('A1', 'B1')
names(dimnames(Ctx2)) <- c('A1', 'B2')
names(dimnames(Ctx3)) <- c('A2', 'B1')
names(dimnames(Ctx4)) <- c('A2', 'B2')

AliceBob <- RVSystem(Ctx1, Ctx2, Ctx3, Ctx4)


test_that('AliceBob passes dichotomous structure test',
          expect_true(isDRVSystem(x = AliceBob)))

test_that('goodData fails dichotomous structure test',
          expect_false(isDRVSystem(x = goodData)))





# # cntMatrices

# # Cyclic Rank 4

aliceCnt <- cntMatrices(AliceBob)

test_that('AliceBob cntMatrices generated object is of class cntMatrices',
          expect_true(isCntMatrices(aliceCnt)))


# # Single content  three contexts

ctxt1 <- prop.table(table(V1 = c(rep(1, 3), rep(0, 7)), V2 = c(rep(1, 3), rep(0, 7))))
ctxt2 <- prop.table(table(V1 = c(rep(1, 5), rep(0, 5)), V2 = c(rep(1, 2), rep(0, 8))))
ctxt3 <- prop.table(table(V1 = c(rep(1, 7), rep(0, 3)), V2 = c(rep(1, 7), rep(0, 3))))

severalCtxts <- RVSystem(ctxt1, ctxt2, ctxt3)

severalMatrices <- cntMatrices(severalCtxts)

test_that('Example of connection with more than two variables - p_l',
          expect_true(all(severalMatrices$lowMargins == c(1, .3, .3, .5, .2, .7, .7))))
test_that('Example of connection with more than two variables - p_c',
          expect_true(all(severalMatrices$connectionProbabilities == c(.3, .3, .5, .2, .3, .2))))
test_that('Example of connection with more than two variables - p_b',
          expect_true(all(severalMatrices$bunchProbabilities == c(.3, .2, .7))))


test_that('Example of connection with more than two variables - det M_l',
          expect_true(det(matprod_simple_triplet_matrix(severalMatrices$cntMatrixL,
                                                        t(severalMatrices$cntMatrixL))) != 0))
test_that('Example of connection with more than two variables - det M_c',
          expect_true(det(matprod_simple_triplet_matrix(severalMatrices$cntMatrixC,
                                                        t(severalMatrices$cntMatrixC))) != 0))
test_that('Example of connection with more than two variables - det M_d',
          expect_true(det(matprod_simple_triplet_matrix(severalMatrices$cntMatrixB,
                                                        t(severalMatrices$cntMatrixB))) != 0))


test_that('Example of connection with more than two variables - rowSums M_l',
          expect_true(all(c('32', '64') %in% names(table(row_sums(severalMatrices$cntMatrixL))))))

test_that('Example of connection with more than two variables - rowSums M_c',
          expect_true(all('16' == names(table(row_sums(severalMatrices$cntMatrixC))))))

test_that('Example of connection with more than two variables - rowSums M_b',
          expect_true(all('16' == names(table(row_sums(severalMatrices$cntMatrixB))))))



test_that('Example of connection with more than two variables - distribution rowsums M_l',
          expect_true(all(table(col_sums(severalMatrices$cntMatrixL)) == choose(6, 0:6))))


lowMatrix <- t(expand.grid(c(1L, 0L), c(1L, 0L), c(1L, 0L), c(1L, 0L), c(1L, 0L), c(1L, 0L), 1L))[7:1, ]

test_that('Example of connection with more than two variables - distribution rowsums M_l',
          expect_true(all(lowMatrix == as.matrix(severalMatrices$cntMatrixL))))


# #

# # System with a variable that appears on a single context

ctxt1 <- prop.table(table(V1 = c(rep(1, 3), rep(0, 7))))
ctxt2 <- prop.table(table(V1 = c(rep(1, 4), rep(0, 6))))
ctxt3 <- prop.table(table(V1 = c(rep(1, 5), rep(0, 5))))

singleCtxt <- RVSystem(ctxt1, ctxt2, ctxt3)

singleMatrices <- cntMatrices(singleCtxt)

lowMatrix <- t(expand.grid(c(1L, 0L), c(1L, 0L), c(1L, 0L), 1L))[4:1, ]

test_that('Example of connection with more than two variables - distribution rowsums M_l',
          expect_true(all(lowMatrix == as.matrix(singleMatrices$cntMatrixL))))


# # System with just a bunch

ctxt1 <- prop.table(table(V1 = c(rep(1, 3), rep(0, 7))))
ctxt2 <- prop.table(table(V1 = c(rep(1, 5), rep(0, 5)), V2 = c(rep(1, 2), rep(0, 8))))

aloneCtxt <- RVSystem(ctxt1, ctxt2)
aloneMatrices <- cntMatrices(aloneCtxt)

