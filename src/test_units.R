################################################################################
# # Measures of contextuality and noncontextuality
# #
# # Example data and test units
# #
# # Víctor H Cervantes
# # 2019
################################################################################

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

maxAllowed <- as.table(matrix(c(cos(pi / 4), 1 - cos(pi / 4), 1 - cos(pi / 4), cos(pi / 4)) / 2, ncol = 2))
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

aliceCnt <- cntMatrices(AliceBob)

test_that('AliceBob cntMatrices generated object is of class cntMatrices',
          expect_true(isCntMatrices(aliceCnt)))
