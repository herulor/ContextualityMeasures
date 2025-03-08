
consistifyRank2 <- function (rank2System) {
    if (!isDRVSystem(rank2System)) {
        stop("rank2System should be a system of dichotomous random variables.")
    }

    if (length(rank2System) != 2L) {
        stop("rank2System should be of rank 2")
    }

    # Setup consistified system

    namesProperties <- unique(unlist(lapply(rank2System,
                                            function (y) names(dimnames(y)))))

    Ctxt1 <- rank2System[[1]]
    Ctxt3 <- rank2System[[2]]

    names(dimnames(Ctxt1)) <- paste0(names(dimnames(Ctxt1)), "C1")
    names(dimnames(Ctxt3)) <- paste0(names(dimnames(Ctxt3)), "C2")

    Ctxt2 <- Ctxt4 <- Ctxt1

    newPropertiesNames <- unique(unlist(lapply(list(
        Ctxt1, Ctxt2, Ctxt3, Ctxt4
    ),
                                               function (y) names(dimnames(y)))))


    names(dimnames(Ctxt2)) <- subset(newPropertiesNames, grepl(namesProperties[2], newPropertiesNames))
    names(dimnames(Ctxt4)) <- subset(newPropertiesNames, grepl(namesProperties[1], newPropertiesNames))


    refValues        <- character(length = length(newPropertiesNames))
    ferValues        <- character(length = length(newPropertiesNames))
    names(ferValues) <- names(refValues) <- newPropertiesNames

    for (iiName in newPropertiesNames) {

        for (jjContext in paste0("Ctxt", 1:4)) {

            if (iiName %in% names(dimnames(get(jjContext)))) {

                refValues[iiName] <-  dimnames(get(jjContext))[[iiName]][2]
                ferValues[iiName] <-  dimnames(get(jjContext))[[iiName]][1]

                break
            }
        }
    }


    # Create distributions for contexts 2 and 4 from multimaximal couplings

    varsCtxt2 <- names(dimnames(Ctxt2))

    refRowIndex <- (which(dimnames(Ctxt2)[[varsCtxt2[1]]] == refValues[varsCtxt2[1]]))
    refColIndex <- (which(dimnames(Ctxt2)[[varsCtxt2[2]]] == refValues[varsCtxt2[2]]))

    if (grepl("C1", varsCtxt2[1])) {
        marRow <- apply(Ctxt1, which(names(dimnames(Ctxt1)) == varsCtxt2[1]), sum)
    } else {
        marRow <- apply(Ctxt1, which(names(dimnames(Ctxt3)) == varsCtxt2[1]), sum)
    }

    if (grepl("C1", varsCtxt2[2])) {
        marCol <- apply(Ctxt1, which(names(dimnames(Ctxt1)) == varsCtxt2[2]), sum)
    } else {
        marCol <- apply(Ctxt1, which(names(dimnames(Ctxt3)) == varsCtxt2[2]), sum)
    }


    Ctxt2[refRowIndex, refColIndex] <- min(
        marRow[refValues[varsCtxt2[1]]],
        marCol[refValues[varsCtxt2[2]]]
    )

    Ctxt2[refRowIndex, 3 - refColIndex] <- (
        marRow[refValues[varsCtxt2[1]]] -
            Ctxt2[refRowIndex, refColIndex]
    )

    Ctxt2[3 - refRowIndex, refColIndex] <- (
        marCol[refValues[varsCtxt2[2]]] -
            Ctxt2[refRowIndex, refColIndex]
    )

    Ctxt2[3 - refRowIndex, 3 - refColIndex] <- (
        1 -
        marRow[refValues[varsCtxt2[1]]] -
        marCol[refValues[varsCtxt2[2]]] +
            Ctxt2[refRowIndex, refColIndex]
    )


    varsCtxt4 <- names(dimnames(Ctxt4))

    refRowIndex <- (which(dimnames(Ctxt4)[[varsCtxt4[1]]] == refValues[varsCtxt4[1]]))
    refColIndex <- (which(dimnames(Ctxt4)[[varsCtxt4[2]]] == refValues[varsCtxt4[2]]))

    if (grepl("C1", varsCtxt4[1])) {
        marRow <- apply(Ctxt1, which(names(dimnames(Ctxt1)) == varsCtxt4[1]), sum)
    } else {
        marRow <- apply(Ctxt1, which(names(dimnames(Ctxt3)) == varsCtxt4[1]), sum)
    }

    if (grepl("C1", varsCtxt4[2])) {
        marCol <- apply(Ctxt1, which(names(dimnames(Ctxt1)) == varsCtxt4[2]), sum)
    } else {
        marCol <- apply(Ctxt1, which(names(dimnames(Ctxt3)) == varsCtxt4[2]), sum)
    }


    Ctxt4[refRowIndex, refColIndex] <- min(
        marRow[refValues[varsCtxt4[1]]],
        marCol[refValues[varsCtxt4[2]]]
    )

    Ctxt4[refRowIndex, 3 - refColIndex] <- (
        marRow[refValues[varsCtxt4[1]]] -
            Ctxt4[refRowIndex, refColIndex]
    )

    Ctxt4[3 - refRowIndex, refColIndex] <- (
        marCol[refValues[varsCtxt4[2]]] -
            Ctxt4[refRowIndex, refColIndex]
    )

    Ctxt4[3 - refRowIndex, 3 - refColIndex] <- (
        1 -
        marRow[refValues[varsCtxt4[1]]] -
        marCol[refValues[varsCtxt4[2]]] +
            Ctxt4[refRowIndex, refColIndex]
    )



    return(list(Ctxt1 = Ctxt1, Ctxt2 = Ctxt2, Ctxt3 = Ctxt3, Ctxt4 = Ctxt4))

}