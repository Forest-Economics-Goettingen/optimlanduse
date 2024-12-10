##--#####################################################--##
#### Transform the input table in an optimLanduse object ####
##--#####################################################--##

# Tue Jul  5 17:18:58 2022 ------------------------------

# Main developer: Kai Husmann

#' Initialize the robust optimization
#'
#' The function initializes an \emph{optimLanduse} S3 object on the
#' basis of a coefficients table. Please note that the coefficients table must follow
#' the expected \emph{optimLanduse} format. The expected format is explained in the example on the
#' \href{https://github.com/Forest-Economics-Goettingen/optimLanduse/}{GitHub project page} and in the publication in Methods in Ecology and Evolution (Husmann et al. ,2022)
#'
#'  Separating the initialization from the optimization is to save
#'  computation time in batch analysis. The separated function calls allow the
#'  user to perform multiple
#'  optimizations from one initialized object. This could save time in the scenario or
#'  sensitivity analysis.
#'
#'  A detailed description of the input parameters can be found in Husmann et al. (2022).
#'
#' @param coefTable Coefficient table in the expected \emph{optimLanduse} format.
#' @param uValue \emph{u} Value. The uncertainty value delivered in the coefTable is
#' multiplied with this u value. The value, therefore, enables scenario analyses with differing
#' uncertainties in relation to indicator values. Higher u values can be interpreted as a higher
#' risk aversion of the decision maker.
#' @param optimisticRule Either \emph{expectation} or \emph{uncertaintyAdjustedExpectation}.
#' The rule indicates whether the optimistic outcomes of an indicator are directly
#' reflected by their expectations or if the indicator is calculated as expectation +
#' uncertainty when "more is better" or expectation - uncertainty respectively when "less is better".
#' An optimization based on \emph{expectation} considers only downside risks.
#' @param fixDistance This optional numeric value allows to define distinct uncertainty levels for the
#' calculation of the uncertainty space and the averaged distances of a certain land-cover composition
#' (see Equation 9 in Husmann et al. (2020)). Passing NA disables fixDistance. In this case,
#' the uncertainty space is defined by uValue.
#' @return An initialized optimLanduse S3 object ready for optimization.
#' @references Husmann, K., von Groß, V., Bödeker, K., Fuchs, J. M., Paul, C., & Knoke, T. (2022). optimLanduse: A package for multiobjective land-cover composition optimization under uncertainty. Methods in Ecology and Evolution, 00, 1– 10. https://doi.org/10.1111/2041-210X.14000
#' @examples
#' require(readxl)
#' require(tidyr)
#' require(dplyr)
#' coefTable <- read_xlsx(exampleData("exampleGosling.xlsx"))
#' landUseRestriction <- c("Crops", "Pasture", "Alley Cropping", "Silvopasture")
#'
#' init <- initScenario(dat,
#'                      uValue = 2,
#'                      optimisticRule = "expectation",
#'                      fixDistance = 3)

#' @import dplyr
#' @import tidyr
#' @importFrom stats setNames
#'
#' @export
initScenario <- function(coefTable,  uValue = 1,
                         optimisticRule = "expectation", fixDistance = 3,
                         landUseRestriction = NA) {

  #-----------------------------------------#
  #### Check the format of the coefTable ####
  #-----------------------------------------#

  if (!all(c("indicator", "direction", "landUse", "indicatorValue", "indicatorUncertainty") %in% names(coefTable))) {
    stop ("At least one necessary variable for the optimization is not available. Are the requirements of the data structure met? Check the variable names.")
  }

  ## Drop unnecessary colums ##
  if(any(!(names(coefTable) %in%
           c("indicator", "direction", "landUse", "indicatorValue",
             "indicatorUncertainty", "indicatorGroup")))) {
    warning("Non-necessary columns detected and neglected.")
    coefTable <- coefTable[, names(coefTable) %in% c("indicator", "direction", "landUse", "indicatorValue",
                                                     "indicatorUncertainty", "indicatorGroup")]
  }

  indicatorNames <- as.character(unique(coefTable$indicator))

  testLandUseIndicators <- function (x) {
    all(indicatorNames %in% x)
  }

  checkLanduseTemp <- stats::aggregate(indicator ~ landUse, FUN = testLandUseIndicators, data = coefTable)


  if (!all(checkLanduseTemp$indicator)) {
    stop ("At least one indicator is not available for at least one land-use option.")
  }
  if (!length(indicatorNames) * length(unique(coefTable$landUse)) == nrow(coefTable)) {
    stop ("The indicator names are not unique. Have you assigned an indicator name twice?")
  }

  #----------------------------#
  #### Initialise the table ####
  #----------------------------#

  landUse <- as.character(unique(coefTable$landUse))

  expandList <- list()
  expandList[landUse] <- list(c("High", "Low"))

  expandMatrix1 <- as.matrix(expand.grid(expandList, stringsAsFactors = FALSE))
  expandMatrix2 <- do.call(rbind, replicate(length(indicatorNames), expandMatrix1, simplify = FALSE))
  scenarioTable <- tibble(indicator = rep(indicatorNames, each = dim(expandMatrix1)[1])) %>%
    bind_cols(as_tibble(expandMatrix2))

  names(scenarioTable)[names(scenarioTable) %in% landUse] <- paste0("outcome",names(scenarioTable)[names(scenarioTable) %in% landUse])

  #--------------------#
  ## Attach direction ##
  #--------------------#

  scenarioTableTemp1 <- scenarioTable
  scenarioTable <- merge(scenarioTable, unique(coefTable[, c("indicator", "direction")]), by = "indicator")
  if(!dim(scenarioTableTemp1)[1] == dim(scenarioTable)[1]) {cat("Error: Direction mising or wrong.")}

  #---------------------------------------------#
  ## Attach indicator values and uncertainties ##
  #---------------------------------------------#

  scenarioTableTemp2 <- scenarioTable

  spread1 <- tidyr::spread(data = coefTable[, !names(coefTable) == "indicatorUncertainty"],
                           key = landUse,
                           value = "indicatorValue")
  names(spread1)[names(spread1) %in% eval(landUse)] <- paste0("mean", names(spread1)[names(spread1) %in% eval(landUse)])


  spread2 <- tidyr::spread(data = coefTable[, !names(coefTable) == "indicatorValue"],
                           key = landUse,
                           value = "indicatorUncertainty")
  names(spread2)[names(spread2) %in% eval(landUse)] <- paste0("sem", names(spread2)[names(spread2) %in% eval(landUse)])


  for(i in landUse) {
    byIndicator <- c("indicator")
    names(byIndicator) <- "indicator"
    scenarioTable <- left_join(scenarioTable, spread1[, c("indicator", paste0("mean", i))], by = byIndicator)
    scenarioTable <- left_join(scenarioTable, spread2[, c("indicator", paste0("sem", i))], by = byIndicator)
  }


  scenarioTable <- scenarioTable %>% select(-contains("sem"), everything()) # Alternatively, but slower, a second loop would be suitable

  if(!dim(scenarioTableTemp1)[1] == dim(scenarioTable)[1]) {cat("Error: Attaching expectation or uncertainty failed.")}

  #--------------------------------------------#
  ## Calculate indicator uncertainty adjusted ##
  #--------------------------------------------#

  scenarioTableTemp3 <- scenarioTable

  # add Adjusted SEM to the scenarioTable

  scenarioTable <- addAdjSEM(scenarioTable = scenarioTable, landUse = landUse,
                             uValue = uValue, optimisticRule = optimisticRule)


  if (!(fixDistance >=0 & fixDistance <= 10)  & !is.na(fixDistance)) {
    fixDistance <- NA
    warning("The fixDistance did not meet the requirements and therefore set to NA. Please find the possible values for the fixDistance in the help.")
  }

  if ((fixDistance >=0 & fixDistance <= 10)  & !is.na(fixDistance)) {
    scenarioTableFix <- addAdjSEM(scenarioTable = scenarioTableTemp3,
                                  landUse = landUse,
                                  uValue = fixDistance,
                                  optimisticRule = optimisticRule)
  }


  if(!optimisticRule %in% c("uncertaintyAdjustedExpectation", "expectation")) {cat("optimisticRule must be uncertaintyAdjustedExpectation or expectation")}
  if(!dim(scenarioTableTemp3)[1] == dim(scenarioTable)[1] | any(is.na(scenarioTable))) {cat("Error: Calculation of adjusted uncertainty.")}

  #--------------------------#
  ## calculate Min Max Diff ##
  #--------------------------#


  if (is.na(fixDistance)) {
    scenarioTable[, c("minAdjSem", "maxAdjSem", "diffAdjSem")] <-
      apply(scenarioTable[, startsWith(names(scenarioTable), "adjSem")], 1,
            function(x) {c(min(x), max(x), (max(x) - min(x)))}) %>% t()
  } else {
    scenarioTable[, c("minAdjSem", "maxAdjSem", "diffAdjSem")] <-
      apply(scenarioTableFix[, startsWith(names(scenarioTableFix), "adjSem")], 1,
            function(x) {c(min(x), max(x), (max(x) - min(x)))}) %>% t()

  }

  #------------------------------#
  ## Include restriced data set ##
  #------------------------------#
  if(any(!is.na(landUseRestriction)) & !all(landUseRestriction %in% landUse)) {
    print(landUseRestriction[!landUseRestriction %in% landUse])
    stop("The landUseRestriction argument must be a subset of the landUse options.")
  }

  if(all(!is.na(landUseRestriction))) {

    #----------------------------#
    #### Initialise the table ####
    #----------------------------#

    landUse_restricted <- landUseRestriction
    coefTable_restricted <- coefTable %>%
      filter(landUse %in% landUse_restricted)


    expandList_restricted <- list()
    expandList_restricted[landUse_restricted] <- list(c("High", "Low"))

    expandMatrix1_restricted <- as.matrix(expand.grid(expandList_restricted, stringsAsFactors = FALSE))
    expandMatrix2_restricted <- do.call(rbind, replicate(length(indicatorNames), expandMatrix1_restricted, simplify = FALSE))
    scenarioTable_restricted <- tibble(indicator = rep(indicatorNames, each = dim(expandMatrix1_restricted)[1])) %>%
      bind_cols(as_tibble(expandMatrix2_restricted))

    names(scenarioTable_restricted)[names(scenarioTable_restricted) %in% landUse_restricted] <-
      paste0("outcome",names(scenarioTable_restricted)[names(scenarioTable_restricted) %in% landUse_restricted])

    #--------------------#
    ## Attach direction ##
    #--------------------#

    scenarioTableTemp1_restricted<- scenarioTable_restricted
    scenarioTable_restricted <- merge(scenarioTable_restricted, unique(coefTable_restricted[, c("indicator", "direction")]), by = "indicator")
    if(!dim(scenarioTableTemp1_restricted)[1] == dim(scenarioTable_restricted)[1]) {cat("Error: Direction mising or wrong.")}

    #---------------------------------------------#
    ## Attach indicator values and uncertainties ##
    #---------------------------------------------#

    scenarioTableTemp2_restricted <- scenarioTable_restricted

    spread1_restricted <- tidyr::spread(data = coefTable_restricted[, !names(coefTable_restricted) == "indicatorUncertainty"],
                             key = landUse,
                             value = "indicatorValue")
    names(spread1_restricted)[names(spread1_restricted) %in% eval(landUse_restricted)] <-
      paste0("mean", names(spread1_restricted)[names(spread1_restricted) %in% eval(landUse_restricted)])


    spread2_restricted <- tidyr::spread(data = coefTable_restricted[, !names(coefTable_restricted) == "indicatorValue"],
                             key = landUse,
                             value = "indicatorUncertainty")
    names(spread2_restricted)[names(spread2_restricted) %in% eval(landUse_restricted)] <-
      paste0("sem", names(spread2_restricted)[names(spread2_restricted) %in% eval(landUse_restricted)])


    for(j in landUse_restricted) {
      byIndicator <- c("indicator")
      names(byIndicator) <- "indicator"
      scenarioTable_restricted <- left_join(scenarioTable_restricted, spread1_restricted[, c("indicator", paste0("mean", j))], by = byIndicator)
      scenarioTable_restricted <- left_join(scenarioTable_restricted, spread2_restricted[, c("indicator", paste0("sem", j))], by = byIndicator)
    }

    scenarioTable_restricted <- scenarioTable_restricted %>% select(-contains("sem"), everything()) # Alternatively, but slower, a second loop would be suitable

    if(!dim(scenarioTableTemp1_restricted)[1] == dim(scenarioTable_restricted)[1]) {cat("Error: Attaching expectation or uncertainty failed.")}

    #--------------------------------------------#
    ## Calculate indicator uncertainty adjusted ##
    #--------------------------------------------#

    scenarioTableTemp3_restricted <- scenarioTable_restricted

    # add Adjusted SEM to the scenarioTable

    scenarioTable_restricted <- addAdjSEM(scenarioTable = scenarioTable_restricted,
                                         landUse = landUse_restricted, uValue = uValue,
                                         optimisticRule = optimisticRule)


    if (!(fixDistance >=0 & fixDistance <= 10)  & !is.na(fixDistance)) {
      fixDistance <- NA
      warning("The fixDistance did not meet the requirements and therefore set to NA. Please find the possible values for the fixDistance in the help.")
    }

    if ((fixDistance >=0 & fixDistance <= 10)  & !is.na(fixDistance)) {
      scenarioTableFix_restriced <- addAdjSEM(scenarioTable = scenarioTableTemp3_restricted,
                                    landUse = landUse_restricted,
                                    uValue = fixDistance,
                                    optimisticRule = optimisticRule)
    }


    if(!optimisticRule %in% c("uncertaintyAdjustedExpectation", "expectation")) {cat("optimisticRule must be uncertaintyAdjustedExpectation or expectation")}
    if(!dim(scenarioTableTemp3_restricted)[1] == dim(scenarioTable_restricted)[1] | any(is.na(scenarioTable_restricted))) {cat("Error: Calculation of adjusted uncertainty.")}

    #--------------------------#
    ## calculate Min Max Diff ##
    #--------------------------#


    if (is.na(fixDistance)) {
      scenarioTable_restricted[, c("minAdjSem", "maxAdjSem", "diffAdjSem")] <-
        apply(scenarioTable_restricted[, startsWith(names(scenarioTable_restricted), "adjSem")], 1,
              function(x) {c(min(x), max(x), (max(x) - min(x)))}) %>% t()
    } else {
      scenarioTable_restricted[, c("minAdjSem", "maxAdjSem", "diffAdjSem")] <-
        apply(scenarioTableFix_restriced[, startsWith(names(scenarioTableFix_restriced), "adjSem")], 1,
              function(x) {c(min(x), max(x), (max(x) - min(x)))}) %>% t()

    }


    for(k in unique(scenarioTable_restricted$indicator)){

      scenarioTable_restricted_temp <- scenarioTable_restricted %>%
        filter(indicator == k)

      for(l in seq_len(nrow(scenarioTable_restricted_temp))){

        scenario_temp_restricted <- scenarioTable_restricted_temp[l, startsWith(names(scenarioTable_restricted_temp), "outcome")]

        col_names_restricted <- names(scenario_temp_restricted)
        values_restricted <- as.character(unlist(scenario_temp_restricted))

        filtered_df <- scenarioTable[scenarioTable$indicator == k,]
        for (m in col_names_restricted) {
          filtered_df <- filtered_df %>% filter(!!sym(m) == values_restricted[match(m, col_names_restricted)])
        }

        values_restricted_mean <- mean(as.numeric(filtered_df[1, startsWith(names(filtered_df), "mean")]))
        filtered_df$test_beta <- ifelse(filtered_df$direction == "more is better",
                                        (filtered_df$maxAdjSem - values_restricted_mean)/filtered_df$diffAdjSem,
                                        (values_restricted_mean - filtered_df$minAdjSem)/filtered_df$diffAdjSem)

        filtered_df <- filtered_df %>%
          filter(test_beta == max(test_beta)) %>%
          filter(!duplicated(test_beta))

        scenarioTable_restricted[scenarioTable_restricted$indicator == k, ][l,"maxAdjSem"] <- filtered_df$maxAdjSem
        scenarioTable_restricted[scenarioTable_restricted$indicator == k, ][l,"minAdjSem"] <- filtered_df$minAdjSem
        scenarioTable_restricted[scenarioTable_restricted$indicator == k, ][l,"diffAdjSem"] <- filtered_df$diffAdjSem

      }
    }

    scenarioTable <- scenarioTable_restricted
    landUse <- landUse_restricted

  }

  # Check if all min / max / diff adj values are the same in restricted and unrestricted scenario
  # scenarioTable_check <- scenarioTable[scenarioTable$diffAdjSem %in% scenarioTable_restricted_check2,]
  #
  # scenarioTable_check2 <- scenarioTable_check %>%
  #   select(indicator, col_names_restricted, minAdjSem, maxAdjSem, diffAdjSem) %>%
  #   filter(!duplicated(.))
  #
  # scenarioTable_restricted_check <- scenarioTable_restricted %>%
  #   select(indicator, col_names_restricted, minAdjSem, maxAdjSem, diffAdjSem)
  #
  # scenarioTable_check2 <- scenarioTable_check2 %>% arrange(!!!syms(col_names_restricted))
  # scenarioTable_restricted_check <- scenarioTable_restricted_check %>% arrange(!!!syms(col_names_restricted))
  #
  # all(scenarioTable_restricted_check==scenarioTable_check2)

  #-------------------------------------------------------------#
  ## Define the coefficients for the linear objective function ##
  #-------------------------------------------------------------#

  #and the restrictions. (Simplify the scenario to a row problem)
  coefObjective <- defineObjectiveCoefficients(scenarioTable)

  #-------------------------------------#
  #### Define the constraints matrix ####
  #-------------------------------------#

  constraintCoefficients <- defineConstraintCoefficients(scenarioTable)

  retList <- list(scenarioSettings = data.frame(uValue = uValue,
                                                optimisticRule = optimisticRule, stringsAsFactors = FALSE),
                  scenarioTable = scenarioTable,
                  coefObjective = coefObjective,
                  coefConstraint = constraintCoefficients,
                  distance = scenarioTable[, c("minAdjSem", "maxAdjSem")],
                  status = "initialized",
                  beta = NA,
                  landUse = setNames(data.frame(matrix(rep(NA, length(landUse)), ncol = length(landUse), nrow = 1)), landUse),
                  optimDetails = list()
  )
  class(retList) <- "optimLanduse"
  return(retList)
}


