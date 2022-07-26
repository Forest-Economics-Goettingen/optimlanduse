
## optimLanduse

<h3>
<a name="menu"></a>
</h3>
<ul>
<li>
<a href="#1. Einleitung">1. Introduction</a>
</li>
<li>
<a href="#3. Input und Output">2. Detailed description of the functions’
in- and outputs</a>
</li>
<li>
<a href="#5. Beispielhafte Anwendung">3. Example application</a>
</li>
<li>
<a href="#6. Erweiterte Anwendung">4. Sophisticated application</a>
</li>
<li>
<a href="#7. Suggested">5. Suggested citation</a>
</li>
<li>
<a href="#8. Literatur">6. Literature</a>
</li>
</ul>
<h3>
<a name="1. Einleitung">1. Introduction</a>
</h3>

When managing ecosystems the question is how to simultaneously combat
biodiversity loss and maintain ecosystem functioning while increasing
human welfare. Designing multifunctional landscapes means deciding upon
the composition and configuration of land-cover types, given that
landscape patterns drive the landscape\`s ecological value and its
ability to satisfy societal demands. Land-cover allocation models,
particularly optimization approaches, have proven helpful in revealing
trade-offs between multiple objectives and services provided by
different land-cover compositions.

The R package **optimLanduse** provides functions for easy and
systematic applications of the robust multiobjective land-cover
composition optimization approach of Knoke et al. (2016). It includes
tools to determine the land-cover composition that best balances the
multiple functions and services a landscape can provide and tools for
deeper understanding and visualizing the contributions of the distinct
indicators. The method has been developed and applied in a series of
studies before. You find some examples in the
<a href="#8. Literatur">literature</a> section. The methodological
background of the approach can e.g. be found in Knoke et al. (2016) and
Husmann et al. (n.d). We here refer to the nomenclature of Husmann et
al. (n.d.). The package opens the approach of Knoke et al. (2016) to the
community of landscape planners and provides opportunities for
systematic or batch applications. To further enhance this, we have
designed a shiny dashboard for the package to get a quick idea of its
functionalities and the modeling philosophy of the overall approach, see
<http://134.76.17.50/optimlanduse_shiny/>.

<h3>
<a name="3. Input und Output">2. Package structure</a>
</h3>

This chapter provides a brief overview over the package functions (Fig.
1). For detailed information about methodological background, functions,
and workflow please also refer to Husmann et al. (n.d.). We further
refer the reader to the respective help pages of the package for more
information.

The stable version of the package can be installed using the CRAN
server. The development version can be found on the GitHub project page.

``` r
# If not already installed
install.packages("optimLanduse")
```

<p align="center">
<img width="781.6" height="452" src="./man/figures/flowchart.png">
</p>

*Fig. 1: Overview of the functions of the* ***optimLanduse*** *package.
Green diamonds: input and output data; blue rectangles: functions; gray
parallelograms: optional function settings.*

#### Initialization and Input

The *initScenario()* function combines the user settings with the data
into an *optimLanduse*-object ready for solving. The following input
data are required:

-   *coefTable*: The package is only capable of processing a
    long-oriented type of data structure (Table 1). All combinations of
    land-cover (landUse) alternatives and indicators have to be listed
    vertically. Each row must contain the average expectation, the
    uncertainty, and the direction of the respective land-cover and
    indicator combination. The column names of the table must follow the
    expected nomenclature displayed below. You also find this format in
    the built-in example tables **exampleGosling.xlsx** or
    **exampleEmpty.xlsx**. The help files of the *exampleData()* and
    *initScenario()* functions contain more details regarding the
    required data format. Also an empty template incl. the predefined
    headings can be accessed via *exampleData(“exampleEmpty.xlsx”)*. All
    further columns will be dropped if passed.

*Table 1: Example of the data set from Gosling et al. (2020) to
illustrate the required data structure.*
<p align="center">
<img width="673.4" height="298.2" src="./man/figures/exampleGraphic.png">
</p>

-   *uValue*: The argument for the uncertainty level
    (![f_u](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;f_u "f_u"),
    equation 4 in Husmann et al. (n.d.)). A higher uValue reflects a
    higher risk aversion of the decision maker. See the help file of the
    *initScenario* function for more details.

-   *optimisticRule*: Specifies whether the optimistic contributions of
    each indicator should be defined either directly by their average,
    or by their average plus their uncertainty (if more is better) or
    minus its uncertainty (if less is better). The former option is most
    frequently used in recent literature and therefore builds the
    default.

-   *fixDistance*: Optional numeric value that defines distinct
    uncertainty levels for the calculation of the uncertainty space and
    the averaged distances of a certain land-cover composition (see
    Equation 9 in Husmann et al. (n. d.)). Passing NA disables
    fixDistance. The uncertainty space is then defined by the uValue.

#### Solver and List with results

The *solveScenario()* function takes the initialized
*optimLanduse*-object and only a few optional solver-specific arguments.

-   *digitsPrecision*: Provides the only possibility for the user to
    influence the calculation time. As the solving process has no
    stochastic element, the calculation times depend practically on the
    number of digits calculated only.

-   *lowerBound* & *upperBound*: Optional bounds for the land-cover
    alternatives. The lower bound must be 0 or a vector with lower
    bounds in the dimension of the number of land-cover alternatives.
    The upper bound, respectively, 1 or a vector with upper bounds.
    Choosing 0 and 1 (the defaults) as boundaries for all decision
    variables, means that no land-cover alternative is forced into the
    portfolio and that no land-cover alternative is assigned a maximum.

The returned *list with results* contains different Information of the
optimization model. First the information of the *initScenario()*
function are displayed again in this list. These include:

-   *scenarioSettings*: Data frame with *uValue* and *optimisticRule*
    used.
-   *scenarioTable*: Data frame with one row for each combination of
    worst-case and best-case outcomes of all indicators (thus the number
    of rows equals the number of scenarios
    ![N_S](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;N_S "N_S")
    in Husmann et al (n.d.)). The columns contain relevant calculation
    steps of the optimization program. *adjSem\** are the uncertainty
    adjusted indicators
    (![R\_{liu}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;R_%7Bliu%7D "R_{liu}")
    in Husmann et al. (n.d.)). *minAdjSem* are the minimal uncertainty
    adjusted indicators
    ![min(R\_{liu})](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;min%28R_%7Bliu%7D%29 "min(R_{liu})")
    and
    ![max(R\_{liu})](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;max%28R_%7Bliu%7D%29 "max(R_{liu})")
    the maximal uncertainty adjusted indicators. diffAdjSem are the
    range between these uncertainties adjusted indicators
    ![\delta\_{\text{min,max}\_{iu}}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cdelta_%7B%5Ctext%7Bmin%2Cmax%7D_%7Biu%7D%7D "\delta_{\text{min,max}_{iu}}").
-   *coefObjective*: The optimization program is translated into a
    linear program with
    ![N_L](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;N_L "N_L")
    (number of land-cover options) coefficients for the inner solution.
    The inner solution is solved using the *lpSolveAPI()* function of
    the **lpSolveApi** package. *coefObjective* is the data frame that
    contains these coefficients.
-   *coefConstraing*: A data frame with the respective constraints for
    the inner solution solved by *lpsolveAPI()*.
-   *distances*: The distance of each scenario to its own theoretically
    best-achievable contribution (reference). See equation 3 in Husmann
    et al. (n.d.).

This is followed by a summary of the results of the optimization:

-   *β*: The maximum distance of the worst performing scenario (equation
    1 in Husmann et al. (n.d.)).
-   *landUse*: The resulting land-cover composition in the optimum.

#### Post-Processing

-   *calcPerformance()*: Attaches the portfolio performances of all
    indicators and scenarios as data frame. The data can be used for
    straightforward visualization of the performance (e.g. Fig. 3). The
    performance is defined as the distance to the maximum achievable
    level for each indicator and uncertainty scenario.

<h3>
<a name="5. Beispielhafte Anwendung">3. Example application</a>
</h3>

We here present the basic workflow on a literature example. Aim of this
chapter is to introduce the functionality of the packages’ functions and
to explain the relevant in- and output on the example of a use-case in
Eastern Panama. The data of this study is accessible in the *Appendix A*
of Gosling et al. (2020) and is also firmly integrated into the
**optimLanduse** package. It can be accessed via
*exampleData(“exampleGosling.xlsx”)*. The data integrated in the package
comes already in the required *optimLanduse* format, such that it can be
used without any data processing.

Enriching agricultural farms with agroforestry has been promoted as a
means to enhance ecosystem functioning of farms in Panama, while
maintaining important economic functions. Gosling et al. (2020)
therefore used the optimization model presented here to understand
smallholder farmers’ perceptions and values of agroforestry systems.
They identified 10 relevant indicators for a predefined set of
land-cover alternatives, which represent the farmers’ goals (such as
long and short-term income or labor demand, as well as carbon and water
regulation). A survey with local farmers provided the empirical basis in
the form of the farmer’s expectations on the indicator performance of
each land-cover (arithmetic mean) and its uncertainties (using the
standard error of the mean across the survey’s respondents).
Descriptions of the land-cover alternatives and indicators can be found
in tables 1 and 2 in Gosling et al. (2020).

### Installing **optimLanduse**, Loading Required Packages and Importing the Data

``` r
library(optimLanduse)
library(readxl)
library(ggplot2)
library(tidyverse)
library(ggsci)

# Loading the example data
path <- exampleData("exampleGosling.xlsx")
dat <- read_excel(path)
```

*dat* is in the required format. Refer to the help of the
*initScenario()* function or to the
<a href="#3. Input und Output">detailed description of the functions’
in- and outputs</a> chapter for more details.

### Initializing an *optimLanduse* Object

``` r
# Initializing an optimLanduse-object
init <- initScenario(coefTable = dat,
                     uValue = 2,
                     optimisticRule = "expectation", 
                     # optimistic contribution of each indicator directly defined by their average 
                     fixDistance = NA) 
# 3 is the default
```

In line with Gosling et al. (2020), we chose the expected value of the
indicator as optimistic outcomes *(optimisticRule = “expectation”)* and
the same uncertainty level for the calculation of the averaged distances
and the uncertainty space (*fixDistance = NA*, see equations 4 and 9 in
Husmann et al. (n.d.) for more details).

### Solving the Initialized *optimLanduse* Object

``` r
# Solve the initialized optimLanduse object using the solveScenario() function                     
result <- solveScenario(x = init)

# Visualize the farm composition
result$landUse %>% gather(key = landCoverOption, value = landCoverShare, 1 : 6) %>% 
  mutate(portfolio = "Optimal farm composition",
         landCoverShare = landCoverShare * 100) %>% 
  ggplot(aes(y = landCoverShare, x = portfolio, fill = landCoverOption)) + 
  geom_bar(position = "stack", stat = "identity") + 
  theme_classic() +
  theme(text = element_text(size = 14)) +
  scale_fill_startrek() +
  labs(y = "Allocated share (%)") +
  scale_y_continuous(breaks = seq(0, 100, 10), 
                     limits = c(0, 100)) +
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank()) + 
  guides(fill=guide_legend(title = ""))
```

<img src="README_files/figure-gfm/unnamed-chunk-4-1.png" style="display: block; margin: auto;" />

*Fig. 2: Composition of the optimized farm (based on data of Gosling et
al. (2020)), including all indicators. Each land-cover option is shown
in an allocated share (%).*

The resulting optimized farm composition (Fig. 2) corresponds to Fig. 3
(![f_u=2](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;f_u%3D2 "f_u=2"))
in Gosling et al. (2020). It can be seen that the farm composition that
best fulfills all 10 indicators (i.e. the multi-functional portfolio) is
dominated by silvopasture and forest. According to Gosling et
al. (2020), this reveals the potential of agroforestry to serve as a
compromise solution to fulfill multiple ecological and economic
functions. Recently, however, the average farm portfolio of the surveyed
farms was mainly composed of pasture and cropland with only a small
share of forest (14 %). This reveals that not all of the selected
objectives (and their weights) currently drive farmers’ land-cover
decisions. The optimization approach can then be used to dive deeper
into the effect of different goals on the resulting optimized land-cover
composition and the effects of uncertainty.

### Calculating the Portfolio Performances of the Optimized *optimLanduse* Object

``` r
# Performance calculations
performance <- calcPerformance(result)

performance$scenarioTable$performance <- performance$scenarioTable$performance * 100

ggplot(performance$scenarioTable,
       aes(x = indicator,
           y = performance,
           color = indicator)) +
  geom_point() +
  geom_hline(yintercept =
               min(performance$scenarioTable$performance),
             linetype = "dashed", color = "red") +
  guides(color = guide_legend(title = "",
                              nrow = 10)) +
  theme_classic() + 
  theme(text = element_text(size = 18),
        legend.position="right",
        axis.ticks.x = element_blank()) +
  scale_x_discrete(labels = seq(1, 10)) +
  labs(y = "Min-max normalized indicator value (%)",
       x = "Indicators") + 
  scale_y_continuous(breaks = seq(0, 101, 10), 
                     limits = c(0, 101)) +
  geom_hline(aes(yintercept=100), size = 1) + 
  annotate(geom = "Text", x = 6, y = 100, label = "Maximum achievable indicator level",
           vjust = -1)
```

<img src="README_files/figure-gfm/unnamed-chunk-5-1.png" style="display: block; margin: auto;" />

*Fig. 3: The performances of each of the 10 indicators for the ideal
farm composition. The colored points are the achieved levels of the
indicators of all scenarios. The dotted, horizontal red line illustrates
the guaranteed performances*
![(1-\beta)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%281-%5Cbeta%29 "(1-\beta)")*,
thus the robust feasible solution of the program (Equation 1 in Husmann
et al (n.d.)).*

<!-- *CP: Im Folgenden Absatz ist nicht ganz klar auf was Ihr Euch da bezieht. Ich denke ich würde für alle Grafiken eine Fig. caption einfügen, das macht es für den Leser deutlich leichter und Ihr könnt dann besser darauf verweisen.* -->

Fig. 3 can be used to further explore the effects of the indicators on
the modeled land-use decisions. Looking at the performances of this
multi-functional farm reveals, which indicator equals
![\beta](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cbeta "\beta")
and therefore defines the result (Equation 1 in Husmann et al. (n.d.)).
<!-- *(VVG: vielleicht: ... which indicators perform poorly and therefore restrict the solution. - das ist der Satz aus unserem TeX-Dokument)* (*CP: Wirklich restrict the result oder drive the largest distance (beta...) under specific uncertainty scenarios. Es stimmt natürlich, aber hier vielleicht noch technischer und mit mehr Bezug zur Formel arbeiten, das macht es denke ich leichter für den Leser*) -->
Here, the worst performing scenarios of indicators 1 (financial
stability), 3 (investment costs) and 8 (meeting household needs) have
equally the largest distances. It can be seen that the portfolio is
apparently driven by these 3 indicators. In the worst-performing
uncertainty scenarios, these 3 indicators show the maximum distances
across all indicators. In other words, the guaranteed performance
![1-\beta](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;1-%5Cbeta "1-\beta")
(lowestPerformance) is defined by these indicators.
<!-- (*CP: Die Sätze hier sind ein bisschen wiederholend, vielleicht nochmal straffen, aber im Zweifel lieber doppelt drin lassen, ist ja ein readme und hat keine Beschränkung*) (*VVG: Wenn wir im ersten Satz nicht direkt auf die Distanz eingehen, finde ich kann man es einfach so lassen, dann eist es nicht wiederholend sondern man steigert das Detail*); KH Ich finde, das ist keine Nachteil im Readme  -->
A full list with performance of all individual scenarios is provided by
the output *scenarioTable* after using the *calcPerformance()* function
(Table 2).

It can be followed that these 3 indicators are crucial when discussing
future land-cover alternatives and concepts. According to Gosling et
al. (2020), this result is in line with current observed behavior, since
the need for short-term liquidity mainly drives decisions of smallholder
farmers in the study region. Intermediate-term economic success is not
relevant until consumption of the household is secured. While the
performances of indicator 1 differ relatively strong among the
scenarios, the performances of indicators 3 and 8 are similar in the
scenarios. This is attributed to the larger standard errors of this
indicator. In contrast, Indicators 3 and 8 still show a low guaranteed
performance level. It thus may be worth investigating the particular
reasons for this high uncertainty of indicator 1.

<!-- *(VVG: Wie Carola direkt hier drunter beschreibt, drüfte es ja vor allem an den verschiedenen Unsicherheiten liegen oder? Wenn man u = 0 setzt gibt es ja gar keine Unterschiede zwischen den einzelnen Szenarien, dann liegt jeder dieser kreise übereinander, daher finde ich den ersten Satz von Carola ganz gut... Auch wenn ich darüber nachdenken muss, warum es nur "in parts" sein soll - wenn ich mich nicht täusche, dürfte es fast ausschließlich durch die verschiedenen standardfehler/unsicherheiten getrieben werden - und dann könnte es beispielweise interessant sein, dass zu untersuchen, also warum diese standardfehler zwischen den Indikatoren so unterschiedlich sind)* -->
<!-- *CP: Hab nicht ganz verstanden was Euch hier wichtig ist, vielleicht eine alternative Formulierung: Fig. 3 can be used to further explore the effects of uncertainty on the modeled land-use decisions. Indicator 1 shows a larger deviation in indicator performance pij (oder wie ihr es nennt; VVG: dürfte 1-diu sein oder? 1-beta ist guaranteed performance, 1-diu dürfte einfach die performance jedes szenarios sein) between uncertainty scenarios. This might be in parts be attributed to the larger standard errors of this indicator. In contrast, Indicators 3 and 8 still show a low guaranteed performance level (VVG: " a low guaranteed performance finde ich irreführend, da es ja die exakt gleiche guaranteed performance ist im worst-case szenario mit geringen differenzen zwischen den einzelnen szenarien), while differences between uncertainty scenarios are less prominent. ... oder so irgendwie: Such in-depth analyses help to identify trade-offs among different indicators and their likely effects on land-use decisions.....Soll hier nochmal gesagt werden wann die Grafik gemacht werden sollte und wann nicht? Vielleicht auch nicht unbedingt notwendig? VVG: Glaube wann die Grafik gemacht werden sollte, finde ich nicht so relevant hier im Text* -->

``` r
performance$beta
```

    ## [1] 0.6132

``` r
# lowestPerformance <- performance$scenarioTable[order(performance$scenarioTable$performance,
#                                                      decreasing = FALSE),]
# performanceExample <- head(lowestPerformance[,c(1 : 8, 31)], n = 8)
# 
# knitr::kable(performanceExample, row.names = F)

performanceExample <- head(performance$scenarioTable[,c(1 : 8, 31)], n = 4)

knitr::kable(performanceExample, row.names = F)
```

*Table 2: An extract of the scenario table of all indicators created
through the calcPerformance() function with the worst performing
scenarios*

| indicator           | outcomeCrops | outcomePasture | outcomeAlley Cropping | outcomeSilvopasture | outcomePlantation | outcomeForest | direction      | performance |
|:--------------------|:-------------|:---------------|:----------------------|:--------------------|:------------------|:--------------|:---------------|------------:|
| Financial stability | High         | High           | High                  | High                | High              | High          | more is better |    61.31992 |
| Financial stability | Low          | High           | High                  | High                | High              | High          | more is better |    72.31477 |
| Financial stability | High         | Low            | High                  | High                | High              | High          | more is better |    61.31992 |
| Financial stability | Low          | Low            | High                  | High                | High              | High          | more is better |    72.31477 |

### Comparison of the performance of the currently observed land-cover composition to the optimized composition

``` r
result_current <- solveScenario(x = init,
                        lowerBound = c(0.26, 0.59, 0, 0, 0.01, 0.14),
                        upperBound = c(0.26, 0.59, 0, 0, 0.01, 0.14))

performance_current <- calcPerformance(result_current)

performance_current$scenarioTable$performance <- performance_current$scenarioTable$performance * 100
```

``` r
performance_current$beta
```

    ## [1] 0.9114

``` r
ggplot(performance_current$scenarioTable,
       aes(x = indicator,
           y = performance,
           color = indicator)) +
  geom_point() +
  geom_hline(yintercept =
               min(performance_current$scenarioTable$performance),
             linetype = "dashed", color = "red") +
  guides(color = guide_legend(title = "",
                              nrow = 10)) +
  theme_classic() + 
  theme(text = element_text(size = 18),
        legend.position="right",
        axis.ticks.x = element_blank()) +
  scale_x_discrete(labels = seq(1, 10)) +
  labs(y = "Min-max normalized indicator value (%)",
       x = "Indicators") + 
  scale_y_continuous(breaks = seq(0, 101, 10), 
                     limits = c(0, 101)) +
  geom_hline(aes(yintercept=100), size = 1) + 
  annotate(geom = "Text", x = 6, y = 100, label = "Maximum achievable indicator level",
           vjust = -1)
```

<img src="README_files/figure-gfm/unnamed-chunk-11-1.png" style="display: block; margin: auto;" />

*Fig. 4: The performance of each of the 10 indicators for the result of
the currently observed land-cover composition. The colored points are
the achieved levels of the indicators of all scenarios s. The dotted,
horizontal red line illustrates the guaranteed performance*
![(1-\beta)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%281-%5Cbeta%29 "(1-\beta)")*,
thus the robust feasible solution of the program (Equation 1 in (Husmann
et al, n.d.)).*

The lower and upper bound arguments allow to model the current
land-cover composition (Fig. 4). It allows to compare the difference of
the optimized land-cover composition to the currently observed
composition. Comparing e.g. the guaranteed performances
![1-\beta](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;1-%5Cbeta "1-\beta")
of the scenarios provides objective measure how an optimization enhances
the achievements of all indicator. A deeper look at the performance of
the scenarios shows which indicators particularly benefit from
optimization. The guaranteed performance is way higher when optimized.
It increases by 0.3 …

The comparison of the performances of the currently observed land-cover
composition (Fig. 4) with the performances of the multi-functional
portfolio (Fig. 3) shows that the worst-case scenario of, e.g.,
liquidity and meeting household needs in the multi-functional portfolio
is way beyond the the worst-case scenario of the same indicators in the
currently observed land-cover composition.
<!-- Ich sehe da keinen Kontrast. Das müssen wir morgen noch Mal besprechen -->
In contrast, the performance of, e.g., the worst-case scenarios of
financial stability, protecting soil resources, and protecting water
supply shows a significant increase. Overall, the maximum distance of
the worst performing scenario of the currently observed land-cover
composition
(![\beta](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cbeta "\beta"))
is reduced by about 30 percentage points, from
![\beta](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cbeta "\beta")
= 0.9114 to
![\beta](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cbeta "\beta")
= 0.6132 (comparing the horizontal red line of Fig. 3 and Fig. 4). This
interaction shows that the optimization process can help to reduce
trade-offs between the different indicators.

<h3>
<a name="6. Erweiterte Anwendung">4. Sophisticated application</a>
</h3>

### Batch Analysis: Solving Multiple Uncertainty Values

``` r
applyDf <- data.frame(u = seq(0, 3, .5))

applyFun <- function(x) {
  init <- initScenario(dat, uValue = x, optimisticRule = "expectation", fixDistance = NA)
  result <- solveScenario(x = init)
  return(c(result$beta, as.matrix(result$landUse)))
}

applyDf <- cbind(applyDf,
                 t(apply(applyDf, 1, applyFun)))

names(applyDf) <- c("u", "beta", names(result$landUse))

applyDf[, c(3 : 8)] <- applyDf[, c(3 : 8)] * 100

applyDf %>% gather(key = "land-cover option", value = "land-cover share", -u, -beta) %>%
  ggplot(aes(y = `land-cover share`, x = u, fill = `land-cover option`)) + 
  geom_area(alpha = .8, color = "white") + theme_minimal()+
  labs(x = "Uncertainty level", y = "Allocated share (%)") + 
  guides(fill=guide_legend(title="")) + 
  scale_y_continuous(breaks = seq(0, 100, 10), 
                     limits = c(0, 100.01)) +
  scale_x_continuous(breaks = seq(0, 3, 0.5),
                     limits = c(0, 3)) + 
  scale_fill_startrek() +
  theme(text = element_text(size = 18),
        legend.position = "bottom")
```

![](README_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

*Fig. 5: Theoretically ideal farm compositions under increasing levels
of uncertainty.*

Solving the portfolio (Fig. 5) under increasing assumptions for the
uncertainty levels (uValue, respectively
![f_u](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;f_u "f_u")
in equation 4 (Husmann et al, n.d.)) gives the sensitivity of the
land-cover compositions towards increasing risk aversion of the farmer.
Fig. 5 corresponds to Fig. 3 in Gosling et al. (2020). The higher the
uncertainty level, the higher the uncertainty spaces of the indicators.
Here, the composition of land-cover alternatives is quite stable across
different uncertainty levels (Fig. 5). Comparing portfolios of u-value 0
with u-value 3, the share of forest decreases slightly from 41.2% to
34.3% and silvopasture from 45.9% to 44.9%. The share of crops decreases
from 12.9% to 1%. At the same time, the shares of pasture increased from
0% to 9.6% and of plantation from 0% to 10.2%.

Alley cropping does not appear in any portfolio at any uncertainty
level. Alley cropping does in average not contribute best to any
indicator. Each indicator is better provided by at least one another
land-cover type. Also contributes worst (highest) to management
complexity. This overall negative contribution does not change with
increasing uncertainty levels. By trend, higher uncertainty levels lead
to more diverse portfolios. The uncertainty spaces of all indicators
increase with increasing uncertainty levels. These broadened individual
uncertainty spaces then lead to a broader state space with a higher
number of possible candidates for lowest-performing scenarios
(i.e. scenarios that can under lower uncertainty not become part of the
solution, as their distances could not be the maximum distance of any
land-cover composition). Plantation, for example, is not part of the
portfolio till an uncertainty level of 1.5. It only provides best to the
long-term income, while providing by far worst to the general
preferences and also providing generally only minor contributions to the
indicators. Under an uncertainty level of 1, for example, plantation
provides worst to the general preferences even if all other indicators
are considered as worst-possible contributions. This ranking, however,
changes after uncertainty levels above 1.5. The worst-possible
contribution of forests then contributes worst among all land-cover
types.

<!-- *VVG: Habe hier im Absatz drüber mal genau die prozente angegeben, da Carola im vorherigen Absatz (direkt hier drunter) "Genau angeben" angemrkt hatte um von dem circa wegzukommen. Es viel mir jedoch schwer das ganze genau anzugeben ohne mehrere Landnutzungsformen aufzuzählen. Vielleicht aber zu genau.* -->
<!-- *Only the small proportion of crops (circa 10 %* ***CP:Genau angeben****) is interchanged with small proportions of pasture and plantation under uncertainty levels larger than 2. The share of forest and silvopasture is hardly affected by the uncertainty level.*  -->

### Selecting Specific Indicator Bundles - Investigate the Indicators’ Sensitivity

<!-- *Soll es bei "leave-indicators-out-analysis" bleiben oder ginge auch etwas wie .. selection of specific indicator (bundles) oder "in and exclusion of indicators"* -->

The sensitivity of the land-cover compositions towards indicators or
groups of indicators can be analysed by either excluding od adding
indicators of interest and interpreting the differences in the results
of the distinct optimization. To do so individual and independent
optimization runs are carried out in and excluding different (sets of)
indicators. The set of indicators considered is representative for the
stakeholders’ preferences and perceptions. Comparison of optimal
land-cover compositions under differing indicator combinations may help
to understand how preferences of stakeholders design the land-cover
compositions. The following code exemplifies such optimizations for 3
subsets of indicators presented in Gosling et al. (2020). The shiny app
of *optimLanduse* (<http://134.76.17.50/optimlanduse_shiny/>) provides
functionality to straightforwardly define sets of indicators by a single
click. Further explanation and instructions are given in the app.

``` r
dat_socioeconomic <- dat[dat$indicator != "Protecting soil resources" & dat$indicator !="Protecting water supply",]

init_socioeconomic <- initScenario(dat_socioeconomic,
                                   uValue = 2,
                                   optimisticRule = "expectation", 
                                   fixDistance = NA) 

result_socioeconomic <- solveScenario(x = init_socioeconomic)

result_socioeconomic$landUse %>% gather(key = landCoverOption, value = landCoverShare, 1 : 6) %>% 
  mutate(portfolio = "Socio-economic",
         landCoverShare = landCoverShare * 100) %>% 
  ggplot(aes(y = landCoverShare, x = portfolio, fill = landCoverOption)) + 
  geom_bar(position = "stack", stat = "identity") + 
  theme_classic() +
  theme(text = element_text(size = 14)) +
  scale_fill_startrek() +
  labs(y = "Allocated share (%)") +
  scale_y_continuous(breaks = seq(0, 100, 10), 
                     limits = c(0, 100)) +
  theme(axis.title.x=element_blank(),
        axis.ticks.x=element_blank()) + 
  guides(fill=guide_legend(title=""))
```

![](README_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

*Fig. 6: Composition of the optimized farm (based on data of Gosling et
al. (2020)), including only socio-economic indicators. Each land-cover
option is shown in an allocated share (%).*

The first example considers socio-economic indicators only (Fig. 6; see
also Fig. 5 of Gosling et al. (2020)). The result corresponds to the
above shown multi-functional portfolio (Fig. 2). This is expected, as
all indicators relevant for the solution of the multi-functional
portfolio (financial stability, investment costs and meeting households
needs) are captured also in the socio-economic bundle.

``` r
performance_socioeconomic <- calcPerformance(result_socioeconomic)

performance_socioeconomic$scenarioTable$performance <-
  performance_socioeconomic$scenarioTable$performance * 100

performance_socioeconomic$beta
```

    ## [1] 0.6132

``` r
ggplot(performance_socioeconomic$scenarioTable,
       aes(x = indicator,
           y = performance,
           color = indicator)) +
  geom_point() +
  geom_hline(yintercept =
               min(performance$scenarioTable$performance),
             linetype = "dashed", color = "red") +
  guides(color = guide_legend(title = "",
                              nrow = 10)) +
  theme_classic() + 
  theme(text = element_text(size = 18),
        legend.position="right",
        axis.ticks.x = element_blank()) +
  scale_x_discrete(labels = seq(1, 10)) +
  labs(y = "Min-max normalized indicator value (%)",
       x = "Indicators") + 
  scale_y_continuous(breaks = seq(0, 101, 10), 
                     limits = c(0, 101)) +
  geom_hline(aes(yintercept=100), size = 1) + 
  annotate(geom = "Text", x = 6, y = 100, label = "Maximum achievable indicator level",
           vjust = -1)
```

<img src="README_files/figure-gfm/unnamed-chunk-15-1.png" style="display: block; margin: auto;" />

*Fig. 7: The performance of each of the socio-economic indicators. The
colored points are the achieved levels of the indicators of all
scenarios s. The dotted, horizontal red line illustrates the guaranteed
performance*
![(1-\beta)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%281-%5Cbeta%29 "(1-\beta)")*,
thus the robust feasible solution of the program (Equation 1 in (Husmann
et al, n.d.)).*

An analysis of the performance of the socio-economic indicators shows
that the performances of the 3 relevant indicators equal the
multi-functional portfolio (Fig. 7). The result is still defined by
financial stability, investment costs and meeting household needs.
Consequently, the guaranteed performance
![(1-\beta)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%281-%5Cbeta%29 "(1-\beta)")
also equals the multi-functional portfolio. Therefore, also this
socio-economic portfolio still does not perfectly reflect the currently
observed land-cover composition. This means that further indicators
appear to be relevant for the actual farmers’ decisions.
<!-- (*CP: Or altnertave weighings are likely...*). KH: Weighting Fass würde ich ungern öffnen, da diese Funktion nicht nicht implementiert ist.-->

``` r
dat_ecologic <- dat[dat$indicator %in% c("Protecting soil resources",
                                         "Protecting water supply"),]

init_ecologic <- initScenario(dat_ecologic,
                              uValue = 2,
                              optimisticRule = "expectation", 
                              fixDistance = NA) 

result_ecologic <- solveScenario(x = init_ecologic)

result_ecologic$landUse %>% gather(key = landCoverOption, value = landCoverShare, 1 : 6) %>% 
  mutate(portfolio = "Ecologic",
         landCoverShare = landCoverShare * 100) %>% 
  ggplot(aes(y = landCoverShare, x = portfolio, fill = landCoverOption)) + 
  geom_bar(position = "stack", stat = "identity") + 
  theme_classic() +
  theme(text = element_text(size = 14)) +
  scale_fill_startrek() +
  labs(y = "Allocated share (%)") +
  scale_y_continuous(breaks = seq(0, 100, 10), 
                     limits = c(0, 100)) +
  theme(axis.title.x=element_blank(),
        axis.ticks.x=element_blank()) + 
  guides(fill=guide_legend(title=""))
```

![](README_files/figure-gfm/unnamed-chunk-16-1.png)<!-- -->

*Fig. 8: Composition of the optimized farm (based on data of Gosling et
al. (2020)), including only ecologic indicators. Each land-cover option
is shown in an allocated share (%).*

The ecological indicator group, as the second example, leads to a
land-cover portfolio comprising forests only (Fig. 8; corresponds to
Fig. 5 of Gosling et al. (2020)). It can be concluded that all
contributions of all other land-cover alternatives in all scenarios
(even the optimistic ones) to the ecological indicators are lower than
these of forests. The land-cover composition of the ecologic bundle
differs fundamentally from the currently observed portfolio. Also the
ecological indicators are therefore apparently not sufficient to
approximate the farmers’ current perceptions.

### Immediate Economic Success

``` r
dat_short <- dat[dat$indicator %in% c("Meeting household needs",
                                      "Liquidity"),]

init_short<- initScenario(dat_short,
                          uValue = 2,
                          optimisticRule = "expectation", 
                          fixDistance = NA) 

result_short <- solveScenario(x = init_short)

result_short$landUse %>% gather(key = landCoverOption, value = landCoverShare, 1 : 6) %>% 
  mutate(portfolio = "Immediate",
         landCoverShare = landCoverShare * 100) %>% 
  ggplot(aes(y = landCoverShare, x = portfolio, fill = landCoverOption)) + 
  geom_bar(position = "stack", stat = "identity") + 
  theme_classic() +
  theme(text = element_text(size = 14)) +
  scale_fill_startrek() +
  labs(y = "Allocated share (%)") +
  scale_y_continuous(breaks = seq(0, 100, 10), 
                     limits = c(0, 100)) +
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank()) + 
  guides(fill = guide_legend(title = ""))
```

![](README_files/figure-gfm/unnamed-chunk-17-1.png)<!-- -->

*Fig. 9: Composition of the optimized farm (based on data of Gosling et
al. (2020)), including the prospective relevant indicators of the
farmers only. Each land-cover option is shown in an allocated share
(%).*

The third example is composed of a bundle of indicators that
prospectively reflect the farmers’ needs and perceptions (Fig. 9).
Corresponding to Fig. 5 of Gosling et al. (2020), this scenario only
consists of indicators that approximate the immediate economic success.
Indeed, the land-cover composition of this portfolio reflects the
portfolio observed in Eastern Panama best. Hence, these indicators best
reflect the farmers’ goals and perceptions in Eastern Panama. The
difference between this portfolio and the desired multi-functional
portfolio (Fig. 2) highlights the requirements a land-cover alternative
must fulfill to meet the farmers’ needs and goals. The silvopasture as
defined in Gosling et al. (2020) may not serve the requirements of the
farmers sufficiently. Since farmers rate liquidity and meeting household
needs higher than long-term profit and economic stability, pasture
outperforms silvopasture in the realistic scenario. Policies or
development plans may consider these indicators as key elements when
promoting a landscape development towards multifunctional landscapes.

### The use of fixDistance

It can be advantageous to define distinct uncertainty levels for the
calculation of the distances to the maximum achievable level (the
reference)
![d\_{iu}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;d_%7Biu%7D "d_{iu}")
(equation 3 in Husmann et al. (n.d.)) and the actual distances under a
certain land-cover composition
![R\_{iu}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;R_%7Biu%7D "R_{iu}")
(equation 5 in Husmann et al. (n.d.), see also equation 9). *CP: Mit dem
folgenden Satz tue ich mich sehr scchwer… Vielleicht sowas wie: For
instance by using an uncertainty level of 3 the uncertainty space
becomes larger, while the uncertainty space reamins constant…. and thus
leads to… VVG: Hier im Satz von Carola stimmt was nicht oder? der
uncertainty space kann ja nicht “lager” und “remains constant”
gleichzeitig abdecken* The broadened uncertainty space allows for a
broader state space of the distances
![d\_{ij}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;d_%7Bij%7D "d_{ij}"),
which usually leads to smoother transitions (*CP: between what?*) when
the uncertainty values are raised. The disadvantage is that the
distances cannot be interpreted straightforwardly any more. The
different uncertainty spaces considered in the denominator and in the
counter of equation 3 (Husmann et al., n.d.) to (*CP:do? An dem Satz ist
was komisch*) not normalize the distance to {0, 1} any more.

``` r
#### uValue 3 ####
path <- exampleData("exampleGosling.xlsx")
dat <- read_excel(path)

applyDf <- data.frame(u = seq(0, 3, .5))

applyFun <- function(x) {
  init <- initScenario(dat, uValue = x, optimisticRule = "expectation", fixDistance = 3)
  result <- solveScenario(x = init)
  return(c(result$beta, as.matrix(result$landUse)))
}

applyDf <- cbind(applyDf,
                 t(apply(applyDf, 1, applyFun)))

names(applyDf) <- c("u", "beta", names(result$landUse))

applyDf[, c(3 : 8)] <- applyDf[, c(3 : 8)] * 100

applyDf %>% gather(key = "land-cover option", value = "land-cover share", -u, -beta) %>%
  ggplot(aes(y = `land-cover share`, x = u, fill = `land-cover option`)) + 
  geom_area(alpha = .8, color = "white") + theme_minimal()+
  labs(x = "Uncertainty level", y = "Allocated share (%)") + 
  guides(fill=guide_legend(title="")) + 
  scale_y_continuous(breaks = seq(0, 100, 10), 
                     limits = c(0, 100.01)) +
  scale_x_continuous(breaks = seq(0, 3, 0.5),
                     limits = c(0, 3)) + 
  scale_fill_startrek() +
  theme(text = element_text(size = 18),
        legend.position = "bottom")
```

![](README_files/figure-gfm/unnamed-chunk-18-1.png)<!-- -->

*Fig. 10: Theoretically ideal farm compositions using the fixDistance
argument and increasing levels of uncertainty.*

*CP: Hier ist auch was komisch* It can be seen that the land-cover
allocation transition under raising uncertainty (Fig. 10) differs
slightly from the multifunctional scenario shown above (Fig. 2). The
here broadened state space leads an earlier increase of pasture (*CP:
Ist hier gemeint: leads to selection of a larger share of pasture under
low uncertainty levels as compared to…*).

<h3>
<a name="7. Suggested">5. Suggested citation </a>
</h3>

Husmann, K., von Groß, V., Fuchs J.M., Bödeker, K. (2022): optimLanduse:
Robust Land-Use Optimization. R package version 1.1.0.
<https://CRAN.R-project.org/package=optimLanduse>.

<h3>
<a name="8. Literatur">6. Literature</a>
</h3>

Gosling, E., Reith, E., Knoke T., Paul, C. (2020): A goal programming
approach to evaluate agroforestry systems in Eastern Panama<em>Journal
of Environmental Management</em> **261**.
<https://doi.org/10.1016/j.jenvman.2020.110248>

Husmann, K.; von Groß, V.; Bödeker, K.; Fuchs, J.; Paul, C.; Knoke, T.
(no date): *optimLanduse*: A Package for Multiobjective Land-cover
Composition Optimization under Uncertainty. In: *Methods Ecol Evol.*
Under revision.

Knoke, T., Paul, C., Hildebrandt, P. et al. (2016): Compositional
diversity of rehabilitated tropical lands supports multiple ecosystem
services and buffers uncertainties. <em>Nat Commun</EM> **7**, 11877.
<https://doi.org/10.1038/ncomms11877>

Paul, C., Weber, M., Knoke, T. (2017): Agroforestry versus farm mosaic
systems – Comparing land-use efficiency, economic returns and risks
under climate change effects. <em>Sci. Total Environ.</em> **587-588**.
<https://doi.org/10.1016/j.scitotenv.2017.02.037>.

Knoke, T., Paul, C., et al. (2020). Accounting for multiple ecosystem
services in a simulation of land‐use decisions: Does it reduce tropical
deforestation?. <em>Global change biology</em> **26(4)**.
<https://doi.org/10.1111/gcb.15003>
