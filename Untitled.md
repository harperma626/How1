Data Analysis on world's richest people
================
Zhuoqin Ma
2/10/2018

*In this project we look at dataset containing information on the world’s richest people from the World Top Incomes Database (WTID) hosted by the Paris School of Economics \[<http://wid.world>\]. This is derived from income tax reports, and compiles information about the very highest incomes in various countries over time, trying as hard as possible to produce numbers that are comparable across time and space.*

*For most countries in most time periods, the upper end of the income distribution roughly follows a Pareto distribution, with probability density function*

$$
f(x)=\\frac{a−1}{x\_{min}}(\\frac{x}{x\_{min}})^{−a}
$$

*for incomes* *X* ≥ *x*<sub>*m**i**n*</sub>. *(Typically*, *x*<sub>*m**i**n*</sub> *is large enough that only the richest 3%-4% of the population falls above it.) As the **Pareto exponent**, *a*, gets smaller, the distribution of income becomes more unequal, that is, more of the population’s total income is concentrated among the very richest people.*

*The proportion of people whose income is at least* *x*<sub>*m**i**n*</sub> *and whose income is also at or above any level* *w* ≥ *x*<sub>*m**i**n*</sub> *is thus*

$$
Pr(X≥w)=\\int\_{w}^{∞}f(x)dx=\\int\_{w}^{∞}\\frac{a−1}{x\_{min}}(\\frac{x}{x\_{min}})^{−a}dx=(\\frac{w}{x\_{min}})^{−a+1}
$$

*We will use this to estimate how income inequality changed in the US over the last hundred years or so. (Whether the trends are good or bad or a mix is beyond our scope here.) WTID exports its data sets as .xlsx spreadsheets. For this project, we have extracted the relevant data and saved it as* `wtid-report.csv`.

$$
(\\frac{P99}{P99.9})^{-a+1}=10 
$$

Estimating a on US data
-----------------------

We will use the fact that we can estimate the exponent using the following formula, which we refer to as Equation (1):

$$
a=1−\\frac{log10}{log(\\frac{P99}{P99.9})}.
$$

``` r
wtid <- read.csv("wtid-report.csv", as.is = TRUE)
wtid <- wtid[, c("Year", "P99.income.threshold", "P99.5.income.threshold", "P99.9.income.threshold")]
names(wtid) <- c("Year", "P99", "P99.5", "P99.9")

exponent.est_ratio <- function(p99, p999) {
  return(1 - log(10)/(log(p99/p999)))
}

ahat <- exponent.est_ratio(wtid$P99, wtid$P99.9)
```

At the same time, the logic leading to Equation (1) also implies that
$$
(\\frac{P99.5}{P99.9})^{−a+1}=5
$$

Goodness of fit
---------------

#### Step 1

We'll calculates the left-hand side of that equation. Plot the values for each year using `ggplot`, and discover how good is the fit using the data and estimates of the exponent.

``` r
library(ggplot2)
exponent.est_ratio2 <- function(p995, p999, a) {
  return((p995/p999)^(-a+1))
}
ests <- exponent.est_ratio2(wtid$P99.5, wtid$P99.9, ahat)
ggplot(data = wtid) +
  geom_point(mapping = aes(x = Year, y = ests)) +
  geom_abline(intercept = 5, slope = 0, color = "orange") +
  labs(x = "Year", y = "a Estimate")
```

![](Untitled_files/figure-markdown_github/unnamed-chunk-2-1.png)

#### Step 2

By parallel reasoning, we should have (*P*99/*P*99.5)<sup>−*a* + 1</sup> = 2. Repeat the previous step with this formula.

``` r
exponent.est_ratio3 <- function(p99, p995, a) {
  return((p99/p995)^(-a+1))
}
ests2 <- exponent.est_ratio3(wtid$P99, wtid$P99.5, ahat)
ggplot(data = wtid) +
  geom_point(mapping = aes(x = Year, y = ests2)) +
  geom_abline(intercept = 2, slope = 0, color = "orange") +
  labs(x = "Year", y = "a Estimate")
```

![](Untitled_files/figure-markdown_github/unnamed-chunk-3-1.png)

#### Step 3

We have shown that if the upper tail of the income distribution followed a perfect Pareto distribution, then

$$
\\begin{eqnarray}
(\\frac{P99}{P99.9})^{−a+1}=10\\\\
(\\frac{P99.5}{P99.9})^{−a+1}=5\\\\
(\\frac{P99}{P99.5})^{−a+1}=2
\\end{eqnarray}
$$

We could estimate the Pareto exponent by solving any one of these equations for *a*. Because of measurement error and sampling noise, we can’t find one value of *a* which will work for all three equations. Generally, trying to make all three equations come close to balancing gives a better estimate of *a* than just solving one of them. (This is analogous to finding the slope and intercept of a regression line by trying to come close to all the points in a scatterplot, and not just running a line through two of them.)

We will therefore estimate a by minimizing

$$
((\\frac{P99}{P99.9})^{−a+1}−10)^2+((\\frac{P99.5}{P99.9})^{−a+1}−5)^2+((\\frac{P99}{P99.5})^{−a+1}−2)^2.
$$

We'll Write a function, `percentile_ratio_discrepancies`, which takes as inputs `P99`, `P99.5`, `P99.9` and *a* , and returns the value of the expression above.

``` r
percentile_ratio_discrepancies <- function(a, p99, p995, p999) {
  return(((p99/p999)^(-a+1) - 10)^2 + ((p995/p999)^(-a+1) - 5)^2 + ((p99/p995)^(-a+1) - 2)^2)
}
```

We’d like to write a function, `exponent.multi_ratios_est`, which takes as inputs the vectors `P99`, `P99.5`, `P99.9`, and estimates *a*. It should minimize the function `percentile_ratio_discrepancies` wrote above.

There are several ways to achieve that result. We could used gradient descent to minimize the mean squared error (MSE) of a model fit depending on a parameter *β* to data. For the gradient descent algorithm, we'll approximate the derivative of the MSE, and adjust our estimate of *β* by an amount proportional (and opposite) to that approximation. We'll stop the algorithm when the derivative became small (assuming, then, that we will be near a minimum). But for this project we will use a built-in R optimization function to do the minimization.

``` r
exponent.multi_ratios_est <- function(p99, p995, p999) {
  a0  <- 1 - (log(10))/(log(p99/p999))
  min <- nlm(percentile_ratio_discrepancies, a0, p99, p995, p999)$estimate
  return(min)
}
```

Conduct estimation for each year
--------------------------------

#### Step 1

We will follow the logic and write a function which uses `exponent.multi_ratios_est` to estimate a for the US for every year from 1913 to 2012.

``` r
multi_ratios_allyears <- function(country_data) {
  
  # country_data should be a dataframe with the variables "Year", "P99", "P99.5", and "P99.9" where each row corresponds to one year. 
  
  n <- nrow(country_data)
  estimates <- rep(NA, n)
  names(estimates) <- country_data$Year
  for (i in 1:n) {
    if (all(is.na(country_data[i, c("P99", "P99.5", "P99.9")]))) {
      estimates[i] <- NA
    } else {
      estimates[i] <- exponent.multi_ratios_est(country_data$P99[i], country_data$P99.5[i], country_data$P99.9[i])
    }
  }  
  return(estimates)
}

USestimates <- multi_ratios_allyears(wtid)

ggplot(data = wtid) +
  geom_point(mapping = aes(x = Year, y = USestimates)) +
  labs(title = "Pareto Exponent Estimates Over the Years", x = "Year", y = "Estimated Pareto Exponent")
```

![](Untitled_files/figure-markdown_github/unnamed-chunk-6-1.png)

#### Step 2

To use Equation (1) to estimate *a* for the US for every year we will make a scatter-plot of these estimates against those from (v.) using `ggplot`.

``` r
single_ratio  <- 1 - (log(10))/(log(wtid$P99/wtid$P99.9))
names(single_ratio) <- wtid$Year

ggplot(data = wtid) +
  geom_point(mapping = aes(x = USestimates, y = single_ratio)) +
  labs(title = "Pareto Exponent Estimates Two Ways", x = "Using All Ratios", y = "Using One Ratio")
```

![](Untitled_files/figure-markdown_github/unnamed-chunk-7-1.png)

Estimate a for other countries
------------------------------

We’re now going to look at this same data for some other countries: Canada, China, Colombia, India, Italy, Japan, and Sweden. This data is in the file `wtid-report.csv`. The WTID website also has data on the average income per `tax unit` (roughly, household) for the US and the other countries. This info is in stored in theAverageIncome column.

Use function from (v) to estimate *a* over time for each of them. Note that the size of the dataset is different for each of these countries, and there may be some NA values.

``` r
wtid2 <- read.csv("wtid-homework.csv", as.is = TRUE)
country_num <- length(unique(wtid2$Country))

for(i in 1:country_num) {
  these_rows <- wtid2$Country == unique(wtid2$Country)[i]
  this_data  <- wtid2[these_rows, c("Year", "P99", "P99.5", "P99.9")]
  estimates  <- multi_ratios_allyears(this_data)
  wtid2$Estimate[these_rows] <- estimates
}
```

We will plot estimates of *a* over time for all the countries and the series of average income per “tax unit” for the US and the countries against time using `ggplot`. Note that the years covered by the data are different for each country.

``` r
ggplot(data = wtid2) +
  geom_point(mapping = aes(x = Year, y = Estimate, color = Country)) +
  labs(title = "Pareto Exponent Estimate", y = "Estimate", x = "Year")
```

    ## Warning: Removed 104 rows containing missing values (geom_point).

![](Untitled_files/figure-markdown_github/unnamed-chunk-9-1.png)

``` r
ggplot(data = wtid2) +
  geom_point(mapping = aes(x = Year, y = AverageIncome, color = Country)) +
  labs(title = "Average Income Over Time", y = "Average Income", x = "Year")
```

    ## Warning: Removed 10 rows containing missing values (geom_point).

![](Untitled_files/figure-markdown_github/unnamed-chunk-9-2.png)

Testing hypothesis on US data
-----------------------------

The most influential hypothesis about how inequality is linked to economic growth is the “U-curve” hypothesis proposed by the great economist Simon Kuznets in the 1950s. According tho this idea, inequality rises during the early, industrializing phases of economic growth, but then declines as growth continues.

We will make a scatter-plot of estimated exponents for the US against the average income for the US in ggplot.

``` r
ggplot(data = wtid2[wtid2$Country == "United States", ]) +
  geom_point(mapping = aes(x = AverageIncome, y = Estimate))
```

![](Untitled_files/figure-markdown_github/unnamed-chunk-10-1.png)

However, our scatterplot doesn’t seem to support this hypothosis. Smaller *a* coefficients mean more inequality, so our scatterplot shows inequality high, then decreasing for a while but now rising again – the opposite of the “U-curve” hypothesis.

If we want a more quantitative check on the Kuznets hypothesis, we could use `lm()` to regress our estimated exponents on the average income for the US, including a quadratic term for income.

``` r
lm0 <- lm(Estimate ~ AverageIncome + I(AverageIncome^2), data = wtid2[wtid2$Country == "United States", ])
summary(lm0)
```

    ## 
    ## Call:
    ## lm(formula = Estimate ~ AverageIncome + I(AverageIncome^2), data = wtid2[wtid2$Country == 
    ##     "United States", ])
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.50724 -0.18364 -0.02531  0.18689  0.54918 
    ## 
    ## Coefficients:
    ##                      Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)         8.230e-01  1.515e-01   5.432 3.93e-07 ***
    ## AverageIncome       1.394e-04  1.015e-05  13.740  < 2e-16 ***
    ## I(AverageIncome^2) -1.891e-09  1.451e-10 -13.027  < 2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.2466 on 100 degrees of freedom
    ## Multiple R-squared:  0.6679, Adjusted R-squared:  0.6612 
    ## F-statistic: 100.5 on 2 and 100 DF,  p-value: < 2.2e-16

``` r
ggplot(data = wtid2[wtid2$Country == "United States", ]) +
  geom_point(mapping = aes(x = AverageIncome, y = fitted(lm0))) +
  labs(main = "United States", x = "Average Income", y = "Fitted Values")
```

![](Untitled_files/figure-markdown_github/unnamed-chunk-12-1.png)

From the plot we conclude that our coefficients are not consistent with the hypothesis. Since the coefficient in front of the quadratic term is negative, the model is a parabola opening downwards, not upwards, like a ‘U’ would.

Does the hypothesis hold true for other countries?
--------------------------------------------------

Next we will do a separate quadratic regression for each country to find out which ones have estimates compatible with the hypothesis.

``` r
reg_fit <- function(data, country) {
  # data is a data frame with the variables Country, Estimate, and Average Income. This function returns the linear regression coefficients and plots the fitted values against the average income.
  
  if (all(is.na(data$Estimate[data$Country == country]))) {
    return(rep(NA, 3))
  } else {  
    lm0 <- lm(Estimate ~ AverageIncome + I(AverageIncome^2), data = data[data$Country == country, ])
    summary(lm0)
    return(lm0$coef)
  }  
}  

country_num <- length(unique(wtid2$Country)) 

coefs <- matrix(NA, ncol = 3, nrow = country_num)
colnames(coefs) <- c("Intercept", "AverageIncome", "AverageIncome2")
rownames(coefs) <- unique(wtid2$Country)

for (i in 1:country_num) {
  coefs[i, ] <- reg_fit(wtid2, unique(wtid2$Country)[i])
}  
coefs
```

    ##                Intercept AverageIncome AverageIncome2
    ## Canada         2.2660536  1.240966e-04  -3.360837e-09
    ## China         10.3978092 -1.126763e-03   5.257536e-08
    ## Colombia      34.6124009 -6.095234e-03   2.867133e-07
    ## Italy          2.5824163  1.594300e-04  -6.591048e-09
    ## Japan          3.7291069 -5.136191e-04   1.889447e-07
    ## Sweden         1.0123533  4.414454e-05  -1.496762e-10
    ## United States  0.8230049  1.394435e-04  -1.890556e-09

Conclusion
----------

From the coefficients output, the relationship between the coefficient estimates and the average income is modeled by a “U” shape for China, Colombia, and Japan, but not for the others.

(If we were doing a more rigorous check of the Kuznet hypothesis, we would want to control for other factors, and not just assume that a quadratic was the right functional form for the curve.)
