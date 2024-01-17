---
title: "Analysis of duration of response in oncology studies"
author: Lu Mao, Lu Tian, Bo Huang, Lee-Jen Wei
execute: 
  eval: false
bibliography: references.bib
---


This R-code implements the analysis of duration-of-response data described in @mao:2024 for oncology studies .

### Background

Duration-of-response (DOR) is often utilized to assess efficacy of oncology therapies. While valid procedures are available for estimating the mean [@huang2018], it is unclear how to estimate the distribution of DOR, which provides a fuller picture of patient experience. Using simple Kaplan–Meier estimates based on the responders can result in distorted conclusions due to the existence of non-responders (whose DOR is zero).

### Notation and Methodology

For a responder, let $T_1$ denote the time to response. Let $T_2$ denote the time to disease progression or death. If a patient experiences progression/death first, then $T_1=\infty$.
The DOR can be defined as $[T_2-T_1]_+$, where $x_+=xI(x>0).$ Thus, DOR is $T_2-T_1$ for responders, i.e., those with $T_1<\infty$, and 0 for nonresponders.
For notational convenience, set $T_1=T_2$ for nonresponders. 

### Usage and Example

#### R-program

To use the program, download `DORfunctions.R` from <https://github.com/lmaowisc/dor>. The file contains all functions needed to perform the analysis.


::: {.cell}

```{.r .cell-code}
source("DORfunctions.R")
```
:::


The main function is `dorfit(x1, delta1, x2, delta2, tau, med_inf = FALSE)`.

**Input**

-   `x1` time to earliest of response, outcome event, or censoring
-   `delta1` indicator of response
-   `x2` time to earlier of outcome event and censoring 
-   `delta2` indicator of outcome event
-   `tau` restriction time
-   `med_inf` whether to make inference on the median DOR

**output**

-   `x1` time to earliest of response, outcome event, or censoring
-   `delta1` indicator of response
-   `x2` time to earlier of outcome event and censoring 


#### A simulated example

ss
