# Full text

The full dissertation is the file master.pdf.

# Abstract

The chapters of this dissertation are intended to be three independent, publishable papers, but they nevertheless share the theme of predictive inferences for explanatory item models.
Chapter 1 describes the differences between the Bayesian and frequentist statistical frameworks in the context of explanatory item response models. The particular model of focus, the "doubly explanatory model", is a model for dichotomous item responses that includes covariates for person ability and covariates for item difficulty. It includes many Rasch-family models as special cases. Differences in how the model is understood and specified within the two frameworks are discussed. The various predictive inferences available from the model are defined for the two frameworks. 

Chapter 2 is situated in the frequentist framework and focuses on approaches for explaining or predicting the difficulties of items. Within the frequentist framework, the linear logistic test model (LLTM) is likely to be used for this purpose, which in essence regresses item difficulty on covariates for characteristics of the items. However, this regression does not include an error term, and so the model is in general misspecified. Meanwhile, adding an error term to the LLTM makes maximum likelihood estimation infeasible. To address this problem, a two-stage modeling strategy (LLTM-E2S) is proposed: in the first stage Rasch model maximum likelihood estimates for item difficulties and standard errors are obtained, and in the second stage a random effects meta-analysis regression of the Rasch difficulties on covariates is performed that incorporates the uncertainty in the item difficulty estimates.
In addition, holdout validation, cross-validation, and Akaike information criteria (AIC) are discussed as means of comparing models that have different sets of item predictors.
I argue that AIC used with the LLTM estimates the expected deviance of the fitted model when applied to new observations from the *same* sample of items and persons, which is unsuitable for assessing the ability of the model to predict item difficulties.
On the other hand, AIC applied to the LLTM-E2S provides the expected deviance related to new observations arising from *new* items, which is what is needed.
A simulation study compares parameter recovery and model comparison results for the two modeling strategies.

Chapter 3 takes a Bayesian outlook and focuses on models that explain or predict person abilities. I argue that the usual application of Bayesian forms of information criteria to these models yields the wrong inference.
Specifically, when using likelihoods that are conditional on person ability, information criteria estimate the expected fit of the model to new data arising from the *same* persons.
What are needed are likelihoods that are marginal over the distribution for ability, which may be used with information criteria to estimate the expected fit to new data from a *new* sample of persons. 
The widely applicable information criterion (WAIC), Pareto-smoothed importance sampling approximation to leave-one-out cross-validation, and deviance information criterion (DIC) are discussed in the context of these conditional and marginal likelihoods. 
An adaptive quadrature scheme for use within Markov chain Monte Carlo estimation is proposed to obtain the marginal likelihoods.
Also, the moving block bootstrap is investigated as a means to estimate the Monte Carlo error for Bayesian information criteria estimates. 
A simulation study using a linear random intercept model is conducted to assess the accuracy of the adaptive quadrature scheme and the bootstrap estimates of Monte Carlo error.
These methods are then applied to an real item response dataset, demonstrating the practical difference between conditional and marginal forms of information criteria.

# Compilation

The dissertation was compiled by executing these files in the order shown below:

* chapter_2/simulation part 1.do
* chapter_2/simulation part 1.R
* chapter_2/chapter_2.Rnw
* chapter_3/simulation/simulation.R
* chapter_3/simulation/bootstrap.R
* chapter_3/application/application.R
* chapter_3/chapter_3.Rnw
* appendix/appendix.Rnw
* compile_chapters.R
* master.tex
