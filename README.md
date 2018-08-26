# Statistical-Modelling-Project

# Simulation study about post-selection inference

We set up a simulation study to show the effect of model selection on the construction of a confidence interval. More precisely, we have a candidate set of models and we use an information criterion (AIC, BIC) to select a model. Once the model is selected, we use it to construct a 95% confidence interval for the model parameters. Via simulation we investigate the empirical coverage probabilities. Generate a large number of datasets from a known model. If there were no model selection, with each dataset one would construct a confidence interval and afterwards, count how many of those intervals contain the true parameter value. With model selection this procedures is different since each generated dataset might lead to a different selected model. For this reason, we focus on a single model M and generate datasets such that in all the datasets the selected model is M. This is called a conditional data generating process. 

We generate data from a linear regression model of the form Y = X+" where Y is a vector with length n = 200, X is a 200 * 5 matrix. The first column of X is for the intercept, so it is a vector of ones. The other four columns of X are generated from a multivariate normal distribution with mean vector 0 and a covariance matrix with diagonal elements equal to 1 and off-diagonal elements equal to ρ. 


You
have to answer the following questions for two values for : 0.25 and 0.75. The parameter
vector  = (1; 2; 3; 0; 0)t is dened by you. The last two elements of  should be set to zero and the
rst three elements can be any number in the interval [2; 6] or [􀀀2;􀀀6]. You choose the values and keep
them xed for this question. The error vector " with length 200 contains i.i.d. elements from a standard
normal distribution. The model M that you focus on contains the parameters  = (1; 2; 3; 4; 0)t.
