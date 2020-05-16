# Functional Space
The code in this repository has been developed in support to the assymptotic formulae presented in Chapter 15 of the book **Perturbation Methods in Credit Derivatives** by Dr. Colin Turfus.

The class `Func` defines a general function together with all the algebraic operations.


The function `Functional` converts the function into the [functional](https://en.wikipedia.org/wiki/Functional_(mathematics)). We use the `quadrature` from `scipy.integrate` package for the integration.

Please see this [jupyter notebook](https://nbviewer.jupyter.org/github/ashubertt/functional_space/blob/uat/jupyter/Two-Factor%20Black-Karasinski%20Pricing%20Kernel.ipynb) for the application to **Black-Karasinski Rates Asymptotic Expansion** formulae.

For the test purposes we've created a benchmark MathCad implementation. The results of such an implementation are stored in `unittest/data` folder. To run tests execute:
`python unittests/test_bk2factor_model.py -v`.

