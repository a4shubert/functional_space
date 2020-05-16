# Functional Space
The code in this repository has been developed in support to the assymptotic formulae presented in Chapter 15 of the book **Perturbation Methods in Credit Derivatives** by Dr. Colin Turfus.

The class `Func` defines a general function together with all the algebraic operations.


The funciton `Functional` converts the function into the [functional](https://en.wikipedia.org/wiki/Functional_(mathematics)). We use the `quadrature` from `scipy.integrate` package for the integration.

Please see this jupyter notebook for the application to **Black-Karasinski  Rates Asymptotic Expansion** formulae.

For the test purposes we've created a benchmark MathCad implementation. The results of such an implementation are stored in `unittest/data` folder. To run tests execute:
`python unittests/test_bk2factor_model.py -v`. You should see:

`
test_D (__main__.BlackKarasinskiTwoFactorModelTest) ... ok
test_F1 (__main__.BlackKarasinskiTwoFactorModelTest) ... ok
test_F2 (__main__.BlackKarasinskiTwoFactorModelTest) ... ok
test_F_xy (__main__.BlackKarasinskiTwoFactorModelTest) ... ok
test_R1 (__main__.BlackKarasinskiTwoFactorModelTest) ... ok
test_R2 (__main__.BlackKarasinskiTwoFactorModelTest) ... ok
test_delta (__main__.BlackKarasinskiTwoFactorModelTest) ... ok
test_f1 (__main__.BlackKarasinskiTwoFactorModelTest) ... ok
test_f2 (__main__.BlackKarasinskiTwoFactorModelTest) ... ok
test_f3 (__main__.BlackKarasinskiTwoFactorModelTest) ... ok
test_f3_prime (__main__.BlackKarasinskiTwoFactorModelTest) ... ok
test_phi_x (__main__.BlackKarasinskiTwoFactorModelTest) ... ok
test_phi_y (__main__.BlackKarasinskiTwoFactorModelTest) ... ok
test_r1_tilde (__main__.BlackKarasinskiTwoFactorModelTest) ... ok
test_r2_tilde (__main__.BlackKarasinskiTwoFactorModelTest) ... ok
test_sigma_star (__main__.BlackKarasinskiTwoFactorModelTest) ... ok
test_sigma_xx (__main__.BlackKarasinskiTwoFactorModelTest) ... ok
test_sigma_xy (__main__.BlackKarasinskiTwoFactorModelTest) ... ok
test_sigma_yy (__main__.BlackKarasinskiTwoFactorModelTest) ... ok

----------------------------------------------------------------------
Ran 19 tests in 3.657s

OK
`



