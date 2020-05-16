import os, sys

import pandas as pd
from pandas.util.testing import assert_frame_equal

import unittest
from functional_space.lib.fspace import Func, Functional, funcs_1d
import numpy as np


class BlackKarasinskiTwoFactorModelTest(unittest.TestCase):
    def setUp(self):
        core_path = os.path.dirname(__file__)
        sub_path = os.path.join(core_path, 'data')
        self.path = os.path.abspath(sub_path)
        self.r_bar = 0.02
        self.alpha_x = 0.15
        self.alpha_y = 0.02
        self.sigma_x = 0.20
        self.sigma_y = 0.10
        self.rho = -0.30
        self.s = 0.0
        self.t0 = 1.0
        self.Time = np.linspace(start=1, stop=5, num=9)
        self.x0 = 0.4
        self.X = np.linspace(start=-1, stop=1., num=11)
        self.y0 = 0.5
        self.Y = np.linspace(start=-1, stop=1., num=11)
        self.pth = lambda x: os.path.join(self.path, x + '.txt')

    @staticmethod
    def assert_mathcad_python_equal(benchmark_file_path, result_list, tol=8):
        with open(benchmark_file_path, 'r') as f:
            r = f.read()
        bench_df = pd.DataFrame([float(k) for k in r.split('\n')])
        bench_df = bench_df.round(tol)
        result_df = pd.DataFrame(result_list)
        df = pd.concat([bench_df, result_df], axis=1)
        df.columns = ['benchmark', 'result']
        assert_frame_equal(bench_df.round(tol), result_df.round(tol))

    def test_delta(self):
        self.delta = Func(lambda t, v: Functional(self.r_bar, x0=t, x1=v, du=True)())
        self.assert_mathcad_python_equal(self.pth('2delta_t_T'), self.delta(self.Time, t=self.t0))

    def test_D(self):
        self.test_delta()
        self.D = Func(lambda t, v: funcs_1d.exp(-self.delta(t, v=v)))
        self.assert_mathcad_python_equal(self.pth('2D0'), self.D(self.Time, t=self.t0))

    def test_phi_x(self):
        self.phi_x = Func(lambda t, v: funcs_1d.exp(-self.alpha_x * (v - t)))
        self.assert_mathcad_python_equal(self.pth('2_phi_x'), self.phi_x(self.Time, t=self.t0))

    def test_phi_y(self):
        self.phi_y = Func(lambda t, v: funcs_1d.exp(-self.alpha_y * (v - t)))
        self.assert_mathcad_python_equal(self.pth('2_phi_y'), self.phi_y(self.Time, t=self.t0))

    def test_sigma_xx(self):
        self.sigma_xx = Func(
            lambda t, v: (pow(self.sigma_x, 2)) / (2 * self.alpha_x) * (1 - funcs_1d.exp(-2 * self.alpha_x * (v - t))))
        self.assert_mathcad_python_equal(self.pth('2_sigma_xx'), self.sigma_xx(self.Time, t=self.t0))

    def test_sigma_xx(self):
        self.sigma_xx = Func(
            lambda t, v: (pow(self.sigma_x, 2)) / (2 * self.alpha_x) * (1 - funcs_1d.exp(-2 * self.alpha_x * (v - t))))
        self.assert_mathcad_python_equal(self.pth('2_sigma_xx'), self.sigma_xx(self.Time, t=self.t0))

    def test_sigma_xy(self):
        self.sigma_xy = Func(lambda t, v:
                             ((self.rho * self.sigma_x * self.sigma_y) / (self.alpha_x + self.alpha_y)) * (
                                     1 - funcs_1d.exp(-(self.alpha_x + self.alpha_y) * (v - t)))
                             )
        self.assert_mathcad_python_equal(self.pth('2_sigma_xy'), self.sigma_xy(self.Time, t=self.t0))

    def test_sigma_yy(self):
        self.sigma_yy = Func(lambda t, v: (pow(self.sigma_y, 2)) / (2 * self.alpha_y) * (
                1 - funcs_1d.exp(-2 * self.alpha_y * (v - t))))
        self.assert_mathcad_python_equal(self.pth('2_sigma_yy'), self.sigma_yy(self.Time, t=self.t0))

    def test_sigma_star(self):
        self.test_phi_x()
        self.test_phi_y()
        self.test_sigma_xx()
        self.test_sigma_yy()
        self.test_sigma_xy()
        self.sigma_star = Func(lambda t, u:
                               pow(self.phi_x(t, v=u), 2) * self.sigma_xx(0, v=t) +
                               2 * self.phi_x(t, v=u) * self.phi_y(t, v=u) * self.sigma_xy(0, v=t) +
                               pow(self.phi_y(t, v=u), 2) * self.sigma_yy(0, v=t)
                               )
        self.assert_mathcad_python_equal(self.pth('2_sigma_star'), self.sigma_star(self.Time, t=self.t0))

    def test_r1_tilde(self):
        self.r1_tilde = Func(lambda t: self.r_bar - self.s)
        self.assert_mathcad_python_equal(self.pth('2_r1_tilde'), self.r1_tilde(self.Time))

    def test_R1(self):
        self.test_r1_tilde()
        self.test_phi_x()
        self.test_phi_y()
        self.test_sigma_star()
        self.R1 = Func(lambda x, y, t, t1:
                       self.r1_tilde(t1) *
                       funcs_1d.exp(self.phi_x(t, v=t1) * x +
                                    self.phi_y(t, v=t1) * y -
                                    0.5 * self.sigma_star(t, u=t1))
                       )
        self.assert_mathcad_python_equal(self.pth('2_R1'), self.R1(self.Time, t=self.t0, x=self.x0, y=self.y0))

    def test_F_xy(self):
        self.test_phi_x()
        self.test_phi_y()
        self.test_sigma_xx()
        self.test_sigma_xy()
        self.test_sigma_yy()
        self.test_sigma_yy()
        self.test_sigma_star()
        self.F_xy = Func(lambda t, t1, t2:
                         self.phi_x(t1, v=t2) * self.sigma_xx(t, v=t1) +
                         (self.phi_x(t1, v=t2) + self.phi_y(t1, v=t2)) * self.sigma_xy(t, v=t1)
                         + self.phi_y(t1, v=t2) * self.sigma_yy(t, v=t1))

        self.assert_mathcad_python_equal(self.pth('2_Fxy'), self.F_xy(self.Time, t=self.t0, t2=1.5))

    def test_r2_tilde(self):
        self.test_r1_tilde()
        self.test_F_xy()

        r2_integrand = Func(lambda t, t1:
                            self.r1_tilde(t1) * (
                                    funcs_1d.exp(self.F_xy(0, t1=t1, t2=t)) - 1))
        self.r2_tilde = Func(lambda t: self.r1_tilde(t) * Functional(r2_integrand(t=t), x0=0, x1=t, du=True)())
        self.assert_mathcad_python_equal(self.pth('2_r2_tilde'), [1e3 * a for a in self.r2_tilde(self.Time)])

    def test_R2(self):
        self.test_r2_tilde()
        self.test_sigma_star()
        self.test_phi_y()
        self.test_phi_x()
        self.R2 = Func(lambda x, y, t, t1:
                       self.r2_tilde(t1) *
                       funcs_1d.exp(self.phi_x(t, v=t1) * x +
                                    self.phi_y(t, v=t1) * y -
                                    0.5 * self.sigma_star(t, u=t1)))
        self.assert_mathcad_python_equal(self.pth('2_R2'), self.R2(self.Time, t=self.t0, x=self.x0, y=self.y0))

    def test_f1(self):
        self.test_R1()
        self.test_r1_tilde()

        self.f1 = Func(lambda x, y, t, T: Functional(
            Func(lambda x, y, t, t1: self.R1(t1, x=x, y=y, t=t) - self.r1_tilde(t1))(x=x, y=y, t=t),
            x0=t, x1=T,
            du=True)())
        self.assert_mathcad_python_equal(self.pth('2_f1'), self.f1(self.Time, x=self.x0, y=self.y0, t=self.t0))

    def test_f2(self):
        self.test_R1()
        self.test_r1_tilde()
        self.f2 = Func(lambda x, y, t, T:
                       0.5 * pow(Functional(
                           Func(lambda x, y, t, t1:
                                self.R1(t1, x=x, y=y, t=t) - self.r1_tilde(t1))
                           (x=x, y=y, t=t),
                           x0=t, x1=T, du=True)(), 2))
        self.assert_mathcad_python_equal(self.pth('2_f2'),
                                         [1e3 * a for a in self.f2(self.Time, x=self.x0, y=self.y0, t=self.t0)])

    def test_f3_prime(self):
        self.test_R1()
        self.test_F_xy()
        self.f3_prime = Func(lambda x, y, t, t2:
                             Functional(
                                 Func(
                                     lambda x, y, t, t1, t2:
                                     self.R1(t1, x=x, y=y, t=t) * (
                                             funcs_1d.exp(self.F_xy(t1, t=t, t2=t2)) - 1)
                                 )
                                 (x=x, y=y, t=t, t2=t2),
                                 x0=t, x1=t2, du=True)())

        self.assert_mathcad_python_equal(self.pth('2_f3_prime'),
                                         [a * 1e3 for a in self.f3_prime(x=self.x0, y=self.y0, t=self.t0)(self.Time)])

    def test_f3(self):
        self.test_R1()
        self.test_R2()
        self.test_f3_prime()
        f3_integrand = Func(lambda x, y, t, t2:
                            self.R1(t2, x=x, y=y, t=t) * self.f3_prime(t2, x=x, y=y, t=t) -
                            self.R2(t2, x=x, y=y, t=t))
        self.f3 = Func(lambda x, y, t, T:
                       Functional(
                           f3_integrand(x=x, y=y, t=t),
                           x0=t, x1=T, du=True)())
        self.assert_mathcad_python_equal(self.pth('2_f3'),
                                         [1e6 * a for a in self.f3(x=self.x0, y=self.y0, t=self.t0)(self.Time)])

    def test_F1(self):
        self.test_f1()
        self.test_D()
        self.F1 = Func(lambda x, y, t, T: self.D(t, v=T) * (1 - self.f1(t, x=x, y=y, T=T)))
        self.assert_mathcad_python_equal(self.pth('2_FCapital1'), self.F1(self.Time, x=self.x0, y=self.y0, t=self.t0))

    def test_F2(self):
        self.test_D()
        self.test_f1()
        self.test_f2()
        self.test_f3()
        self.F2 = Func(lambda x, y, t, T:
                       self.D(t, v=T) * (1
                                         - self.f1(t, x=x, y=y, T=T)
                                         + self.f2(t, x=x, y=y, T=T)
                                         + self.f3(t, x=x, y=y, T=T)))

        self.assert_mathcad_python_equal(self.pth('2_FCapital2'), self.F2(self.Time, x=self.x0, y=self.y0, t=self.t0))


if __name__ == '__main__':
    unittest.main()
