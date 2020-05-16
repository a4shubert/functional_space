from collections import namedtuple
import inspect
from inspect import getfullargspec, ismethod
from functools import partial

from scipy.integrate import quadrature
from copy import deepcopy
import numpy as np


def keywords_first(f):
    def wrapper(*a, **k):
        a = list(a)
        for idx, arg in enumerate(getfullargspec(f).args, -ismethod(f)):
            if arg in k:
                if idx < len(a):
                    a.insert(idx, k.pop(arg))
                else:
                    break
        return f(*a, **k)
    return wrapper


def resolve_all_args(f, **kwargs):
    return partial(keywords_first(f), **kwargs)


def to_namedtuple(d, recursive=True):
    if isinstance(d, dict):
        d = d.copy()
        if recursive:
            for k, v in d.items():
                d[k] = to_namedtuple(v, recursive)
        d = namedtuple('_', d.keys())(**d)
    return d


def Integral(_func, x0, x1):
    result = quadrature(_func, a=x0, b=x1, maxiter=100000)[0]
    return result


def Functional(_func, x0, x1=None, du=False):
    if not isinstance(_func, Func):
        _func = funcs_1d.c(_func)
    if not _func.resolved:
        raise ValueError("function has to be resolved!")
    if du:
        if x1 is None:
            return _func * (Func(lambda x1: x1) - x0)
        else:
            return Func(f=lambda _: Integral(_func, x0=x0, x1=x1))
    else:
        if x1 is None:
            return Func(f=lambda x1: Integral(_func, x0=x0, x1=x1))
        else:
            return Func(f=lambda _: Integral(_func, x0=x0, x1=x1))


class Func(object):

    def __init__(self, f, alias='General'):
        self.f = f
        self.res_f = deepcopy(f)
        self.f_args = inspect.getfullargspec(f).args
        self.dim = len(self.f_args)
        self.alias = alias
        self.resolved = True if self.dim == 1 else False

    def resolve(self, **kwargs):
        if self.dim > 1:
            self.res_f = resolve_all_args(self.f, **kwargs)
            self.resolved = True
        else:
            raise ValueError("one-argument function cannot be resolved")
        return self

    def unresolve(self):
        if self.dim > 1:
            self.res_f = self.f
            self.resolved = False

    def __neg__(self):
        return -1 * self

    def __call__(self, other=None, **kwargs):
        if isinstance(other, Func):
            other_obj = other
            _new_obj = self.compose(other_obj)
            new_obj = Func(f=_new_obj)
            return new_obj
        if len(kwargs) > 0:
            res = self.resolve(**kwargs)
            if other is None: return res
        else:
            if not self.resolved:
                raise ValueError("Please resolve function for all but one argument!")
        try:
            if len(other) > 0:
                if isinstance(other, list) and len(other) == 1:
                    other = other[0]
                res = deepcopy([self.res_f(v) for v in other])
                return res
        except Exception as e:
            if not self.resolved: raise ValueError("Cannot call an unresolved function!")
            res = deepcopy(self.res_f(other))
            return res

    def __radd__(self, other):
        new_func = None
        if not isinstance(other, Func):
            new_func = Func(lambda _: other)
        return self + new_func

    def __rdiv__(self, other):
        new_func = None
        if not isinstance(other, Func):
            new_func = Func(lambda x: other)
        return self.__truediv__(new_func)

    def __add__(self, other):
        if not isinstance(other, Func):
            _new_obj = Func(lambda x: other)
        else:
            _new_obj = deepcopy(other)
        f = lambda _: deepcopy(self)(_) + _new_obj(_)
        new_func = Func(f=f, alias='from add')
        return new_func

    def __rsub__(self, other):
        new_obj = Func(f=lambda _: other, alias='from r_sub')
        return new_obj - self

    def __sub__(self, other):
        if not isinstance(other, Func):
            _new_obj = Func(lambda _: other)
        else:
            _new_obj = deepcopy(other)
        f = lambda _: self(_) - _new_obj(_)
        new_func = Func(f=f, alias='from subtract')
        return new_func

    def __rmul__(self, other):
        if not isinstance(other, Func):
            new_obj = Func(lambda _: other)
        else:
            new_obj = deepcopy(other)
        return self * new_obj

    def __mul__(self, other):
        if not isinstance(other, Func):
            _new_obj = Func(lambda _: other)
        else:
            _new_obj = deepcopy(other)
        f = lambda _: deepcopy(self)(_) * _new_obj(_)
        new_func = Func(f=f, alias='from mult')
        return new_func

    def __truediv__(self, other):
        if not isinstance(other, Func):
            _new_obj = Func(lambda _: other)
        else:
            _new_obj = deepcopy(other)
        f = lambda _: self(_) / _new_obj(_)
        new_func = Func(f=f, alias='from div')
        return new_func

    def compose(self, other):
        if not isinstance(other, Func):
            raise TypeError('Composition with strange functions')
        from functools import reduce

        def _compose(*fs):
            def compose2(k, l):
                return lambda _: k(l(_))
            return reduce(compose2, fs)
        f = _compose(deepcopy(self), other)
        return f


funcs_1d = to_namedtuple({'c': lambda c, alias='const': Func(f=lambda _: c, alias=alias),
                          'lin': Func(f=lambda _: _),
                          'exp': Func(f=lambda _: np.exp(_))})
