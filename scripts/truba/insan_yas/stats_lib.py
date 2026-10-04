"""Small statistics helpers for the human KIF13A age analysis (numpy/scipy only).

Everything here works on donor-level values: one number per donor and cell type.
No test in this project ever uses cells as replicates.
"""
import re

import numpy as np
import pandas as pd
from scipy import stats


def parse_age(label):
    """Years from a CELLxGENE development_stage label. Prenatal -> negative, unknown -> NaN.

    '80 year-old and over stage' is open ended and returns NaN on purpose: it cannot
    be placed on a continuous age axis.
    """
    if label is None or (isinstance(label, float) and np.isnan(label)):
        return np.nan
    s = str(label)
    m = re.match(r"^(\d+)-year-old", s)
    if m:
        return float(m.group(1))
    m = re.match(r"^(\d+)-month-old", s)
    if m:
        return float(m.group(1)) / 12.0
    m = re.match(r"^(\d+)-week-old", s)
    if m:
        return float(m.group(1)) / 52.0
    if s.startswith("newborn"):
        return 0.04
    if "LMP month" in s or "post-fertilization" in s or "Carnegie" in s or "fetal" in s:
        return -0.3
    return np.nan


def ols(Y, X):
    """OLS for many responses sharing one design. Y: n x k, X: n x p (with intercept).

    Returns beta (p x k), se (p x k), t, p-values, residual df, residual sd (k).
    """
    Y = np.asarray(Y, float)
    X = np.asarray(X, float)
    if Y.ndim == 1:
        Y = Y[:, None]
    n, p = X.shape
    XtX_inv = np.linalg.pinv(X.T @ X)
    B = XtX_inv @ X.T @ Y
    R = Y - X @ B
    df = n - p
    s2 = (R ** 2).sum(axis=0) / df
    se = np.sqrt(np.outer(np.diag(XtX_inv), s2))
    t = B / se
    pv = 2 * stats.t.sf(np.abs(t), df)
    return B, se, t, pv, df, np.sqrt(s2)


def cliffs_delta(a, b):
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    gt = (a[:, None] > b[None, :]).sum()
    lt = (a[:, None] < b[None, :]).sum()
    return (gt - lt) / (len(a) * len(b))


def two_group(old, young):
    old = np.asarray(old, float)
    young = np.asarray(young, float)
    old, young = old[~np.isnan(old)], young[~np.isnan(young)]
    if len(old) < 3 or len(young) < 3:
        return dict(n_young=len(young), n_old=len(old))
    _, p_mw = stats.mannwhitneyu(old, young, alternative="two-sided")
    _, p_t = stats.ttest_ind(old, young, equal_var=False)
    return dict(n_young=len(young), n_old=len(old), mean_young=young.mean(), mean_old=old.mean(),
                delta=old.mean() - young.mean(), fold=2 ** (old.mean() - young.mean()),
                mannwhitney_p=p_mw, welch_p=p_t, cliffs_delta=cliffs_delta(old, young),
                sd_young=young.std(ddof=1), sd_old=old.std(ddof=1))


def power_two_group(delta, sd1, sd2, n1, n2, alpha=0.05):
    """Power of a two-sided Welch t test for a true mean difference delta."""
    se = np.sqrt(sd1 ** 2 / n1 + sd2 ** 2 / n2)
    df = (sd1 ** 2 / n1 + sd2 ** 2 / n2) ** 2 / ((sd1 ** 2 / n1) ** 2 / (n1 - 1) + (sd2 ** 2 / n2) ** 2 / (n2 - 1))
    nc = delta / se
    tc = stats.t.ppf(1 - alpha / 2, df)
    return float(stats.nct.sf(tc, df, nc) + stats.nct.cdf(-tc, df, nc))


def power_slope(slope, resid_sd, x, alpha=0.05, n_cov=0):
    """Power of the two-sided t test of a regression slope, given the observed x values."""
    x = np.asarray(x, float)
    n = len(x)
    se = resid_sd / np.sqrt(((x - x.mean()) ** 2).sum())
    df = n - 2 - n_cov
    nc = slope / se
    tc = stats.t.ppf(1 - alpha / 2, df)
    return float(stats.nct.sf(tc, df, nc) + stats.nct.cdf(-tc, df, nc))


def meta_analysis(est, se):
    """Fixed effect and DerSimonian-Laird random effects pooling."""
    est = np.asarray(est, float)
    se = np.asarray(se, float)
    w = 1 / se ** 2
    fe = (w * est).sum() / w.sum()
    fe_se = np.sqrt(1 / w.sum())
    q = (w * (est - fe) ** 2).sum()
    k = len(est)
    tau2 = max(0.0, (q - (k - 1)) / (w.sum() - (w ** 2).sum() / w.sum())) if k > 1 else 0.0
    wr = 1 / (se ** 2 + tau2)
    re_ = (wr * est).sum() / wr.sum()
    re_se = np.sqrt(1 / wr.sum())
    i2 = max(0.0, (q - (k - 1)) / q) if q > 0 and k > 1 else 0.0
    return dict(k=k, fe=fe, fe_se=fe_se, fe_p=2 * stats.norm.sf(abs(fe / fe_se)),
                re=re_, re_se=re_se, re_p=2 * stats.norm.sf(abs(re_ / re_se)),
                tau2=tau2, Q=q, Q_p=stats.chi2.sf(q, k - 1) if k > 1 else np.nan, I2=i2)


def lowess_curve(x, y, frac=0.35, grid=None):
    """Tiny LOWESS (tricube, local linear) for descriptive lifespan curves."""
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    ok = ~(np.isnan(x) | np.isnan(y))
    x, y = x[ok], y[ok]
    if grid is None:
        grid = np.linspace(x.min(), x.max(), 60)
    k = max(5, int(np.ceil(frac * len(x))))
    out = []
    for g in grid:
        d = np.abs(x - g)
        h = np.sort(d)[k - 1] + 1e-12
        w = np.clip(1 - (d / h) ** 3, 0, None) ** 3
        W = w.sum()
        if W <= 0:
            out.append(np.nan)
            continue
        xm = (w * x).sum() / W
        ym = (w * y).sum() / W
        sxx = (w * (x - xm) ** 2).sum()
        b = (w * (x - xm) * (y - ym)).sum() / sxx if sxx > 0 else 0.0
        out.append(ym + b * (g - xm))
    return np.asarray(grid), np.asarray(out)
