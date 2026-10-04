"""Python re-implementation of coloc.abf (Giambartolomei 2014; Wallace coloc v5 defaults).

Mirrors coloc's R code (claudia.R / approx.bf.*, combine.abf, logsum, logdiff):
  * prior SD of the true effect W: 0.15 * sdY for quantitative traits, 0.2 for case-control
  * ABF (Wakefield 2009):  lABF = 0.5*log(1-r) + 0.5*r*z^2,  r = W^2 / (W^2 + V)
  * from beta/varbeta: z = beta/sqrt(V);  from p-values: z^2 = qchisq(p, 1, upper)
    and V = sdY^2 / (2 N f (1-f))  (quant) or 1 / (2 N f (1-f) s (1-s))  (case-control)
  * default priors p1 = 1e-4, p2 = 1e-4, p12 = 1e-5
It is validated against the CRAN coloc package on identical inputs (s05_validate_coloc.*)
and against the PP values Bryois et al. 2022 published (positive control).

Also: a minimal SuSiE-RSS (Wang 2020 / Zou 2022, IBSS on summary statistics, fixed prior
variance estimated per effect by EM) and coloc.bf_bf (pairwise coloc over credible sets).
"""
from __future__ import annotations

import numpy as np
from scipy import stats
from scipy.special import logsumexp


def logdiff(a: float, b: float) -> float:
    """log(exp(a) - exp(b)) for a > b, as in coloc."""
    mx = max(a, b)
    return mx + np.log(np.exp(a - mx) - np.exp(b - mx))


def labf_beta(beta, varbeta, sdY=1.0, kind="quant"):
    beta = np.asarray(beta, float)
    V = np.asarray(varbeta, float)
    W = 0.15 * sdY if kind == "quant" else 0.2
    z = beta / np.sqrt(V)
    r = W ** 2 / (W ** 2 + V)
    return 0.5 * np.log(1 - r) + 0.5 * r * z ** 2


def labf_p(p, maf, N, sdY=1.0, kind="quant", s=None):
    p = np.asarray(p, float)
    f = np.asarray(maf, float)
    f = np.minimum(f, 1 - f)
    N = np.asarray(N, float)
    if kind == "quant":
        V = sdY ** 2 / (2 * N * f * (1 - f))
        W = 0.15 * sdY
    else:
        V = 1.0 / (2 * N * f * (1 - f) * s * (1 - s))
        W = 0.2
    z2 = stats.chi2.isf(p, 1)
    r = W ** 2 / (W ** 2 + V)
    return 0.5 * np.log(1 - r) + 0.5 * r * z2


def combine_abf(l1, l2, p1=1e-4, p2=1e-4, p12=1e-5):
    l1 = np.asarray(l1, float)
    l2 = np.asarray(l2, float)
    lsum = l1 + l2
    lH0 = 0.0
    lH1 = np.log(p1) + logsumexp(l1)
    lH2 = np.log(p2) + logsumexp(l2)
    lH3 = np.log(p1) + np.log(p2) + logdiff(logsumexp(l1) + logsumexp(l2), logsumexp(lsum))
    lH4 = np.log(p12) + logsumexp(lsum)
    all_ = np.array([lH0, lH1, lH2, lH3, lH4])
    pp = np.exp(all_ - logsumexp(all_))
    snp_pp = np.exp(lsum - logsumexp(lsum))
    return {"nsnps": int(len(l1)), "PP.H0": pp[0], "PP.H1": pp[1], "PP.H2": pp[2],
            "PP.H3": pp[3], "PP.H4": pp[4], "snp_pp_h4": snp_pp}


def coloc_abf(l1, l2, p1=1e-4, p2=1e-4, p12=1e-5):
    return combine_abf(l1, l2, p1, p2, p12)


def p12_sensitivity(l1, l2, p1=1e-4, p2=1e-4, grid=None):
    """PP.H4 / PP.H3 across p12 (coloc::sensitivity analogue)."""
    grid = grid if grid is not None else [1e-8, 1e-7, 5e-7, 1e-6, 5e-6, 1e-5, 5e-5, 1e-4]
    out = []
    for g in grid:
        r = combine_abf(l1, l2, p1, p2, min(g, min(p1, p2)))
        out.append({"p12": g, "PP.H3": r["PP.H3"], "PP.H4": r["PP.H4"]})
    return out


def p12_for_h4(l1, l2, target=0.5, p1=1e-4, p2=1e-4):
    """Smallest p12 at which PP.H4 reaches target (log-grid search); None if never."""
    for g in np.logspace(-10, -4, 241):
        if combine_abf(l1, l2, p1, p2, g)["PP.H4"] >= target:
            return float(g)
    return None


# ---------------------------------------------------------------------------------------------
# SuSiE-RSS (summary statistics + LD), single-trait, used to split a locus into signals
# ---------------------------------------------------------------------------------------------

def susie_rss(z, R, n, L=5, max_iter=200, tol=1e-4, prior_var_init=None, coverage=0.95,
              min_abs_corr=0.5, lam=1e-3):
    """IBSS for SuSiE-RSS on z-scores with LD matrix R (regularised R + lam*I).

    Works on the standardized scale: XtX = (n-1) R, Xty = sqrt(n-1) z (approximately, with
    var(y) = 1). Returns per-effect alpha (PIP vectors), lbf per SNP per effect, credible sets.
    """
    z = np.asarray(z, float)
    p = len(z)
    R = (1 - lam) * np.asarray(R, float) + lam * np.eye(p)
    nm1 = n - 1.0
    XtX = nm1 * R
    Xty = np.sqrt(nm1) * z / np.sqrt(1 + z ** 2 / nm1)  # PVE-adjusted (susie_rss default)
    d = np.diag(XtX)
    sigma2 = 1.0  # var(y) fixed at 1 for z-score input (susie_rss with estimate_residual_variance=False)
    alpha = np.full((L, p), 1.0 / p)
    mu = np.zeros((L, p))
    mu2 = np.zeros((L, p))
    V = np.full(L, 0.2 if prior_var_init is None else prior_var_init)
    lbf = np.zeros((L, p))
    lbf_model = np.zeros(L)
    XtXr = XtX @ (alpha * mu).sum(0)
    for it in range(max_iter):
        old = alpha.copy()
        for l in range(L):
            bl = alpha[l] * mu[l]
            XtXr = XtXr - XtX @ bl
            XtR = Xty - XtXr
            betahat = XtR / d
            shat2 = sigma2 / d
            # EM-free optimisation of prior variance (susie "optim" ~ grid search on log scale)
            def lbf_of(v):
                return 0.5 * np.log(shat2 / (v + shat2)) + 0.5 * betahat ** 2 / shat2 * v / (v + shat2)
            grid = np.concatenate([[0.0], np.logspace(-4, 1, 60)])
            best, bestv = -np.inf, 0.0
            for v in grid:
                if v == 0:
                    val = 0.0
                else:
                    lb = lbf_of(v)
                    val = logsumexp(lb) - np.log(p)
                if val > best:
                    best, bestv = val, v
            V[l] = bestv
            if V[l] == 0:
                lb = np.zeros(p)
            else:
                lb = lbf_of(V[l])
            lbf[l] = lb
            lbf_model[l] = logsumexp(lb) - np.log(p)
            w = lb - lb.max()
            a = np.exp(w)
            a /= a.sum()
            post_var = 1.0 / (1.0 / V[l] + d / sigma2) if V[l] > 0 else np.zeros(p)
            post_mean = post_var * XtR / sigma2 if V[l] > 0 else np.zeros(p)
            alpha[l] = a
            mu[l] = post_mean
            mu2[l] = post_var + post_mean ** 2
            XtXr = XtXr + XtX @ (alpha[l] * mu[l])
        if np.max(np.abs(alpha - old)) < tol:
            break
    pip = 1 - np.prod(1 - alpha, axis=0)
    cs = []
    Rabs = np.abs(R)
    for l in range(L):
        if V[l] <= 0:
            continue
        order = np.argsort(-alpha[l])
        cum = np.cumsum(alpha[l][order])
        k = int(np.searchsorted(cum, coverage) + 1)
        idx = order[:k]
        purity = Rabs[np.ix_(idx, idx)].min() if len(idx) > 1 else 1.0
        if purity >= min_abs_corr:
            cs.append({"effect": l, "idx": idx.tolist(), "coverage": float(cum[k - 1]),
                       "min_abs_corr": float(purity), "prior_var": float(V[l]),
                       "lbf_model": float(lbf_model[l])})
    return {"alpha": alpha, "lbf": lbf, "pip": pip, "cs": cs, "V": V, "iters": it + 1}


def coloc_bf_bf(lbf1, lbf2, p1=1e-4, p2=1e-4, p12=5e-6):
    """coloc.bf_bf for one pair of SNP-wise log BF vectors (one SuSiE effect from each trait).

    coloc.susie uses p12 = 5e-6 by default.
    """
    return combine_abf(lbf1, lbf2, p1, p2, p12)
