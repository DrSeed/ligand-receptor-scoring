# Ligand-receptor scoring three ways, with label permutation and an agreement-weighted consensus.
import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr

LR_PAIRS = [("CD226", "NECTIN2"), ("CD226", "PVR"), ("TIGIT", "PVR"), ("CD28", "CD80"), ("CTLA4", "CD80"),
            ("PDCD1", "CD274"), ("CXCL12", "CXCR4"), ("CCL19", "CCR7"), ("TNF", "TNFRSF1A"),
            ("SPP1", "CD44"), ("MIF", "CD74"), ("APOE", "LRP1")]


def type_means(X, labels, types):
    onehot = (labels[:, None] == np.array(types)[None, :]).astype(float)
    return (onehot.T @ X) / onehot.sum(0)[:, None]


def candidates(types):
    return [(s, r, l, rc) for s in types for r in types for l, rc in LR_PAIRS]


def score_all(X, labels, genes, types, n_perm=500, seed=0):
    """For every sender, receiver and ligand-receptor pair:
    product      mean ligand in sender x mean receptor in receiver (favours abundant genes)
    cpdb_p       CellPhoneDB-style: mean of the two means vs a cell-label permutation null (p-value)
    perm_fold    the same statistic divided by its permutation average (an effect size)
    specificity  z-score of the ligand across senders x z-score of the receptor across receivers,
                 counted only if both are detected in >= 10% of the cells in question."""
    gi = {g: i for i, g in enumerate(genes)}
    M = type_means(X, labels, types)
    F = type_means((X > 0).astype(float), labels, types)
    Z = (M - M.mean(0)) / (M.std(0) + 1e-9)
    ti = {t: i for i, t in enumerate(types)}
    cand = candidates(types)
    s_idx = np.array([ti[s] for s, _, _, _ in cand]); r_idx = np.array([ti[r] for _, r, _, _ in cand])
    l_idx = np.array([gi[l] for _, _, l, _ in cand]); rc_idx = np.array([gi[rc] for _, _, _, rc in cand])
    lm, rm = M[s_idx, l_idx], M[r_idx, rc_idx]
    obs = (lm + rm) / 2
    rng = np.random.default_rng(seed)
    exceed, null_sum = np.zeros(len(cand)), np.zeros(len(cand))
    for _ in range(n_perm):
        Mp = type_means(X, rng.permutation(labels), types)
        null = (Mp[s_idx, l_idx] + Mp[r_idx, rc_idx]) / 2
        exceed += null >= obs
        null_sum += null
    detected = (F[s_idx, l_idx] >= 0.1) & (F[r_idx, rc_idx] >= 0.1)
    df = pd.DataFrame(cand, columns=["sender", "receiver", "ligand", "receptor"])
    df["product"] = lm * rm
    df["cpdb_p"] = (exceed + 1) / (n_perm + 1)
    df["perm_fold"] = obs / (null_sum / n_perm + 1e-9)  # effect size: how far above the shuffled-label average
    df["specificity"] = np.where(detected, np.clip(Z[s_idx, l_idx], 0, None) * np.clip(Z[r_idx, rc_idx], 0, None), 0)
    return df


def consensus(df, methods=("product", "cpdb", "perm_fold", "specificity")):
    """CellCom2-style integration: turn each method into ranks, weight each method by how well it agrees
    with the others (mean Spearman correlation, floored at 0), and average the weighted ranks."""
    ranks = pd.DataFrame({
        "product": rankdata(-df["product"]),
        "cpdb": rankdata(df["cpdb_p"] - 1e-9 * df["product"]),  # ties in p broken by expression
        "perm_fold": rankdata(-df["perm_fold"]),
        "specificity": rankdata(-df["specificity"]),
    })[list(methods)]
    corr = ranks.corr(method="spearman").values
    agree = (corr.sum(1) - 1) / (len(methods) - 1)
    w = np.clip(agree, 0, None)
    w = w / w.sum() if w.sum() > 0 else np.full(len(methods), 1 / len(methods))
    out = df.copy()
    for m in methods:
        out[f"rank_{m}"] = ranks[m].values
    out["consensus_rank"] = rankdata(ranks.values @ w)
    return out, dict(zip(methods, np.round(w, 3)))
