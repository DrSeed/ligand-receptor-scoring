# Self-contained demo: inflamed bile-duct niche with planted T cell-cholangiocyte signalling.
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from lr import LR_PAIRS, score_all, consensus

TYPES = ["T cell", "cholangiocyte", "hepatocyte", "macrophage", "endothelial"]
TRUE = {("T cell", "cholangiocyte", "CD226", "NECTIN2"), ("endothelial", "T cell", "CXCL12", "CXCR4"),
        ("macrophage", "T cell", "SPP1", "CD44"), ("T cell", "macrophage", "PDCD1", "CD274")}


def simulate(seed=21, n_bg=150):
    rng = np.random.default_rng(seed)
    lr_genes = sorted({g for pair in LR_PAIRS for g in pair})
    genes = lr_genes + [f"bg{i}" for i in range(n_bg)]
    gi = {g: i for i, g in enumerate(genes)}
    base = rng.lognormal(-1.5, 0.8, len(genes))
    prof = np.tile(base, (len(TYPES), 1))
    for g, level in [("MIF", 6.0), ("CD74", 4.0), ("APOE", 5.0), ("LRP1", 2.0)]:  # abundant, everywhere
        prof[:, gi[g]] = level * rng.uniform(0.7, 1.3, len(TYPES))
    up = {("T cell", "CD226"): 1.2, ("cholangiocyte", "NECTIN2"): 1.0, ("endothelial", "CXCL12"): 1.5,
          ("T cell", "CXCR4"): 1.0, ("macrophage", "SPP1"): 1.5, ("T cell", "CD44"): 0.8,
          ("T cell", "PDCD1"): 0.6, ("macrophage", "CD274"): 0.8}
    for (t, g), level in up.items():
        prof[TYPES.index(t), gi[g]] = level
    n_per = [300, 150, 400, 200, 120]
    labels = np.repeat(np.array(TYPES), n_per)
    size = rng.lognormal(0, 0.3, len(labels))[:, None]
    X = np.log1p(rng.poisson(prof[[TYPES.index(t) for t in labels]] * size * 3))
    return X, labels, genes


def precision_at(df, col, k):
    top = df.nsmallest(k, col)
    return np.mean([(s, r, l, rc) in TRUE for s, r, l, rc in top[["sender", "receiver", "ligand", "receptor"]].values])


def main():
    os.makedirs("figures", exist_ok=True)
    os.makedirs("results", exist_ok=True)
    X, labels, genes = simulate()
    df = score_all(X, labels, genes, TYPES)
    df, weights = consensus(df)
    df["true_interaction"] = [(s, r, l, rc) in TRUE for s, r, l, rc in df[["sender", "receiver", "ligand", "receptor"]].values]
    df.sort_values("consensus_rank").to_csv("results/all_interactions.csv", index=False)
    rows = []
    for m in ["product", "cpdb", "perm_fold", "specificity", "consensus"]:
        col = "consensus_rank" if m == "consensus" else f"rank_{m}"
        rows.append({"method": m, "weight_in_consensus": weights.get(m, np.nan),
                     "precision_at_5": precision_at(df, col, 5), "precision_at_10": precision_at(df, col, 10),
                     "median_rank_of_true": float(df.loc[df.true_interaction, col].median())})
    table = pd.DataFrame(rows).round(3)
    table.to_csv("results/summary.csv", index=False)

    fig, axes = plt.subplots(1, 3, figsize=(20, 5.6), gridspec_kw={"width_ratios": [1, 1, 1.15]})
    for ax, col, title in [(axes[0], "rank_product", "Top 10 by expression product"),
                           (axes[1], "consensus_rank", "Top 10 by agreement-weighted consensus")]:
        top = df.nsmallest(10, col)
        names = [f"{s} -> {r}: {l}-{rc}" for s, r, l, rc in top[["sender", "receiver", "ligand", "receptor"]].values]
        ax.barh(names, 11 - np.arange(1, 11), color=["#55a868" if t else "#bbbbbb" for t in top.true_interaction])
        ax.invert_yaxis()
        ax.set_xticks([])
        ax.tick_params(axis="y", labelsize=8)
        ax.set_title(f"{title}\n(green = planted, cell-type-specific signalling)")
    x = np.arange(len(table))
    axes[2].bar(x - 0.2, table.precision_at_5, 0.4, color="#4c72b0", label="precision@5")
    axes[2].bar(x + 0.2, table.precision_at_10, 0.4, color="#9db7de", label="precision@10")
    labels_ = [f"{m}\n(weight {w:.2f})" if not np.isnan(w) else m for m, w in zip(table.method, table.weight_in_consensus)]
    axes[2].set_xticks(x, labels_, fontsize=8, rotation=25)
    axes[2].set_ylim(0, 1.05)
    axes[2].set_ylabel("fraction of top hits that are planted interactions")
    axes[2].axhline(len(TRUE) / 5, color="k", ls=":", lw=0.8)
    axes[2].set_title("Abundant genes (MIF-CD74) swamp product and p-value\nrankings; effect sizes find the niche (dotted: best @5)")
    axes[2].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig("figures/demo.png", dpi=120)
    print(table.to_string(index=False))
    print(df.nsmallest(6, "consensus_rank")[["sender", "receiver", "ligand", "receptor", "cpdb_p", "true_interaction"]].to_string(index=False))


if __name__ == "__main__":
    main()
