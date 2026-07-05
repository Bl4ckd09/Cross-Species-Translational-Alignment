#!/usr/bin/env python3
"""
run_mlp.py — the deferred neural net, finally run.  (needs torch -> .venv-chembert)

The main experiment used per-assay logistic regression because the MLP was blocked on the
Intel-Mac torch env. Now we run it, as a REGULARISED MULTITASK net: one shared trunk feeding
12 assay heads. The multitask sharing is the scientific bet — low-prevalence assays
(NR-PPAR-gamma n=12, SR-ATAD5 n=14) may borrow signal from the others, which a set of 12
independent logistic heads cannot.

Everything else is held identical to run_experiment.py so the ONLY change is the head:
  - same ECFP4-binary -> PCA-128 (structure) and logFC -> PCA-100 (GE), z-scored, fit in-fold
  - same repeated stratified CV, same folds/repeats/seeds, same masked-label handling
  - masked BCE loss (loss only on measured assays); per-assay pos_weight from the train fold
  - dropout + weight decay + a train-internal early-stopping split -> guard against N=177 overfit

Heads compared: structure_only (128), expr_only (100), fusion (228).
Writes: data/results/mlp_results.csv (+ .txt).  Run: .venv-chembert/bin/python scripts/run_mlp.py
"""
import os, json
import numpy as np
import pandas as pd
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.metrics import roc_auc_score, average_precision_score
import torch
import torch.nn as nn

import sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from structure_embed import featurize

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SIG  = os.path.join(ROOT, "data", "signatures")
RES  = os.path.join(ROOT, "data", "results"); os.makedirs(RES, exist_ok=True)

ASSAYS = ["NR-AR","NR-AR-LBD","NR-AhR","NR-Aromatase","NR-ER","NR-ER-LBD","NR-PPAR-gamma",
          "SR-ARE","SR-ATAD5","SR-HSE","SR-MMP","SR-p53"]
SR = ["SR-ARE","SR-ATAD5","SR-HSE","SR-MMP","SR-p53"]
NR = ["NR-AR","NR-AR-LBD","NR-AhR","NR-Aromatase","NR-ER","NR-ER-LBD","NR-PPAR-gamma"]

CFG = dict(n_pca_expr=100, n_pca_struct=128, folds=5, repeats=10, seed0=1000,
           hidden=64, dropout=0.4, weight_decay=1e-3, lr=1e-3, max_epochs=300, patience=25,
           signatures="drugmatrix_liver_logfc.parquet", labels="labels.csv")
APPROACHES = ["structure_only", "expr_only", "fusion"]


def load(cfg=CFG):
    sig = pd.read_parquet(os.path.join(SIG, cfg["signatures"]))
    lab = pd.read_csv(os.path.join(SIG, cfg["labels"])).set_index("connectivity").loc[sig.index]
    st  = featurize(list(sig.index), kind="ecfp4").loc[sig.index]     # SAME features as logistic run
    Y = lab[ASSAYS].apply(pd.to_numeric, errors="coerce").values.astype(float)
    return sig.index.to_numpy(), sig.values.astype("float32"), st.values.astype("float32"), Y


class MultiTaskMLP(nn.Module):
    """Shared trunk -> 12 assay logits. Small + heavily regularised for N~140/fold."""
    def __init__(self, in_dim, n_tasks=12, hidden=64, dropout=0.4):
        super().__init__()
        self.trunk = nn.Sequential(
            nn.Linear(in_dim, hidden), nn.ReLU(), nn.Dropout(dropout),
            nn.Linear(hidden, hidden // 2), nn.ReLU(), nn.Dropout(dropout),
        )
        self.head = nn.Linear(hidden // 2, n_tasks)

    def forward(self, x):
        return self.head(self.trunk(x))


def masked_bce(logits, y, mask, pos_weight):
    loss = nn.functional.binary_cross_entropy_with_logits(
        logits, torch.nan_to_num(y), weight=None, pos_weight=pos_weight, reduction="none")
    return (loss * mask).sum() / mask.sum().clamp(min=1)


def train_predict(Ftr, Ytr, Fte, cfg, seed):
    """Train a multitask MLP on train fold (with an internal val split for early stopping),
    return sigmoid probabilities for the test fold. Ytr may contain NaN (unmeasured)."""
    torch.manual_seed(seed)
    Ntr = Ftr.shape[0]
    # internal 85/15 stratified-ish split on active-count for early stopping
    rng = np.random.default_rng(seed)
    perm = rng.permutation(Ntr)
    nval = max(8, int(0.15 * Ntr))
    vidx, tidx = perm[:nval], perm[nval:]

    Xt = torch.tensor(Ftr[tidx]); Yt = torch.tensor(Ytr[tidx])
    Xv = torch.tensor(Ftr[vidx]); Yv = torch.tensor(Ytr[vidx])
    Mt = (~torch.isnan(Yt)).float(); Mv = (~torch.isnan(Yv)).float()

    # per-assay pos_weight from measured TRAIN labels (neg/pos), clamped
    pw = []
    for j in range(len(ASSAYS)):
        col = Ytr[tidx, j]; col = col[~np.isnan(col)]
        npos = float((col == 1).sum()); nneg = float((col == 0).sum())
        pw.append(np.clip(nneg / npos, 0.2, 20.0) if npos > 0 else 1.0)
    pos_weight = torch.tensor(pw, dtype=torch.float32)

    model = MultiTaskMLP(Ftr.shape[1], len(ASSAYS), cfg["hidden"], cfg["dropout"])
    opt = torch.optim.AdamW(model.parameters(), lr=cfg["lr"], weight_decay=cfg["weight_decay"])

    best_val, best_state, bad = float("inf"), None, 0
    for ep in range(cfg["max_epochs"]):
        model.train(); opt.zero_grad()
        loss = masked_bce(model(Xt), Yt, Mt, pos_weight)
        loss.backward(); opt.step()
        model.eval()
        with torch.no_grad():
            vloss = masked_bce(model(Xv), Yv, Mv, pos_weight).item()
        if vloss < best_val - 1e-4:
            best_val, best_state, bad = vloss, {k: v.clone() for k, v in model.state_dict().items()}, 0
        else:
            bad += 1
            if bad >= cfg["patience"]:
                break
    if best_state is not None:
        model.load_state_dict(best_state)
    model.eval()
    with torch.no_grad():
        p = torch.sigmoid(model(torch.tensor(Fte))).numpy()
    return p


def build_blocks(GE, ST, tr, te, cfg):
    """In-fold standardise+PCA+z-score for GE and structure; return z-scored blocks."""
    def block(mat, npc):
        sc = StandardScaler().fit(mat[tr])
        pca = PCA(n_components=min(npc, len(tr) - 1), random_state=0)
        Ztr = pca.fit_transform(sc.transform(mat[tr])); Zte = pca.transform(sc.transform(mat[te]))
        z = StandardScaler().fit(Ztr)
        return z.transform(Ztr).astype("float32"), z.transform(Zte).astype("float32")
    Zge_tr, Zge_te = block(GE, cfg["n_pca_expr"])
    Zst_tr, Zst_te = block(ST, cfg["n_pca_struct"])
    return {
        "structure_only": (Zst_tr, Zst_te),
        "expr_only":      (Zge_tr, Zge_te),
        "fusion":         (np.hstack([Zst_tr, Zge_tr]), np.hstack([Zst_te, Zge_te])),
    }


def evaluate(cfg, data):
    conns, GE, ST, Y = data
    N = len(conns)
    strat = pd.qcut(np.nansum(Y == 1, axis=1), 4, labels=False, duplicates="drop")
    rec = {ap: {a: {"auc": [], "ap": []} for a in ASSAYS} for ap in APPROACHES}
    for r in range(cfg["repeats"]):
        skf = StratifiedKFold(cfg["folds"], shuffle=True, random_state=cfg["seed0"] + r)
        oof = {ap: np.full((N, len(ASSAYS)), np.nan) for ap in APPROACHES}
        for fold, (tr, te) in enumerate(skf.split(np.zeros(N), strat)):
            blocks = build_blocks(GE, ST, tr, te, cfg)
            for ap, (Ftr, Fte) in blocks.items():
                p = train_predict(Ftr, Y[tr], Fte, cfg, seed=cfg["seed0"] + r * 100 + fold)
                for j in range(len(ASSAYS)):
                    mte = ~np.isnan(Y[te, j])
                    oof[ap][te[np.where(mte)[0]], j] = p[mte, j]
        for ap in APPROACHES:
            for j, a in enumerate(ASSAYS):
                yj, pj = Y[:, j], oof[ap][:, j]
                m = ~np.isnan(yj) & ~np.isnan(pj)
                if yj[m].sum() >= 2 and (1 - yj[m]).sum() >= 2:
                    rec[ap][a]["auc"].append(roc_auc_score(yj[m], pj[m]))
                    rec[ap][a]["ap"].append(average_precision_score(yj[m], pj[m]))
        print(f"  repeat {r+1}/{cfg['repeats']} done", flush=True)
    return rec


def macro(rec, ap, assays=ASSAYS):
    return float(np.mean([np.mean(rec[ap][a]["auc"]) for a in assays]))


def run():
    torch.set_num_threads(max(1, (os.cpu_count() or 2) - 1))
    data = load(CFG)
    print("loaded:", len(data[0]), "compounds | torch", torch.__version__, "\n")
    rec = evaluate(CFG, data)

    def ci(x):
        x = np.array(x); return (x.mean(), 1.96 * x.std(ddof=1) / np.sqrt(len(x))) if len(x) > 1 else (np.nan, np.nan)

    rows = []
    for j, a in enumerate(ASSAYS):
        row = {"assay": a, "panel": "SR" if a in SR else "NR",
               "n_active": int(np.nansum(data[3][:, j] == 1))}
        for ap in APPROACHES:
            m, e = ci(rec[ap][a]["auc"]); row[f"{ap}_AUC"] = round(m, 3)
        d = np.array(rec["fusion"][a]["auc"]) - np.array(rec["structure_only"][a]["auc"])
        row["dAUC_fusion_minus_struct"] = round(d.mean(), 3)
        rows.append(row)
    tab = pd.DataFrame(rows)

    def add_summary(name, assays):
        r = {"assay": name, "panel": "", "n_active": "", "dAUC_fusion_minus_struct": ""}
        for ap in APPROACHES:
            r[f"{ap}_AUC"] = round(macro(rec, ap, assays), 3)
        return r
    tab = pd.concat([tab, pd.DataFrame([add_summary("MACRO", ASSAYS),
                                        add_summary("SR_mean", SR),
                                        add_summary("NR_mean", NR)])], ignore_index=True)
    tab.to_csv(os.path.join(RES, "mlp_results.csv"), index=False)

    def dstats(assays):
        ds = []
        for a in assays:
            x = np.array(rec["fusion"][a]["auc"]); y = np.array(rec["structure_only"][a]["auc"])
            n = min(len(x), len(y)); ds.append((x[:n] - y[:n]).mean())
        return float(np.mean(ds))

    lines = [
        "MULTITASK MLP  (N=177 DrugMatrix liver, %dx%d CV)" % (CFG["repeats"], CFG["folds"]),
        f"trunk={CFG['hidden']}->{CFG['hidden']//2}  dropout={CFG['dropout']}  wd={CFG['weight_decay']}  "
        f"early-stop patience={CFG['patience']}",
        "",
        "MACRO ROC-AUC:",
        *[f"  {ap:16s} {macro(rec, ap):.3f}" for ap in APPROACHES],
        "",
        "reference — per-assay LOGISTIC (same features):  structure 0.757  expr 0.679  fusion 0.766",
        "",
        f"GE-adds delta (fusion - structure):  macro={dstats(ASSAYS):+.3f}  "
        f"SR={dstats(SR):+.3f}  NR={dstats(NR):+.3f}",
    ]
    txt = "\n".join(lines)
    with open(os.path.join(RES, "mlp_results.txt"), "w") as f:
        f.write(txt + "\n")
    with open(os.path.join(RES, "mlp_config.json"), "w") as f:
        json.dump({"config": CFG, "n_compounds": int(len(data[0])), "assays": ASSAYS}, f, indent=2)

    pd.set_option("display.width", 220, "display.max_columns", 30)
    print(tab.to_string(index=False))
    print("\n" + txt)
    print(f"\nwrote {RES}/mlp_results.csv  and  mlp_results.txt")


if __name__ == "__main__":
    run()
