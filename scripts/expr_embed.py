#!/usr/bin/env python3
"""
expr_embed.py — the SWAPPABLE expression-embedding interface.

    embedder = make_embedder("pca", n_components=100)
    Z_train  = embedder.fit_transform(X_train)   # fit uses TRAIN FOLD ONLY
    Z_test   = embedder.transform(X_test)         # test never seen during fit

Today: PCA(100). Tomorrow: drop in a rat->human translation model (e.g. PLOS 2023, 64-D
latent) with the SAME fit_transform / transform contract — no downstream change.
Every embedder standardises genes on the train fold before reducing (leakage-safe).
"""
import numpy as np
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler


class PCAEmbedder:
    """Standardise (train-fit) -> PCA (train-fit). The leakage-safe default."""
    def __init__(self, n_components=100, random_state=0):
        self.n_components = n_components
        self.scaler = StandardScaler()
        self.pca = None
        self.random_state = random_state

    def fit_transform(self, X):
        Xs = self.scaler.fit_transform(X)
        n = min(self.n_components, min(Xs.shape) - 1)
        self.pca = PCA(n_components=n, random_state=self.random_state)
        return self.pca.fit_transform(Xs).astype("float32")

    def transform(self, X):
        return self.pca.transform(self.scaler.transform(X)).astype("float32")

    @property
    def out_dim(self):
        return self.pca.n_components_


class IdentityEmbedder:
    """Passthrough (standardise only) — for ablating the reduction step."""
    def __init__(self): self.scaler = StandardScaler()
    def fit_transform(self, X): return self.scaler.fit_transform(X).astype("float32")
    def transform(self, X):     return self.scaler.transform(X).astype("float32")


# --- future: class RatToHumanEmbedder with the identical contract -------------
# class RatToHumanEmbedder:
#     def fit_transform(self, X): ...   # learn/apply rat->human mapping on train fold
#     def transform(self, X): ...

def make_embedder(kind="pca", **kw):
    return {"pca": PCAEmbedder, "identity": IdentityEmbedder}[kind](**kw)
