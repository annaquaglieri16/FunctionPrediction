"""A simple frozen-embedding probe: logistic regression on top of PLM
embeddings. Deliberately simple -- the point of this benchmark is the
embeddings, not the probe."""

from __future__ import annotations

import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, f1_score
from sklearn.preprocessing import StandardScaler


def train_probe(X: np.ndarray, y, max_iter: int = 2000):
    scaler = StandardScaler().fit(X)
    clf = LogisticRegression(max_iter=max_iter)
    clf.fit(scaler.transform(X), y)
    return scaler, clf


def evaluate_probe(scaler, clf, X: np.ndarray, y) -> dict:
    preds = clf.predict(scaler.transform(X))
    return {
        "accuracy": accuracy_score(y, preds),
        "f1_macro": f1_score(y, preds, average="macro"),
        "n": len(y),
    }
