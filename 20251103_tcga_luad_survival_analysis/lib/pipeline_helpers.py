import numpy as np
from sksurv.linear_model import CoxPHSurvivalAnalysis


def fit_and_score_features(X, y):
    n_features = X.shape[1]
    scores = np.empty(n_features)
    m = CoxPHSurvivalAnalysis(alpha=0, tol=1e-3, n_iter=50)
    for j in range(n_features):
        Xj = X[:, j : j + 1]
        m.fit(Xj, y)
        scores[j] = m.score(Xj, y)
    return scores


def p_to_var(p):
    return p * (1-p)