from copy import deepcopy

import matplotlib.pyplot as plt
import pandas as pd
from sksurv.nonparametric import kaplan_meier_estimator


def plot_bar(df, col, dropna=False, **kwargs):
    ax = df[col].value_counts(dropna=dropna).sort_index().plot(kind='bar', title=col, **kwargs)
    ax.set_xlabel('')
    return ax


def plot_hist(df, col, **kwargs):
    ax = df[col].plot(kind='hist', title=col, **kwargs)
    return ax


def side_by_side_bar_plots(df, col1, col2, figsize=(10, 4)):
    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=figsize)
    df[col1].value_counts(dropna=False).rename(col1).plot(kind='bar', ax=axes[0], subplots=True)
    df[col2].value_counts(dropna=False).rename(col2).plot(kind='bar', ax=axes[1], subplots=True)
    axes[0].set_xlabel('')
    axes[1].set_xlabel('')


def plot_kaplan_meier(df, feature_col, event_col='OS', time_col='OS.time'):
    _ = deepcopy(df[[event_col, time_col, feature_col]].dropna())
    if pd.api.types.is_numeric_dtype(_[feature_col]):
        if _[feature_col].nunique() > 2:
            _[feature_col] = pd.qcut(_[feature_col], 2, labels=["low", "high"], duplicates='drop')
    _ = _.sort_values(feature_col)
    classes = _[feature_col].unique()
    for cls in classes:
        mask = _[feature_col] == cls
        time, survival_prob, conf_int = kaplan_meier_estimator(
            _[event_col][mask] == 1,
            _[time_col][mask],
            conf_type="log-log",
        )
        plt.step(time, survival_prob, where="post", label=f"{feature_col} = {cls}")
        plt.fill_between(time, conf_int[0], conf_int[1], alpha=0.25, step="post")
    plt.ylim(0, 1)
    plt.ylabel(r"est. probability of survival $\hat{S}(t)$")
    plt.xlabel("time $t$")
    plt.legend(loc="best")
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.15))
    return plt