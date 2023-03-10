# HMM Experiments [![DOI](https://zenodo.org/badge/260986338.svg)](https://zenodo.org/badge/latestdoi/260986338)


This repository holds data, experimental setups, and code closely related to the paper ___"Exploiting Relations, Sojourn-Times and Joint Conditional Probabilities for Automated Commit Classification"___.

While the name implies hidden Markov Models, those were just one type of model tested here.
We attempt also to fit _dependent mixture models_, _joint conditional density models_, as well as attempting traditional machine learning.
We [**reverse-engineer**](https://github.com/MrShoenel/hmm-paper-2020-R-experiments/blob/5cdc1e5adfcae7e98d9d7ed5d57e53aa00d822e6/reverse-engineered-rules.md) rules for manually labeling commits and create labels for a few hundred commits.
Then, we produce datasets of **adjacent**, labeled commits that can be used with the mentioned models:

- [**`commits_t-0.csv`**](data/commits_t-0.csv): Approx. 300 newly labeled commits. Some of these were contained previously in Levin's dataset and we labeled them here again to see if we would reach the same consensus. All commits have size properties from `Git-Density` attached and come also with the usual information, like author, committer, email, timestamps, messages, hashes, etc.
- [**`commits_t-1.csv`**](data/commits_t-1.csv), [**`commits_t-2.csv`**](data/commits_t-2.csv), [**`commits_t-3.csv`**](data/commits_t-3.csv): Those are the "interesting" datasets. These have the same feature names as `commits_t-0.csv`, but also come with 1/2/3 directly predecessing commits, so they 2/3/4 times the features as `commits_t-0.csv`. The names are the same, but suffixed by `_t_{1,2,3}`.


The latter type of **commit chains** can be exploited for sequential learning.
For example, is there any value in knowing the *activity* that was carried out in the previous commit(s)?
