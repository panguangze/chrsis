#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
    tenxtools.clustering
    ~~~~~~~~~~~~~~~~~~~~



    @Copyright: (c) 2018-09 by Lingxi Chen (chanlingxi@gmail.com).
    @License: LICENSE_NAME, see LICENSE for more details.
"""

from sklearn import decomposition
import seaborn as sns
import matplotlib.pyplot as plt
sns.set(color_codes=True)
plt.switch_backend('agg')


def pca(X, y=None, n_components=2):
    pca = decomposition.PCA(n_components=n_components)
    pca.fit(X)
    X = pca.transform(X)
    return X


def draw_hierarchical_clustering(out_prefix, X, y=None,
                                 method='average', metric='euclidean',
                                 standard_scale=None, z_score=None):
    if y is not None:
        lut = dict(zip(y.unique(), "rb"))
        row_colors = y.map(lut)
        row_colors.index = X.index
        sns.clustermap(X, row_colors=row_colors,
                       method=method, metric=metric,
                       standard_scale=standard_scale, z_score=z_score,
                       )
    else:
        sns.clustermap(X,
                       method=method, metric=metric,
                       standard_scale=standard_scale, z_score=z_score,
                       )

    plt.suptitle('Hierarchical Clustering:' + method + '.' + metric)
    plt.savefig(out_prefix + '_' + method + '_' + metric + '.png', format='png')


def draw_PC2(out_prefix, X, y):
    plt.figure(figsize=(8, 6))
    tumor = X[0:8]
    normal = X[8:16]
    plt.plot(tumor['PC1'], tumor['PC2'], '*r', label='Tumor')
    plt.plot(normal['PC1'], normal['PC2'], '*b', label='Normal')

    plt.legend(loc='upper right')
    plt.xlabel('PC1')
    plt.ylabel('PC2')
    plt.suptitle('PC2 Cluster')
    plt.savefig(out_prefix + '_PC2.png', format='png')
