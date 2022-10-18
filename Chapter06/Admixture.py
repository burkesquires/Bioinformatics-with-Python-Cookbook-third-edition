# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.13.3
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# +
from collections import defaultdict
import os

import matplotlib.pyplot as plt

from genomics.popgen.admix import cluster, plot

# %matplotlib notebook
# -

k_range = range(2, 10)  # 2..9

with open('hapmap10_auto_noofs_ld.fam') as f:
    ind_order = []
    for l in f:
        toks = l.rstrip().replace(' ', '\t').split('\t')
        fam_id = toks[0]
        ind_id = toks[1]
        ind_order.append((fam_id, ind_id))
# ## CV-plot

CVs = []
for k in k_range:
    with open('admix.%d' % k) as f:
        for l in f:
            if l.find('CV error') > -1:
                CVs.append(float(l.rstrip().split(' ')[-1]))
                break
fig = plt.figure(figsize=(16, 9))
ax = fig.add_subplot(111)
ax.plot(k_range, CVs)
ax.set_title('Cross-Validation error')
ax.set_xlabel('K')

with open('relationships_w_pops_121708.txt') as f:
    pop_ind = defaultdict(list)
    f.readline()  # header
    for l in f:
        toks = l.rstrip().split('\t')
        fam_id = toks[0]
        ind_id = toks[1]
        if (fam_id, ind_id) not in ind_order:
            continue
        mom = toks[2]
        dad = toks[3]
        if mom != '0' or dad != '0':
            continue
        pop = toks[-1]
        pop_ind[pop].append((fam_id, ind_id))


def load_Q(fname, ind_order):
    ind_comps = {}
    with open(fname) as f:
        for i, l in enumerate(f):
            comps = [float(x) for x in l.rstrip().split(' ')]
            ind_comps[ind_order[i]] = comps
    return ind_comps


comps = {
    k: load_Q('hapmap10_auto_noofs_ld.%d.Q' % k, ind_order) for k in k_range
}

ordering = {k: cluster(comps[k], pop_ind) for k in k_range}
fig = plt.figure(figsize=(9, 9))
plot.single(comps[4], ordering[4], fig)
None

fig = plt.figure(figsize=(16, 9))
plot.stacked(comps, ordering[7], fig)

# ## Q files?

# ## Log-likelihood


