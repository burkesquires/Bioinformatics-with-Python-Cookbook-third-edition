# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.14.0
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# + jupyter={"outputs_hidden": false}
import os

from sklearn.decomposition import PCA
import numpy as np

from genomics.popgen.pca import plot
with open('../Chapter06/relationships_w_pops_041510.txt') as f:
    ind_pop = {}
    f.readline()  # header
    for l in f:
        toks = l.rstrip().split('\t')
        fam_id = toks[0]
        ind_id = toks[1]
        pop = toks[-1]
        ind_pop['/'.join([fam_id, ind_id])] = pop
with open('../Chapter06/hapmap10_auto_noofs_ld_12.ped') as f:
    ninds = 0
    ind_order = []
    for line in f:
        ninds += 1
        toks = line[:100].replace(' ', '\t').split('\t') #  for speed
        fam_id = toks[0]
        ind_id = toks[1]
        ind_order.append(f'{fam_id}/{ind_id}')
    nsnps = (len(line.replace(' ', '\t').split('\t')) - 6) // 2
# + jupyter={"outputs_hidden": false}
pca_array = np.empty((ninds, nsnps), dtype=int)
print(pca_array.shape)
with open('../Chapter06/hapmap10_auto_noofs_ld_12.ped') as f:
    for ind, line in enumerate(f):
        snps = line.replace(' ', '\t').split('\t')[6:]
        for pos in range(len(snps) // 2):
            a1 = int(snps[2 * pos])
            a2 = int(snps[2 * pos])
            my_code = a1 + a2 - 2
            pca_array[ind, pos] = my_code
# + jupyter={"outputs_hidden": false}
my_pca = PCA(n_components=8)
my_pca.fit(pca_array)
trans = my_pca.transform(pca_array)
#Memory required

# + jupyter={"outputs_hidden": false}
sc_ind_comp = {ind_order[i]: ind_pca for i, ind_pca in enumerate(trans)}
plot.render_pca_eight(sc_ind_comp, cluster=ind_pop)

# + jupyter={"outputs_hidden": false}

