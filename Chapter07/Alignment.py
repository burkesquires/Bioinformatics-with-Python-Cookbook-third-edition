# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.13.6
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# +
import os

import dendropy
# -

# ## Genome alignment

from Bio.Align.Applications import MafftCommandline
mafft_cline = MafftCommandline(input='sample.fasta', ep=0.123, reorder=True, maxiterate=1000, localpair=True)
print(mafft_cline)
stdout, stderr = mafft_cline()
with open('align.fasta', 'w') as w:
    w.write(stdout)

os.system('trimal -automated1 -in align.fasta -out trim.fasta -fasta')


# ## Protein alignment

# +
from Bio.Align.Applications import MuscleCommandline

my_genes = ['NP', 'L', 'VP35', 'VP40']

for gene in my_genes:
    muscle_cline = MuscleCommandline(input=f'{gene}_P.fasta')
    print(muscle_cline)
    stdout, stderr = muscle_cline()
    with open(f'{gene}_P_align.fasta', 'w') as w:
        w.write(stdout)

# +
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
# XXX vvv
# from Bio.Alphabet import generic_protein

for gene in my_genes:
    unal_gene = SeqIO.parse(f'{gene}.fasta', 'fasta')
    gene_seqs = {rec.id: rec.seq for rec in unal_gene}
    al_prot = SeqIO.parse(f'{gene}_P_align.fasta', 'fasta')
    al_genes = []
    for protein in al_prot:
        my_id = protein.id
        seq = ''
        pos = 0
        for c in protein.seq:
            if c == '-':
                seq += '---'
            else:
                seq += str(gene_seqs[my_id][pos:pos + 3])
                pos += 3
        al_genes.append(SeqRecord(Seq(seq), id=my_id))


    SeqIO.write(al_genes, f'{gene}_align.fasta', 'fasta')
# -


