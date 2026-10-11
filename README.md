# INSL3 fetal-brain strict sex-state screen

This public, reproducible workflow tests a deliberately strict descriptive criterion in human fetal-brain transcriptomics:

- `min(XY) > max(XX)`, or
- `max(XY) < min(XX)`

The unit of replication is the **donor**, not the individual RNA-seq sample or nucleus.

## Stage 1
Downloads public `E-MTAB-4840` data from EMBL-EBI Expression Atlas, excludes sex-chromosome genes, normalizes RNA-seq counts with edgeR TMM/log-CPM, and screens region × developmental-window donor distributions.

It also looks for **temporal emergence**: sex distributions overlap in an earlier window, then become completely separated later.

## Stage 2
If Stage 1 finds temporal-emergence candidates, downloads the public Velmeshev et al. 2023 CELLxGENE H5AD and asks whether the same genes show complete donor-level separation within a fetal cortical cell type.

## Interpretation
A strict survivor is a candidate fetal neural molecular state. It does **not** establish INSL3 causality. Causality requires an INSL3/RXFP2 perturbation experiment.

The workflow also reports the null probability of complete separation given the donor counts. With very small cohorts, complete separation can occur surprisingly often by chance, so independent replication is mandatory.
