import os, re, sys, gc, requests, math
from pathlib import Path
import numpy as np
import pandas as pd
import anndata as ad
from scipy import sparse

Path("results").mkdir(exist_ok=True)
temporal_path = Path("results/stage1_temporal_emergence.tsv")
if not temporal_path.exists() or temporal_path.stat().st_size == 0:
    print("No Stage-1 temporal result file; skipping Stage 2.")
    sys.exit(0)

try:
    temporal = pd.read_csv(temporal_path, sep="\t")
except Exception:
    print("Temporal file is empty/unreadable; skipping Stage 2.")
    sys.exit(0)

if temporal.empty or "gene_id" not in temporal:
    print("Zero Stage-1 temporal candidates; skipping Stage 2.")
    pd.DataFrame().to_csv("results/stage2_strict_celltypes.tsv", sep="\t", index=False)
    sys.exit(0)

candidates = sorted(set(temporal["gene_id"].dropna().astype(str)))
print("Stage-2 candidates:", len(candidates))

# Public CELLxGENE asset for Velmeshev et al. 2023, supplied by the dataset release.
url = os.environ.get(
    "VELMESHEV_H5AD_URL",
    "https://datasets.cellxgene.cziscience.com/b126d3e9-a1a2-419a-a049-b4a27d3ce3ab.h5ad"
)
path = Path("velmeshev_2023.h5ad")
if not path.exists():
    with requests.get(url, stream=True, timeout=180) as r:
        r.raise_for_status()
        with path.open("wb") as f:
            for chunk in r.iter_content(8*1024*1024):
                if chunk:
                    f.write(chunk)
print("Downloaded", path, "bytes", path.stat().st_size)

adata = ad.read_h5ad(path, backed="r")
obs = adata.obs.copy()
var = adata.var.copy()

def first_col(patterns):
    for p in patterns:
        hits = [c for c in obs.columns if re.search(p, str(c), re.I)]
        if hits: return hits[0]
    return None

sex_col = first_col([r"^sex$",r"karyotype",r"sex"])
age_col = first_col([r"development_stage",r"age",r"gest"])
donor_col = first_col([r"donor_id",r"donor",r"individual",r"subject"])
cell_col = first_col([r"cell_type",r"subclass",r"cell.*class",r"cluster"])
if None in [sex_col,age_col,donor_col,cell_col]:
    pd.DataFrame({"obs_columns":list(obs.columns)}).to_csv("results/stage2_metadata_columns.tsv",sep="\t",index=False)
    raise RuntimeError(f"Could not resolve metadata: sex={sex_col} age={age_col} donor={donor_col} cell={cell_col}")

def norm_sex(x):
    s=str(x).strip().lower().replace(" ","")
    if s in {"male","m","xy","46,xy","46xy"}: return "XY"
    if s in {"female","f","xx","46,xx","46xx"}: return "XX"
    return None

def parse_wpc(x):
    if pd.isna(x): return np.nan
    s=str(x).lower().replace("–","-").replace("—","-")
    if any(k in s for k in ["year","month","postnatal","adult"]): return np.nan
    nums=[float(z) for z in re.findall(r"\d+(?:\.\d+)?",s)]
    if not nums: return np.nan
    return float(np.mean(nums[:2])) if "-" in s and len(nums)>=2 else nums[0]

obs["sex_std"]=obs[sex_col].map(norm_sex)
obs["age_wpc"]=obs[age_col].map(parse_wpc)
obs["donor_std"]=obs[donor_col].astype(str)
obs["cell_std"]=obs[cell_col].astype(str)

# Find candidate genes by Ensembl ID in index or common var columns.
gene_pos={}
sources=[("index",var.index.astype(str).values)]
for c in ["feature_id","gene_ids","gene_id"]:
    if c in var.columns:
        sources.append((c,var[c].astype(str).values))
for _,vals in sources:
    for i,g in enumerate(vals):
        if g in candidates and g not in gene_pos:
            gene_pos[g]=i

missing=[g for g in candidates if g not in gene_pos]
pd.DataFrame({"missing_gene_id":missing}).to_csv("results/stage2_missing_candidates.tsv",sep="\t",index=False)
if not gene_pos:
    print("No candidates mapped into H5AD.")
    sys.exit(0)

# Candidate developmental interval. Use 11-20 wpc to capture the reported circulating INSL3 window.
mask=(obs["age_wpc"].between(11,20,inclusive="both") &
      obs["sex_std"].isin(["XX","XY"]) &
      obs["donor_std"].notna() &
      obs["cell_std"].notna())
idx=np.flatnonzero(mask.to_numpy())
o=obs.iloc[idx][["donor_std","sex_std","age_wpc","cell_std"]].copy()

# Require 25 nuclei per donor×celltype and >=3 donors per sex.
counts=o.groupby(["donor_std","sex_std","cell_std"],observed=True).size().rename("n").reset_index()
valid=counts[counts["n"]>=25]
valid_keys=set(map(tuple,valid[["donor_std","sex_std","cell_std"]].astype(str).to_numpy()))

gene_order=list(gene_pos)
pos=[gene_pos[g] for g in gene_order]
X=adata.X[idx,pos]
if sparse.issparse(X): X=X.toarray()
else: X=np.asarray(X)

# If X looks integer-like/raw, compute library size from full matrix by chunk.
xnorm=str(adata.uns.get("X_normalization","")).lower()
if not xnorm or xnorm in {"raw","none"}:
    lib=np.zeros(len(idx),float)
    chunk=10000
    for a in range(0,len(idx),chunk):
        sub_idx=idx[a:a+chunk]
        block=adata.X[sub_idx,:]
        lib[a:a+chunk]=np.asarray(block.sum(axis=1)).ravel() if sparse.issparse(block) else np.asarray(block).sum(axis=1)
    lib[lib<=0]=np.nan
    X=np.log2((X/lib[:,None])*1e6+1)
else:
    X=X.astype(float)

frame=o.reset_index(drop=True)
for j,g in enumerate(gene_order): frame[g]=X[:,j]
frame["_valid"]=[(str(d),str(s),str(c)) in valid_keys for d,s,c in zip(frame.donor_std,frame.sex_std,frame.cell_std)]
frame=frame[frame["_valid"]].drop(columns="_valid")

pb=frame.groupby(["donor_std","sex_std","cell_std"],observed=True)[gene_order].mean().reset_index()

def strict(xy,xx):
    xy=np.asarray(xy,float); xx=np.asarray(xx,float)
    xy=xy[np.isfinite(xy)]; xx=xx[np.isfinite(xx)]
    if len(xy)<3 or len(xx)<3: return None
    if xy.min()>xx.max(): return ("XY>XX",xy.min()-xx.max(),xy.min(),xy.max(),xx.min(),xx.max())
    if xy.max()<xx.min(): return ("XY<XX",xx.min()-xy.max(),xy.min(),xy.max(),xx.min(),xx.max())
    return None

rows=[]
for ct,sub in pb.groupby("cell_std",observed=True):
    for g in gene_order:
        xy=sub.loc[sub.sex_std=="XY",g].dropna().to_numpy()
        xx=sub.loc[sub.sex_std=="XX",g].dropna().to_numpy()
        st=strict(xy,xx)
        if st is None: continue
        rows.append(dict(gene_id=g,cell_type=ct,n_xy=len(xy),n_xx=len(xx),
                         direction=st[0],margin=st[1],xy_min=st[2],xy_max=st[3],
                         xx_min=st[4],xx_max=st[5],
                         null_complete_sep_prob=2/math.comb(len(xy)+len(xx),len(xy))))
out=pd.DataFrame(rows)
out.to_csv("results/stage2_strict_celltypes.tsv",sep="\t",index=False)
print("Stage 2 strict gene×celltype survivors:",len(out))
print(out.head(50).to_string(index=False) if len(out) else "ZERO")
