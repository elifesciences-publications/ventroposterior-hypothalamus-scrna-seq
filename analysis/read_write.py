#!/usr/bin/env python
# -*- coding: utf-8 -*-

import scanpy as sc
import numpy as np
import pandas as pd

from datetime import datetime

from config import *


def save_for_harmony(adata, output_dir, covariates=["sampleid"], thetas=[2], suffix=""):
    harmony_input_dir = os.path.join(output_dir, f"harmony_input{suffix}")
    harmony_output_dir = os.path.join(output_dir, f"harmony_output{suffix}")    
    if not os.path.exists(harmony_input_dir):
        os.makedirs(harmony_input_dir)
        print(f"Created input dir: {harmony_input_dir}")
    if not os.path.exists(harmony_output_dir):
        os.makedirs(harmony_output_dir)
        print(f"Created output dir: {harmony_output_dir}")

    embeddings_file = os.path.join(harmony_input_dir, "pc_embeddings.csv")
    metadata_file = os.path.join(harmony_input_dir, "metadata.csv")
    covariates_file = os.path.join(harmony_input_dir, "covariates.csv")
    thetas_file = os.path.join(harmony_input_dir, "thetas.csv")

    assert len(covariates) == len(thetas)
    assert pd.Index(covariates).isin(adata.obs_keys()).all()

    np.savetxt(embeddings_file, adata.obsm["X_pca"], delimiter=",")
    adata.obs[covariates].to_csv(metadata_file)
    with open(covariates_file, "w") as fout1:
        fout1.write(",".join(covariates) + "\n")
    with open(thetas_file, "w") as fout2:
        fout2.write(",".join(map(str, thetas)) + "\n")

def load_from_harmony(adata, output_dir, suffix=""):
    harmony_output_dir = os.path.join(output_dir, f"harmony_output{suffix}")
    harmony_pca = pd.read_csv(os.path.join(harmony_output_dir, "harmony_pc_embeddings.csv"), header=0, index_col=0)
    adata.obsm["X_pca"] = harmony_pca.values
    adata.obsm["X_pca_harmony"] = harmony_pca.values
    return adata

def timestamp():
    return datetime.now().strftime("%Y-%m-%dT%H-%M-%S")

def save_adata(adata, suffix="", subdir=""):
    filename = f"{adata.uns['sampleid']}{'-' + suffix if suffix else ''}-{timestamp()}.h5ad"
    sc.write(Path(adata.uns["output_dir"]) / subdir / filename, adata)

def load_sample(record):
    possible_input_paths = (DATA_DIR / "processed" / record.library_id).glob("*/filtered*.h5")
    input_path = list(filter(lambda p: p.exists(), possible_input_paths))[0]
    input_dir = input_path.parent

    adata = sc.read_10x_h5(input_path, record.genome)
    adata.var_names_make_unique()

    adata.obs["sequencing_saturation"] = np.nan
    seqsat_file = input_dir / "sequencing_saturation.csv"
    if seqsat_file.exists():
        seqsat = pd.read_csv(seqsat_file, index_col=0, header=0)
        adata.obs.loc[seqsat.index, "sequencing_saturation"] = seqsat["saturation"].values

    for key, value in record.items():
        adata.uns[key] = value

    adata.uns["analyst"] = ANALYST
    adata.uns["analyst_email"] = ANALYST_EMAIL

    adata.uns["datetime_created"] = timestamp()

    output_dir = DATA_DIR / "h5ad" / record.library_id
    adata.uns["input_file"] = input_path
    adata.uns["input_dir"] = input_dir
    adata.uns["output_dir"] = output_dir

    hemo_genes = pd.read_csv(DATA_DIR / "databases" / "hemoglobin_genes.GRCh38.csv", header=None, index_col=1).index
    mito_genes = pd.read_csv(DATA_DIR / "databases" / "mito_genes.GRCh38.csv", header=None, index_col=1).index
    ribo_genes = adata.var_names.str.startswith("RPS") |\
                 adata.var_names.str.startswith("RPL") |\
                 adata.var_names.str.startswith("MRPS") |\
                 adata.var_names.str.startswith("MRPL")
    # adapted from https://academic.oup.com/jmcb/article/11/8/703/5188008#164389921
    # exclude genes with annotations "Other", "Function known but...", and "Uncharacterized"
    cc_genes = pd.read_csv(DATA_DIR / "databases" / "cell_cycle.csv", header=0, index_col=1).index[:484]
    sr_genes = pd.read_csv(DATA_DIR / "databases" / "coregene_stress-response.csv", header=0, index_col=6).index

    adata.var["mitochondrial"] = adata.var_names.isin(mito_genes)
    adata.var["hemoglobin"] = adata.var_names.isin(hemo_genes)
    adata.var["ribosomal"] = ribo_genes
    adata.var["cell_cycle"] = adata.var_names.isin(cc_genes)
    adata.var["stress_response"] = adata.var_names.isin(sr_genes)

    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=("mitochondrial", "hemoglobin", "ribosomal", "cell_cycle", "stress_response"),
        inplace=True
    )

    if HAS_SCRUBLET:
        scrub = scr.Scrublet(adata.X)
        doublet_scores, predicted_doublets = scrub.scrub_doublets()
        adata.obs["scrublet_putative_doublet"] = predicted_doublets
        adata.obs["scrublet_score"] = doublet_scores

    save_adata(adata, suffix="raw")

    return record.library_id, adata
    pass

def load_samplesheet():
    assert SAMPLESHEET.exists()

    library_list = pd.read_csv(SAMPLESHEET, index_col=None, header=0)
    assert len(library_list) == N_SAMPLES

    return library_list
