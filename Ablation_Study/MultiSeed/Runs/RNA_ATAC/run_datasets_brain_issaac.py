import random
import os
import numpy as np
import torch
import scanpy as sc
import scvi
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed


sys.path.append("../../multiHIVE/src")
from multiHIVE.model import multiHIVE


def set_all_seeds(seed):
    random.seed(seed)
    os.environ['PYTHONHASHSEED'] = str(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed(seed)
        torch.cuda.manual_seed_all(seed)
        # Guarantees deterministic operations on GPU
        torch.backends.cudnn.deterministic = True 



def run_multiHIVE(seed):

    set_all_seeds(seed)

    adata = sc.read_h5ad("/Data/RNA_ATAC/Brain-ISSAAC/Brain-ISSAAC-seq.h5ad")
    del adata.obsm
    del adata.obsp


    adata = scvi.data.organize_multiome_anndatas(adata)
    adata = adata[:, adata.var["modality"].argsort()].copy()
    sc.pp.filter_genes(adata, min_cells=int(adata.shape[0] * 0.01))
    multiHIVE.setup_anndata(adata, batch_key="modality")

    vae = multiHIVE(adata, latent_distribution="normal",
                    n_genes=(adata.var["modality"] == "Gene Expression").sum(),
                    n_regions=(adata.var["modality"] == "Peaks").sum(),
                    n_proteins=0,
                    kl_dot_product  = True,
                    deep_network = False,
                )

    vae.train(max_epochs=500, early_stopping=False)
    vae.get_latent_representation()

    del adata.X
    del adata.obs
    os.makedirs("./outputs/Brain-ISSAAC", exist_ok=True)
    embeddings = adata.obsm['Z_multiHIVE']
    np.save(f"./outputs/Brain-ISSAAC/multiHIVE_seed_{seed}_embeddings.npy", embeddings)
    # adata.write(f"./outputs/Brain-ISSAAC/multiHIVE_seed_{seed}.h5ad")

if __name__ == "__main__":
    # for seed in [234,45346,123,567,5678]:
    #     run_multiHIVE(seed)

    with ProcessPoolExecutor(max_workers=5) as executor:
        futures = {executor.submit(run_multiHIVE, seed): seed for seed in [4055, 2120, 4519, 3868, 8785, 234, 45346, 123, 567, 5678]}
        for future in as_completed(futures):
            try:
                print(future.result())
                print(f"Completed seed: {futures[future]}")
            except Exception as e:
                print(f"Error occurred: {e}")