# celltypist

# !pip install celltypist

#cd '/data/objects'

python

import scanpy as sc
import celltypist
from celltypist import models

# load dataset

adata = sc.read('Tcells.h5ad')

adata.shape
adata.obs

# Download the latest CellTypist models.
# Enabling `force_update = True` will overwrite existing (old) models.
#models.download_models(force_update = True)
#models.models_path
#models.models_description()

# choose model
model = models.Model.load(model = 'Immune_All_Low.pkl')
model
model.cell_types

# transfer labels
# Not run; predict cell identities using this loaded model.
#predictions = celltypist.annotate(adata_2000, model = model, majority_voting = True)
# Alternatively, just specify the model name (recommended as this ensures the model is intact every time it is loaded).

predictions = celltypist.annotate(adata, model = 'Immune_All_Low.pkl', majority_voting = True)

predictions.predicted_labels

# Get an `AnnData` with predicted labels embedded into the cell metadata columns.
adata = predictions.to_adata()
adata.obs

# If the UMAP or any cell embeddings are already available in the `AnnData`, skip this command.
# sc.tl.umap(adata)

# visualize
sc.pl.umap(adata, color = ['cell_type', 'predicted_labels', 'majority_voting'], legend_loc = 'on data')

# compare predictions with preliminary labels
celltypist.dotplot(predictions, use_as_reference = 'cell_type', use_as_prediction = 'majority_voting')


# saving adata is not working so just keeping the metadata

import pandas as pd
pd.DataFrame(adata.obs).to_csv("metadata.Tcells.celltypist.csv")

# this saves adata, but not opening correctly in R
# required this to fix an index error before saving
adata.__dict__['_raw'].__dict__['_var'] = adata.__dict__['_raw'].__dict__['_var'].rename(columns={'_index': 'features'})
#del(adata.var['_index'])

adata.write_h5ad('Tcells_celltypist.h5ad')

exit()


