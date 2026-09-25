import time, tiledbsoma, cellxgene_census

uri = "/mnt/bulk/dl_bio/cellxgene_census/2025-11-08/soma"
t0 = time.monotonic()
census = tiledbsoma.Collection.open(uri)
print("open:", time.monotonic() - t0, "s")

t0 = time.monotonic()
adata = cellxgene_census.get_anndata(
    census,
    organism="Homo sapiens",
    var_value_filter="feature_name in ['ADORA1','ADORA2A','ADORA2B','ADORA3']",
    obs_value_filter="tissue_general == 'brain' and is_primary_data == True",
)
print("brain ADORA fetch:", time.monotonic() - t0, "s")
print(adata)
