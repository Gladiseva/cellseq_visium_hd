import scanpy as sc


def perform_pca(adata, n_components=50):
    """Perform PCA on the AnnData object."""
    sc.pp.pca(adata, n_comps=n_components)
    return adata

def build_neighbors_graph(adata, n_neighbors=15, n_pcs=40):
    """Construct a neighborhood graph using PCA components."""
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)
    return adata


def perform_leiden_clustering(adata, resolution=0.5, key_added='clusters'):
    """Apply Leiden clustering on the neighborhood graph."""
    sc.tl.leiden(adata, resolution=resolution, key_added=key_added)
    return adata
