import os
import glob
import time
import requests
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from gprofiler import GProfiler
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform

# ========== Parameter Settings ==========
base_dir = '/content/gdrive/MyDrive/bio'  # Path to the root directory containing data
splits = ['train', 'test', 'validation']  # Dataset splits
background_size = 20000  # Whole genome background size for enrichment analysis

# Sources used for enrichment analysis
enrichment_sources = ["KEGG", "GO:BP", "REAC", "CP", "CORUM", "WP"]

# ========== Initialize g:Profiler ==========
# The parameter no_evidences=False ensures the results include matched gene lists in the 'intersection' column
gp = GProfiler(return_dataframe=True)

# ========== Define Backup SNP Online Mapping Function (MyVariant.info) ==========
def fetch_gene_from_snp_fallback(rsid, retries=3):
    """Queries MyVariant.info API for gene symbols associated with a given SNP (rsID)."""
    url = f"https://myvariant.info/v1/variant/{rsid}?fields=gene.symbol"
    attempt = 0
    while attempt < retries:
        try:
            response = requests.get(url, timeout=10)
            if response.status_code != 200:
                attempt += 1
                time.sleep(2)
                continue
            data = response.json()
            genes = set()
            if "gene" in data:
                gene_info = data["gene"]
                if isinstance(gene_info, dict) and "symbol" in gene_info:
                    genes.add(gene_info["symbol"])
                elif isinstance(gene_info, list):
                    for item in gene_info:
                        if "symbol" in item:
                            genes.add(item["symbol"])
            return list(genes)
        except Exception:
            attempt += 1
            time.sleep(2)
    return []

# ========== Define Main SNP Online Mapping Function (Using Ensembl GRCh37) ==========
def fetch_gene_from_snp(rsid, retries=3):
    """Queries Ensembl GRCh37 API to retrieve gene symbols for a given SNP (rsID)."""
    url = f"https://grch37.rest.ensembl.org/vep/human/id/{rsid}?content-type=application/json"
    headers = {"Content-Type": "application/json"}
    attempt = 0
    while attempt < retries:
        try:
            response = requests.get(url, headers=headers, timeout=10)
            if response.status_code != 200:
                attempt += 1
                time.sleep(2)
                continue
            data = response.json()
            if not data:
                return []
            genes = set()
            for entry in data:
                if "transcript_consequences" in entry:
                    for t in entry["transcript_consequences"]:
                        if "gene_symbol" in t:
                            genes.add(t["gene_symbol"])
                        elif "gene_id" in t:
                            genes.add(t["gene_id"])
            if genes:
                return list(genes)
            else:
                break
        except Exception:
            attempt += 1
            time.sleep(2)
    return fetch_gene_from_snp_fallback(rsid, retries)

# Caching SNP-to-gene mapping results to avoid redundant queries
snp_gene_cache = {}
def get_genes_from_snps(snps):
    """Retrieves associated genes for a list of SNPs, utilizing caching to minimize API calls."""
    genes = set()
    for snp in snps:
        if snp in snp_gene_cache:
            result = snp_gene_cache[snp]
        else:
            result = fetch_gene_from_snp(snp)
            snp_gene_cache[snp] = result
        if result:
            genes.update(result)
    return list(genes)

# ========== Define Plotting Function ==========
def plot_enrichment_results(df, title, filename):
    """Generates a bar plot of the top enrichment terms based on p-values."""
    if df.empty:
        return
    term_col = 'term_name' if 'term_name' in df.columns else ('name' if 'name' in df.columns else None)
    if term_col is None:
        return
    df_top = df.sort_values('p_value').head(10)
    plt.figure(figsize=(10, 6))
    plt.barh(df_top[term_col], -np.log10(df_top['p_value']))
    plt.xlabel('-log10(p_value)')
    plt.title(title)
    plt.gca().invert_yaxis()
    plt.tight_layout()
    plt.savefig(filename)
    plt.close()

# ========== Define Clustering Function: Compute Kappa Score Based on Intersection ==========
def calculate_kappa(setA, setB, total):
    """Calculates the kappa score based on the overlap of two gene sets."""
    n11 = len(setA & setB)
    n10 = len(setA) - n11
    n01 = len(setB) - n11
    n00 = total - len(setA | setB)
    Po = (n11 + n00) / total
    Pe = (len(setA)/total * len(setB)/total) + ((total - len(setA))/total * (total - len(setB))/total)
    return 0 if (1 - Pe) == 0 else (Po - Pe) / (1 - Pe)

# ========== Define Helper Function: Parse Intersection Column ==========
def parse_intersection(x):
    """Parses the 'intersection' column into a set of genes."""
    if pd.isnull(x):
        return set()
    if isinstance(x, str):
        return set(x.split(','))
    if isinstance(x, list):
        return set(x)
    return set()

# ========== Summary List ==========
summary_list = []
