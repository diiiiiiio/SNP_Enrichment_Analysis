# Patient-Specific Enrichment Analysis Workflow

This README outlines the entire process of our patient-specific analysis, from data organization to enrichment analysis, filtering, clustering, and final result output.

---

## 1. Data Organization

- **Dataset Partitioning:**  
  The patient SNP data are organized into three main folders: **train**, **test**, and **validation**.

- **File Format:**  
  Each file (e.g., `train_snps.xlsx`) contains SNP data for multiple patients. For each patient, global and local SNP selections have already been merged.

---

## 2. Global and Local SNP Selection

- **Global Selection:**  
  Includes 106 globally important SNPs identified by ProtoGate.

- **Local Selection:**  
  Each patient has their own set of locally selected SNPs based on their unique profile.

- **Merging:**  
  For each patient, the global and local SNP lists are merged to form the final SNP input for further analysis.

---

## 3. SNP-to-Gene Mapping

- **Primary Method:**  
  Use the Ensembl GRCh37 VEP API (based on hg19) to map each SNP to its corresponding gene.
  - If the API returns an error (e.g., a 400 error), a retry mechanism is employed.

- **Fallback:**  
  If Ensembl mapping fails after retries, the script automatically queries the MyVariant.info API as a backup.

- **Caching:**  
  Mapping results are cached to avoid repeated queries for the same SNP, thereby improving efficiency.

---

## 4. Enrichment Analysis

- **Tool:**  
  g:Profiler is used to perform enrichment analysis on the gene list.

- **Enrichment Sources:**  
  The following ontology sources are used:
  - **KEGG Pathway**
  - **GO Biological Processes (GO:BP)**
  - **Reactome Gene Sets (REAC)**
  - **Canonical Pathways (CP)**
  - **CORUM**
  - **WikiPathways (WP)**

- **Background Setting:**  
  The entire genome is used as the enrichment background (set to a background size of 20,000).

---

## 5. Filtering of Enrichment Results

- **Filtering Criteria:**  
  - p-value < 0.01  
  - At least 3 observed genes (intersection_size ≥ 3)  
  - Enrichment Factor > 1.5  
    - **Calculation:**  
      - Expected Count = (term_size / background_size) * query_size  
      - Enrichment Factor = intersection_size / Expected Count

---

## 6. Clustering and Redundancy Reduction

- **Objective:**  
  To reduce redundancy by clustering similar enrichment terms and selecting a representative term from each cluster.

- **Method:**
  - **Similarity Calculation:**  
    For each enriched term, convert the "intersection" (comma-separated list of query genes) into a gene set.  
    Calculate the Kappa score between every pair of terms based on these gene sets, with the total number of query genes as the background.
    
  - **Hierarchical Clustering:**  
    Construct a distance matrix using (1 – Kappa) and perform average-linkage hierarchical clustering.
    
  - **Cluster Definition:**  
    Use a distance threshold of 0.7 (corresponding to Kappa > 0.3) to group similar terms.
    
  - **Representative Selection:**  
    Within each cluster, select the term with the smallest p-value as the representative.

---

## 7. Result Output

- **Per-Patient Results:**  
  For each patient, an output folder is created within the respective partition (train/test/validation). The folder contains:
  - A CSV file with the clustered enrichment results.
  - A bar chart (PNG format) visualizing the top enriched terms.

- **Aggregated Summary:**  
  A master CSV file (`aggregated_summary.csv`) is generated to summarize key information across all patients, facilitating cross-patient comparisons.

---

## 8. Dependencies and Execution

- **Dependencies:**  
  - Python 3.x  
  - pandas  
  - numpy  
  - matplotlib  
  - gprofiler-official  
  - requests  
  - scipy

- **Installation:**  
  Install the required packages via:
  ```bash
  pip install pandas numpy matplotlib gprofiler-official requests scipy
