# Connectome-seq Code

## Overview
The core function of this codebase is to generate single-cell level connectome information based on gene and barcode information from single-nucleus and single-synaptosome sequencing.

---

## Code Usage
1. **`BC_collapse`**: Processing raw barcode through collapsing cell barcodes.  
2. **`cseq_NucUMIClean`**: Collapsing nucleus barcode information based on UMI.  
3. **`cseq_synaptosome`**: Preprocessing synaptosomal barcode, applying filters to extract cleaned data.  
4. **`cseq_HarmonyIntegration`**: Integrate different connectome-seq datasets to eliminate batch effect.  
5. **`cseq_Annotation`**: Annotate the integrated connectome-seq dataset with Allen Institute mouse brain cell type atlas.  
6. **`cseq_Cluster`**: Cluster the integrated dataset and annotate each cluster from marker gene expression.  
7. **`cseq_Match`**: Conduct barcode matching based on synapse barcode and cell barcode.  
8. **`cseq_PurkinjeAnalysis`**: Analyze possible Pons–Purkinje cell connection.  
9. **`cseq_Rfunction`**: Wrapper for R functions used in `cseq_PurkinjeAnalysis`.  
10. **`cseq_virusLib`**: Analyze virus library sequencing results, extracting barcode information.  
11. **`cseq_preseq`**: Predict virus library barcode diversity.

---

## Barcode Processing
1. Raw single-nucleus and single-synaptosome sequencing data were processed with CellRanger to extract barcode information. Barcode data were further transferred into BAM files, which serve as input for **`BC_collapse.py`**.  
2. Running **`BC_collapse.py`** on barcode BAM files generates text files such as:  
   - `synp.post.BC.table.hamming1.txt`  
   - `synp.pre.BC.table.hamming1.txt`  
   These contain barcode sequence, UMI, CB, and read count information. `BC_collapse.py` can also collapse barcodes at a given Hamming distance.  
3. Synaptosome barcode information (TXT format) is processed with **`cseq_synaptosome.Rmd`**, where extensive quality control steps are taken. Output:  
   - `synp.post.clean.csv`  
   - `synp.pre.clean.csv`  
   These data serve as input for **`cseq_Match.Rmd`**.  
4. Cell barcode information (e.g., `cere.post.BC.hamm1.txt`, `pons.pre.BC.hamm1.txt`) are UMI collapsed with **`cseq_NucUMIClean.ipynb`**, producing:  
   - `CN.umi.csv`  
   - `PN.umi.csv`  
   These serve as input for **`cseq_Match.Rmd`**, where cell barcodes undergo disambiguous assignment and cleaned barcode information is generated.

---

## Transcriptome Processing
1. Raw single-nucleus sequencing data were processed with CellRanger to extract feature_bc_matrix files:  
   - `barcodes.tsv.gz`  
   - `features.tsv.gz`  
   - `matrix.mtx.gz`  
   These are used as input for **`cseq_HarmonyIntegration.Rmd`**, where Seurat objects are created and biological replicates integrated with Harmony. Output Seurat objects:  
   - `cn.af`  
   - `pn.af`  
2. Integrated data go through cell type prediction by **MapMyCell** in **`cseq_Annotation.Rmd`**, generating:  
   - `scdata.pn`  
   - `scdata.cn`  
   Further QC, dimension reduction, and annotation are performed in **`cseq_Cluster.Rmd`**, producing:  
   - `celltype.pons.new`  
   - `celltype.cere.new`  
   - `pn_barcodes`  
   - `cn_barcodes`  
   These provide cell type information for different CBs and corresponding barcodes, serving as input for **`cseq_Match.Rmd`**.

---

## Connectome Processing
All outputs above serve as input for **`cseq_Match.Rmd`**, which generates the final connectome. Outputs include:  
- `PN.neuron.umiClean.csv`  
- `CN.neuron.umiClean.csv`  
- `synp.pre.clean.csv`  
- `synp.post.clean.csv`  
- `connectome.match.count.csv`

---

## Purkinje Cell Analysis
Purkinje cells are extracted from annotated snRNA-seq data and processed with **`cseq_PurkinjeAnalysis.Rmd`** to analyze potential connectivity-related gene markers.

---

## Virus Library Analysis
1. Raw sequencing data of Barcode AAV genomic DNA are processed with **`cseq_VirusLib.ipynb`**. Barcode sequence and count information are extracted, and distribution calculated.  
2. **`cseq_preseq.sh`** is applied to project barcode sequence diversity at future sequencing depth.

---

