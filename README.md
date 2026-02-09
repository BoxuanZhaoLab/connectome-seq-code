# Connectome-seq Code

## Overview
The connectome-seq workflow consists of **two main parts**:

### a. Single-cell mRNA Sequencing Analysis

#### i. Standardized Workflow
1. Dataset integration  
2. Clustering and quality control  
3. Cell type annotation (automatic and manual)

#### ii. Output
- Cell type information

---

### b. Barcode Matching Analysis

#### i. Customized Workflow

1. **Cell Ranger processing**  
   Raw sequencing files (`.fastq`) are transformed into BAM files.  
   - mRNA features and Pre/PostRNA features come from separate FASTQ files.  
   - UMIs that differ within **Hamming distance = 1** are automatically collapsed by Cell Ranger.

2. **Custom Python collapsing of Pre/PostRNAs**  
   Pre/PostRNAs with the same cell barcode (CB) that differ by **Hamming distance = 1** are further collapsed.

3. **(Configurable) Disambiguation of Pre/PostRNA assignment**  
   If a Pre/PostRNA is linked to more than two CBs, it is assigned to the CB with the highest Pre/PostRNA UMI count after passing quality checks.  
   
   Current criteria:
   - CB #1 PreRNA > 5 × CB #2 PreRNA  
   - CB #1 PostRNA > 2 × CB #2 PostRNA  
   - CBs ranked by Pre/PostRNA UMI count

4. **Matching between nuclei and synaptosome Pre/PostRNAs**  
   Matching is performed at a specified Hamming distance threshold.

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

Demonstrating with cseq25 dataset.

### Step 1. `cellranger count`
Run one-line commands to quantify features from raw sequencing data.

```bash
cellranger count   --id=cseq25-PN-new   --transcriptome=/home/boxuan/data/hdd2/mouse_2022   --libraries=library_PN.csv   --feature-ref=feature_ref_2.csv   --force-cells=21000   --localmem 128   --localcores=32   --check-library-compatibility=false
```

```bash
cellranger count   --id=cseq25-CN2   --transcriptome=/home/boxuan/data/hdd2/mouse_cseq6   --libraries=library_CN.csv   --feature-ref=feature_ref_2.csv   --force-cells=4000   --localmem 128   --localcores=32
```

```bash
cellranger count   --id=cseq25-CS-pre2   --transcriptome=/home/boxuan/data/hdd2/mouse_cseq6   --libraries=library_CSpre2.csv   --feature-ref=feature_ref_pre.csv   --expect-cells=8000   --localmem 128   --localcores=32
```

```bash
cellranger count   --id=cseq25-CS-post2   --transcriptome=/home/boxuan/data/hdd2/mouse_cseq6   --libraries=library_CSpost2.csv   --feature-ref=feature_ref_post.csv   --expect-cells=8000   --localmem 128   --localcores=32
```

#### Instructions

1. Change `id`, `libraries`, and `feature-ref` as needed.
2. Use `--force-cells` for **nucleus libraries**.
3. Use `--expect-cells` for **synaptosome libraries**.

---

### Step 2. `cellranger subset-bam_linux`

```bash
gunzip barcodes.tsv.gz
```

```bash
~/cellranger/subset-bam_linux   --bam ../../cseq1-PN2/outs/possorted_genome_bam.bam   --cell-barcodes barcodes.tsv   --cores 32   --out-bam pons.pre.barcodes.bam
```

#### Instructions

1. Navigate to the `outs/` directory of each `cellranger count` result.
2. Unzip barcode files:
   ```bash
   gunzip barcodes.tsv.gz
   ```
3. Create a new folder for each sample.
4. Copy:
   - **Filtered** `barcodes.tsv` from nucleus samples
   - **Unfiltered (raw)** `barcodes.tsv` from synaptosome samples

---

### Step 3. `connectome.py`

```bash
../connectome.py pons.pre.barcodes.bam pons.pre.barcodes.hamm1.txt 1
```

#### Instructions

1. Extract reads and collapse barcodes at a given Hamming distance.

---

### Step 4.1. `cseq_synaptosome.Rmd`

#### Instructions

1. Synaptosome barcode information (TXT format) is processed with **`cseq_synaptosome.Rmd`**, where extensive quality control steps are taken. Output:  
   - `synp.post.clean.csv`  
   - `synp.pre.clean.csv`  
   These data serve as input for **`cseq_Match.Rmd`**.

### Step 4.2. `cseq_NucUMIClean.ipynb`

#### Instructions

1. Cell barcode information (e.g., `cere.post.BC.hamm1.txt`, `pons.pre.BC.hamm1.txt`) are UMI collapsed with **`cseq_NucUMIClean.ipynb`**, producing:  
   - `CN.umi.csv`  
   - `PN.umi.csv`  
   These serve as input for **`cseq_Match.Rmd`**, where cell barcodes undergo disambiguous assignment and cleaned barcode information is generated.

---

### Step 5. `cseq_Match.Rmd`

#### Instructions

1. Follow code annotation to run this RMarkdown.

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

