# connectome-seq-code
Code usage:
1. BC_collapse: Processing raw barcode through collapsing cell barcodes.
2. cseq_NucUMIClean: Collapsing nucleus barcode information base on UMI.
3. cseq_synaptosome: Preprocessing synaptosomal barcode, apply filters to extract cleaned data.
4. cseq_HarmonyIntegration: Integrate different connectome-seq datasets to eliminate batch effect.
5. cseq_Annotation: Annotate the integrated connectome-seq dataset with Allen Institue mouse brain cell type atlas.
6. cseq_Cluster: Cluster the integrated dataset and annotate each cluster from marker gene expression.
7. cseq_Match: Conduct barcode matching base on synapse barcode and cell barcode.
8. cseq_PurkinjeAnalysis: Analyze possible Pons-Purkinje cell connection.
9. cseq_Rfunction: Wrapper for R functions used in cseq_PurkinjeAnalysis.
10. cseq_virusLib: Code for analyzing virus libraray sequencing results, extracting barcode information.
11. cseq_preseq: Code for predicting virus library barcode diversity.


The core function of our code is to generate single-cell level connectome information based on gene and barcode information from single-nucleus and single-synaptosome sequencing.

Barcode processing:
1. Raw single-nuclues and single-synaptosome sequencing data were processed with cellranger to extract barcode information. Barcode data were further transfered into bam files which serves as input for BC_collapse.py.
2. Through running BC_collaspe.py on barcode bam files, a text file ("synp.post.BC.table.hamming1.txt", "synp.pre.BC.table.hamming1.txt") containing barcode sequence, UMI, CB and read count information was generated. BC_collapse.py can also collapse barcode at given hamming distance.
3. Synaptosome barcode information (in txt format) are than processed with cseq_synaptosome.Rmd, where extensive quality control steps were taken. Output files contain cleaned barcode information (in csv format, "synp.post.clean.csv", "synp.pre.clean.csv"). These data serve as one of the input data for cseq_Match.Rmd.
4. Cell barcode information ("cere.post.BC.hamm1.txt", "pons.pre.BC.hamm1.txt") were UMI collapsed with cseq_NucUMIClean.ipynb ("CN.umi.csv", "PN.umi.csv") and serves as one of the input data for cseq_Match.Rmd, where cell barcode will undergo disambiguous assignment. Cleaned barcode information were than generated.

Transcriptome processing:
1. Raw single-nucleus sequencing data were processed with cellranger to extract feature_bc_matrix ("/filtered_feature_bc_matrix": "barcodes.tsv.gz", "features.tsv.gz", "matrix.mtx.gz") and then used as input for cseq_HarmonyIntegration.Rmd, where seurat object of snRNA-seq data will be created and data from different biological replicates will be integrated with Harmony (output seurat object: "cn.af", "pn.af").
2. Integrated data will then go through cell type prediction by MapMyCell in cseq_Annotation.Rmd (output seurat object: "scdata.pn", "scdata.cn") and extensive quality control, dimension reduction and annotation in cseq_Cluster.Rmd (output seurat object: "celltype.pons.new", "celltype.cere.new", "pn_barcodes", "cn_barcodes") to recover cell type information. The resulting seurat object, serves as input for cseq_Match.Rmd, which provide cell type information for different CB and corresponding barcode .

Connectome processing:
1. Output data mentioned above all serves as input for cseq_Match.Rmd to generate final connectome, which generate final connectome (output: "PN.neuron.umiClean.csv", "CN.neuron.umiClean.csv", "synp.pre.clean.csv", "synp.post.clean.csv", "connectome.match.count.csv").

Purkinje cell analysis:
1. Purkinje cell were extracted from annotated snRNA-seq data and processed with cseq_PurkinjeAnalysis.Rmd to analyze potential connectivity-related gene markers.

Virus library analysis:
1. Raw sequencing data of Barcode AAV genomic DNA were served as input for cseq_VirusLib.ipynb. Barcode sequence and count information were extracted and calculated distribution. cseq_preseq.sh is then applied to project barcode sequence diversity at future sequencing depth.

