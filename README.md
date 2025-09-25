# connectome-seq-code
Code usage:
1. BC_collapse: Processing raw barcode through collapsing cell barcodes.
2. cseq_synaptosome: Preprocessing synaptosomal barcode, apply filters to extract cleaned data.
3. cseq_HarmonyIntegration: Integrate different connectome-seq datasets to eliminate batch effect.
4. cseq_Annotation: Annotate the integrated connectome-seq dataset with Allen Institue mouse brain cell type atlas.
5. cseq_Cluster: Cluster the integrated dataset and annotate each cluster from marker gene expression.
6. cseq_Match: Conduct barcode matching base on synapse barcode and cell barcode.
7. cseq_PurkinjeAnalysis: Analyze possible Pons-Purkinje cell connection.
8. cseq_Rfunction: Wrapper for R functions used in cseq_PurkinjeAnalysis.
9. cseq_virusLib: Code for analyzing virus libraray sequencing results, extracting barcode information.
10. cseq_preseq: Code for predicting virus library barcode diversity.


Barcode processing:
1. Raw sequencing fastq files: raw single-nuclues and single-synaptosome sequencing data were processed with cellranger to extract both gene and barcode information. Barcode data were further transfered into bam files which serves as input for BC_collapse.py.
2. Through running BC_collaspe.py on barcode bam files, a text file containing barcode sequence, UMI, CB and read count information was generated. BC_collapse.py can also collapse barcode at given hamming distance.
3. Synaptosome barcode information (in txt format) are than processed with cseq_synaptosome.Rmd, where extensive quality control steps were taken. Output files contain cleaned barcode information (in csv format). These data serve as one of the input data for cseq_Match.
4. Cell barcode 
