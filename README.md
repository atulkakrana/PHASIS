## *PHASIS* suite for phased clusters discovery, comparison, annotation and to identify miRNA triggers

**Please see *PHASIS* [wiki](https://github.com/atulkakrana/PHASIS/wiki) for details on installation, usage, output files and more.**

version: v3     
release: [https://github.com/atulkakrana/PHASIS/releases](https://github.com/atulkakrana/PHASIS/releases)     
updated: 03/10/2017    
citation: [Preprint-Copy](http://www.biorxiv.org/content/early/2017/07/03/158832.full.pdf+html)   

***PHASIS*** (*PHAS* Inspection Suite) is a parallelized suite that provides an integrated solution for the large-scale survey of tens to hundreds of sRNA libraries for the following applications: a) *de novo* discovery of *PHAS* loci and precursor transcripts, b) a summarization of *PHAS* loci from specific groups of sRNA libraries, c) a comparison of *PHAS* summaries between groups corresponding to samples from different stages, tissues and treatments, d) quantification and annotations of *PHAS* loci, and e) discovery of their miRNA triggers. *PHASIS* generates easily parsed output files for downstream bioinformatics analysis, formatted result files for immediate consumption and organized ancillary data to facilitate optimizations like a re-summarization to exclude or include libraries. The three *PHASIS* components are:

**1.** ***phasdetect***: Identifies library-specific phased-loci

**2.** ***phasmerge***: Summarizes the libraries-specific results to a non-overlapping genome- or transcriptome-level set, comapares *PHAS* summaries between groups and annotates *PHAS* loci

**3.** ***phastrigs***: Predicts miRNAs triggers for *PHAS* locus using sPARTA miRNA-target predictions. If supplied with PARE/Degradome data then data supported triggers are reported



**Authors:**     
Atul kakrana     
kakrana@udel.edu       

Pingchuan Li   
lipingchuan@gmail.com
