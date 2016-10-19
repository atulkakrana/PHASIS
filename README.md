### phasTER tool-set for phased clusters discovery
**version**: beta-8  
**updated**: 10/13/2015  

phasTER has four user accesible components and these will be released in particular

####1. phaser: library-specific phased-loci using the phasTER-core-scripts
  * Massively parallel – both at genome and library-level
  * Optimizes for minimum runtime
  * Works at both genome- and transcriptome-level

####2. collapser: Collapses loci from libraries-specific results to a non-overlapping genome-level set
  * Optimizes for best locus across different libraries
  * Provides most-accurate 5’ and 3’-ends

##### Usage
  1. Copy the collapser binary to same folder where phaser analysis was done
  2. Run collapser `python3 collapser -pval <p-value cutoff>`
  3. Results will be placed in *collapsed* folder

##### Note
Authors recommend using a robust p-value cutoff like 1e-05 or less. If you do not use the `-pval` switch, and execute 
`collapser` then a default pvalue of 1e-05 is used. `collapser` can be executed without `-pval` switch by running command:
`python3 collapser`

####3. summarizer: coming later this week

####4. revFerno: coming soon


Author: kakrana@udel.edu
