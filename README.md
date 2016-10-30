##phasTER tool-set for phased clusters discovery

version: rc1
updated: 10/24/2015

**phasTER** is a parallelized toolset for phased siRNAs clusters discovery. It has three user accessible components that enables discovery of phased siRNA loci from tens to hundreds of sRNA libraries, prepares summary of these loci containing library-specific abundances and other quality-filtering parameters and finally identifies miRNA triggers of these PHAS loci. The three **phasTER** components are:

**1.** **phaser**: Identifies library-specific phased-loci

**2.** **collapser**: Summarizes the libraries-specific results to a non-overlapping genome- or transcriptome-level set

**3.** **revferno**: Predicts miRNAs triggers for PHAS locus using sPARTA miRNA-target predictions

Please see phasTER [wiki] (https://github.com/atulkakrana/phasTER/wiki) for more details on installation, usage, output files and more.

author: kakrana@udel.edu
