##*PHASworks* suite for phased clusters discovery, comparison, annotation and to identify miRNA triggers

version: release canditate v3    
updated: 03/05/2017

***PHASworks*** is a parallelized suite for phased siRNAs clusters characterization. It has three user accessible components that enables discovery of phased siRNA loci from tens to hundreds of sRNA libraries, prepares summary of these loci containing library-specific abundances and other quality-filtering parameters and finally identifies miRNA triggers of these PHAS loci. The three ***PHASworks*** components are:

**1.** **phasdetect**: Identifies library-specific phased-loci

**2.** **phasmerge**: Summarizes the libraries-specific results to a non-overlapping genome- or transcriptome-level set, comapares PHAS summaries between groups and annotates PHAS loci

**3.** **phastrigs**: Predicts miRNAs triggers for PHAS locus using sPARTA miRNA-target predictions. If supplied with PARE/Degradome data then data supported triggers are reported

Please see *PHASworks* [wiki] (https://github.com/atulkakrana/PHASworks/wiki) for details on installation, usage, output files and more.

Author: Atul kakrana
E-mail: kakrana@udel.edu or kakrana@gmail.com
