
## phasTER tool-set for phased clusters discovery
**version**: rc1  
**updated**: 10/24/2015  

phasTER has three user accessible components and these explained below in order of their use.

### 0. Installation
Basic knowledge of UNIX/Linux command-line interface is expected to install and use phasTER tool set. If you are not comfortable using the command-line then take help from your Linux administrator or geeky colleague. Follow these steps to generate phasTER binaries for your machine:

1. Download the archive for latest phasTER version from our release page: https://github.com/atulkakrana/phasTER/releases
2.  Copy the archive to the machine, which you intend to use for your analysis and decompress:

    For *.zip archive:
``` unzip phasTER*.zip```
and for *.tar.gz archive:
```tar -xvzf phasTER*.tar.gz```
3. Install or compile the scripts by using bundled `install.sh` script:
    ```sh install.sh```
4. At this point you should be able to see binaries for `phaser`, `collapser` and `revferno`
5. Copy these binaries to your working folder and use them for you analyses.

**Note**: These binaries are good to use of the machine on which these were compiled. Please repeat the above process, if you intend to use phasTER on a different machine at some point in your study. 



### 1. phaser: Library-specific phased-loci using the phasTER-core-scripts
  * Massively parallel – both at genome and library-level
  * Optimizes for minimum runtime
  * Works at both genome- and transcriptome-level

### 2. collapser: Collapses loci from libraries-specific results to a non-overlapping genome-level set
  * Optimizes for best locus across different libraries
  * Provides most-accurate 5’ and 3’-ends

**Usage**
  1. Copy the `collapser` binary to same folder where phaser analysis was done
  2. Run collapser `python3 collapser -pval <p-value cutoff>`
  3. Results will be placed in *collapsed* folder

**Note**
Authors recommend using a robust p-value cutoff like 1e-05 or less. If you do not use the `-pval` switch, and execute 
`collapser` then a default pvalue of 1e-05 is used. `collapser` can be executed without `-pval` switch by running command:
`python3 collapser`

####3. summarizer: coming later this week

####4. revFerno: coming soon


Author: kakrana@udel.edu
