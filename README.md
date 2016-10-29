


#phasTER tool-set for phased clusters discovery
**version**: rc1  
**updated**: 10/24/2015  

phasTER has three user accessible components and these explained below in order of their use.

## 0. Installation
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

**Note**: These binaries are good to use on the machine on which these were compiled. Please repeat the above process, if you intend to use phasTER on a different machine. 



## 1. phaser: Identifies library-specific phased-loci 

`phaser` identifies the phased sRNA clusters using the phasTER-core scripts.  It uses a parallelized algorithm to reduce run time and works at both genome and transcriptome level.

###**Usage**

1. Copy the `phaser` binary and settings file `phaser.set` to a new folder or some existing analysis folder along with your sRNA libraries
2. Provide required information via the settings file  `phaser.set`. These settings are described below
3. Execute `phaser` from command line: `python3 phaser`
4. A successful run will create a new results folder with its name encoding information of time the analysis was started - `phased_month_day_hour_seconds`

###**Output**
The result folder contains library-specific files of two different types:
a) `*.list`  - List of phased coordinates identified at different *p-value* cutoffs and 
b) `*.clust` - Cluster files corresponding to different *p-values* and containing information on siRNA for each of these loci. 

###**phaser settings**
Here is an example of settings file:

```
<<< Mandatory Settings, value in text>>>
@runType    	= Y
@reference  	= /home/kakrana/99.genomes/rice.msu.7/test.fa
@index 			= 
@userLibs 		= 2599_chopped.txt,2600_chopped.txt
@libFormat 		= T					
@phase          = 24
```

**Description**

`@runType` - This is to differentiate and to specify the type of reference sequence provided by the user i.e. whole genome or transcriptome/scaffolds. Accepted values are `Y` when reference FASTA is whole genome or  `N` if reference is transcriptome or scaffold FASTA file.

`@reference` - A full path to the reference FASTA file. This will be used to make Bowtie index, if not provided in `@index` parameter, and will be used to rename some output files.

`@userLibs` - sRNA library filenames, provided as a comma separated list. Since, library names are used by `phaser` to name the output files, and other downstream scripts to pick up the result files, users are advised to provide simple names to their sRNA libraries i.e. no special characters.

`@libFormat` - File format for sRNA sequencing library. In present version on tag-count format is supported. Accepted value is `T`. Please see *tally* (http://wwwdev.ebi.ac.uk/enright-dev/kraken/reaper/src/reaper-latest/doc/tally.html) to convert your FASTA file to tag-count. Instructions on how to install *tally* could be found here: https://github.com/atulkakrana/helper.github#install-prerequisites

`@phase` - Periodicity of phased siRNA clusters. Accepted values are `21` or `22` or `24`.

`@index` - A full path to the Bowtie index of your reference FASTA. If not provided i.e. left blank, then index will be prepared in first run, and will be reused by the subsequent `phaser` runs, until value of `@reference` parameter is changed or reference FASTA file is modified by the user. If index is provided then the `phaser` run will take much less time. 



##2. collapser: Summarizes the libraries-specific results to a non-overlapping genome- or transcriptome-level set

`phaser` can be used to generate output files for tens- to hundreds of libraries, but a user might be interested in phased-loci from a specific sub-set of libraries. `collapser` does three jobs:

**a.)** summarizes the results for set of libraries, and provides additional information required to filter out false positives  along with the library-specific siRNAs abundances

**b.)** collapses the redundant loci from multiple libraries to a non-redundant genome or transcriptome level set

**c.)** computes precise start and stop coordinates of phasing for each loci using data from multiple libraries which is imperative for detection of their miRNA triggers.

###**Usage**
  **1**. Copy the `collapser` binary to the folder which contains phaser results folder
  **2**. Run collapser `python3 collapser -pval <p-value cutoff> -dir <phased results folder>`
  **3**. Results will be placed in *collapsed* folder

###**Output**
**1**. `*.collapsed.txt`: A non-redundant list of PHAS loci for given genome or transcriptome.

**2**. `*.summary.txt`: Summary of non-redundant PHAS, containing library-specific abundances of phasIRNAs and other information for quality filtering.

**3.**`*.phasi.csv*`: A comma-separated file containing phasiRNAs sequences in first column and additional PHAS-specific information for quality filtering. The CSV file can be filtered on quality parameters to get a confident set of phasiRNAs.

**Quality Parameters**

`bestLib`:  sRNA library name which has longest PHAS length i.e. maximum cycles of phasiRNAs detected
`bestKVal`: The maximum and uninterrupted cycle of phasiRNAs detected for specific PHAS locus.  PHAS locus with >= 6 cycles can be considered as a true PHAS locus.
`phasiDominance`: Ratio of phasiRNAs abundances against total siRNAs from specific PHAS locus. A high ratio (~0.75) is expected from a true PHAS locus.
`MaxTagRatio`: Ratio of most abundant siRNA against all the phasiRNAs from specific PHAS locus. A higher ratio would indicate that PHAS locus could be a lineage or species-specific miRNA.
`totalPHASAbundance`: Sum of the abundances of phasiRNAs. 

**Note**
Authors recommend using a robust p-value cutoff like 1e-05 or less. If user executes `collapser` without the `-pval` switch, then a default *p-value* of 1e-05 is used. `collapser` can be executed without `-pval` switch by running command: `python3 collapser -dir <phased results folder>`



##3. revferno: Identifies miRNAs triggers for PHAS locus using sPARTA miRNA-target predictions

`revferno` performs an exhaustive check for miRNA triggers of PHAS loci using the predicted or PARE-validated target sites from sPARTA. It works downstream of `miRferno`, which is an exhaustive algorithm to predict miRNA targets, thereby increasing the likelihood of capturing the miRNA trigger of PHAS loci. `revferno` can test hundred thousands `mirferno` predicted targets as triggers of PHAS loci, in a reasonable time ranging from few seconds to minutes. In addition, it scans both 5'-head and 3'-tail of PHAS loci, accounts for missing phasiRNAs, dicer-and strand-offsets for high sensitivity. 

###**Arguments**
| Switches		|Description															|
|:--------------|:-----------------------------------------------------------------------------------------|
|-coord 		| Concatenated `all_coords` file from 'genic' and 'intergenic' target prediction or validation analyses using sPARTA																											|
|-predfile		| Concatenated  `All.targs.csv` or `All.libs.validated.uniq.csv` file from 'genic' and 'intergenic' target prediction or PARE-validation output of sPARTA						|
|-predtype	| Type of miRNA target results provided as input to revferno. P for predicted. D for degradome validated																													|
|-score			| Cutoff for complementarity based miRNA-target score. Default is 6. See [sPARTA] (http://nar.oxfordjournals.org/content/42/18/e139) paper for more details			|
| -phas			| Path to 'collapsed' file generated by the `collapser`					|


###**Usage**
**1**. Create a new folder for `revferno` analysis. For this README and as an example, we are using the name REVFERNO for folder.

**2**. Prepare for `revferno` analysis by first identifying the genome-wide or transcriptome-specific target sites of candidate miRNAs through sPARTA. For this use the same reference FASTA that was used for `phaser` analysis.

	**a.)** In case your reference is whole genome, predict or PARE-validate targets from both 'genic' and 'intergenic' regions. If your reference is transcriptome then predict targets by using it with `-featureFile` option and `genomeFeature` value of 0. See sPARTA README for more details: https://github.com/atulkakrana/sPARTA.github/tree/master/sparta.
	
	**b.)** For reference genome, you will have `allcoords` file from analysis of 'genic' and 'intergenic' regions. Concatenate these two files . If the reference was transcriptome then there is just one `allcoords` file. Copy this concatenated or single `allcoords` file to REVFERNO folder which you indent to use for `revferno` analysis.
	
	**c.)** In case of genome reference, concatenate the predicted targets file generated by sPARTA from  'genic' and 'intergenic' analysis. These files could be found in `predicted` folder and named as `All.targs.csv`. If your  reference was transcriptome then no concatenation is required as there is just one of these.  Copy this concatenated or single `All.targs.csv` file to REVFERNO folder. Alternatively, you can also use the PARE-validated targets file `All.libs.validated.uniq.csv` from `output` folder instead of predicted targets `All.targs.csv` file.
	
**3.** Execute `revferno` from command line: ```python3 revFerno.py -coord <allcoords file> -predtype <P/D> -predfile <PATH/TO/miRNA target file> -phas <PATH/TO/COLLAPSED/PHAS LOCI FILE>```

**4.** The analysis will generate a new folder `triggers` with results inside it.

###**Output**
`*.trigger.csv` File contains information from miRNA target prediction/validation file along with PHAS loci, PHASindex and matFlag. PHASindex states the PHAS cycle where trigger most-likely cleaved. These are most likely 0 which indicates that predicted cleave site matches with first phasiRNA (P1) 5-start site. The +/- sigh reflects left (outside) and right (inside) of first phasiRNA i.e. +1 index would mean that cleave site matches one cycle before the actual 5'-phasiRNA (P1). matFlag indicates dicer offset (+1/-1 nt) or strand offset that might have been observed by revferno.


Author: kakrana@udel.edu
