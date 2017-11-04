## Welcome to VirGenA home page

VirGenA is a reference guided assembler of highly variable viral genomes, based on iterative mapping and de novo reassembling of highly variable regions, which can handle with distant reference sequence due to specially designed read mapper. VirGenA can separate mixtures of strains of different intraspecies genetic groups (genotypes, subtypes, clades, etc.) and assemble a separate consensus sequence for each group in a mixture.

If provided with multiple sequence alignment (MSA) of target references VirGenA selects optimal reference set, sorts reads to selected references and outputs consensus sequences corresponding to these references. For each consensus sequence the multiple sequence alignment of its constituent reads is printed in BAM format.

If no MSA provided, VirGenA works in single-reference mode and use user-provided reference.

Multi-fragment references are supported in single-reference mode.

You can use VirGenA for full genome assembly or just to find optimal reference set for given fastq files with Illumina paired end reads.

### Documentation

Complete documentation is provided in [wiki](https://github.com/gFedonin/VirGenA/wiki) format.

### Installation

VirGenA is a java application: it runs on any platform supporting JVM. Simply download [release](https://github.com/gFedonin/VirGenA/releases) file and run according to [usage instructions](https://github.com/gFedonin/VirGenA/wiki).


### Required dependencies

The following are required to run VirGenA:

-Java version 8 or higher

-[USEARCH](http://www.drive5.com/usearch/) binary in any location. Path to the binary is set in configuration file. Recomended version is included in the distribution.

-[Blast](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download) installed locally

-Some external libraries:

  1. [JDOM](http://www.jdom.org/): Java-based solution for accessing, manipulating, and outputting XML data from Java code.
  
  2. [Picard Tools](https://broadinstitute.github.io/picard/): A set of command line tools (in Java) for manipulating high-throughput sequencing (HTS) data and formats such as SAM/BAM/CRAM and VCF.
  
  3. [The Trove library](http://trove.starlight-systems.com/) provides high speed regular and primitive collections for Java.

Corresponding jars should be in the same directory with VirGenA.jar. The recommended versions are: jdom.jar version 2.0.5, picard.jar version 2.6.0 and trove.jar version 3.0.3. These jars are included in the the [release](https://github.com/gFedonin/VirGenA/releases).

### Toy example

To run VirGenA with test data download and unzip [release](https://github.com/gFedonin/VirGenA/releases) files.

on Windows:

You can set number of threads in config_test_win.xml by changing value of **ThreadNumber** element.

Using Windows command promt change dir to unzipped folder and type:

**java -cp ./VirGenA.jar RefBasedAssembler config_test_win.xml**   

on Linux:  

You can set number of threads in config_test_linux.xml by changing value of **ThreadNumber** element.

Change permissions of ./tools/usearch8.1.1861_i86linux32 to make it executable. After that using shell change dir to unzipped folder and type:

**java -cp ./VirGenA.jar RefBasedAssembler config_test_linux.xml**

Test data is an artificial mixture containing 100000 HIV paired reads of three different subtypes (01_AE, B and C) in equal proportions. VirGenA should detect these components and assemble genome-length consensus sequences for all components.

Results will be stored in ./res/ folder. Expected output is:  
1. Files (fasta) with assemblies of three mixture components named after the selected references:
01_AE.TH.90.CM240.U54771_assembly.fasta, B.FR.83.HXB2_LAI_IIIB_BRU.K03455_assembly.fasta, C.BW.96.96BW0502.AF110967_assembly.fasta   
2. Sorted bam files with read alignments and corresponding index files (bai): 'reference_name'_mapped_reads.bam and 'reference_name'_mapped_reads.bai  
3. Log file.

### How to cite:

Fedonin GG, Fantin YS, Favorov AV, Shipulin GA, Neverov AD. _VirGenA: a reference-based assembler for variable viral genomes._ Brief Bioinform, 2017 Jul 28. doi: 10.1093/bib/bbx079.
