## Welcome to ViRGA home page

ViRGA is a reference guided assembler of highly variable viral genomes, based on iterative mapping and de novo reassembling of highly variable regions, which can handle with distant reference sequence due to specially designed read mapper. ViRGA can separate mixtures of strains of different intraspecies genetic groups (genotypes, subtypes, clades, etc.) and assemble a separate consensus sequence for each group in a mixture.

If provided with multiple sequence alignment (MSA) of target references ViRGA selects optimal reference set, sorts reads to selected references and outputs consensus sequences corresponding to these references. For each consensus sequence the multiple sequence alignment of its constituent reads is printed in BAM format.

If no MSA provided, ViRGA works in single-reference mode and use user-provided reference.

Multi-fragment references are supported in single-reference mode.

You can use ViRGA for full genome assembly or just to find optimal reference set for given fastq files with Illumina paired end reads.

### Installation

ViRGA is a pure java application: it runs on any platform supporting JVM. Simply download jar file and run according to [usage instructions](https://github.com/gFedonin/ViRGA/wiki) or integrate in your own java application.


### Required dependencies

The following are required to run ViRGA:

-Java version 8 or higher

-[USEARCH](http://www.drive5.com/usearch/) binary in any location. Path to the binary is set in configuration file. Recomended version is included in the distribution.

-Some external libraries:

  1. [JDOM](http://www.jdom.org/): Java-based solution for accessing, manipulating, and outputting XML data from Java code.
  
  2. [Picard Tools](https://broadinstitute.github.io/picard/): A set of command line tools (in Java) for manipulating high-throughput sequencing (HTS) data and formats such as SAM/BAM/CRAM and VCF.
  
  3. [The Trove library](http://trove.starlight-systems.com/) provides high speed regular and primitive collections for Java.

Corresponding jars should be in the same directory with ViRGA.jar. The recommended versions are: jdom.jar version 2.0.5, picard.jar version 2.6.0 and trove.jar version 3.0.3. These jars are included in the the distribution.

For usage help and examples, see the [ViRGA wiki page](https://github.com/gFedonin/ViRGA/wiki).
