# GSA_Mouse-SNP-Panel

Overview

This repository hosts the computational scripts developed for a study focused on screening minimal single nucleotide polymorphism (SNP) panels from whole-genome data. The core goal is to identify SNP subsets that can accurately replicate the phylogenetic relationships among mouse strains, providing a robust foundation for high-throughput genetic monitoring.

Dataset & Methodology

Data Source: Whole-genome SNP data of 51 mouse strains were obtained from the Mouse Genome Informatics (MGI) database, encompassing a total of 315,669 SNP loci.



Benchmark Phylogeny Construction: Hamming distance was used to quantify genetic divergence between strains. A benchmark phylogenetic tree was then constructed using the Neighbor-Joining (NJ) method based on the full SNP dataset.



GSA-Based SNP Screening: The Generalized Simulated Annealing (GSA) algorithm was implemented for iterative dimensionality reduction of SNP markers. The fitness function was defined as achieving a Robinson-Foulds (RF) distance of zero, which ensures the phylogenetic topology derived from the screened SNP panel is identical to the whole-genome benchmark tree.



Quality Control: Stringent quality control filters were applied to exclude low-quality variants, ensuring the final SNP panels meet high genetic data standards.

Key Results

Through algorithmic screening and quality control, two minimal SNP panels were successfully extracted:
A 512-locus panel containing strain-specific haplotypes (hSNPs)
A 297-locus core panel

Phylogenetic trees constructed using both panels showed complete topological consistency with the whole-genome benchmark tree, with no branching discrepancies observed. This confirms that the selected SNP subsets retain all critical phylogenetic information while eliminating redundant genetic signals.

Significance

This study introduces a novel paradigm for genetic marker feature selection, using phylogenetic topological invariance as a strict screening criterion. The methodology not only drastically reduces the number of required genotyping loci but also lays a theoretical and technical framework for developing automated, high-throughput SNP monitoring chips for mice and other model organisms.

Usage

Detailed instructions for script execution, data preprocessing, and result validation are provided in the docs directory. Users can adapt the scripts to their own SNP datasets by modifying configuration parameters specified in the config files.

Citation

If you use the resources in this repository in your research, please cite our corresponding publication (citation information will be updated upon publication). 
