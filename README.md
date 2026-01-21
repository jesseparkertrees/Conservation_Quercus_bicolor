# Quercus_bicolor_conservation
PROJECT DESCRIPTION This repository contains data and analysis code for the paper "Genomically informed seed orchard design for trailing-edge tree populations: a perspective from Quercus bicolor conservation"
Jesse B Parker, Sean Hoban, Laura M Thompson, Scott E Schlarbaum DOI: ______

ABSTRACT The loss of forest genetic resources due to land use change is both substantial and ongoing. Climate change, along with other anthropogenic environmental alterations, accelerates this loss and underscores the need to conserve diverse genetic material with potential for adaptation. Range-marginal populations, particularly those at trailing edges, are critical conservation targets, as they may harbor alleles important for resilience, yet they face heightened risk of extirpation under environmental change. However, effectively conserving these populations is challenging: their marginal position often results in small, fragmented populations that are difficult to locate and adequately sample, and may have reduced genetic diversity. Additionally, elevated rates of interspecific introgression at range edges can complicate conservation decision-making. In this work, we discuss considerations broadly relevant to conservation of range-marginal germplasm and illustrate how these considerations can be practically addressed, using a widespread and important tree for restoration, Quercus bicolor Willd., as a case study. We outline considerations and approaches to (1) locating source populations, using a combination of biodiversity dataset review, species distribution modeling, and field surveying; (2) genomically assessing hybridization and selecting tolerances for introgression in ex situ conservation; (3) evaluating patterns of genetic variation using population genomics techniques; and (4) using these data to optimize the selection of genotypes for the establishment of conservation seed orchards, with attention to introgression, population structure, genetic diversity, minimum sample size estimates, and inbreeding potential.    We find that our approach offers a substantial improvement over an “uninformed” sampling strategy with potential long-term benefits for conservation and restoration.

Sequence reads available at SRA BioProject: PRJNA1260989.
Data processing and analysis for RADseq data of Quercus bicolor available at https://github.com/jesseparkertrees/Quercus_bicolor

DESCRIPTION OF FILES IN main directory

File: "rcode.R" Description: R file containing code for all analyses.

File: "qubi.vcf" Description: VCF file of the 135 sampled individuals filtered to a call rate of 1 (7955 biallelic SNPs). 

File: "qubi95.vcf" Description: VCF file of the 135 sampled individuals filtered to a call rate of 0.95 (100968 biallelic SNPs). 

File: "ch33.csv" Description: Summary of individuals selected for final orchard.

File: "Qbicolor_data.csv" Description: Summary data of all sampled individuals.
