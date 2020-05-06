# Mutation-Simulation
## Introduction:
A small project using bioinformatic tools in R to analyse point mutation in oncogenes. Oncogenes are define as genes that are prone to transform a cell into a cancerous cell via the slightest mutation. This script analyses the probability that a random point mutation in the three selected oncogenes (KRas, OR1A1, and PTPN11) will result in a mutation that can potentially transform the cell into a cancerous one. Below gives a brief summary of the genes.

The KRas gene gives functionality to the GTPase K-Ras which plays a major role in the RAS/MAPK pathway of a cell. This pathway play important roles in apoptosis, cell division, and cell differentiation. OR1A1 is a protein coding gene that plays an important role in pathways for peptide ligand-binding receptors and signaling by GPCR. The PTPN11 gene codes for a protein called SHP-2, and this protein helps regulate the RAS/MAPK signaling pathway. 

## Visualization of Results
![](https://github.com/a-rigido/Mutation-Simulation/blob/master/MutationSummary.png)
## Interpretation of Results:

Observing the general trend of the mutation probability of my R code output, the most likely type of mutation is a missense mutation, following by a silent mutation, followed by a nonsense mutation. For the 3 given genes, this general trend is consistent with data observed from the IntOGen website. This, of course, does not paint the entire picture because upon examination the chances of each type of mutation generated from the R script can vary wildly from the experimental data from the site. An obvious example of this is seen between the KRas probabilities. For missense 0.99 recorded, 0.73144 generated, for silent 0.01 recorded, 0.2059 generated, and for nonsense 0 recorded, 0.05277 generated. 

The difference between the two is quite clear, which indicates there must be a selection bias favoring missense mutations over the other types. I believe further study about the driving forces behind this extreme bias is needed, as this behavior appears to be severely inconsistent with a random process. Of all the genes compared, the generated PTPN11 mutation probabilities share the closest resemblance to that of IntOGen's with all values differing by less than 7%, which all things considered is a nice result, as one could make an argument for random processes being the driving force behind these probability distributions. My overall conclusion is that a mutation probability that is generated randomly cannot be consistently accurate enough to be looked upon as a viable model for real experimental observation, it is simply appropriate too show general trends regarding mutation types and their probabilities.

<p align="center">
Copyright (c) 2019 Alex Rigido
</p> 
