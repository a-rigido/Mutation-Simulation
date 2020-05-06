# Alex Rigido
# Mutation simulaiton

# Load appropriate libraries
library(Biostrings)

# Load data file
load(file = "MutationSimulationCodons.RData")

# The following codons that have been prepared from this file
KRascodons
OR1A1codons
PTPN11codons

#====================================================================
# Task: Write code that executes a loop N times 
# (for N <- 100000) to create a point mutation 
# randomly in each of the three genes. Keep track 
# of the number of missense, silent ("synonymous"), 
# and nonsense ("truncating")" mutations you find.

N <- 100000
NTs <- c("A","G","T","C") # DNA Nucleotides

# Let's start by creating a storage for the each
# gene and type of mutation

# for KRas (K)
K_silMut <- 0 
K_misMut <- 0
K_nonMut <- 0  

# for OR1A1 (O)
O_silMut <- 0
O_misMut <- 0 
O_nonMut <- 0 
# for PTPN11 (P)
P_silMut <- 0
P_misMut <- 0
P_nonMut <- 0

# Now that we have somewhere to store data about
# the generated mutations, let's build a for-loop
# to accomplish this.

# set seed to duplicate randomness
set.seed(1002743944)
for (i in 1:N){
  # Generate a Random mutation via random sampling of each of the genes
  # which keeps assigning a random index number to the variable with
  # each iteration of the for-loop. This will feed forward into the 
  # code next to select mutation postion. 
  
  randK <- sample(grep(sample(KRascodons,1),KRascodons),1) 
  randO <- sample(grep(sample(OR1A1codons,1),OR1A1codons),1)  
  randP <- sample(grep(sample(PTPN11codons,1),PTPN11codons),1)
  
  # Simulate the point mutations
  
  K_code <- GENETIC_CODE[KRascodons[randK]]
  mtK <- KRascodons[randK]                    
  triplet <- unlist(strsplit(mtK, ""))           # split into three nucl.
  iNuc <- sample(1:3, 1)                         # choose one of the three
  mutNuc <- sample(NTs[NTs != triplet[iNuc]], 1) # chose a mutated nucleotide
  triplet[iNuc] <- mutNuc                        # replace the original 
  mtK <- paste0(triplet, collapse = "")          
  mtK_code <- GENETIC_CODE[mtK]
  
  O_code <- GENETIC_CODE[OR1A1codons[randO]]
  mtO <- OR1A1codons[randO]                    
  triplet <- unlist(strsplit(mtO, ""))           # split into three nucl.
  iNuc <- sample(1:3, 1)                         # choose one of the three
  mutNuc <- sample(NTs[NTs != triplet[iNuc]], 1) # chose a mutated nucleotide
  triplet[iNuc] <- mutNuc                        # replace the original 
  mtO <- paste0(triplet, collapse = "")          
  mtO_code <- GENETIC_CODE[mtO]
  
  P_code <- GENETIC_CODE[PTPN11codons[randP]]
  mtP <- PTPN11codons[randP]                    
  triplet <- unlist(strsplit(mtP, ""))           # split into three nucl.
  iNuc <- sample(1:3, 1)                         # choose one of the three
  mutNuc <- sample(NTs[NTs != triplet[iNuc]], 1) # chose a mutated nucleotide
  triplet[iNuc] <- mutNuc                        # replace the original 
  mtP <- paste0(triplet, collapse = "")          
  mtP_code <- GENETIC_CODE[mtP]
  
  # Now each iteration of the for-loop will have a mutation, and the following
  # part of the for-loop will differentiate between the aforementioned types.
  
  # KRas Mutations:
  if(K_code != mtK_code){
    if(mtK_code == "*" | K_code == "*"){ # add/remove stop codon
      K_nonMut <- K_nonMut + 1
    } else{ # all other non-silent mutations
      K_misMut <- K_misMut + 1 
    }
  } else{ # silent mutations
    K_silMut <- K_silMut + 1
  }
  
  # OR1A1 Mutations:
  if(O_code != mtO_code){
    if(mtO_code == "*" | O_code == "*"){ # add/remove stop codon
      O_nonMut <- O_nonMut + 1
    } else{ # all other non-silent mutations
      O_misMut <- O_misMut + 1 
    }
  } else{ # silent mutations
    O_silMut <- O_silMut + 1
  }
  
  # PTPN11 Mutations:
  if(P_code != mtP_code){
    if(mtP_code == "*" | P_code == "*"){ # add/remove stop codon
      P_nonMut <- P_nonMut + 1
    } else{
      P_misMut <- P_misMut + 1 
    }
  } else{ # silent mutations
    P_silMut <- P_silMut + 1
  }
}

# Give the for-loop time to run. Once done the Mutations have been generated.

# for KRas (K)
(K_silMutFREQ <- (K_silMut)/(N))
(K_misMutFREQ <- (K_misMut)/(N))  
(K_nonMutFREQ <- (K_nonMut)/(N))
# Check if normalized
K <- c(K_silMutFREQ, K_misMutFREQ, K_nonMutFREQ)
all.equal(sum(K), 1) # True

# for OR1A1 (O)
(O_silMutFREQ <- (O_silMut)/(N))
(O_misMutFREQ <- (O_misMut)/(N))
(O_nonMutFREQ <- (O_nonMut)/(N))
# Check if normalized
O <- c(O_silMutFREQ, O_misMutFREQ, O_nonMutFREQ)
all.equal(sum(O), 1) # True

# for PTPN11 (P)
(P_silMutFREQ <- (P_silMut)/(N)) 
(P_misMutFREQ <- (P_misMut)/(N)) 
(P_nonMutFREQ <- (P_nonMut)/(N))
# Check if normalized
P <- c(P_silMutFREQ, P_misMutFREQ, P_nonMutFREQ)
all.equal(sum(P), 1) # True

# visualize the results
tab <- matrix(c(K,O,P), ncol=3, byrow = T)
colnames(tab) <- c("Silent", "Mismatch", "Nonsense")
rownames(tab) <- c("KRas", "ORA1A", "PTPN11")
tab <- as.table(tab)

tab   # summary table of mutation frequencies

{barplot(tab, 
        main = "Summary of Mutation Types from Simulation", 
        ylab="Mutation Frequency", xlab = "Type of Mutation",
        ylim=c(0,1),
        col = c("red","blue", "green"),
        beside = T)
legend("topright", title = "Oncogenes",
       c("KRas", "ORA1A", "PTPN11"),
       fill = c("red","blue", "green"))
}

#
# [END]

#============================================================================
