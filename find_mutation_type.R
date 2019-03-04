library("stringr")

setClass("Codon",
         representation = representation(sequence = "character",
                                         start = "numeric",
                                         end = "numeric",
                                         protein = "character"),
         prototype = prototype(sequence = NA_character_,
                               start = NA_real_,
                               end = NA_real_,
                               protein = NA_character_))




orf3NT <- "ATGGCTGGCGCACTCTTTGGTGCGATTGGAGGTGGCCTGATGGGCATAATTGGCAATTCCATCTCAACAGTCCAGAATCTTCAGGCAAATAAACAATTGGCTGCACAGCAATTTGGCTATAATTCCTCTCTGCTTGCAACGCAAATTCAGGCCCAGAAGGATCTCACACTGATGGGGCAGCAGTTCAACCAGCAGCTCCAAGCCAACTCTTTCAAGCATGACCTTGAGATGCTTGGCGCCCAGGTGCAAGCCCAGGCGCAGGCCCAGGAGAACGCTATCAACATCAGGTCGGCGCAGCTCCAGGCCGCAGGCTTTTCAAAGTCCGACGCCATTCGCTTGGCCTCGGGGCAGCAACCGACGAGGGCCGTTGACTGGTCTGGGACGCGGTATTACGCCGCTAACCAGCCGGTTACGGGCTTCTCGGGTGGCTTCACCCCAAGTTACACTCCAGGTAGGCAAATGGCAGTCCGCCCTGTGGACACATCCCCTCTACCGGTCTCGGGTGGACGCATGCCGTCCCTTCGTGGAGGTTCCTGGTCTCCGCGTGATTACACGCCGCAGACCCAAGGCACCTACACGAACGGGCGGTTTGTGTCCTTCCCAAAGATCGGGAGTAGCAGGGCATAG"

orf3NTPositions <- seq(from = 1, to = nchar(orf3NT), by = 1)

# . matches any character
# .{3} match any character 3 x's
# note in the backreference there is a white space char "\\1<ws>"
codons <- gsub("(.{3})", "\\1-", orf3NT)
head(codons)

codonsSeparate <- str_split(codons, pattern = "-") #contains empty trailing "" that should be removed

orf3AA <- "MAGALFGAIGGGLMGIIGNSISTVQNLQANKQLAAQQFGYNSSLLATQIQAQKDLTLMGQQFNQQLQANSFKHDLEMLGAQVQAQAQAQENAINIRSAQLQAAGFSKSDAIRLASGQQPTRAVDWSGTRYYAANQPVTGFSGGFTPSYTPGRQMAVRPVDTSPLPVSGGRMPSLRGGSWSPRDYTPQTQGTYTNGRFVSFPKIGSSRA"
aa <- gsub("(.{1})", "\\1-", orf3AA)
aaSeparate <- str_split(aa, pattern = "-")

codonVec <- codonsSeparate[[1]]
# remove that last weird element
codonVec <- codonVec[1:length(codonVec)-1]
aaVec <- aaSeparate[[1]]
aaVec[aaVec == ""] <- "STOP"

# name the codon vector
names(codonVec) <- aaVec
head(codonVec)




codonClassList <- vector(mode = "list", length = length(codonVec))

for (i in 1:length(codonClassList)) {
  if (i == 1) {
    newCodonClass <- new("Codon",
                         sequence = codonVec[i],
                         start = i,
                         end = (i + 2),
                         protein = names(codonVec)[i])
    codonClassList[[i]] <- newCodonClass
  } else {
    newCodonClass <- new("Codon",
                         sequence = codonVec[i],
                         start = ((i*3) - 2),
                         end = (i*3),
                         protein = names(codonVec)[i])
    codonClassList[[i]] <- newCodonClass
  }
}


# These positions will need to be adjusted for their actual position within the genome....
#  They don't all start at codon 1....for example this is codon 3 and it starts at position 6680
#  Need to get this correct to be able to correctly place variants from the VCF.
