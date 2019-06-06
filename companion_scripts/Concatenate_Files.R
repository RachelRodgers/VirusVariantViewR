# Concatenate_Files.R

# This script is used to concatenate certain EH MNV sequencing files with
#   certain Baldridge_rumspringa MNV sequencing files to see if the combining
#   helps to improve variant detection.  This is part of the VariantViewR project.

library("ShortRead")
library("stringr")

allFilesList <- list.files("./raw_data/", pattern = "*.fastq.gz")

# First generate a map to link patterns (regex) from the EH files to the rumspringa files.
#testPattern <- grepl("^EH_Baldridge_4_.*_R1_.*", x = allFilesList)
#any(testPattern) == TRUE
#whosTrue <- allFilesList[testPattern2]
#whosTrue

patternMap <- list("^EH_Baldridge_2_.*_R1_.*" = "^Baldridge_1_.*_R1_.*",
                   "^EH_Baldridge_2_.*_R2_.*" = "^Baldridge_1_.*_R2_.*",
                   "^EH_Baldridge_3_.*_R1_.*" = "^Baldridge_5_.*_R1_.*",
                   "^EH_Baldridge_3_.*_R2_.*" = "^Baldridge_5_.*_R2_.*",
                   #"^EH_Baldridge_4_.*_R1_.*" = "^Baldridge_6_.*_R1_.*",
                   #"^EH_Baldridge_4_.*_R2_.*" = "^Baldridge_6_.*_R2_.*",
                   "^EH_Baldridge_5_.*_R1_.*" = "^Baldridge_7_.*_R1_.*",
                   "^EH_Baldridge_5_.*_R2_.*" = "^Baldridge_7_.*_R2_.*",
                   "^EH_Baldridge_6_.*_R1_.*" = "^Baldridge_8_.*_R1_.*",
                   "^EH_Baldridge_6_.*_R2_.*" = "^Baldridge_8_.*_R2_.*",
                   "^EH_Baldridge_7_.*_R1_.*" = "^Baldridge_9_.*_R1_.*",
                   "^EH_Baldridge_7_.*_R2_.*" = "^Baldridge_9_.*_R2_.*",
                   "^EH_Baldridge_8_.*_R1_.*" = "^Baldridge_15_.*_R1_.*",
                   "^EH_Baldridge_8_.*_R2_.*" = "^Baldridge_15_.*_R2_.*",
                   "^EH_Baldridge_9_.*_R1_.*" = "^Baldridge_16_.*_R1_.*",
                   "^EH_Baldridge_9_.*_R2_.*" = "^Baldridge_16_.*_R2_.*",
                   "^EH_Baldridge_10_.*_R1_.*" = "^Baldridge_17_.*_R1_.*",
                   "^EH_Baldridge_10_.*_R2_.*" = "^Baldridge_17_.*_R2_.*",
                   "^EH_Baldridge_13_.*_R1_.*" = "^Baldridge_21_.*_R1_.*",
                   "^EH_Baldridge_13_.*_R2_.*" = "^Baldridge_21_.*_R2_.*",
                   "^EH_Baldridge_14_.*_R1_.*" = "^Baldridge_22_.*_R1_.*",
                   "^EH_Baldridge_14_.*_R2_.*" = "^Baldridge_22_.*_R2_.*",
                   "^EH_Baldridge_16_.*_R1_.*" = "^Baldridge_28_.*_R1_.*",
                   "^EH_Baldridge_16_.*_R2_.*" = "^Baldridge_28_.*_R2_.*",
                   "^EH_Baldridge_17_.*_R1_.*" = "^Baldridge_29_.*_R1_.*",
                   "^EH_Baldridge_17_.*_R2_.*" = "^Baldridge_29_.*_R2_.*",
                   "^EH_Baldridge_18_.*_R1_.*" = "^Baldridge_30_.*_R1_.*",
                   "^EH_Baldridge_18_.*_R2_.*" = "^Baldridge_30_.*_R2_.*",
                   "^EH_Baldridge_19_.*_R1_.*" = "^Baldridge_31_.*_R1_.*",
                   "^EH_Baldridge_19_.*_R2_.*" = "^Baldridge_31_.*_R2_.*",
                   "^EH_Baldridge_20_.*_R1_.*" = "^Baldridge_33_.*_R1_.*",
                   "^EH_Baldridge_20_.*_R2_.*" = "^Baldridge_33_.*_R2_.*")

ehFilePatterns <- names(patternMap)

for (i in 1:length(ehFilePatterns)) {
  currentEHPattern <- ehFilePatterns[i]
  ehFile <- allFilesList[grepl(currentEHPattern, allFilesList)]
  
  currentRumPattern <- patternMap[[currentEHPattern]]
  rumFile <- allFilesList[grepl(currentRumPattern, allFilesList)]
  
  # Get name for output file from rum file
  outFilePrefix <- str_extract(rumFile, "^Baldridge_[:digit:]{1,2}_")
  outFileSuffix <- str_extract(rumFile, "_R[1,2]")
  outFileName <- paste0(outFilePrefix, "concatenated", outFileSuffix, ".fastq.gz")
  
  filesToConcat <- c(ehFile, rumFile)
  fileOut <- paste0("./concatenated_files/", outFileName)
  
  # Talk to me
  print(paste0("Combining ", ehFile, " and ", rumFile, " into ", fileOut))
  
  for(j in 1:length(filesToConcat)) {
    fileIn <- readFastq(dirPath = "./raw_data/", pattern = filesToConcat[j])
    writeFastq(fileIn, fileOut, mode = "a", compress = TRUE)
  }
}



#firstEHFilePattern <- ehFilePatterns[1]
#firstEHFile <- allFilesList[grepl(firstEHFilePattern, allFilesList)]

#firstRumFilePattern <- patternMap[[firstEHFilePattern]]
#firstRumFile <- allFilesList[grepl(firstRumFilePattern, allFilesList)]

#outFilePrefix <- str_extract("Baldridge_17_SIC_index_0496_SIC_index_0545_TCACTACAG_TATCAGCAG_S17_L001_R2_001.fastq.gz",
                             "^Baldridge_[:digit:]{1,2}_")
#outFileSuffix <- str_extract(firstRumFile, "_R[1,2]")
#outFileName <- paste0(outFilePrefix, "concatenated", outFileSuffix, ".fastq.gz")

#filesToConcat <- c(firstEHFile, firstRumFile)
#fout <- "./concatenated_files/testConcat.fastq.gz"

#for (i in 1:length(filesToConcat)) {
#  fq <- readFastq(dirPath = "./raw_data/", pattern = filesToConcat[i])
#  writeFastq(fq, fout, mode = "a", compress = TRUE)
#}

#print(paste0("Combining ", firstEHFile, "and ", firstRumFile))

