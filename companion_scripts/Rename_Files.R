# Rename_Files.R

# This script is used to rename files in the raw_data directory.  These raw data
#   files consist of both original baldridge_rumspringa files and 
#   baldridge_rumspringa files that have been concatenated with re-sequenced
#   EH files.

library("stringr")
library("ShortRead")

allFiles <- list.files(path = "./raw_data", pattern = "*.fastq.gz")

for (i in 1:length(allFiles)) {
  currentFile <- allFiles[i]
  
  originalPrefix <- str_extract(string = currentFile,
                                pattern = "^Baldridge_[:digit:]{1,2}_")
  
  newPrefix <- str_remove_all(string = originalPrefix, pattern = "[Baldridge_, _]")
  
  # Do we need to add a 0 to the new prefix?
  if (nchar(newPrefix) == 1) {
    newPrefix <- paste0("0", newPrefix)
  }
  
  originalSuffix <- str_extract(string = currentFile,
                                pattern = "_R[1,2][_, .]")
  newSuffix <- str_remove_all(string = originalSuffix, pattern = "[_, .]")
  
  newFileName <- paste0(newPrefix, "_", newSuffix, ".fastq.gz")
  
  print(paste0("Renaming ", currentFile, " to ", newFileName))
  
  # Write out the current file with its new name
  fileIn <- readFastq(dirPath = "./raw_data/", pattern = currentFile)
  fileOut <- paste0("./raw_data/renamed/", newFileName)
  writeFastq(fileIn, fileOut, mode = "w", compress = TRUE)
}



