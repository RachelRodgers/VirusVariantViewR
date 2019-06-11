# Rename_Files_Again.R

# Ebrahim has a different re-numbering scheme he wants me to use.  In this script
#   I'll be re-naming the concatenated and unconcatenated baldridge_rumspringa
#   files I previously re-named in the script Rename_Files.R

library("stringr")
library("ShortRead")

allFiles <- list.files(path = "../combined_dataset/raw_data/", 
                       pattern = "*.fastq.gz")

numberingMap <- list("01" = "02", 
                     "02" = "03", 
                     "03" = "04", 
                     "04" = "05", 
                     "05" = "06", 
                     "06" = "07", 
                     "07" = "08",
                     "08" = "09", 
                     "09" = "10", 
                     "10" = "11",
                     "11" = "12", 
                     "12" = "13", 
                     "13" = "14", 
                     "14" = "15", 
                     "15" = "16",
                     "16" = "17", 
                     "17" = "18",
                     "18" = "20",
                     "19" = "21",
                     "20" = "23",
                     "21" = "24",
                     "22" = "25",
                     "23" = "26",
                     "24" = "27",
                     "25" = "28",
                     "26" = "29",
                     "27" = "30",
                     "28" = "32",
                     "29" = "33",
                     "30" = "34",
                     "31" = "35",
                     "32" = "36",
                     "33" = "37")

for (i in 1:length(allFiles)) {
  oldFileName <- allFiles[[i]]
  currentNumber <- str_extract(oldFileName, pattern = "^[:digit:]{2}")
  newNumber <- numberingMap[[currentNumber]]
  newFileName <- str_replace(oldFileName, pattern = currentNumber, 
                             replacement = newNumber)
  # Write out the current file with its new name
  fileIn <- readFastq(dirPath = "../combined_dataset/raw_data/", pattern = oldFileName)
  fileOut <- paste0("../combined_dataset/raw_data/renamed/", newFileName)
  writeFastq(fileIn, fileOut, mode = "w", compress = TRUE)
}

# Manually moving over the old EH samples to re-name.
