# ScanCoverageValues.R

# Check coverage files generated with bedtools genomecov -d to see if any
#   baldridge_rumspringa files have coverage of at least 200 at each position

library("data.table")
library("dplyr")
library("stringr")

genomeCovFiles <- list.files("../baldridge_rumspringa_coverageInfo/",
                             pattern = "_genomeCov.txt",
                             full.names = TRUE)

fileCoverage <- vector(mode = "character", length = length(genomeCovFiles))

#----- For positions 500 to 6500 -----#

for (i in 1:length(genomeCovFiles)) {
  currentFile <- genomeCovFiles[[i]]
  currentFileContents <- fread(currentFile,
                               col.names = c("genome", "position", "depth"),
                               skip = 499,
                               nrows = 6001)
  goodCoverage <- all(currentFileContents$depth >= 200)
  fileCoverage[i] <- goodCoverage
  names(fileCoverage)[i] <- basename(currentFile)
}

#-------------------------------------#

resultsDF <- data.frame(fileCoverage)
resultsDF <- resultsDF %>%
  tibble::rownames_to_column(var = "Coverage_File") %>%
  mutate("Sample" = str_remove(string = Coverage_File, pattern = "_genomeCov.txt")) %>%
  select("Sample", "fileCoverage") %>%
  setnames("fileCoverage", "Meets_Criteria")

write.table(resultsDF, file = "../baldridge_rumspringa_coverageInfo/CoverageResults.txt",
            quote = FALSE, row.names = FALSE, sep = "\t")


#----- For all positions -----#
for (i in 1:length(genomeCovFiles)) {
  currentFile <- genomeCovFiles[[i]]
  currentFileContents <- fread(currentFile,
                               col.names = c("genome", "position", "depth"))
  goodCoverage <- all(currentFileContents$depth >= 200)
  fileCoverage[i] <- goodCoverage
  names(fileCoverage)[i] <- basename(currentFile)
}
#---------------------------#