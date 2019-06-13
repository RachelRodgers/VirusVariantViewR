# Info_6595.R

# In this script I am compiling the following information for each sample in 
#   the Combined_Data set:
# 1. sample name
# 2. reference allele @ position 6595 (will always be the same)
# 3. alternative allele @ position 6595
# 4. raw read depth @ position 6595
# 5. allelic depth @ position 6595

library("stringr")
library("tidyverse")

vcfFiles <- list.files("../Combined_Data/full_VCF_files/",
                       pattern = ".vcf", full.names = TRUE)

vcfHeaders <- c("Reference Genome", "Position", 
                "ID", "Reference", "Alternative",
                "Quality", "Filter", "Info", "Format", "Values")

outputDF <- data.frame(matrix(ncol = 5, nrow = length(vcfFiles)))
outputColNames <- c("sample", "reference", "alternative", "raw_read_depth", "allelic_depth")
colnames(outputDF) <- outputColNames

for (i in 1:length(vcfFiles)) {
  print(paste0("i = ", i))
  # set flag to determine if file is empty
  currentVCFFile <- vcfFiles[[i]]
  currentSample <- str_extract(basename(currentVCFFile),
                               "^[:digit:]{2}")
  outputDF[i, "sample"] <- currentSample
  
  # Read file:
  vcfContents <- tryCatch({
    read.delim(currentVCFFile, comment.char = "#",
               header = FALSE, colClasses = "character")
  },
  error = function(e) {
    print(paste0("No data in sample ", currentSample))
    vcfContents <- NULL
  })
  
  if(is.null(vcfContents)) {
    outputDF[i, 2:5] <- "empty file"
    next()
  }
  
  names(vcfContents) <- vcfHeaders
  
  infoPos6595 <- filter(vcfContents, Position == "6595") # keep just the row of interest 
  
  # Do something the row is not found
  if(nrow(infoPos6595) == 0) {
    print(paste0("sample ", currentSample, " is empty at position 6595"))
    outputDF[i, 2:5] <- "no coverage at 6595"
    next()
  }
  
  reference <- infoPos6595$Reference
  outputDF[i, "reference"] <- reference
  
  rawDepth <- infoPos6595$Info %>% str_extract("DP=[:digit:]+") %>% str_remove("DP=")
  outputDF[i, "raw_read_depth"] <- rawDepth
  
  alternative <- infoPos6595$Alternative
  outputDF[i, "alternative"] <- alternative
  
  if (alternative == ".") {
    allelicDepth <- NA
  } else {
    allelicDepth <- infoPos6595$Values %>%
      map_chr(.x = str_split(string = ., pattern = ","), .f = tail, n = 1)
  }
  outputDF[i, "allelic_depth"] <- allelicDepth
}

write.table(outputDF, file = "../Combined_Data/Info_6595.txt",
            quote = FALSE, sep = "\t", row.names = FALSE,
            col.names = TRUE)







