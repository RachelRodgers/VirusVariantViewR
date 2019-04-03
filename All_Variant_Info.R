# All_Variant_Info.R

# For each sample in the baldridge_rumspringa data set, write large table:
#   Sample Name
#   Average Genome Coverage
#   Mutations: Position, Ref, Alt, Qual, AD, DP

library("tidyverse")

# Get all sample names:
sampleNamesLong <- read.delim(file = "../baldridge_rumspringa_lookup.txt",
                              header = FALSE, stringsAsFactors = FALSE)
sampleNames <- str_extract(sampleNamesLong$V1,
                           pattern = "^Baldridge_[:digit:]{1,2}")

dfList <- list()

for (i in 1:length(sampleNames)) {
  
  currentSample <- sampleNames[i]
  
  #----- Average Coverage -----#
  columnNames <- c("chromosome", "depth", "number_of_bases",
                   "chromosome_size", "fraction_of_bases")
  coverageFile <- read.delim(paste0("../baldridge_rumspringa/sample_data/genome_coverage/",
                                    currentSample, "_coverage.txt"), 
                             header = FALSE,
                             col.names = columnNames,
                             stringsAsFactors = FALSE)
  coverageFile <- coverageFile %>%
    dplyr::select(depth, number_of_bases) %>%
    dplyr::mutate(product = depth * number_of_bases)
  
  averageCoverage <- round(sum(coverageFile$product)/7383, digits = 0)
  
  #----- Mutation Info -----#
  mutationFile <- read.delim(file = paste0("../baldridge_rumspringa/variants/annotated_variants/",
                                           currentSample, "_variants_annotated.txt"),
                             header = TRUE, stringsAsFactors = FALSE,
                             colClasses = "character")
  mutationFileModified <- mutationFile %>%
    select(Position, Reference, Alternative, Quality, Info, Values) %>%
    mutate("Allelic Depth" = map_chr(.x = str_split(Values, pattern = ","),
                                     .f = tail, n = 1),
           "Total Depth" = str_remove(str_extract(Info, "DP=[:digit:]+"),
                                      "DP="),
           "Allelic Frequency (%)" = 
             round((100 * as.numeric(`Allelic Depth`)/as.numeric(`Total Depth`)), 
                   digits = 2)) %>%
    select(-c("Info", "Values"))
  
  #----- Combine Info -----#
  allInfo <- data.frame("Sample" = currentSample,
                        "Average Coverage" = averageCoverage,
                        mutationFileModified)
  
  dfList[[i]] <- allInfo
  
}

# Merge all the lists in dfList

fullTable <- Reduce(f = function(df1, df2) 
  {rbind(x = df1, y = df2, make.row.names = FALSE)}, x = dfList)

# Save
write.table(fullTable, file = "../sample_info_baldridge_rumspringa.txt",
            quote = FALSE, sep = "\t", row.names = FALSE)

