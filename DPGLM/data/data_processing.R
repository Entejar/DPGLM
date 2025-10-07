rm(list = ls())
# Load Libraries
source("load.R")

## Read in raw data
hustad <- read.csv("data/intelligibility-data-for-paul.csv")
summary(hustad)
names(hustad)
for (i in c(1,3,6,7)) {
  hustad[[i]] <- factor(hustad[[i]])
}
summary(hustad)

## Subset TD data
hustadTD <- data.frame(hustad[hustad$group=="TD",])
ooo <- order(hustadTD$age)
hustadTD <- hustadTD[ooo,]
rm(ooo)
summary(hustadTD)
hustadTD <- hustadTD[, c(5,6,7,8)]

## Subset TD Data: Multi Word, ordered y --------------------------------------------------

# Convert to a tibble using tidyr
hustadTDMW <- as_tibble(hustadTD)

# Filter rows where 'intelligibility_type' equals "multiword"
hustadTDMW <- hustadTDMW[hustadTDMW$intelligibility_type == "multiword", ]

# Remove the 'intelligibility_type' column
hustadTDMW$intelligibility_type <- NULL

# Arrange the data by 'mean_intelligibility'
hustadTDMW <- hustadTDMW[order(hustadTDMW$mean_intelligibility), ]

saveRDS(hustadTDMW, "data/hustadTDMW.rds")

# Convert to a tibble using tidyr
hustadTDSW <- as_tibble(hustadTD)

# Filter rows where 'intelligibility_type' equals "multiword"
hustadTDSW <- hustadTDSW[hustadTDSW$intelligibility_type == "single-word", ]

# Remove the 'intelligibility_type' column
hustadTDSW$intelligibility_type <- NULL

# Arrange the data by 'mean_intelligibility'
hustadTDSW <- hustadTDSW[order(hustadTDSW$mean_intelligibility), ]

saveRDS(hustadTDSW, "data/hustadTDSW.rds")