install.packages(c("sf", "dplyr", "dbscan", "dtw", "ggplot2", "tidyr", "zoo", "lubridate", "geosphere","mlogit", "devtools"))
library(sf)
library(dplyr)
library(dbscan)
library(dtw)
library(ggplot2)
library(mlogit)
library(InformationValue)
library(tidyr)
library(zoo)
library(lubridate)
library(stringr)
library(geosphere)
library(purrr)
#------------------------------------------------------------------------------------------------------------------------------
#This is the well data cleanup section:


df <- read.csv("C:/Users/epmos/OneDrive - Washington State University (email.wsu.edu)/z_gradschool/fall 2024/CPT_575/proposal/data/GWL_data.csv")

# Pivot the table from wide format to long format
long_df <- df %>%
  pivot_longer(cols = -Dates, # all columns except 'year' are pivoted
               names_to = "well_name", # new column for well names
               values_to = "measurement") # new column for measurement values
#keep only the year
long_df <- long_df %>%
  mutate(year = str_sub(Dates, -4, -1)) # Extract last 4 characters
#keep only measurement rows
unique_long_df <- long_df %>%
  filter(!is.na(measurement) & measurement != "")

# Reorder columns for clarity (optional)
final_long_df <- unique_long_df %>%
  select(well_name, year, measurement)

#---------------------------------------------------------------------------------------------------------------

#This is creating a dataframe of each well within 5 miles of every well on the list. 
#The end product has two columns, one with the well of interest and 
#the second with the wells within 5 miles of it. all rows should be a unique combination.

wells_df <- read.csv("C:/Users/epmos/OneDrive - Washington State University (email.wsu.edu)/z_gradschool/fall 2024/CPT_575/proposal/data/Lat-Lon.csv")

# Calculate pairwise distances between wells using the Haversine formula
# distm() calculates distance in meters between two lat/lon points
distance_matrix <- distm(
  wells_df[, c("Longitude", "Latitude")], 
  wells_df[, c("Longitude", "Latitude")], 
  fun = distHaversine
)

# Convert distance matrix into a long format data frame
distance_df <- as.data.frame(as.table(distance_matrix))
colnames(distance_df) <- c("well_of_interest", "nearby_well", "distance_meters")


# Add well names to the data frame (convert row/column indices to well names)
distance_df <- distance_df %>%
  mutate(
    well_of_interest = wells_df$Location_ID[as.numeric(well_of_interest)],
    nearby_well = wells_df$Location_ID[as.numeric(nearby_well)]
  )

# Filter for wells within 5 miles (1 mile = 1609.34 meters)

filtered_distance_df <- distance_df %>%
  filter(distance_meters > 0 & distance_meters <= 5 * 1609.34) # Exclude self-distances and keep distances <= 5 miles

#Remove duplicate pairs to keep unique well-to-well combinations
filtered_duplicate_distance_df <- filtered_distance_df %>%
  rowwise() %>%
  mutate(pair_id = paste(sort(c(well_of_interest, nearby_well)), collapse = "_")) %>%
  distinct(pair_id, .keep_all = TRUE) %>%
  select(well_of_interest, nearby_well)

#-------------------------------------------------------------------------------------------------------------------------
#This step takes each pair of wells that are within 5 miles of each other and compares their 
#hydrographs to see if they have similar behavior:

#this needs adjusted to account for slight name diff
wells_of_interest <- unique(c(distance_df$well_of_interest, distance_df$nearby_well))

# Filter water level measurements for wells that are within 5 miles
filtered_df <- final_long_df %>%
  filter(well_name %in% wells_of_interest)
#clean up well names to be consistant
final_long_df$well_name <- trimws(toupper(final_long_df$well_name))
distance_df$well_of_interest <- trimws(toupper(distance_df$well_of_interest))
distance_df$nearby_well <- trimws(toupper(distance_df$nearby_well))
#substitue . for -
final_long_df$well_name <- gsub("\\.", "-", final_long_df$well_name)
#interpolate missing values with adjacent values
final_long_df$measurement <- zoo::na.approx(final_long_df$measurement, rule = 2)
#remove double values
final_long_df <- final_long_df %>%
  group_by(well_name, year) %>%
  summarize(measurement = mean(measurement, na.rm = TRUE), .groups = 'drop')



#Step 3: Identify unique well combinations (to avoid duplicates like Well1-Well2 and Well2-Well1)
# --------------------------------------------------


# Create a "pair_id" that allows sorting of well pairs

unique_pairs <- filtered_duplicate_distance_df %>%
  rowwise() %>%
  mutate(pair_id = paste(sort(c(well_of_interest, nearby_well)), collapse = "_")) %>%
  distinct(pair_id, .keep_all = TRUE) %>%
  select(well_of_interest, nearby_well)


#Create a loop to compute First-Difference and Pearson Correlation Coefficient for each pair
results_list <- list()


#Filter unique_pairs to include only wells in final_long_df

filtered_pairs <- unique_pairs %>%
  filter(
    well_of_interest %in% final_long_df$well_name &
      nearby_well %in% final_long_df$well_name
  )


results_list <- list()  # Create an empty list to store results

#loop through each pair

for (i in 1:nrow(filtered_pairs)) {
  well_1 <- filtered_pairs$well_of_interest[i]
  well_2 <- filtered_pairs$nearby_well[i]
  
  # Filter and pivot the data for the two wells
  pair_data <- final_long_df %>%
    filter(well_name %in% c(well_1, well_2)) %>%
    pivot_wider(names_from = well_name, values_from = measurement)
  
  # Check if the required columns exist
  if (!(well_1 %in% colnames(pair_data)) || !(well_2 %in% colnames(pair_data))) {
    next  # Skip this pair if either well column is missing
  }
  
  # Flatten list-columns to numeric if necessary
  pair_data <- pair_data %>%
    mutate(across(-year, ~ as.numeric(unlist(.))))
  
  # Interpolate missing values
  pair_data <- pair_data %>%
    mutate(across(-year, ~ zoo::na.approx(., rule = 2)))
  
  # Apply first-difference transformation
  diff_data <- pair_data %>%
    mutate(across(-year, ~ c(NA, diff(.)))) %>%
    drop_na()  # Drop rows with NA values from the first difference
  
  # Check if there are at least 2 finite (non-NA) values in both columns
  if (sum(is.finite(diff_data[[well_1]])) < 2 || sum(is.finite(diff_data[[well_2]])) < 2) {
    results_list[[i]] <- data.frame(
      well_1 = well_1,
      well_2 = well_2,
      pearson_correlation = NA,
      p_value = NA
    )
    next  # Skip this pair
  }
  
  # Check if either well has zero variance
  if (sd(diff_data[[well_1]], na.rm = TRUE) == 0 || sd(diff_data[[well_2]], na.rm = TRUE) == 0) {
    results_list[[i]] <- data.frame(
      well_1 = well_1,
      well_2 = well_2,
      pearson_correlation = NA,
      p_value = NA
    )
    next  # Skip this pair
  }
  
  # Check if all differences are zero for either well
  if (all(diff_data[[well_1]] == 0, na.rm = TRUE) || all(diff_data[[well_2]] == 0, na.rm = TRUE)) {
    results_list[[i]] <- data.frame(
      well_1 = well_1,
      well_2 = well_2,
      pearson_correlation = NA,
      p_value = NA
    )
    next  # Skip this pair
  }
  
  # Calculate Pearson Correlation Coefficient and P-Value
  cor_test_result <- cor.test(diff_data[[well_1]], diff_data[[well_2]], use = "pairwise.complete.obs")
  
  # Extract the Pearson correlation and p-value
  correlation_value <- cor_test_result$estimate  # Pearson correlation coefficient
  p_value <- cor_test_result$p.value             # p-value for the test
  
  # Add results to the list
  results_list[[i]] <- data.frame(
    well_1 = well_1,
    well_2 = well_2,
    pearson_correlation = correlation_value,
    p_value = p_value
  )
}

#Combine all results into a single DataFrame
final_results <- bind_rows(results_list)


anom <- final_results %>%
  filter(final_results$p < 0.05)



#CREATING TEST AND TRAIN DATA


# Read the data
anom = read.csv("C:/Users/epmos/OneDrive - Washington State University (email.wsu.edu)/z_gradschool/fall 2024/CPT_575/proposal/data/anom.csv")

# Prepare data
missing.anom = anom[is.na(anom$dist_geofeature),]
nonmissing.anom = anom[!(is.na(anom$dist_geofeature)),]

# Check representativeness
describe(anom)
describe(nonmissing.anom)

#make this reproducible
set.seed(1)

#Use 70% of dataset as training set and remaining 30% as testing set
sample = sample(c(TRUE, FALSE), nrow(nonmissing.anom), replace=TRUE, prob=c(0.7,0.3))
train = nonmissing.anom[sample, ]
test = nonmissing.anom[!sample, ]  
#FITTING THE MODEL
# Logistic Regression
n = nrow(nonmissing.anom)
model1 = glm(anom ~ dist_geofeature, family = "binomial", data = train,)
summary(model1)

# Odds ratios and confidence intervals
exp(model1$coefficients)
exp(confint(model1))

# Model fit
pscl::pR2(model1)["McFadden"]

# Check assumptions
car::vif(model1)
1 / car::vif(model1)
print("only one variable")

#predict probability of anomolies
test_predicted = predict(model1, test, type="response")
test_predicted
# Linearity of the logit
logit(test_predicted)
nonmissing.anom$logxDist_geofeature = log(nonmissing.anom$dist_geofeature) * nonmissing.anom$dist_geofeature

optimal <- optimalCutoff(test$anom, test_predicted)[1]
optimal
confusionMatrix(test$anom, test_predicted)
#MODEL DIAGNOSTICS

#calculate true positive rate
sensitivity(test$anom, test_predicted)

#calculate specificity
specificity(test$anom, test_predicted)

#calculate total misclassification error rate
misClassError(test$anom, test_predicted, threshold = optimal)

plotROC(test$anom, test_predicted)


