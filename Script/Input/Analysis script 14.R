start_time <- Sys.time()

library(tidyverse)
library(data.table)
library(plotrix)
library(zoo)
library(reshape2)
library(broom)
library(rstudioapi)
library(gtools)
library(ggplot2)
library(ggpubr)
library(mgcv)
library(emmeans)
library(gamlss)
library(MASS)
library(quantreg)
library(openxlsx)
library(PMCMRplus)
library(car)
library(rstatix)
library(parallel)
library(Rmpfr)


####!!!! Step-By-Step Instructions !!!!####

#!!It is strongly recommended that you make a copy of this file before
#making any edits!!

#### Step 1: Make sure that this R file is in the same directory as the folders that contain each replicate of your experiment
#Note: !This step is required!

#### Step 2: Input how your plate is organized ####
#Note: !This section is required!

#Note: If there are wells that have no fish, put "empty" in the for the wells description.
#Make sure that it is all lower case.

conds <- data.frame(
  c(	"WT+DMSO",	"WT+DMSO",	"WT+DMSO",	"WT+DMSO",	"WT+DMSO",	"WT+DMSO",	"WT+DMSO",	"WT+DMSO"	), #column 1
  c(	"WT+DMSO",	"WT+DMSO",	"WT+DMSO",	"WT+DMSO",	"WT+DMSO",	"WT+DMSO",	"WT+DMSO",	"WT+DMSO"	), #column 2
  c(	"Mut1",	"Mut1",	"Mut1",	"Mut1",	"Mut1",	"Mut1",	"Mut1",	"Mut1"	), #column 3
  c(	"Mut1",	"Mut1",	"Mut1",	"Mut1",	"Mut1",	"Mut1",	"Mut1",	"Mut1"	), #column 4
  c(	"Mut2",	"Mut2",	"Mut2",	"Mut2",	"Mut2",	"Mut2",	"Mut2",	"Mut2"	), #column 5
  c(	"Mut2",	"Mut2",	"Mut2",	"Mut2",	"Mut2",	"Mut2",	"Mut2",	"Mut2"	), #column 6
  c(	"Mut3",	"Mut3",	"Mut3",	"Mut3",	"Mut3",	"Mut3",	"Mut3",	"Mut3"	), #column 7
  c(	"Mut3",	"Mut3",	"Mut3",	"Mut3",	"Mut3",	"Mut3",	"Mut3",	"Mut3"	), #column 8
  c(	"Mut4",	"Mut4",	"Mut4",	"Mut4",	"Mut4",	"Mut4",	"Mut4",	"Mut4"	), #column 9
  c(	"Mut4",	"Mut4",	"Mut4",	"Mut4",	"Mut4",	"Mut4",	"Mut4",	"Mut4"	), #column 10
  c(	"MUT+DMSO",	"MUT+DMSO",	"MUT+DMSO",	"MUT+DMSO",	"MUT+DMSO",	"MUT+DMSO",	"MUT+DMSO",	"MUT+DMSO"	), #column 11
  c(	"MUT+DMSO",	"MUT+DMSO",	"MUT+DMSO",	"MUT+DMSO",	"MUT+DMSO",	"MUT+DMSO",	"MUT+DMSO",	"MUT+DMSO"	) #column 12
)
colnames(conds) <- c(1:12)
rownames(conds) <- LETTERS[1:8]
#This plate reads from row A to row H, left to right, with sample 1 being A1 and sample 2 being A2 and so on.
#each c() is a column in the 96 well plate, and the first column here is column 1
#so you have to mentally flip and this array along its horizontal axis and rotate it 90 degrees clockwise for it to have 
#the same orientation as a 96-well plate.
#highlight and run the 'conds' object to see the structure and well-condition assignments of your plate

#-------------------------------------------------------------------------------------------------------------------------------------

#### Step 3: Declare the times in seconds for each phase of your experiment ####
#Note: !This section is required!

#!!!Note: All of the phases labeled similarly, e.g. "Dark" will be concatenated and treated as one phase, and compared against phases 
#with other names. The 'Group' column refers to the grouping of steady-state phases that will be compared 
#directly, e.g. 'Light' and 'Dark' is group 1, '600Hz' and '200Hz' are group 2, etc. 

#The transition phases are meant to allow for analysis of the parts of the phases where activity spikes during the phase transitions
#cause the distributions of pixel differences to become skewed. The conditions will be compared within the transition phases, but 
#the transition phases will not be compared to the steady-state phases or to other transition phases. 
#The transition phases are optional and usually overlap the junctions of adjacent phases from the beginning to the ends of any 
#movement spikes. This may take some trial and error to find the best breaking point between a steady-phase and
#an adjacent transition phase.
#The steady-state state phases occur after movement has settled after a spike, within a phase. These are referred to as the standard phase
#descriptions (e.g. "Dark", "Light", "600Hz", "200Hz", etc.)

#Example:
#It is necessary to fill all columns. Case for the Trans column is important.
#time_frame <- tribble(
#    ~Group,  ~Trans,   ~Phase,    ~Start,   ~End,
#  #|-------|--------|-----------|--------|--------|
#       1,    FALSE,     "Dark",      1,       60,
#       1,     TRUE,  "Trans 1",     45,       85, 
#       1,    FALSE,    "Light",     61,      120,
#       1,     TRUE,  "Trans 2",    105,      145,
#       1,    FALSE,     "Dark",    121,      180,
#       2,    FALSE,    "600Hz",    181,      210,
#       2,    FALSE,    "200Hz",    211,      240
#)

time_frame <- tribble(
  ~Group,  ~Trans,   ~Phase,    ~Start,   ~End,
  #|-------|--------|-----------|--------|--------|
  1,    FALSE,     "Dark",     1,       56,
  1,     TRUE,  "Trans 1",    57,       72, 
  1,    FALSE,    "Light",    73,      115,
  1,     TRUE,  "Trans 2",   116,      145,
  1,    FALSE,     "Dark",   146,      180,
  2,    FALSE,    "600Hz",   181,      210,
  2,    FALSE,    "200Hz",   211,      240
)

#-------------------------------------------------------------------------------------------------------------------------------------

#### Step 4: Determine the p-value adjustment methods that will be used for multiple hypotheses testing. Descriptions of p-value adjustment methods used in this script:
#Note: !This section is required!

#"BH" (Benjamini-Hochberg): This controls the False Discovery Rate, to control the proportion of Type I errors (false positive, 
#incorrect rejections of the null hypothesis) among all the hypotheses that are identified as significant. In other words, 
#if the experimental design is exploratory, the results can be confirmed in a later study, or a small number of 
#false positives (percentage less than your alpha value, or critical p-value threshold) are expected and can be tolerated. This is more
#common in studies where there are a large number of groups or phases being compared. 
#p-values will likely be smaller with this adjustment and will scale to the number of adjustments made.

#"holm": This controls the Family-Wise Error Rate, or the probability of making even one Type I error (false positive) across all tests. 
#This method is usually used when the cost of having even one false positive is high and is less tolerable. This is more common in 
#studies with less groups or phases being compared. 

#For p-value adjustment in the figures comparing movement between groups of animals within one phase
#(Differences in movement of Mut 1 vs Control within the Dark phase), decide the appropriate method
#Example:
#phase_adjustment_method <- "BH"

phase_adjustment_method <- "holm"

#For p-value adjustment in the figures comparing the movement of each condition between multiple phases of the same group
#(Difference in Mut 1 movement between Dark and Light phases), decide the appropriate method
#Example:
#group_adjustment_method <- "holm"

group_adjustment_method <- "holm"

#-------------------------------------------------------------------------------------------------------------------------------------

#### Step 5: Declare all of the controls being used in these experiments ####
#Note: !This section is required!
#Example: controls <- c("WT+DMSO", "MUT+DMSO")

controls <- c("WT+DMSO", "MUT+DMSO")

#-------------------------------------------------------------------------------------------------------------------------------------

#### Step 6: Declare which Condition all of your other Conditions will be normalized against, e.g. WT controls ####
#Note: This value needs to be updated with whatever Condition your experimental Conditions are being normalized against, e.g. 'WT' ####
#Note: !This section is required!
#Example: normalizing_control <- c("MUT+DMSO")

normalizing_control <- c("MUT+DMSO")

#-------------------------------------------------------------------------------------------------------------------------------------

#### Step 7: Optional: Choose which colors, if any, that you would like to represent your Conditions in figures 1-3
#Un-comment the lines below if you want to set your own colors for your figures.The variable must be named 'Condition_colours'.
#The left side argument is the Condition names, exactly as it appears in the 'conds' table, and
#the right side argument is the color of choice.
#Example: 
#Condition_colours = c("WT+DMSO"="red",
#                      "MUT+DMSO"="blue",
#                      "Mut1" = "green",
#                      "Mut2" = "orange")


######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
#Below here should only be edited if you are comfortable with programming and have read the documentation for this script

#This sets the current working directory to be whatever directory the file is
#located in
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
getwd()

#The Zebrabox will give several csv files per replicate. This will sort throught each folder, which contains the csv files of each 
#replicate, and add one column per replicate the Recursive function is to add the data from the different csv files, per replicate, 
#into a data column. So each replicate must be divided into each own folder.
#Get a vector for paths to all the sub-directories
observations_dirs <- list.dirs(path = dirname("."), full.names = TRUE, recursive = FALSE)

# Check if "./Output" is in the vector and remove it if present
if ("./Output" %in% observations_dirs) {
  observations_dirs <- observations_dirs[observations_dirs != "./Output"]
}

#Search recursive through each sub-directory to get a vector of paths of all .csv files
for (d in 1:length(observations_dirs)) {
  file_vec <- list.files(observations_dirs[d], pattern=".csv$", full.names=TRUE, recursive=TRUE)
  
  # Create empty list in which to add data.tables
  res_list <- list()
  
  #### Create a list of all the .csv files in the selected directory ####
  # In a loop, load and remove unwanted rows and columns from each .csv
  # and combine them into one variable, which will be combined into one data.table
  for (i in seq_along(file_vec)) {
    file_path <- file_vec[i]
    tmp <- fread(file_path, na.strings = c("NA", "")) # Load count file in as data.table named tmp.
    set(tmp, j=c("abstime", "type", "area", "data2", # Remove unnecessary columns
                 "data3", "data4", "data5", "data6",
                 "data7", "data8", "data9", "data10",
                 "data11", "data12", "data13", "data14",
                 "data15", "data16"), value=NULL)
    tmp <- tmp[!(data1 == "SOUND")]   # Remove non data rows, check your data for others
    tmp <- tmp[!(data1 == "BACK_LIGHT")] # Remove non data rows
    tmp <- tmp[!(data1 == "TOP_LIGHT")] # Remove non data rows
    tmp <- tmp[!(data1 == "Not identify")] # Remove non data rows
    tmp$data1 <- as.integer(tmp$data1) # Set values to integer
    
    res_list[[i]] <- tmp
  }
  
  #Create the full_tab with all columns, including data1
  if (d == 1) {
    full_tab <- rbindlist(res_list)
  }
  
  #Add only the data1 column from the other replicates
  else if (d >= 2) {
    consec_data <- rbindlist(res_list)
    
    row_num_diff <- abs(nrow(full_tab) - nrow(consec_data))
    
    if (nrow(full_tab) > nrow(consec_data)) {
      full_tab <- slice(full_tab, 1:(n() - row_num_diff)) 
    }
    
    else if(nrow(full_tab) < nrow(consec_data)) {
      consec_data <- slice(consec_data, 1:(n() - row_num_diff)) 
    }
    
    full_tab <- full_tab[, paste0('data', d) := consec_data$data1]
  }
}

#Create a folder for output, if they do not already exist
# Define the path of the new folder
figures_folder_path <- file.path(".", "Output", "Whole experiment", "Figures")

# Check if the folder already exists
if (!dir.exists(figures_folder_path)) {
  # Create the folder if it does not exist
  dir.create(figures_folder_path, recursive = TRUE)
}

data_folder_path <- file.path(".", "Output", "Whole experiment", "Data analysis")

if (!dir.exists(data_folder_path)) {
  # Create the folder if it does not exist
  dir.create(data_folder_path, recursive = TRUE)
}

#### Write the long form data into one data.table ####
fwrite(full_tab, file = file.path(data_folder_path, "long format of pixel differences grouped by time and location.txt"), sep = "\t")

# Get the column names that contain "data"
data_cols <- grep("data", names(full_tab), value = TRUE)

# Loop through each data column
for (c in data_cols) {
  mtab <- dcast(full_tab, formula = time ~ location, value.var = c)
  file_name <- paste0("wide format of pixel differences by time and location, rep ", c, ".txt")
  fwrite(mtab, file = file.path(data_folder_path, file_name), sep = "\t")
}

print("Building the main data.table")
#### Determine number of bins ####
nframe_per_sec <- 30

exp.length = nrow(mtab)/nframe_per_sec
full_tab[, time_sec:=time/1E6] #changing the microseconds in the data to seconds


#### Important to keep all the experiment to the same length ####
#### Loop through table creating bins per location ####
loc_vec = unique(full_tab$location)

location_list = list()

for (i in seq_along(loc_vec)) {
  tmp_tab = full_tab[location == loc_vec[i]]
  tmp_tab[, bin:=cut(time_sec, nbins = nframe_per_sec, 
                     breaks = exp.length)]
  location_list[[i]] = tmp_tab
}

bin_tab = rbindlist(location_list)


dim(bin_tab)

dim(full_tab)

#### Summarize data per bin and per location ####

#Adding all of the movement together that is within a bin (~1/30 of a second), by location so that
#sum_table$time is per second
sum_table <- bin_tab[, list(time=round(time_sec[nframe_per_sec],  0)), 
                     by=list(location, bin)]

int_tab <- bin_tab[, lapply(.SD, sum), .SDcols = patterns('data'), by = list(location, bin)]

sum_table <- merge(sum_table, int_tab, by=c('location', 'bin'))

#### Replace na values ####
#number of na are summed, and any na in the time column are replaced by interpolations
print("printing the inital number of na after the pixel-differeces-per-frame are summed into second-long bins. Next: Interpolation.")
print(sum(is.na(sum_table)))

sum_table <- sum_table[, time := na.approx(time, na.rm=FALSE), by=list(location)]

#Number of na remaining after interpolation. Any remaining rows with na, which should just be the last rows, are replaced by extrapolation
print("printing the number of na after any were interpolated in the time column. Next: Extrapolation")
print(sum(is.na(sum_table)))

sum_table <- sum_table[, time := na.spline(time, na.rm=FALSE), by=list(location)]

#Number of na remaining after locf replacing na with previous data points
print("printing the number of na after any were extrapolated in the time column. Next: Removing rows with na in other columns")
print(sum(is.na(sum_table)))

sum_table <- na.omit(sum_table, na.rm=TRUE)

print("printing the number of na after rows containing na were removed, if any remained in the data.table")

print(sum(is.na(sum_table)))

#### Add Condition information specific to this experiment ####

sum_table[, loc_id:=str_extract(location, "\\d{2}")]

sum_table$loc_id <- as.double(sum_table$loc_id)

sum_table[, Condition:="empty"]

#Assigning Conditions and coordinates to location ID's
print("Assigning Conditions and coordinates to location ID's")
for(r in 1:8) {
  for(c in 1:12) {
    l = LETTERS[r]
    
    #assigning user-input experimental Conditions to the loc_ids
    sum_table[loc_id %in% (((r-1)*12)+c), Condition := paste(conds[r,c])]
    
    #assigning the plate coordinates (ie A1, B2, etc.) to the loc_ids
    sum_table[loc_id == (((r-1)*12)+c), loc_coord := paste0(l,c)]
  }
}

#### Check for wells with no movement for the entire duration of the experiment ####
sum_table <- sum_table[!Condition == "empty"]

sum_table$Condition <- factor(sum_table$Condition,
                              levels = gtools::mixedsort(unique(sum_table$Condition)))

sum_table_title <- "pixel differences per second, includes Condition and location info per replicate"
write.csv(sum_table, file=file.path(data_folder_path, paste0("sum_table - ", sum_table_title, ".csv")), row.names = FALSE)

#This subsets sum_table, but for only the selected time-span and inputs what the phase is for that time span
sum_tablenew <- data.table()

print("Organizing data by phases and groups")
for (t in 1:nrow(time_frame)) {
  x.sub <- sum_table[time >= time_frame$Start[t] & time <= time_frame$End[t]]
  x.sub[, Phase := time_frame$Phase[t]]
  x.sub[, Group := time_frame$Group[t]]
  sum_tablenew <- rbind(sum_tablenew, x.sub)
}
sum_tablenew$Phase <- factor(sum_tablenew$Phase, levels = unique(sum_tablenew$Phase))
sum_tablenew$Group <- factor(sum_tablenew$Group, levels = unique(sum_tablenew$Group))

#Pivot the sum_tablenew data.table so that it the datax columns in one column
#for easier manipulation
print("Pivoting the data table to long format")
sum_table_pivot <- sum_tablenew %>%
  pivot_longer(
    cols = starts_with("data"),
    names_to = "plate",
    names_prefix = "data",
    values_to = "exp_dat",
    values_drop_na = TRUE
  )
setDT(sum_table_pivot)

sum_table_pivot$plate <- as.integer(sum_table_pivot$plate)

sum_table_pivot_title <- "pixel differences per second, Condition and location info per replicate"
write.csv(sum_table_pivot, file=file.path(data_folder_path, paste0("sum_table_pivot - ", sum_table_pivot_title, ".csv")), row.names = FALSE)

#calculate the sum and mean for each Condition over time, per phase
print("Generating the 'table_conds_time' data.table")
table_conds_time <- sum_table_pivot[, list(n_samples=as.double(.N),
                                           exp_sums=as.double(sum(exp_dat)),
                                           mean=as.double(mean(exp_dat)),
                                           median=as.double(median(exp_dat)),
                                           iqr = as.double(IQR(exp_dat)),
                                           sd=as.double(sd(exp_dat)),
                                           sem=as.double(std.error(exp_dat))),
                                    by=list(Condition, time, Phase, Group)]
table_conds_time[, CI_lower := mean - (1.96 * sem)]
table_conds_time[, CI_upper := mean + (1.96 * sem)]

table_conds_time_title <- "pixel differences per second & condition, includes phase, group info"
write.csv(table_conds_time, file=file.path(data_folder_path, paste0("table_conds_time - ", table_conds_time_title, ".csv")), row.names = FALSE)


#calculate the sum and mean for each location, plate, phase, and Condition
table_location_conds_time <- sum_table_pivot[, list(seconds_total=as.double(.N),
                                                    exp_sums=as.double(sum(exp_dat)),
                                                    mean=as.double(mean(exp_dat)),
                                                    sd=as.double(sd(exp_dat)),
                                                    sem=as.double(std.error(exp_dat)),
                                                    median=as.double(median(exp_dat)),
                                                    iqr = as.double(IQR(exp_dat))),
                                             by=list(Condition, loc_coord, loc_id, plate, Phase, Group)]
table_location_conds_time[, CI_lower := mean - (1.96 * sem)]
table_location_conds_time[, CI_upper := mean + (1.96 * sem)]

phase_con_mean_sd <- table_location_conds_time %>%
  group_by(Condition, Phase, Group) %>%
  summarise(phase_con_mean = mean(mean),
            phase_con_sd = sd(mean),
            phase_con_sem = std.error(mean))

table_location_conds_time <- left_join(table_location_conds_time,
                                       phase_con_mean_sd,
                                       by = c("Condition", "Phase", "Group")) 

table_location_conds_time_title <- "movement data through time"
write.csv(table_location_conds_time, file=file.path(data_folder_path, paste0("table_location_conds_time - ", table_location_conds_time_title, ".csv")), row.names = FALSE)

#Testing the distribution of the data, as well as the variance and outliers within groups
#Levene's test for homogeneity of variance across groups
levene_test_Phase <- table_location_conds_time %>%
  group_by(Condition, Group) %>%
  levene_test(exp_sums ~ Phase)

write.csv(levene_test_Phase, file=file.path(data_folder_path, paste0("Levene's test for homogeneity of variance, within conditions between phases.csv")), row.names = FALSE)

levene_test_Conditions <- table_location_conds_time %>%
  group_by(Phase, Group) %>%
  levene_test(exp_sums ~ Condition)

write.csv(levene_test_Conditions, file=file.path(data_folder_path, paste0("Levene's test for homogeneity of variance, within phases between conditions.csv")), row.names = FALSE)

outlier_test <- table_location_conds_time %>%
  select(Condition, loc_coord, loc_id, Phase, Group, exp_sums) %>%
  group_by(Condition, Phase, Group) %>%
  identify_outliers(exp_sums)

write.csv(outlier_test, file=file.path(data_folder_path, paste0("Outlier test of whole experiment by Condition, Phase, and Group.csv")), row.names = FALSE)

# Shapiro-Wilk test
shapiro_test <- table_location_conds_time %>%
  group_by(Condition, Phase, Group) %>%
  shapiro_test(exp_sums)

write.csv(shapiro_test, file=file.path(data_folder_path, paste0("Shapiro-Wilks test of normality for whole experiment by Condition, Phase, and Group.csv")), row.names = FALSE)

#Create a facet grid of qq plots
facet_grid_qqplot <- ggqqplot(table_location_conds_time, "exp_sums", ggtheme = theme_bw()) +
  facet_grid(Condition + Group ~ Phase, labeller = "label_both")

ggsave(facet_grid_qqplot, filename = file.path(figures_folder_path, paste0("Whole experiment, facet grid qqplot by Condition, Phase, and Group.png")), width = 25, height = 15)

#calculate the sum and mean for each plate by phases
print("Generating the 'table_conds_plates' data.table")
table_conds_plates <- table_location_conds_time[, list(n_animals=as.double(.N),
                                                       exp_sums=as.double(sum(exp_sums)),
                                                       mean=as.double(mean(exp_sums)),
                                                       sd=as.double(sd(exp_sums)),
                                                       sem=as.double(std.error(exp_sums)),
                                                       median=as.double(median(exp_sums))),
                                                by=list(Condition, plate, Phase, Group)]
table_conds_plates[, CI_lower := mean - (1.96 * sem)]
table_conds_plates[, CI_upper := mean + (1.96 * sem)]


#calculate the sum and mean for each Condition, comparing plates
print("Generating the 'table_plates_compare' data.table")
table_plates_compare <- table_conds_plates[, list(n_plates=as.double(.N),
                                                  exp_sums=as.double(sum(exp_sums)),
                                                  mean=as.double(mean(exp_sums)),
                                                  sd=as.double(sd(exp_sums)),
                                                  sem=as.double(std.error(exp_sums)),
                                                  median=as.double(median(exp_sums))),
                                           by=list(Condition, Phase, Group)]
table_plates_compare[, CI_lower := mean - (1.96 * sem)]
table_plates_compare[, CI_upper := mean + (1.96 * sem)]

table_plates_compare_title <- "movement data between plates"
write.csv(table_plates_compare, file=file.path(data_folder_path, paste0("table_plates_compare - " , table_plates_compare_title, ".csv")), row.names = FALSE)

# Assuming 'data' is your dataset (ensure 'data' is a vector)
set.seed(123)  # For reproducibility
bootstrap_medians <- replicate(1000, median(sample(table_conds_time$median, replace = TRUE)))
ci <- quantile(bootstrap_medians, probs = c(0.025, 0.975))  # 95% confidence interval

table_conds_time_title <- "comparing movement data between conditions through time"
write.csv(table_conds_time, file=file.path(data_folder_path, paste0("table_conds_time - ", table_conds_time_title, ".csv")), row.names = FALSE)

#calculate the sum and mean for each Condition, per phase
print("Generating the 'table_conds' data.table")
table_conds <- table_location_conds_time[, list(n_animals=as.double(.N),
                                                exp_sums=as.double(sum(exp_sums)),
                                                mean=as.double(mean(exp_sums)),
                                                sd=as.double(sd(exp_sums)),
                                                sem=as.double(std.error(exp_sums)),
                                                median=as.double(median(exp_sums)),
                                                iqr = as.double(IQR(exp_sums))),
                                         by=list(Condition, Phase, Group)]
table_conds[, CI_lower := mean - (1.96 * sem)]
table_conds[, CI_upper := mean + (1.96 * sem)]

table_conds_title <- "movement data of conditions between phases"
write.csv(table_conds, file=file.path(data_folder_path, paste0("table_conds - ", table_conds_title, ".csv")), row.names = FALSE)

#Below is for comparing only phases of only the same group
#vector of the unique groups
group_vec <- unique(time_frame[['Group']])

#For generating figures of individual phases
for (g in 1:length(group_vec)) {
  current_group <- group_vec[[g]]
  
  #Create a folder for output, if they do not already exist
  # Define the path of the new folder
  group_figures_folder_path <- paste0("./Output/Group ", current_group,"/Comparison of phases/Figures")
  
  # Check if the folder already exists
  if (!dir.exists(group_figures_folder_path)) {
    # Create the folder if it does not exist
    dir.create(group_figures_folder_path, recursive = TRUE)
  }
  
  group_data_folder_path <- paste0("./Output/Group ", current_group,"/Comparison of phases/Data analysis")
  
  if (!dir.exists(group_data_folder_path)) {
    # Create the folder if it does not exist
    dir.create(group_data_folder_path, recursive = TRUE)
  }
  
  #subsetting the time_frame tibble so that the phrases of each group can be compared
  time_frame_sub <- time_frame[time_frame$Group == g, ]
  
  #vector of the unique phases for this experiment, regardless of their group
  phase_vec <- unique(time_frame_sub$Phase)
  
  transition_vec <- time_frame_sub %>%
    filter(Trans == TRUE) %>% pull(Phase) %>% unique()
  
  phase_vec <- phase_vec[!phase_vec %in% transition_vec]
  
  #sub-setting the data to only the groups currently iterated
  print("Subsetting each data.table to compare and graph phases within the same groups")
  
  sum_table_group <- droplevels(sum_table_pivot[sum_table_pivot$Group == group_vec[[g]] & (!sum_table_pivot$Phase %in% transition_vec), ])
  
  table_location_conds_time_group <- droplevels(table_location_conds_time[table_location_conds_time$Group == group_vec[[g]] & (!table_location_conds_time$Phase %in% transition_vec), ])
  
  table_conds_time_group <- droplevels(table_conds_time[table_conds_time$Group == group_vec[[g]] & (!table_conds_time$Phase %in% transition_vec), ])
  
  table_conds_group <- droplevels(table_conds[Group == group_vec[[g]] & (!table_conds$Phase %in% transition_vec), ])
  
  table_plates_compare_group <- droplevels(table_plates_compare[Group == group_vec[[g]] & (!table_plates_compare$Phase %in% transition_vec), ])
  
  #convert the Conditions to factors and order them for the figures
  table_location_conds_time_group$Condition <- factor(table_location_conds_time_group$Condition, 
                                                      levels = unique(table_location_conds_time_group$Condition[order(table_location_conds_time_group$Phase)]))
  
  ##Run a t.test of the Conditions across the different phases
  # Create a vector of unique Conditions
  Conditions <- unique(table_location_conds_time_group$Condition)
  
  #### Statistical analysis ####
  #Perform Shapiro-Wilk test for each Condition
  normalcy_test <- shapiro_test %>%
    filter(Group == current_group)
  
  #Subsetting Levene's test based on the current phase
  levene_test <- levene_test_Phase %>%
    filter(Group == current_group)
  
  if (any(normalcy_test$p < 0.05) || any(levene_test$p < 0.05)) {
    test <- 'Friedman Test'
    print("Performing Friedman Test")
    # Remove the Group column from the table
    
    # Friedman test
    friedman_result <- friedman.test(exp_sums ~ Condition | Phase, data = table_plates_compare_group)
    
    # Extracting necessary information
    m4 <- data.frame(
      statistic = friedman_result$statistic,
      p.value = friedman_result$p.value,
      df = friedman_result$parameter,
      method = friedman_result$method,
      data_name = friedman_result$data.name
    )
    write.csv(m4, file = file.path(data_folder_path, paste0("Phase group ", current_group, ", ", test, " of Conditions between phases.csv")), row.names = FALSE)
    
  } else {
    test <- "ANOVA model"
    
    # Fit the ANOVA model
    model <- aov(exp_sums ~ Phase + Error(Condition/Phase), data = table_location_conds_time_group)
    summary_model <- summary(model)
    
    summary_names <- names(summary_model)
    
    # Create a new workbook
    wb <- createWorkbook()
    
    # Extract the ANOVA table
    for (i in seq_along(summary_model)) {
      current_name <- gsub(":", "-", summary_names[[i]])
      names(summary_model[[i]]) <- NULL
      
      addWorksheet(wb, current_name)
      
      combined_anova_df <- data.frame()
      
      # Loop through each data frame within the list
      for(j in 1:length(summary_model[[i]][[1]]$Df)) {
        # Bind the row to the combined data frame
        anova_df <- data.frame(
          'Deg Freedom' = summary_model[[i]][[1]]$Df[[j]],
          'Sum Sq' = summary_model[[i]][[1]]$'Sum Sq'[[j]],
          'Mean Sq' = summary_model[[i]][[1]]$'Mean Sq'[[j]],
          'F value' = summary_model[[i]][[1]]$'F value'[[j]],
          'Pr(>F)' = summary_model[[i]][[1]]$'Pr(>F)'[[j]]
        )
        combined_anova_df <- rbind(combined_anova_df, anova_df)
      }
      
      # Write data to the sheets
      writeData(wb, sheet = current_name, x = combined_anova_df)
    }
    
    # Save the workbook to a file
    saveWorkbook(wb, file = file.path(group_data_folder_path, paste0("Phase group ", current_group, ", ", test, " of Conditions between phases.xlsx")), overwrite = TRUE)
    
  }
  
  # pairwise test loop for statistical difference between Conditions every second
  print("Performing pairwise tests")
  test <- "Non-parametric test"
  
  # Perform the Pairwise tests
  # Create an empty data.table to store the results
  stat_test <- data.table()
  
  # Perform tests for each combination of Condition and phase
  for (c in Conditions) {
    # Get unique phases for the current Condition
    Condition_phases <- unique(table_location_conds_time_group[Condition == c, Phase])
    
    # Check if there are at least 2 unique phases
    if (length(Condition_phases) >= 2) {
      # Generate pairwise combinations of phases
      phase_combs <- combn(Condition_phases, 2, simplify = FALSE)
      
      # Apply the appropriate test to each phase combination
      for (i in 1:length(phase_combs)) {
        phase1 <- phase_combs[[i]][1]
        phase2 <- phase_combs[[i]][2]
        
        # Subset the data for the current phase combination
        subset1 <- table_location_conds_time_group[Condition == c & Phase == phase1, exp_sums]
        subset2 <- table_location_conds_time_group[Condition == c & Phase == phase2, exp_sums]
        
        if (any(normalcy_test$p < 0.05) || any(levene_test$p < 0.05)) {
          # Perform the Wilcoxon Signed-Rank Test
          test_result <- wilcox.test(subset1, subset2, paired = TRUE, exact = FALSE, correct = FALSE, alternative = "two.sided", conf.int = TRUE, conf.level = 0.95)
          
          # Store the results in the data.table
          stat_test <- rbind(stat_test, data.table(Condition = c,
                                                   phase1 = phase1, 
                                                   phase2 = phase2,
                                                   pseudo_median_difference = test_result$estimate,
                                                   conf_int_lower = test_result$conf.int[1],
                                                   conf_int_upper = test_result$conf.int[2],
                                                   conf_level = test_result$conf_level,
                                                   statistic = test_result$statistic,
                                                   tails = test_result$alternative,
                                                   p.value = test_result$p.value,
                                                   method = test_result$method,
                                                   data_name = test_result$data.name))
        } else {
          test <- "Parametric test"
          # Perform the t-test
          test_result <- t.test(subset1, subset2, paired = TRUE, correct = FALSE, alternative = "two.sided", conf.level = 0.95, conf.int = TRUE)
          
          # Store the results in the data.table
          stat_test <- rbind(stat_test, data.table(Condition = c,
                                                   phase1 = phase1, 
                                                   phase2 = phase2,
                                                   mean_difference = test_result$estimate,
                                                   conf_int_lower = test_result$conf.int[1],
                                                   conf_int_upper = test_result$conf.int[2],
                                                   conf_level = test_result$conf_level,
                                                   stderr = test_result$stderr,
                                                   statistic = test_result$statistic,
                                                   df = test_result$parameter,
                                                   p.value = test_result$p.value,
                                                   method = test_result$method,
                                                   data_name = test_result$data.name))
        }
        
        print(paste('pairwise ', test, ' of Conditions between phases', phase1, 'and', phase2,', for Condition', c, 'of group', current_group))
      }
    }
  }
  
  # Get the unique names of the rows in a column
  Condition_vec_len <- unique(sum_table_group$Condition)
  all_combinations <- combn(Condition_vec_len, 2, simplify = FALSE)
  
  #filter out combinations to reduce unnecessary pairwise comparisons
  filtered_combinations <- Filter(function(pair) {
    any(pair %in% controls)
  }, all_combinations)
  
  stat_test$p.adjusted <- p.adjust(stat_test$p.value, method = group_adjustment_method)
  stat_test$adjustment_method <- group_adjustment_method
  
  write.csv(stat_test, file = file.path(group_data_folder_path, paste0("pairwise tests comparing Conditions in different phases, group ", current_group, ".csv")), row.names = FALSE)
  
  # Define the conditional formatting function
  format_with_condition <- function(x, threshold = 0.001) {
    sapply(x, function(value) {
      if (abs(value) < threshold) {
        format(value, scientific = TRUE, digits = 2)
      } else {
        format(value, scientific = FALSE, digits = 3)
      }
    })
  }
  
  formatted_p_values <- format_with_condition(stat_test$p.adjusted)
  
  #This is determines where the geom_signif annotations sit in relation to the data
  max_row_index <- which.max(table_location_conds_time_group$mean)
  max_sd_value <- table_location_conds_time_group[max_row_index, "sd"]
  top_mean <- max(table_location_conds_time_group$mean) + max_sd_value
  top_median <- max(table_location_conds_time_group$median + table_location_conds_time_group$iqr)
  select_mean_plus <- as.integer(top_mean - top_mean*0.05)
  select_median_plus <- as.integer(top_median - top_median*0.05)
  
  #Roughly the width of a bin in the figures
  full_width <- 0.90
  bar_width <- 0.80 / length(phase_vec)
  
  annotation_df <- data.frame(Condition = character(), id = numeric(), y_position_mean = numeric(), y_position_median = numeric())
  
  max_results <- list()
  min_results <- list()
  
  comparisons_num <- choose(length(phase_vec), 2) #gives the number of unique combinations using n! / (r! * (n-r)!)
  
  cond_vec <- as.character(unique(table_location_conds_time_group$Condition))
  
  conds_length <- length(cond_vec)
  
  for (c in 1:conds_length) {
    cond <- cond_vec[[c]]
    id <- c
    num_increases <- length(phase_vec)-1    # Initial repetition count
    step_value <- bar_width
    start_min <- id - (0.5*bar_width) * (length(phase_vec) - 1)
    
    start_max <- id + (0.5*bar_width) * (length(phase_vec) - 1)  
    num_decreases <- length(phase_vec) - 1 
    
    y_position_mean <- select_mean_plus
    y_position_median <- select_median_plus
    
    for (i in 1:comparisons_num) {
      new_row <- data.frame(Condition = cond,
                            id = id,
                            y_position_mean = y_position_mean,
                            y_position_median = y_position_median,
                            stringsAsFactors = FALSE)
      
      annotation_df <- rbind(annotation_df, new_row)
      
      y_position_mean <- y_position_mean + (select_mean_plus*0.1)
      y_position_median <- y_position_median + (select_median_plus*0.1)
    }
    
    min_list <- list()
    max_list <- list()
    
    current_val <- start_min
    
    for (i in 1:num_increases) {
      # Create a sublist with the repeating value  
      sublist <- rep(current_val, num_increases)  
      min_list <- c(min_list, sublist)
      
      # Update variables for the next iteration
      current_val <- current_val + step_value
      num_increases <- num_increases - 1
    }
    
    min_results[[c]] <- unlist(min_list) # Store results for this id
    
    max_list <- numeric(0) # Initialize an empty list to store all numbers
    
    while (num_decreases > 0) {
      current_list <- seq(start_max, step_value, by = -step_value) 
      current_list <- current_list[1:num_decreases] 
      
      max_list <- c(max_list, current_list) # Append generated numbers
      
      num_decreases <- num_decreases - 1
      
    }
    max_results[[c]] <- unlist(max_list) # Store results for this id 
  }
  
  annotation_df$xmin <- unlist(min_results)
  annotation_df$xmax <- unlist(max_results)
  
  diamond_annotation <- table_conds_group %>% select(mean)
  
  #For generating figures of like phases
  for (p in 1:length(phase_vec)) {
    
    current_phase <- phase_vec[[p]]
    
    #Create a folder for output, if they do not already exist
    # Define the path of the new folder
    phase_figures_folder_path <- paste0("./Output/Group ", current_group,"/Phase ", current_phase, "/Figures")
    
    # Check if the folder already exists
    if (!dir.exists(phase_figures_folder_path)) {
      # Create the folder if it does not exist
      dir.create(phase_figures_folder_path, recursive = TRUE)
    }
    
    phase_data_folder_path <- paste0("./Output/Group ", current_group,"/Phase ", current_phase, "/Data analysis")
    
    if (!dir.exists(phase_data_folder_path)) {
      # Create the folder if it does not exist
      dir.create(phase_data_folder_path, recursive = TRUE)
    }
    
    #sub-setting the data to only the phases currently iterated
    print("Subsetting each data.table to analyze and graph each phase")
    
    sum_table_phase <- droplevels(sum_table_pivot[sum_table_pivot$Phase == phase_vec[[p]]])
    
    table_location_conds_time_phases <- droplevels(table_location_conds_time[table_location_conds_time$Phase == phase_vec[[p]]])
    
    table_conds_time_phases <- droplevels(table_conds_time[table_conds_time$Phase == phase_vec[[p]]])
    
    table_conds_phases <- droplevels(table_conds[table_conds$Phase == phase_vec[[p]]])
    
    table_plates_compare_phases <- droplevels(table_plates_compare[table_plates_compare$Phase == phase_vec[[p]]])
    
    #### Statistical analysis ####
    #Perform Shapiro-Wilk test for each Condition
    normalcy_test <- shapiro_test %>%
      filter(Phase == current_phase)
    
    #Subsetting Levene's test based on the current phase
    levene_test <- levene_test_Conditions %>%
      filter(Phase == current_phase)
    
    if (any(normalcy_test$p < 0.05) || any(levene_test$p < 0.05)) {
      test = 'Kruskal-Wallis Test'
      print("Performing Kruskal-Wallis Test")
      # Kruskal-Wallis test
      kruskal_results <- kruskal.test(exp_sums ~ Condition, data = table_location_conds_time_phases)
      
      # Extracting necessary information
      m3 <- data.frame(
        statistic = kruskal_results$statistic,
        p.value = kruskal_results$p.value,
        method = kruskal_results$method,
        parameter = kruskal_results$parameter,
        data_name = kruskal_results$data.name
      )
      write.csv(m3, file = file.path(phase_data_folder_path, paste0("Group , ", current_group, ", ", current_phase, " phase, ", test, " by Conditions.csv")), row.names = FALSE)
      
    } else {
      test = "ANOVA model"
      
      # Fit the ANOVA model
      model <- aov(exp_sums ~ Condition, data = table_location_conds_time_phases)
      
      # Extract the ANOVA table
      anova_df <- data.frame(summary(model)[[1]])
      rownames(anova_df) <- NULL
      
      # Create a new workbook
      wb <- createWorkbook()
      
      # Add sheet to the workbook
      addWorksheet(wb, "ANOVA Table")
      
      # Write data to the sheets
      writeData(wb, sheet = "ANOVA Table", x = anova_df)
      
      # Save the workbook
      saveWorkbook(wb, file = file.path(phase_data_folder_path, paste0("Group , ", current_group, ", ", current_phase, " phase, ", test, " by Conditions.xlsx")), overwrite = TRUE)
    }
    
    # pairwise test loop for statistical difference between Conditions every second
    print("Performing pairwise tests")
    
    if (any(normalcy_test$p < 0.05) || any(levene_test$p < 0.05)) {
      # Perform the Wilcoxon test
      pair_test <- pairwise.wilcox.test(table_location_conds_time_phases$exp_sums, table_location_conds_time_phases$Condition)
    } else {
      # Perform the t-test
      pair_test <- pairwise.t.test(table_location_conds_time_phases$exp_sums, table_location_conds_time_phases$Condition)
    }
    
    # Convert the matrix to a data frame/tibble
    p_matrix <- as.data.frame(pair_test$p.value)
    
    # Add the row names as a column
    p_matrix$Condition1 <- rownames(p_matrix)
    
    #Turn the test results to long format and save the parameters of the tests performed
    test_results <- melt(p_matrix, id.vars = "Condition1", variable.name = "Condition2", value.name = "p.value", na.rm = TRUE)
    test_results$p.adjusted <- p.adjust(test_results$p.value, method = phase_adjustment_method)
    test_results$adjustment_method <- phase_adjustment_method
    test_results$method <- pair_test$method
    
    write.csv(test_results, file = file.path(phase_data_folder_path, paste0("Group , ", current_group, ", ", current_phase, " phase, ", " pairwise comparison by Condition.csv")), row.names = FALSE)
    
    #This determines which Condition differences will be compared for statistical significance 
    test_results <- test_results %>%
      mutate(pair = as.numeric(factor(1:nrow(test_results))),
             xmin = pair - 0.2,
             xmax = pair + 0.2) %>%
      filter(Condition1 %in% controls | Condition2 %in% controls)
    
    #column names of the Conditions
    comparisons <- grep("^Condition", names(test_results), value = TRUE)
    
    #data.table of the Conditions column, with rows of comparisons with controls
    test_results <- as.data.table(test_results)
    comparisonsdt <- test_results[, ..comparisons, with = FALSE]
    
    comparison_list <- lapply(seq_len(nrow(comparisonsdt)), function(i) {
      unlist(list(as.character(comparisonsdt[i, Condition1]), as.character(comparisonsdt[i, Condition2])))
    })
    
    comparison_list <- comparison_list[order(sapply(comparison_list, `[`, 1))]
    
    
    #get the mean of the sums for the normalizing control as a numeric
    norm_mean <- as.numeric(table_conds_phases$mean[table_conds_phases$Condition == normalizing_control])
    
    #get the sd for the normalizing control as a numeric
    norm_sum_sd <- as.numeric(table_conds_phases$sd[table_conds_phases$Condition == normalizing_control])
    
    #calculate the 'hit_cut off', which is 2*sd above the mean(sum(WT control))
    hit_cutoff_num <- as.numeric(norm_mean + (2*norm_sum_sd)) #converting mean + sd of the control to a numeric
    
    norm_cond_loc_tab <- table_location_conds_time_phases[, list(seconds=as.double(.N),
                                                                 normalized_sums = exp_sums/hit_cutoff_num),
                                                          by=list(Condition, loc_coord, loc_id, plate)]
    
    write.csv(norm_cond_loc_tab, file = file.path(phase_data_folder_path, paste0("Group , ", current_group, ", ", current_phase, " phase, ", "normalized sums of movement by location.csv")), row.names = FALSE)
    
    norm_cond_tab <- norm_cond_loc_tab[, list(n_samples=as.double(.N),
                                              exp_sums=as.double(sum(normalized_sums)),
                                              mean=as.double(mean(normalized_sums)),
                                              median=as.double(median(normalized_sums)),
                                              sd=as.double(sd(normalized_sums)),
                                              sem=as.double(std.error(normalized_sums))),
                                       by=list(Condition)]
    norm_cond_tab[, CI_lower := mean - (1.96 * sem)]
    norm_cond_tab[, CI_upper := mean + (1.96 * sem)]
    
    write.csv(norm_cond_tab, file = file.path(phase_data_folder_path, paste0("Group , ", current_group, ", ", current_phase, " phase, ", "normalized sums of movement by Condition.csv")), row.names = FALSE)
    
    #remove the 'control groups' to prepare for graphing
    norm_cond_tab <- norm_cond_tab[!norm_cond_tab$Condition %in% controls, ]
    
    only_hits_above_cutoff <- function(x) { ifelse(norm_cond_tab[["mean"]] - norm_cond_tab[["sd"]] > 1, as.character(x), "") }
    
    norm_cond_tab$hit_conds <- only_hits_above_cutoff(norm_cond_tab[["Condition"]])
    
    labels <- setNames(norm_cond_tab$hit_conds, norm_cond_tab$Condition)
    
    print(paste0("Creating figure Q, phase ", current_phase, " of group ", current_group))
    #plot sum(movement) of the experimental Conditions / (mean + 2 standard deviations) of the control
    fig_Q <- ggplot(norm_cond_tab, aes(x=Condition, y=mean)) + 
      geom_point(size = 3) +
      labs(color = 'Conditions', x = 'Condition', y = 'Sum of Movement Normalized to Normalizing Control Mean + 2 sd') +
      geom_hline(yintercept=1, linetype="dashed", color = "black", linewidth=0.5) +
      geom_errorbar(aes(ymin = mean - sd, 
                        ymax = mean + sd,
                        color = ifelse((mean - sd) > 1, 'Hit', 'Fail'))) +
      scale_color_manual(values = c('#000000', '#FF0000')) + #red is a hit, black is a fail
      scale_x_discrete(labels = labels) +
      theme_classic(base_size = 20) +
      theme(
        legend.position = 'right', #also controls 'no legend'
        legend.title = element_text(size=30),
        legend.text = element_text(size=25),
        plot.title = element_text(face = 'bold', margin = margin(10, 0, 10, 0), hjust = 0.5, size = 12), #margin(t, r, b, l)
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing.x = unit(0.5, "cm"),
        axis.text.y = element_text(size = 25),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 25),
        axis.title = element_text(size=35, color = 'black', face = 'bold'),
        axis.title.x = element_text(margin = margin(t = 12), size = 24), #margin(t, r, b, l)
        axis.title.y = element_text(margin = margin(r = 10), size = 24))
    
    ggsave(fig_Q, filename = file.path(phase_figures_folder_path, paste0("Group , ", current_group, ", ", "Figure Q, ", current_phase, " phase, ", " hit cut-off experimental Conditions.png")), height=15, width=25, units="in", dpi = 150)
    
  }
}

end_time <- Sys.time()

print(end_time - start_time)
