#format data from plate reader to graph growth rates
library(ggplot2)
library(dplyr)

#Written by Renee Geck, July 2021
#Updated October 2021
#Linear fit written by Anja Ollodart and adapted by Renee Geck

##README
#Required input: .txt tab-separated spreadsheet where:
#The first 3 rows are headers for (1) strain, (2) condition, (3) replicate number (integer)
#The first column is time in hours (or minutes, could be any numerical unit)
#Columns 2-n are OD readings from the plate reader, one column per well
#Change file path in line 20

#______________________________
##PART 1: Format the input data

#Import .txt file, no header
plate_input <- read.delim("file.txt", header=FALSE)

#Separate headers from data
plate_headers <- plate_input %>% slice(1:3) #the first 3 rows are header
plate_data <- plate_input %>% slice(-(1:3)) #everything not the first 3 rows is data

#take column 1 from data, that is the time
plate_time <- plate_data %>% pull(1) #take col 1 as a vector
plate_time <- as.numeric(plate_time) #treat as numbers not strings

#find the length of the data (to know how long to make label vectors)
plate_length <- length(plate_time)

#make a dataframe to store the values
all_data=data.frame()

#Make a column counter to iterate over
col_iter <- c(2:ncol(plate_data))

#Go through each column and do the following:
for (iter in col_iter) {
  col_num <- iter
  #make labels
  well_strain <- plate_headers[1,col_num] #strain is first column
  strain_ls <- rep(well_strain, plate_length) #repeat strain name for the length of the data
  well_cond <- plate_headers[2,col_num] #condition is second column
  cond_ls <- rep(well_cond, plate_length) #repeat condition name for the length of the data
  well_rep <- as.integer(plate_headers[3,col_num]) #replicate is second column
  rep_ls <- rep(well_rep, plate_length) #repeat replicate number for the length of the data
  #get od values and make numbers, not strings
  od_ls <- plate_data[,col_num]
  od_ls <- as.numeric(od_ls)
  
  #assemble a df from the parts
  #time, od, strain, condition, rep
  well_data <- data.frame(plate_time, od_ls, strain_ls, cond_ls, rep_ls)
  #append the whole thing onto the bottom of the ongoing df
  all_data <- rbind(all_data, well_data)
  #increase the iter to repeat
  iter <- iter+1
}

#rename the columns to things that make more sense (and don't have other meanings in R)
all_data <- all_data %>% rename(time=plate_time, dens=od_ls, strain=strain_ls, cond=cond_ls, replicate=rep_ls)

#________________________________
##PART 2: Graph the growth curves

#plot smooth line for the set, using method="gam"
ggplot(data=all_data,
       aes(x=time, y=dens, color=strain)) + #plot od/time by strain
  geom_smooth(method='gam', position = "identity") + #curve fit method, default 95%ci
  scale_y_continuous(trans = "log2",n.breaks = 15,
                     labels = function(x) ifelse(x == 0, "0", x)) + #make log2 yaxis
  facet_wrap(vars(cond)) + #make a different plot for each condition
  xlab("Time (hours)") +
  ylab("OD600") +
  theme_classic() #no gridlines

#_______________________________________________
##PART 3: Calculate and graph the doubling times
#Note: Depending on the starting density, may need to change lower bound in line 94
#Note: You can change the Rsq cutoff in line 111 if your curves are not well-fit

#make vector of possible upper bounds of linear range
up_bound_ls <- seq(from=1, to=0.15, by=-0.01)
#make dataframe to hold all the possible linear fits
all_linfits = data.frame()

#do linear fitting, testing each possible upper bound
#records the upper bounds and Rsq of fit to select later
for (up_bound in up_bound_ls) { #for each of the possible upper bounds of linear range
  this_linfit <- all_data %>% 
    group_by(cond, strain, replicate) %>%  #do calculations for each individual sample
    filter(dens>0.13 & dens<up_bound) %>% #take only OD readings from linear range
    do(mod_lin = summary(lm(log2(dens) ~ time, data = .))) %>% # this calculates slope and outputs a list of all the info associated with the Lm function 
    mutate(intercept = mod_lin$coefficients[1], #make a column for each of the described values
           slope = mod_lin$coefficients[2], #determine the slope of the fit
           Rsq = mod_lin$r.squared, #report the goodness of fit
           doubling_time = log(2)/mod_lin$coefficients[2], #calculate the doubling time from the slope
           upper_bound = up_bound) %>% #also report the upper bound used
    select(-mod_lin)
  #append the whole thing onto the bottom of the ongoing df
  all_linfits <- rbind(all_linfits, this_linfit)
  #increase the count to repeat
  up_bound <- up_bound+1
}

#if Rsq>0.98, keep that linear fit
cutoff_linfit <- all_linfits %>%
  as.data.frame() %>%
  filter(Rsq > 0.98)

#for each group, take the fit with Rsq>0.98 that used the highest upper bound
#this is what will be used to determine the doubling time of the strain
best_linfit <- cutoff_linfit %>%
  as.data.frame() %>%
  group_by(cond, strain, replicate) %>% #do for each sample
  filter(upper_bound == max(upper_bound)) #take the one that fit the most points (highest upper bound)

#_______________________________________________
##PART 4: Graph the doubling times
#make a box plot of the doubling times by strain for each condition
ggplot(data=best_linfit, mapping=aes(x=cond, y=doubling_time, color=strain)) +
  geom_boxplot() +
  labs(x="", y="Doubling time (hours)") +
  theme_classic() #no gridlines