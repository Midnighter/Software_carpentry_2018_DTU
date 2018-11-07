library(ggplot2)
library(tidyverse)
data <- read.csv('2018_0817_Wt_131.csv')
data_reduced <- select(data, Well.Name, All.Events..Events, All.Events.FITC.A.Mean ) #3 columns out of 108 that I am interested in
names(data_reduced) <- c('Well', 'count', 'Mean_fluorescence') # Rename the columns

head(data_reduced)

data_reduced$Mean_fluorescence[data_reduced$count < 5000] <- NA    # Remove wells that have less than 5000 cells

Both_plates <- data_reduced
Both_plates[1:96, 'No_ligand'] <- data_reduced$Mean_fluorescence[1:96]           #Plate_A_No_Ligand
Both_plates[1:96, 'VAC'] <- data_reduced$Mean_fluorescence[193:288]      #Plate_A_VAC
Both_plates[1:96, 'Vanillin'] <- data_reduced$Mean_fluorescence[385:480]     #Plate_A_Vanillin

Both_plates[97:192, 'No_ligand'] <- data_reduced$Mean_fluorescence[97:192]      #Plate_B_No_Ligand
Both_plates[97:192, 'VAC'] <- data_reduced$Mean_fluorescence[289:384]     #Plate_B_VAC
Both_plates[97:192, 'Vanillin'] <- data_reduced$Mean_fluorescence[481:576]     #Plate_B_Vanillin
Both_plates <- Both_plates[1:192,]
Both_plates <- select(Both_plates, -count, -Mean_fluorescence)

Both_plates[1:192, 'Type'] <- 'Mutant'       # Add an extra variable
Both_plates[95:96, 'Type'] <- 'Control'
Both_plates[191:192, 'Type'] <- 'Control'
#Both_plates[c(95:96, 191:192), 'Type'] <- 'Control'

Both_plates$Well <- as.character(Both_plates$Well)
Both_plates$Well[1:96] <- paste0('Plate_1_', Both_plates$Well[1:96])
Both_plates$Well[97:192] <-  paste0('Plate_2_', Both_plates$Well[97:192])


myplot_1 <- gather(Both_plates,  No_ligand, VAC, Vanillin, key = 'ligand', value = 'induction' ) #myplot_1 is all raw data
myplot_1 <- ggplot(myplot_1, aes(x= Well, color = interaction(ligand, Type)))
myplot_1 + geom_point(aes(y= induction), size = 3 ) + facet_grid(. ~ ligand) 



Fold_change <- Both_plates

Fold_change[1:192, 'No_ligand'] <- Both_plates$No_ligand / Both_plates$No_ligand
Fold_change[1:192, 'VAC'] <- Both_plates$VAC / Both_plates$No_ligand
Fold_change[1:192, 'Vanillin'] <- Both_plates$Vanillin / Both_plates$No_ligand

myplot_2 <- gather(Fold_change,  No_ligand, VAC, Vanillin, key = 'ligand', value = 'induction' ) #myplot_2 is all fold-induction data
myplot_2 <- ggplot(myplot_2, aes(x= Well, color = ligand))
myplot_2 + geom_point(aes(y= induction), size = 3) + facet_grid(. ~ ligand)


df_sorted <- arrange(Fold_change, -desc(VAC))
df_sorted$Well <- factor(df_sorted$Well, levels = unique(df_sorted$Well))

myplot_3 <- gather(df_sorted,  No_ligand, VAC, Vanillin, key = 'ligand', value = 'induction' )
myplot_3$Well <- factor(myplot_3$Well, levels = unique(myplot_3$Well))

myplot_3 <- ggplot(myplot_3, aes(x= Well, y= induction, color = interaction(ligand, Type)))
myplot_3 + geom_point() + facet_grid(. ~ ligand)


myplot_3_b <- ggplot(df_sorted,aes(x = Well, y = VAC, color = Type))
myplot_3_b + geom_point( size = 4)

top_variants <- filter(df_sorted, VAC < 1.45)
top_variants[29:30, ] <- Fold_change[95:96,]  
top_variants[31:32, ] <- Fold_change[191:192,] 

top_variants <- arrange(top_variants, -desc(Vanillin))

myplot_4 <- gather(top_variants, VAC, Vanillin, key = 'ligand', value = 'induction' )
myplot_4$Well <- factor(myplot_4$Well, levels = unique(myplot_4$Well))
 
myplot_4 <- ggplot(myplot_4, aes(x= Well, y= induction, color = interaction (ligand, Type )))
myplot_4 + geom_point(size = 4 ) + facet_grid(. ~ ligand) + theme(axis.text.x = element_text(angle = 90, hjust = 1))

top_variants
