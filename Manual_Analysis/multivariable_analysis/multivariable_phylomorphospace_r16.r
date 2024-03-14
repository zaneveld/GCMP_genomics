
#A script for running phytools PIC and phylomorphospace analysis from
#the command line

#See below for usage. 
print("Usage: Rscript phylomorphospace_r14.r <path_to_trait_table> <path_to_tree> <x_trait_column> <y_trait_column> <z_trait_column> <filter_column> <filter_value> <output_dir_suffix>")

library(ggplot2)
library(phytools)

#Get user input and assign to variables
print("Parsing command line arguments ....")
args <- commandArgs(trailingOnly=TRUE)
print(args)
trait_table_fp <-args[1]
tree_fp <- args[2]
x_traits <-args[3]
y_trait <-args[4]
filter_column <- args[6]
filter_value <- args[7]
user_specified_output_dir <- args[8]
suffix <- args[9]

#Extract a vector of trait column names from
#comma separated list
x_trait_columns <- unlist(strsplit(x_traits, ","))
print(paste0("Found x trait columns:",x_trait_columns))

print(paste0("Creating set_output_dir_name function..."))
#create function to automatically name the output directory given 
#command line arguments
set_output_dir_name <- function(output_dir_for_all_pics,x_traits,y_trait,filter_column,
  filter_value,suffix){
    output_dir <- paste0(output_dir_for_all_pics,"PIC_",x_traits,"_vs_",y_trait)

    #Add filter column and value to directory name, if they were provided
    if (!is.na(filter_column) & (filter_column != 'None') & (filter_value != 'None')){
        output_dir <- paste0(output_dir,"_",filter_column,"_is_",gsub("\\;","_",filter_value))
    }
    #Add a special suffix to the output dir, if the user requested one
    if (!is.na(suffix)){
        output_dir <- paste0(output_dir,"_",suffix)
    }
    output_dir <- paste0(output_dir,"/")

    print(paste("Outputting results to:",output_dir))   
    return(output_dir)
}

#Set output directory for all PICs 
#(results of this analysis will be created in a subdirectory)
output_dir_for_all_pics = "./PIC_results/"
if (!is.na(user_specified_output_dir) & (user_specified_output_dir != 'None')){
    output_dir_for_all_pics = paste0(user_specified_output_dir,"/")
}
dir.create(output_dir_for_all_pics, showWarnings = FALSE)
print(paste0("Created output dir for all PIC results:",output_dir_for_all_pics))

#Create output subdirectory for this analysis specifically
output_dir <- set_output_dir_name(output_dir_for_all_pics,x_traits,y_trait,filter_column,
  filter_value,suffix)
dir.create(output_dir,showWarning = FALSE)
print(paste0("Created output dir for this analysis:",output_dir))


#create the formula string we'll use for the PIC regression
#Note that x1 * x2 * x3 is shorthand that tests interaction terms *and* main effects of x1,x2,x3
#formula_str <- paste(y_trait, "~",  paste(x_trait_columns, collapse = " * "),"-1")
formula_str <- paste(y_trait, "~",  paste(x_trait_columns, collapse = " + "),"-1")
print(paste0("Setting formula string to:",formula_str))

#Output results to a log file as well as to the screen
sink(paste0(output_dir,x_traits,"_vs_",y_trait,"_results_log.txt"),append=FALSE,split=TRUE)

#Summarize command line parameters for the log
print("Phylogenetic Independent Contrast Analysis Report")
print("--------------------------------------------------")
print(paste("Analyzing",x_traits,"vs.",y_trait))
print(paste("Trait table filepath:",trait_table_fp))
print(paste("Tree filepath:",tree_fp))
print(paste("Filtering data based on column:",filter_column))
print(paste("Including data only if filter column value is:",filter_value))
print(paste("Output dir:",output_dir))
print(paste("Regression formula:",formula_str))
print("--------------------------------------------------")

#Load trait table
trait_table <- read.table(trait_table_fp,comment.char="",header=T,row.names=1,as.is=T,sep="\t")

#Load tree
tree <- read.tree(tree_fp)

#Generate a vector with all relevant column names
trait_columns = c(y_trait,x_trait_columns)

# Filter the trait table to include only the specified columns
filtered_trait_table <- trait_table[, trait_columns]

# Define a function for removing rows with missing data
remove_taxa_with_missing_data <- function(column, tree, trait_table) {
  cat("Analyzing Trait Column:", column, "\n")

  #Coerce trait column to numeric, possibly inducing NAs
  trait_table[,column]<-as.numeric(trait_table[,column])
  #Drop rows with an NA for that trait
  trait_table<- trait_table[!is.na(trait_table[,column]),]
  return(trait_table)
}

#Optionally filter the trait table based on a specific column value
#prior to harmonizing the tree and trait table
if (!is.na(filter_column) & filter_column != 'None' & filter_value != 'None'){
    trait_table <- subset(trait_table,trait_table[[filter_column]] == filter_value)
}

# Filter specified trait columns for missing data
for (column in trait_columns) {
  filtered_trait_table <- remove_taxa_with_missing_data(column, tree, filtered_trait_table)
  print(filtered_trait_table)
}

#Define a function to drop tree tips or trait table rows that don't match
harmonize_tree_and_trait_table <-function(trait_table,tree){
	#Filter table to match tree tips
	filtered_trait_table <- trait_table[rownames(trait_table) %in% tree$tip.label,]
	#Filter tree to match table
	filtered_tree <- drop.tip(tree,tree$tip.label[!tree$tip.label %in% rownames(trait_table)])
	#Dichotomize tree
	tree <- multi2di(filtered_tree)
	return(list(tree = filtered_tree, trait_table = filtered_trait_table))
}

#Drop mismatched data or tips between trait table and tree
harmonized_data <- harmonize_tree_and_trait_table(filtered_trait_table,tree)
trait_table <- harmonized_data$trait_table
tree <- harmonized_data$tree
print(paste0("tree class:",class(tree)))

#Record the filtered tree and trait table
write.tree(tree,paste0(output_dir,x_traits,"_vs_",y_trait,"_filtered_tree.newick"))
write.table(trait_table,paste0(output_dir,x_traits,"_vs_",y_trait,"_filtered_table.csv"))

#****POSITIVIZE PIC VALUES (negative signs are arbitrary)
#Note: this is not just the absolute value since the difference x=-1 y=2
#has a different interpretation than x=1 y=2

#Code snippet from:https://github.com/bomeara/ComparativeMethodsInR/blob/master/ContinuousTrait_Answers.R
#for contrasts, you should positivize them, since the order doesn't matter. This is NOT taking absolute value.

positivize_contrasts <- function(reference_PIC,all_PICs) {
    #Cache the sign of x so it doesn't
    #change as we reflect points about x-axis!
    print("figuring out sign of x")
    sign_of_x <- sign(all_PICs[[reference_PIC]])
    print(sign_of_x)
    for (pic in names(all_PICs)){
    	all_PICs[[pic]] <- all_PICs[[pic]] * sign_of_x
    }
    
    return(all_PICs)
}


#Run PICs and return a list of PICs by their column name
run_PICs <- function(trait_table,tree,trait_columns,reference_trait_column){
	result = list()
	for (trait in trait_columns) {
  		curr_data = trait_table[,trait]
  		names(curr_data) = rownames(trait_table)
  		print(paste0("tree class:",class(tree)))
  		raw_pic <- pic(curr_data,tree)
  		result[[trait]] <- raw_pic
	}
	
	#Positivize results to the first 
	result <- positivize_contrasts(reference_trait_column,result)
	return(result)	
}

reference_trait_column <- x_trait_columns[1]
pic_result <- run_PICs(trait_table,tree,trait_columns,reference_trait_column)

#Convert the pic results to a dataframe
# *Apparently* we do this with do.call.
# This feels cludgy, and can be updated if 
#there's a better way to do it, but should get the job done.

# Convert the list to a dataframe
pic_df <- do.call(data.frame, pic_result)
pic_df <- as.data.frame(pic_df)

names(pic_df) <- trait_columns
print(paste("PIC dataframe names:",names(pic_df)))
print(paste("PIC dataframe:",pic_df))
print(class(pic_df))



# regress through origin (e.g. expect that 0 change in one trait is on average
# correlated with zero change in the other)

print(paste0("Converting formula ",formula_str, "to formula"))
formula = as.formula(formula_str)
print("Fitting linear model")
fit <- lm(formula,data=pic_df)
print(paste0("Summary lm ",formula_str," for" ,x_traits,"(x) and ",y_trait,"(y)"))
print(summary(fit))



## this is a projection of the tree into morphospace
##This code snippit is adapted from a phytools tutorial (http://www.phytools.org/Cordoba2017/ex/3/PICs.html)
for (x_trait in x_trait_columns){
	print(paste("Plotting ",x_trait, "(x-axis) vs. ",y_trait,"(y-axis)"))
	X <- trait_table[[x_trait]]
	Y <- trait_table[[y_trait]]
	print(X)
	print(Y)
	pdf(paste0(output_dir,x_trait,"_vs_",y_trait,"_phylomorphospace.pdf"))
	phylomorphospace(tree,cbind(X,Y),xlab=x_trait,ylab=y_trait,label="off",node.size=c(0,0))
	points(X,Y,pch=21,bg="firebrick",cex=1.4)
	dev.off()

	#Save raw PIC contrasts as a pdf
	pdf(paste0(output_dir,x_trait,"_vs_",y_trait,"_pic_scatter_YX.pdf"))

	print("PIC df")
	print(names(pic_df))

	#Do some quick calculations to figure out the y-axis limits
	#We want them symmetrical above and below the null 0 slope expectation
    pic.Y <- pic_df[[y_trait]]
	pic.X <- pic_df[[x_trait]]
	largest_limit = max(abs(min(pic.Y)),abs(max(pic.Y)))

	#largest limit must be positive, so set the lower bound
	#to be that far below 0
	lower_bound = largest_limit * -1
	upper_bound = largest_limit 

	#Set the null model label 75% of the way along the axis
	null_model_label_x_pos = max(pic.X)*0.75 
	null_model_label_y_pos = 0 + 0.05*upper_bound

	fontsize = 16

	ggplot(pic_df, aes(pic.X,pic.Y)) + 
    geom_smooth(method = "lm",se=TRUE,level=0.95,col = "black",fill="cyan",formula = y ~ x -1) +
    geom_point(size = 3, shape=21, fill = "lightgrey",colour = "black",alpha=0.5) + 
    labs(x = paste("Contrast in ",x_trait), y = paste("Contrast in ",y_trait))+
    scale_y_continuous(limits=c(lower_bound,upper_bound))+ 
    geom_hline(yintercept=0, color = "firebrick", size=1, linetype="dotted")+
    geom_label(
        label="Null expectation", 
        x=null_model_label_x_pos,
        y=null_model_label_y_pos,
        label.size = NA,
        size = 6, #not a font size, but a point size
        color = "firebrick",
        face = "bold",
        fill = NA
    )+
    
    theme_classic()+
    theme(axis.text=element_text(size=16))+
    theme(axis.title=element_text(size=16))+
    theme(text = element_text(family = "sans"))

	dev.off()
}
#Build contmap for trait X
	
fit<-fastAnc(tree,X,vars=TRUE,CI=TRUE)
#Print model fit to screen
print(paste(c("FastAnc ML modelfit for",x_trait)))
print(fit)
obj <- contMap(tree,X,plot=F)
tree_direction <- "rightwards"
inverse_green_colorscheme <-c('black','springgreen3','yellow','white')
obj <- setMap(obj,colors=inverse_green_colorscheme)

#Write contmap for trait X to file
pdf(paste0(output_dir,x_trait,"_asr_contmap.pdf"))
par(mai=c(12.12,1,1.1,1.1))
plot(obj,direction=tree_direction,legend=0.7*max(nodeHeights(tree)),fsize=c(0.222,0.9))
axis(1)
title(xlab="time from the root (mya)")
dev.off()

#Build contmap for trait y
fit<-fastAnc(tree,Y,vars=TRUE,CI=TRUE)
#Print model fit to screen
print(paste(c("FastAnc ML modelfit for",y_trait)))
print(fit)
obj <- contMap(tree,Y,plot=F)
tree_direction <- "leftwards"
inverse_green_colorscheme <-c('black','springgreen3','yellow','white')
obj <- setMap(obj,colors=inverse_green_colorscheme)

#Write contmap for trait X to file
pdf(paste0(output_dir,y_trait,"_asr_contmap_leftwards.pdf"))
par(mai=c(12.12,1,1.1,1.1))
plot(obj,direction=tree_direction,legend=0.7*max(nodeHeights(tree)),fsize=c(0.222,0.9))
axis(1)
title(xlab="time from the root (mya)")
dev.off()


