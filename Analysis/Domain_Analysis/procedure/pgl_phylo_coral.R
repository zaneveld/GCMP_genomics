#Phylogenetic generalized least squares regression and phylogenetic generalized ANOVA
#load libraries
library(ape)
library(nlme)
library(geiger)
library(caper)
library(MuMIn)
library(rr2)
library(phytools)
library(MASS)
library(robustbase)

#create this file only the first time then we append to it
pgl.output = data.frame("Immune_Trait", "Diversity_Trait", "N_Microbiomes", "N_Genomes", "Model_1", "pVal_Model_1", "Rsquared_Model_1", "AIC_Model_1", "AICc_Model_1", "Model_2", "pVal_Model_2", "Rsquared_Model_2", "AIC_Model_2", "AICc_Model_2", "Model_3", "pVal_Model_3", "Rsqared_Model_3", "AIC_Model_3", "AICc_Model_3", "Model_4", "pVal_Model_4", "Rsquared_Model_4", "AIC_Model_4", "AICc_Model_4", "Minimum_AIC", "Minimum_AICc", "PIC_Model", "pVal_PIC", "Rsquared_PIC")
pgl.output
write.table(pgl.output,file="coral_output/PGL_coral_output_genomes.csv",append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)



trait_table_input <-"input/coral_genome_trait_table.csv"

#Pick the columns to analyze
#y traits can be: obs_asvs, ln_asvs, dominance, gini_index, simpson_e, faith_pd
y_trait_column <-"gini_index_tissue"
#x traits can be: TIR_total, TIR_total_unique, IL1R, LRR_total, LRR_total_unique, TLR, Lectin_total, Lectin_unique
x_trait_column<-"TLR"
trait_table <- read.csv(trait_table_input,header=TRUE,row.names=1)
x_trait <-trait_table[[x_trait_column]]
y_trait <-trait_table[[y_trait_column]]

#load the newick file created with the huang roy tree names
tree <- "input/huang_roy_molecular_r2.newick"
tree <- read.tree(tree)
print(tree)

#make sure that data and tree use same species names
coral_obj <- name.check(tree,trait_table)
coral_obj
tree.prune <- coral_obj

#take out tree tips that do not have data
tree.prune<-drop.tip(tree, coral_obj$tree_not_data)

prune_tree <- name.check(tree.prune,trait_table)
prune_tree
plot(tree.prune)

#run phylogenetic GLS using brownian motion process
spp = trait_table$Species
bm<-corBrownian(1, tree.prune)
#bm<-corBrownian(1, tree.prune, form=~spp)
bm

#first model is a simple model
model.1_formula <-"gls(x_trait~y_trait, data=trait_table, correlation=bm)"
model.1<-gls(x_trait~y_trait, data=trait_table, correlation=bm)
#model.1.rlm<-rlm(model.1)
summary(model.1)
#summary(model.1.rlm)

#extract the AIC and AICc for model 1
model.1_aic <- AIC(model.1)
model.1_aic
model.1_aicc <- AICc(model.1)
model.1_aicc
#extract only the coefficient output so that we can more easily extract the p value
summary_model.1<-coef(summary(model.1))
#get only the p value from the summary table for model 1
model.1_pval<- (summary_model.1)[2,4]
#can't get the r2 directly from the summary so need to go using rr2
model.1_rSquared <- R2.lik(model.1)
model.1_rSquared
print(paste ("The pvalue for model 1 is ", model.1_pval, "the R2 value is", model.1_rSquared, "and the AIC/AICc is", model.1_aic, model.1_aicc))

#model 2 we relax brownian motion and use the model of Pagel
model.2_formula<-"gls(x_trait~y_trait, data=trait_table, correlation=corPagel(1,tree.prune))"
model.2<-gls(x_trait~y_trait, data=trait_table, correlation=corPagel(1,tree.prune))
summary(model.2)
#extract the AIC and AICc for model 2
model.2_aic <- AIC(model.2)
model.2_aic
model.2_aicc <- AICc(model.2)
model.2_aicc
#extract only the coefficient output so that we can more easily extract the p value
summary_model.2<-coef(summary(model.2))
#get only the p value from the summary table for model 2
model.2_pval<- (summary_model.2)[2,4]
#need to use a round about way to get the r2 since its not in the summary.
model.2_rSquared <- R2.lik(model.2)
model.2_rSquared
print(paste ("The pvalue for model 2 is ", model.2_pval, "the R2 value is", model.2_rSquared, "and the AIC/AICc is", model.2_aic, model.2_aicc))

#model 3b takes into account multiple parameters. I'm just adding TLR to the list to see.
#Can uncomment this if we want to run multiple parameters
#model.3b_value<-"gls(x_trait~y_trait+3rd_trait, data=trait_table, correlation=corPagel(1,tree.prune))"
#model.3b <-gls(x_trait~y_trait+TIR_total, data=trait_table, correlation=corPagel(1,tree.prune))
#summary(model.3)

#analyze data using the package caper which combines the data and tree into one file
#trait_table <-"input/symbiodinium_trait_table_its2.csv"
#reuse the trait_table_input from the top so that we don't need to reload the entire table
trait_table <- read.csv(trait_table_input,header=TRUE)

comp.data <- comparative.data(tree.prune, trait_table, names.col="Species", vcv.dim=2, warn.dropped = TRUE)

#run data using the basic model
model.3_formula<-"pgls(x_trait~y_trait, data=comp.data)"
model.3<-pgls(x_trait~y_trait, data=comp.data)
summary(model.3)
#extract the AIC and AICc from Model 3
model.3_aic <- AIC(model.3)
model.3_aic
model.3_aicc <- AICc(model.3)
model.3_aicc
#extract only the coefficient output so that we can more easily extract the p value
summary_model.3<-coef(summary(model.3))
#get only the p value from the summary table for model 3
model.3_pval<- (summary_model.3)[2,4]
#extract the R squared value for model 3
model.3_rSquared <- summary(model.3)$r.squared
print(paste ("The pvalue for model 3 is ", model.3_pval, "R squared is ", model.3_rSquared, "and the AIC/AICc is ", model.3_aic, model.3_aicc))

#run data using the maximum likelihood analysis
model.4_formula<-"pgls(x_trait~y_trait, data=comp.data, lambda='ML')"
model.4<-pgls(x_trait~y_trait, data=comp.data, lambda="ML")
summary(model.4)
#extract the AIC and AICc from Model 4
model.4_aic <- AIC(model.4)
model.4_aic
model.4_aicc <- AICc(model.4)
model.4_aicc
#extract only the coefficient output so that we can more easily extract the p value
summary_model.4<-coef(summary(model.4))
#get only the p value from the summary table for model 4
model.4_pval<- (summary_model.4)[2,4]
#extract the R squared value for model 4
model.4_rSquared <- summary(model.4)$r.squared
print(paste ("The pvalue for model 4 is ", model.4_pval, "the R squared is ", model.4_rSquared, "and the AIC/AICc is ", model.4_aic, model.4_aicc))


#run model with multiple parameters.
#uncomment this model if interested in testing multiple variables.
#model.4b_formula<-"pgls(x_trait~y_trait+ 3rd_trait, data=comp.dat,lambda='ML'"
#model.4b<-pgls(x_trait~y_trait+TIR_total, data=comp.data,lambda="ML")
#summary(model.4b)

#write a data frame with each AIC for each combination and identify the AIC and AICc with 
#the lowest value.
aic.df <- data.frame(Model_1_AIC=model.1_aic, Model_2_AIC=model.2_aic, Model_3_AIC=model.3_aic, Model_4_AIC=model.4_aic)
aic.df
aic.df_columns <- colnames(aic.df)
aic.df_columns
aic.df$min <- apply(aic.df,1,which.min)
aic.df$min
aic.df$min_model <- colnames(aic.df)[aic.df$min]
aic.df$min_model

aicc.df <- data.frame(Model_1_AICc=model.1_aicc, Model_2_AICc=model.2_aicc, Model_3_AICc=model.3_aicc, Model_4_AICc=model.4_aicc)
aicc.df
aicc.df_columns <- colnames(aicc.df)
aicc.df_columns
aicc.df$min <- apply(aicc.df,1,which.min)
aicc.df$min
aicc.df$min_model <- colnames(aicc.df)[aicc.df$min]
aicc.df$min_model

#need to redefine the traits for some reason.
trait_table_input <-"input/coral_genome_trait_table.csv"
#Pick the columns to analyze
#y traits can be: obs_asvs, ln_asvs, dominance, gini_index, simpson_e, faith_pd
y_trait_column <-"gini_index_tissue"
#x traits can be: TIR_total, TIR_total_unique, IL1R, LRR_total, LRR_total_unique, Lectin_total, Lectin_unique
x_trait_column<-"TLR"
trait_table <- read.csv(trait_table_input,header=TRUE,row.names=1)
x_trait <-trait_table[[x_trait_column]]
y_trait <-trait_table[[y_trait_column]]

#make plots of the traits and save as a pdf 
phylomorphospace(tree.prune,trait_table[,c(x_trait_column,y_trait_column)],xlab=x_trait_column,ylab=y_trait_column)
pdf(paste("coral_output/",x_trait_column,"_",y_trait_column,"_phylomorphospace.pdf",sep=""))
phylomorphospace(tree.prune,trait_table[,c(x_trait_column,y_trait_column)],xlab=x_trait_column,ylab=y_trait_column)
points(x_trait,y_trait,pch=21,bg="grey",cex=1.4)
dev.off()

# Extract columns
#host <- anoleData[, "hostility"]
#awe <- anoleData[, "awesomeness"]

# Give them names
names(x_trait) <- rownames(trait_table)
names(y_trait) <- rownames(trait_table)

#for contrasts, you should positivize them, since the order doesn't matter. This is NOT taking absolute value.
PositivizeContrasts <- function(x_trait, y_trait) {
  #Cache the sign of x so it doesn't
  #change as we reflect points about x-axis!
  sign_of_x <- sign(x_trait)
  x_trait.positivized <- x_trait * sign_of_x
  y_trait.positivized <- y_trait * sign_of_x
  return(cbind(x_trait.positivized, y_trait.positivized))
}




# Calculate PICs
x_trait_pic <- pic(x_trait, tree.prune)
#print(x_trait_pic)
y_trait_pic <- pic(y_trait, tree.prune)
#print(y_trait_pic)

positivized.results <- PositivizeContrasts(x_trait_pic, y_trait_pic)
x_trait_positive <-positivized.results[,1]
y_trait_positive <-positivized.results[,2]

# Make a pic model
#pic_model <- lm(x_trait ~ y_trait - 1)
# plot pic results
#plot(x_trait ~ y_trait,xlab=x_trait_column,ylab=y_trait_column,bg='gray',pch=16)
#abline(a = 0, b = coef(pic_model))

# Make a pic model that has only positive values
pic_model <- lm(y_trait_positive ~ x_trait_positive - 1)
pic_rlm <- rlm(y_trait_positive ~ x_trait_positive -1)
# plot pic results
plot(y_trait_positive ~ x_trait_positive,xlab=x_trait_column,ylab=y_trait_column,bg='gray',pch=16)
abline(a = 0, b = coef(pic_model))
#Save raw PIC contrasts as a pdf
pdf(paste("coral_output/",x_trait_column,"_",y_trait_column,"_pic_scatter_YX.pdf",sep=""))
plot(y_trait_positive ~ x_trait_positive,xlab=x_trait_column,ylab=y_trait_column,bg='gray',pch=16)
abline(a = 0, b = coef(pic_model))
dev.off()

#summarize the results
summary(pic_model)
summary(pic_rlm)
rSquared_PIC <-summary(pic_model)$r.squared
pval_PIC <- anova(pic_model)$'Pr(>F)'[1]
pval_PIC.rlm <- anova(pic_rlm)$'Pr(>F)'[1]
model_PIC <- "lm(y_trait_positive ~ x_trait_positive -1)"

#rlm does not produce a p value. Will try doing the robust regression using "robustbase"
robust.model = lmrob(y_trait_positive ~ x_trait_positive - 1)
summary(robust.model)


# plot contMap
x_trait_reconstruction<-contMap(tree.prune,x_trait)
plot(x_trait_reconstruction,direction="rightwards")
#save plot as pdf
pdf(paste("coral_output/",x_trait_column,"_contmap.pdf",sep=""))
plot(x_trait_reconstruction,direction="rightwards")
dev.off()

y_trait_reconstruction<-contMap(tree.prune,y_trait)
plot(y_trait_reconstruction,direction="leftwards")
#save plot as pdf
pdf(paste("coral_output/",y_trait_column,"_contmap.pdf",sep=""))
plot(y_trait_reconstruction,direction="leftwards")
dev.off()

#write results from the PGL models to the created csv file

trait_table.df <- data.frame(trait_table)
#this extracts the code from the tree
n_microbe <- sum(trait_table.df$Num_samples_microbe)
n_genome <- sum(trait_table.df$Num_samples_genome)
output_row <- data.frame(y_trait_column, x_trait_column, n_microbe, n_genome, model.1_formula, model.1_pval, model.1_rSquared, model.1_aic, model.1_aicc, model.2_formula, model.2_pval, model.2_rSquared, model.2_aic, model.2_aicc, model.3_formula, model.3_pval, model.3_rSquared, model.3_aic, model.3_aicc, model.4_formula, model.4_pval, model.4_rSquared, model.4_aic, model.4_aicc, aic.df$min_model, aicc.df$min_model, model_PIC, pval_PIC, rSquared_PIC)
print(output_row)

#Write it to the csv file
write.table(output_row, file="coral_output/PGL_coral_output_genomes.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)


#plot at likelyhood surface of the lambda parameter
#lm.lk<-pgls.profile(model.5, which="lambda")
#plot(lm.lk)

#Here we want to run the correction on only the best fit models.
#Is there a way to group by best model?
#correct p-values from the csv file for multiple comparisons
#extract columns we want
pgl_input_table <- "./symbiodinium_pics/PGL_Sym_output_its2&clade.csv"
pgl_table <- read.csv(pgl_input_table,header=TRUE)
pgl_table.df <- data.frame(pgl_table)
p_value_columns <- pgl_table.df[,c("Immune_Trait", "Diversity_Trait","pVal_Model_1", "pVal_Model_2", "pVal_Model_3", "pVal_Model_4")]
#p_value_columns
#Select Immune Trait Values
#Immune traits can be: TIR_total, TIR_total_unique, LRR_total, LRR_total_unique, Lectin_total, Lectin_unique
Immune_trait_variable <- p_value_columns[p_value_columns$Immune_Trait == "Lectin_total",]
#Immune_trait <- p_value_columns$Immune_trait_variable
Immune_trait_variable
#Select by Diversity Trait Value
#Diversity traits can be: obs_asvs, ln_asvs, dominance, gini_index, simpson_e, faith_pd
Diversity_trait_variable <- p_value_columns[p_value_columns$Diversity_Trait =="simpson_e",]
Diversity_trait_variable
#Select the p value column you will be using
#p value columns can be "p_val_Model_1","p_val_Model_2","p_val_Model_3","p_val_model_4"
p_value_select <- Diversity_trait_variable$pVal_Model_1
p_value_select
#run p adjust for bonferroni
adjusted_bonferroni <- p.adjust(p_value_select,method="bonferroni")
adjusted_bonferroni
#Might not need this
#run p adjust for hochberg
adjusted_hochberg <- p.adjust(p_value_select,method="hochberg")
adjusted_hochberg
#run p adjust for BH (fdr)
adjusted_bh <- p.adjust(p_value_select,method="BH")
adjusted_bh

