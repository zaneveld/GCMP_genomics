#Goal: test if rarefaction depth alters relative richness of coral microbiomes

#Rarefy skeleton samples to super high depth, losing most samples
qiime feature-table rarefy --i-table ../input/feature_table_skeleton.qza --p-sampling-depth 100000 --o-rarefied-table ../output/feature_table_skeleton_100k.qza

#Take *just* the high-sequencing-depth samples, and re-rarefy to 1000. This will separate the effects of rarefaction depth itself vs. just losing samples, since these samples all have high rarefaction depth.
 qiime feature-table rarefy --i-table ../output/feature_table_skeleton_100k.qza --p-sampling-depth 1000 --o-rarefied-table ../output/feature_table_skeleton_1k.qza

#Calculate alpha diversity on the 100k rarefied samples
qiime diversity alpha --i-table ../output/feature_table_skeleton_100k.qza --p-metric observed_features --o-alpha-diversity ../output/adiv_100k.qza

#...and on the 1k rarefied samples
 qiime diversity alpha --i-table ../output/feature_table_skeleton_1k.qza --p-metric observed_features --o-alpha-diversity ../output/adiv_1k.qza

#Now extract the raw data
 qiime tools extract --input-path ../output/adiv_1k.qza --output-path ../output/adiv_1k
 qiime tools extract --input-path ../output/adiv_100k.qza --output-path ../output/adiv_100k
