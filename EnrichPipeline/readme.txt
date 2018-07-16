Content
1. introduction
2. Examples
3. output file format introduction
4. miscellaneous useful programs
5. metadata files introduction

1. introduction

This pipeline is designed to do enrichment analysis for the category data, 
such as GO/KEGG/IPR, etc. These data are characterised as one class containing
many genes and one gene is involved in many categories. So, for every class,
a p value is calculated representing the probability that the observed numbers
of counts could have resulted from randomly distributing this class between
the tested gene list and the reference gene list. Usually, the reference gene
list is the total genes of one organism annotated, but one can use customized
gene list as reference. The p value can be approximated by ChiSquare-Fisher
test or hypergeometric test. The p values are corrected by fdr or other
methods.
For GO enrichment, first the hierarchy feature of GO information is converted
to MetaGO_release.RData object from GO.db package of bioconductor. This
package will be updated every half year. Second, the GOdata_organism.RData are
constructed for specific organisms. This file contain the map relationship
for GOs at level 3 to genes and genes to the lowest level GO annotated,
several functions used in the enrichment analysis. Finally the enrichment
analysis is done based on the GOdata_organism.RData and supplied gene list at
level 3 and all levels. To remove the redundance, if the GOs enriched at
different levels with parent-child relationship and have the same gene list,
the lowest level is choose and other levels are filtered. The results file is
suffix "merg".


Reference:
(if you used this program, any of these two papers should be cited)

1. Chen S, Yang P, Jiang F, Wei Y, Ma Z, Kang L: De Novo Analysis of Transcriptome
Dynamics in the Migratory Locust during the Development of Phase Traits. PLoS One 2010,
5(12):e15633.
2. Huang da, W., B.T. Sherman, and R.A. Lempicki. 2009. Bioinformatics
enrichment tools: paths toward the comprehensive functional analysis of large
gene lists. Nucleic Acids Res 37: 1-13.
3. Beissbarth, T. and T.P. Speed. 2004. GOstat: find statistically
overrepresented Gene Ontologies within a group of genes. Bioinformatics 20:
1464-1465.

2. Example

PLEASE ADD THE FOLLOWING LINE TO YOUR ~/.bash_profile, OR ADD TO YOUR SHELL SCRIPT THAT
RUN THE PROGRAM
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/panfs/ANIMAL/GROUP/group001/yangpch/soft/lib/

1. GO enrichment analysis
1.1 construct GOdata
perl bin/EnrichPipeline.pl --cls GO --MetaGO data/MetaGO_200908.RData --wegoF input/gene.wego  --outDir output/GOdata  --GOdata_outF GOdata_example.RData

1.2 do enrichment
perl bin/EnrichPipeline.pl --cls GO --GOdata output/GOdata_example.RData --supplyF input/geneids.list1 --outDir output/EnrichGO --P_Adjust_Method fdr --TestMethod  FisherChiSquare

2. KEGG enrichment analysis
perl bin/EnrichPipeline.pl --cls KEGG --mapGene  input/KEGG.map.gene.txt  --supplyF input/geneids.list1 --outDir output/EnrichKEGG

3. IPR enrichment analysis
perl bin/EnrichPipeline.pl --cls IPR  --ipr2geneF input/ipr2gene.txt.ancest  --supplyF input/geneids.list1 --outDir output/EnrichIPR 
perl bin/EnrichPipeline_just_for_test.pl --cls IPR --ipr2geneF input/ipr2gene.txt.ancest --supplyF input/geneids.list1 --supplyF2 input/geneids.list2 --outDir output/EnrichIPR/

3. output files 

*dif*: 
the main output of the enrichment analysis

ColumnNumber	Name	Explanation
1	Class ID	Class id, such as GO id, IPR id
2	class Title	Class title, such as GO term, KEGG pathway name
3	Pvalue		Enrichment p value
4	AdjustedPv	Adjusted p value
5	x		Number of genes present in this class of all supplied genes
6	y		Number of genes present in this class of all reference genes
7	n		Number of supplied genes
8	N		Number of reference genes
9	EnrichDirect	Enrich direction, Over or Under.
10	GeneIDs		Gene ids that enriched in this class of all supplied genes

*.filt: 
the filtered results with adjusted p value less than 0.05 and enrichment
direction "Over"

*.R
The R script used in the process

*.Rout
The output of the R script, which can be used to monitor the process.

*.merg
specific for GO enrichment.

4. miscellaneous useful programs

4.1 Extract_gene_by_go.pl
extract all genes that belong to defined gos

4.2 prepare_data_*.pl
prepare the file used for ipr/kegg enrichment analysis

4.3 Plot.GOslim.R
draw the barchart graph of supplied GOs for GO enrichment analysis result files

4.4 add_annotation_to_enrich_result.pl
add the annotation to the enrich result

4.5 cat_enrich_result.pl
cat the enrich result into one excel file

4.6 cat_enrich_anno_result.pl
cat the enrich result with added annotation information together into one excel file

4.7 paste_id_files.pl
paste ids together into one excel file

5. used data sets

ParentChildTreeFile.txt
-----------------------
ftp://ftp.ebi.ac.uk/pub/databases/interpro/ParentChildTreeFile.txt

#################
# Bug Report
#################
Pengcheng Yang
pengchy@gmail.com

perl bin/Extract_gene_by_go.pl output/GOdata/GOdata_example.RData input/go_id.list1 output/Extract_Gene.list

