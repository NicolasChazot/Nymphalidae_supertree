library(ape)
library(picante)
library(geiger)
library(phytools)
library(TreePar)

info_table<-read.table("Nymphalidae_subclades_info_test.txt", header=TRUE)

N=1000 #number of backbone trees grafted

#load the backbone tree posterior distribution
BACKBONE_POSTDIST<-read.nexus("backbone_6g_ref_nucCOI12_COMBINED.trees")

#remove burnin
BACKBONE_POSTDIST_burnin<-BACKBONE_POSTDIST[101:1001]

#sample N backbone tree
BACKBONE_POSTDIST_sampled<-sample(x= BACKBONE_POSTDIST,size=N,replace=F)

#load the posterior distribution of subclade trees
ITHOMIINI_POSTDIST_sampled <-read.tree("Ithomiini_1000trees.trees")

APATURINAE_POSTDIST_sampled <-read.tree("Apaturinae_1000trees.trees")

BIBLIDINAE_POSTDIST_sampled <-read.tree("Biblidinae_1000trees.trees")

CALINAGINAE_POSTDIST_sampled <-read.tree("Calinaginae_1000trees.trees")

CHARAXINAE_POSTDIST_sampled <-read.tree("Charaxinae_1000trees.trees")

CYRESTINAE_POSTDIST_sampled <-read.tree("Cyrestinae_1000trees.trees")

DANAINI_POSTDIST_sampled <-read.tree("Danaini_1000trees.trees")

HELICONIINAE_POSTDIST_sampled <-read.tree("Heliconiinae_1000trees.trees")

LIBYTHAEINAE_POSTDIST_sampled <-read.tree("Libythaeinae_1000trees.trees")

LIMENITIDINAE_POSTDIST_sampled <-read.tree("Limenitidinae_1000trees.trees")

MORPHOGROUP_POSTDIST_sampled <-read.tree("Morpho_group_1000trees.trees")

MYCALCOENO_POSTDIST_sampled <-read.tree("Mycalesina_Coenonympha_1000trees.trees")

NYMPHALINAE_POSTDIST_sampled <-read.tree("Nymphalinae_1000trees.trees")

PSEUDERGOLINAE_POSTDIST_POSTDIST_sampled <-read.tree("Pseudergolinae_1000trees.trees")

SATYRINI_POSTDIST_sampled <-read.tree("Satyrini_1000trees.trees")

subclade_dir<-list(ITHOMIINI_POSTDIST_sampled,APATURINAE_POSTDIST_sampled,BIBLIDINAE_POSTDIST_sampled,CALINAGINAE_POSTDIST_sampled,CHARAXINAE_POSTDIST_sampled
            ,CYRESTINAE_POSTDIST_sampled,DANAINI_POSTDIST_sampled,HELICONIINAE_POSTDIST_sampled
            ,LIBYTHAEINAE_POSTDIST_sampled,LIMENITIDINAE_POSTDIST_sampled,MORPHOGROUP_POSTDIST_sampled
            ,MYCALCOENO_POSTDIST_sampled,NYMPHALINAE_POSTDIST_sampled,PSEUDERGOLINAE_POSTDIST_POSTDIST_sampled,SATYRINI_POSTDIST_sampled)

#output object containing all the grafted trees
GRAFTED_POST_DIST<-list()

for (j in 1:N) #the j backbone tree sampled from the posterior distribuion
{

backbonej_started<-BACKBONE_POSTDIST_sampled[[j]]

	for (i in 1:15) #subclade i (Danaini, Ithomiini, etc)
	{
	  if(i==1) #this initiate a new backbone when the first subclade has to be grafted
	  {
	    backbonej<-backbonej_started
	  }
	  
	  else
	  {
	 	backbonej<-backbonej_grafted
	  }

		if(i==4) #special case for Calinaginae grafted up the stem
		{
			subcladei<-subclade_dir[i][[1]][[j]]
			subcladei_rescaled<-rescaleTree(subcladei,1)

			subclade_node_ages<-getx(subcladei_rescaled)
			backbone_node_ages<-as.matrix(getx(backbonej))

			backbone_stem<-backbone_node_ages[rownames(backbone_node_ages)==getMRCA(backbonej,c(as.character(info_table[i,2]),as.character(info_table[i,4]))),1][[1]]
			subclade_crown<-subclade_node_ages[["11"]]
			subclade_stem<-subclade_node_ages[["10"]]
			
			backbone_crown<-backbone_stem*subclade_crown/subclade_stem

			subcladei_rescaled<-drop.tip(subcladei_rescaled,c(as.character(info_table[i,2]),as.character(info_table[i,3])))
			subcladei_forgraft<-rescaleTree(subcladei_rescaled,backbone_crown)
			subcladei_forgraft$root.edge<-0
			
			backbonej_wtout_sub<-backbonej
			backbonej_wtout_sub$tip.label[which(backbonej_wtout_sub$tip.label==info_table[i,4])]<-"NA"

			backbonej_wtout_sub$edge.length[which.edge(backbonej_wtout_sub,"NA")]<-backbone_stem-backbone_crown

			backbonej_grafted<-paste.tree(backbonej_wtout_sub, subcladei_forgraft)
		}
		
		else #all other subclades, grafted at the tip of the backbone stem
		
		{
      		#extract all the tip names corresponding to the subclade to be grafted (outgroups excluded) IN THE BACKBONE
			subclade_tips<-tips(backbonej,getMRCA(backbonej,c(as.character(info_table[i,4]),as.character(info_table[i,5]))))
     		
     		#drop all the tips corresponding to the subclade to be grafted but one, FROM THE BACKBONE
			backbonej_wtout_sub<-drop.tip(backbonej,setdiff(subclade_tips,info_table[i,4]))


			#pick one subclade tree
			subcladei<-subclade_dir[i][[1]][[j]]

			#rescale the subclade tree to 1
			subcladei_rescaled<-rescaleTree(subcladei,1)

			#node ages
			subclade_node_ages<-getx(subcladei_rescaled)
			backbone_node_ages<-as.matrix(getx(backbonej))

 		   	#find the stem age and the crown age of the subclade to be grafted IN THE BACKBONE
			backbone_stem<-backbone_node_ages[rownames(backbone_node_ages)==getMRCA(backbonej,c(as.character(info_table[i,2]),as.character(info_table[i,4]))),1][[1]]
			backbone_crown<-backbone_node_ages[rownames(backbone_node_ages)==getMRCA(backbonej,c(as.character(info_table[i,4]),as.character(info_table[i,5]))),1][[1]]

      		#rescale the subclade tree to the backbone estimate after deleting the outgroups
			subcladei_rescaled<-drop.tip(subcladei_rescaled,c(as.character(info_table[i,2]),as.character(info_table[i,3])))
			subcladei_forgraft<-rescaleTree(subcladei_rescaled,backbone_crown)
			
			#define the stem length IN THE SUBCLADE to 0
			subcladei_forgraft$root.edge<-0

			#change the name of the subclade branch to "NA" IN THE BACKBONE
			backbonej_wtout_sub$tip.label[which(backbonej_wtout_sub$tip.label==info_table[i,4])]<-"NA"

			#change the branch length of the subclade branch IN THE BACKBONE to the stem length
			backbonej_wtout_sub$edge.length[which.edge(backbonej_wtout_sub,"NA")]<-backbone_stem-backbone_crown

			#graft the rescaled subclade tree to the subclade branch in the backbone
			backbonej_grafted<-paste.tree(backbonej_wtout_sub, subcladei_forgraft)
		}
	}
  print(paste(j*100/N,"%",sep=" "))
  GRAFTED_POST_DIST[[j]]<-backbonej_grafted
}


class(GRAFTED_POST_DIST)<-"multiPhylo"
write.tree(GRAFTED_POST_DIST,"GRAFTED_POST_DIST.nex")
