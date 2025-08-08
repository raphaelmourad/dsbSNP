################################################################################
##11/07/2025
#Sébastien AUBER


#this script is used to easaly compute average profile with deeptools then plot them wit ggplot.
#it is used in figure 3 and 4
#bed file for input dsbSNPs must be ordered as "*input_SNP_name*/bed/dsb.bed"
#and corresponding controls as "*input_SNP_name*/bed/ctrl.bed"
#interest protein or DSB mapping bigwig as *protein*/bw/*.bw
###############################################################################

library(readr)
library(plyranges)
library(BSgenome.Hsapiens.UCSC.hg19)
library(reticulate)
library(tidyr)
setwd(paste0("/media/sauber/Elements/DSBsnp_project_seb/dsbSNP_Av_pr/"))

use_condaenv("base")  
# ctrl<-read_bed(  "Legube_dsbSNP/bed/Legube_SNP_ctrl_dist500kb_GR.bed")
# same_size_ctrl<-sort(sample(ctrl,length(ctrl)/10))
# write_bed(same_size_ctrl,"Legube_dsbSNP/bed/Legube_SNP_ctrl_dist500kb_GR_same_size.bed")

#ENd_seq
pair0<-c("Legube_dsbSNP","new_End_seq_noOHT")
pair1<-c("Legube_dsbSNP","End_seq_merge_juillet23")


pairlist<-list(pair0,pair1)
KB_bp = 1600  #half KB
KB="3.2kb"
bin=10
################Other
#pair0<-c("Legube_dsbSNP","End_seq_merge_juillet23")


#################define SNPs and bw to use

#for dsb SNPs generated in legube team
pair2<-c("Legube_dsbSNP","pATM_seb_all")
pair3<-c("Legube_dsbSNP","Bless")
pair4<-c("Legube_dsbSNP","BRCA1")
pair5<-c("Legube_dsbSNP","NBN")
pair6<-c("Legube_dsbSNP","sbliss_NEU")
pair7<-c("Legube_dsbSNP","sbliss_NES")
pair8<-c("Legube_dsbSNP","sbliss_NPC")

pairlist<-list(pair2,pair3,pair4,pair5,pair8)

#for dsb SNPs downlaoded on adastra
pair2<-c("ADASTRA_dsb_SNP","pATM_seb_all")
pair3<-c("ADASTRA_dsb_SNP","Bless")
pair4<-c("ADASTRA_dsb_SNP","BRCA1")
pair5<-c("ADASTRA_dsb_SNP","NBN")
pair6<-c("ADASTRA_dsb_SNP","sbliss_NEU")
pair7<-c("ADASTRA_dsb_SNP","sbliss_NES")
pair8<-c("ADASTRA_dsb_SNP","sbliss_NPC")

pairlist<-list(pair2,pair3,pair4,pair5,pair8)


#define parameters to use in deeptools
KB_bp = 10000  #half KB
KB="20kb"#
bin=2000
#bin=200




# KB_bp = format(1000000,digits=15,sci=15) #half KB
# KB="2MB"
# bin=format(10000,digits=15,sci=15)
# #bin=format(1000000,digits=15,sci=15)

# 
# KB_bp = 1000  #half KB
# KB="2kb"
# bin=200

setwd(paste0("/media/sauber/Elements/DSBsnp_project_seb/dsbSNP_Av_pr/"))

for( pair in pairlist){
  print(pair)
  #get big wig
   orig_bw<-list.files(paste0(pair[[2]],"/bw/"),pattern = "bw")
  path<-getwd()
  
   orig_bw_file<-paste0(pair[[2]],"/bw/",orig_bw)
  
  
  
  
  
#get bed files
  pred_bed<-paste0(pair[[1]],"/bed/dsb.bed")
  pred_neg<-paste0(pair[[1]],"/bed/ctrl.bed")
  
  
   orig_bw<-list.files("bw/",pattern = "bw")
   path<-getwd()
#write a script to launch deeptools
  d.sh<-(paste0('
## bigwig file path with inverted commas
## you can do more than one file at once on separate plots by doing: "file1" "file"
## the same for the sample and plot labels: "lab1" "lab2"

## If you just want to have more than one bigwigs on the same plot then you can
## put the files in one set of inverted commas separated by a space e.g. "file1 file2"
## You will then need to do the same for the sample labels (LABS) but not for the plot "lab1 lab2"
## labels (NAMES), as there is still only one name per plot

declare -a FILES=( 

"',orig_bw_file,'"


)

## sample label, one per bigwig
declare -a LABS=(
"',paste0(pair[[2]]),'"

)

## plot label, one per plot 
declare -a NAMES=(
"',paste0(pair[[2]]),'"
)



## the same rules apply if you want to bed files to be on the same plot




## bed file label
declare -a BED=(

"',pred_bed,'"
"',pred_neg,'"


)


## bed file label
declare -a BED_LABS=(
"',paste0(pair[[1]],"_dsbSNPs"),'"
"',paste0(pair[[1]],"_ctrlSNPs"),'"



)


## file regions name (to be used in the plot title)
declare -a FILE_EXT=(
"',paste0(pair[[1]],"_dsbSNPs"),'"
"',paste0(pair[[1]],"_ctrlSNPs"),'"


)

## out file path

## total region to plot from center of file
## write the total region length for plot title, and then the half length in bp for the command itself
KB_bp="',KB_bp,'"
KB=',KB,'

## bin size - i usually try to keep the profiles to have 200-400 bins
bs=',bin,'

## out file path
OUT="',paste0(pair[[1]],"_",pair[[2]],"_deeptools/"),'"
mkdir $OUT

###########################################################
# you dont need to change anything after this line
###########################################################


nrow_bed=$((${#BED[@]}-1))
nrow=$((${#NAMES[@]}-1))
echo $nrow


for bed in $(seq 0 $nrow_bed)
do
FILE_EXT=$(echo "${FILE_EXT[$bed]}")
BED_LABS=$(echo "${BED_LABS[$bed]}")
BED=$(echo "${BED[$bed]}")

echo $FILE_EXT
echo $BED_LABS
echo $BED

for i in $(seq 0 $nrow)
do
FILE=$(echo "${FILES[$i]}")
LAB=$(echo "${LABS[$i]}")
NAME=$(echo "${NAMES[$i]}")

echo $FILE
echo $LAB
echo $NAME

## -a and -b refer to base pairs to plot from the center
computeMatrix reference-point -S ${FILE} \\
-R ${BED} \\
--referencePoint center -a $KB_bp -b $KB_bp \\
--samplesLabel $LAB \\
--outFileName ${OUT}${NAME}_${FILE_EXT}_${KB}.ma.gz \\
--binSize $bs \\
--sortRegions descend \\
-p 11


## Only use this if plotting over genes

#computeMatrix scale-regions -S ${FILE} \
#-R $BED \
#-a 1000 -b 1000 \
#--samplesLabel $LAB \
#--regionBodyLength 1000 \
#--outFileName ${OUT}${NAME}_${FILE_EXT}_${KB}.ma.gz \
#--binSize 10 \
#--sortRegions descend \
#-p 20


## check out the different color options in the deeptools manual 
plotHeatmap -m ${OUT}${NAME}_${FILE_EXT}_${KB}.ma.gz \\
-out ${OUT}${NAME}_${FILE_EXT}_${KB}_Heatmap.pdf \\
--colorMap viridis -x " " --heatmapWidth 10 \\
--sortRegions descend \\
--regionsLabel ${BED_LABS} \\
--whatToShow "heatmap and colorbar" \\
--heatmapHeight 30 #\\
#--outFileSortedRegions ${OUT}${NAME}_${FILE_EXT}_${KB}_sortedRegions.txt


#echo "Plotting $NAME Profile ..."
plotProfile -m ${OUT}${NAME}_${FILE_EXT}_${KB}.ma.gz \\
--outFileName ${OUT}${NAME}_${FILE_EXT}_${KB}_profile.png \\
--regionsLabel ${BED_LABS} \\
--outFileNameData ${OUT}${NAME}_${FILE_EXT}_${KB}_out.prof

#rm ${OUT}${NAME}_${FILE_EXT}_${KB}.ma.gz

done
done
wait






'
  )
  )
#save script
write_lines(d.sh,paste0("d_",th,".sh"))

system("./launch_d.sh")
}
###################################

#Plotting
####################################

require("dplyr")
require("ggplot2")
require("reshape2")
require("rtracklayer")
library("BSgenome.Hsapiens.UCSC.hg19")
library("plyranges")
library("tidyverse")
library("RColorBrewer")
seqlens = seqlengths( Hsapiens );

setwd(paste0("/media/sauber/Elements/DSBsnp_project_seb/dsbSNP_Av_pr/"))

path=(paste0("/media/sauber/Elements/DSBsnp_project_seb/dsbSNP_Av_pr/"))


##################
## Functions for plotting deeptools
#################


plot_group <- function(FILES, vars, regions){
  
  dat.plot <- NULL
  dat.heatmap <- NULL
  
  
  
  message("Reading in data...")
  ## load in data to dfs
  for(n in 1:length(FILES)){
    df <- read.csv(FILES[n], sep="\t", header = FALSE, skip=1)
    rownames(df) <- df$V4
    df <- na.omit(df[,-c(1:6)])
    
    
    sub.dat.plot <- data.frame(
      Window=c(1:ncol(df)),
      Value=colMeans((df)),
      variable = vars[n],
      Peaks=regions[n]
    )
    
    if(is.null(dat.plot)){
      dat.plot <- sub.dat.plot
    }else{
      dat.plot <- rbind(dat.plot,sub.dat.plot)
    }
    
    rownames(df)<- 1:nrow(df)
    
    colnames(df) <- 1:ncol(df)
    melted <- melt(df)
    colnames(melted) <- c("window","value")
    melted$variable <- vars[n]
    melted$Peaks <- regions[n]
    melted$site <- rownames(df)
    
    if(is.null(dat.heatmap)){
      dat.heatmap <- melted
    }else{
      dat.heatmap <- rbind(dat.heatmap,melted)
    }   
  }
  
  return(list(dat.heatmap,dat.plot))
  
} 





##################################

plot_prof <- function(FILES){
  dat.plot <- NULL
  
  for(f in 1:length(FILES)){
    
    df <- read.csv(FILES[f], sep="\t", header = FALSE, skip=1)
    
    for(c in 3:length(df)){
      (print(paste0("c",c)))
      for(i in 2:nrow(df)){
        sub.dat.plot <- data.frame(
          Window=df[1,c],
          Value=df[i,c],
          variable=df[i,1],
          Peaks=df[i,2]
        )
        if(is.null(dat.plot)){
          dat.plot <- sub.dat.plot
        }else{
          dat.plot <- rbind(dat.plot,sub.dat.plot)
        }
      }
    }
  }
  return(dat.plot)
}

KB_bp=as.numeric(KB_bp)
bin=as.numeric(bin)

for(pair in pairlist){
  
  
  plot_title = paste0("Average profil of ", pair[[1]]," Positions \non ",pair[[2]]," signal") # title at top of plot
  y_axis_title=  "Average Normalized read count" # y axis label for average profile
  
  width_hm = 9 # width in inches for heatmap file
  height_hm = 7  # height in inches for heatmap file
  width_avgP = 12  # width in inches for average profile file
  height_avgP = 3.5 # height in inches for average profile file
  
  var_order =  c(pair[[2]])# order to show samples (needs to be exact sample as "vars" input in function)
  regions_order = c("dsbSNPs","ctrlSNPs") # Order to show regions if more than one used. If not, just name the one
  hm_order_var = "dsbsSNPs" # sample used to order heatmap rows (needs to be exact sample as one of the "vars" input in function)
  basepairs = as.numeric(KB_bp)*2 #Total number of basepairs that is used 
 # brewer_palette=palette_blues#YlBu"#"RdPu"#"YlOrRd"#"RdPu" # brewer color palette used for heatmap
  #filename= "/media/scollins/Sarah_scRNA/figures_thomas/bwCompare_ATACA_80DSB_withRPA_separate_sameScale_80kb" # output file name title - do not add suffix 
  hm_cutoff = 0.2# NULL # heatmap cutoff value (sometimes you will need to scale the heatmaps properly by removing large outlier values)
  
  
  size <- paste0((basepairs / 1000) / 2,"kb")
  
  filess=paste0(path,"/",pair[[1]],"_",pair[[2]],"_deeptools/" ,list.files(paste0(pair[[1]],"_",pair[[2]],"_deeptools/"),pattern=paste0(KB,".ma.gz")))
  getwd()
  df<-plot_group(filess,vars = c(pair[[2]],pair[[2]]), 
                 regions=rep(c("ctrlSNPs","dsbSNPs")  ,1))
  dat.heatmap = df[[1]]
  dat.plot = df[[2]]
  dat.plot2<-dat.plot
  
  print(p)
  new_name=pair[[1]]
  
  cols=c("red","blue",'grey60') #set colors
regions_order = c(paste0(new_name,'\npeaks'),paste0(new_name,'\nprediction>0.5'),paste0(new_name,'\nprediction<0.5')) 
  
  dat.plot2$Peaks <- factor(dat.plot2$Peaks,  levels=c(regions_order))
  dat.plot2$variable <- factor(dat.plot2$variable, levels=c(var_order))

  p <-ggplot( na.omit( dat.plot2), aes( Window, Value) ) +
    labs( list( title = "", x = "", y = " " )) +
    geom_line(aes(color=Peaks), size=2, alpha=0.8) +
    theme_classic(base_size = 15) +
    theme(axis.text.y = element_text(size=13),
          axis.text.x = element_text(size=13),
          axis.title.x = element_text(size=17),
          plot.title = element_text(size=20, face="bold", hjust=0.5),
          legend.position = "right") +ylab("Signal")+ #ylab(y_axis_title) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    
    scale_x_continuous(name = paste0("Distance from center (bp)"),
                       breaks = c(1, length(levels(as.factor(dat.plot$Window))) / 2 ,length(levels(as.factor(dat.plot$Window)))),
                       labels = c(paste0("-",size[1]), 'Center', paste0("+",size[1]))) +
    scale_color_manual(values=cols)+
    #geom_hline(yintercept = 0, linetype="dashed") + 
    #scale_color_manual(values=c("sandybrown","tomato3", "cornflowerblue",  "cornflowerblue" )) +
    facet_wrap(~variable) + 
    ggtitle(paste0(plot_title))
  print(p)
  
  
  
  filepdf=paste0("pdf/average_profil_",pair[[2]],"_signal_",pair[[1]],"_",KB,"dsbSNP_ctrlSNP.pdf")
  pdf(filepdf,5,4)
  #pdf("/media/sauber/Elements/DSBsnp_project_seb/breaks_bw/averageprofil_Crosseto_sbliss_big_txt.pdf",10,4)
  print(p)
 dev.off()

}
