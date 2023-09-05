source("[PATHWAY]/UK_biobank_setup_scripts_for_pub.R")
library("RColorBrewer")


#### Correlating environmental variables to PRS and CAD ####
w<-eid %in% wb

model<-glm(meta_grs_norm[w]~scale(age[w])+as.factor(sex_f31[w])+scale(pca$PC1[w])+scale(pca$PC2[w])+scale(pca$PC3[w])+as.factor(smoking_status_f20116[w])+scale(townsend_deprivation_index_at_recruitment_f189[w])+scale(income[w])+scale(age_completed_education[w])+scale(home_location_at_assessment_north_coordinate_rounded_f20075[w])+scale(home_location_at_assessment_east_coordinate_rounded_f20074[w])+scale(place_of_birth_in_uk_north_coordinate_f129[w])+scale(place_of_birth_in_uk_east_coordinate_f130[w]),family="gaussian")
s<-summary(model)
s$coeff[grep("pca",rownames(s$coeff)),c(1,4)]

model<-glm(k2018_grs_norm[w]~scale(age[w])+as.factor(sex_f31[w])+scale(pca$PC1[w])+scale(pca$PC2[w])+scale(pca$PC3[w])+as.factor(smoking_status_f20116[w])+scale(townsend_deprivation_index_at_recruitment_f189[w])+scale(income[w])+scale(age_completed_education[w])+scale(home_location_at_assessment_north_coordinate_rounded_f20075[w])+scale(home_location_at_assessment_east_coordinate_rounded_f20074[w])+scale(place_of_birth_in_uk_north_coordinate_f129[w])+scale(place_of_birth_in_uk_east_coordinate_f130[w]),family="gaussian")
s<-summary(model)
s$coeff[grep("pca",rownames(s$coeff)),c(1,4)]

model<-glm(e2020_grs_norm[w]~scale(age[w])+as.factor(sex_f31[w])+scale(pca$PC1[w])+scale(pca$PC2[w])+scale(pca$PC3[w])+as.factor(smoking_status_f20116[w])+scale(townsend_deprivation_index_at_recruitment_f189[w])+scale(income[w])+scale(age_completed_education[w])+scale(home_location_at_assessment_north_coordinate_rounded_f20075[w])+scale(home_location_at_assessment_east_coordinate_rounded_f20074[w])+scale(place_of_birth_in_uk_north_coordinate_f129[w])+scale(place_of_birth_in_uk_east_coordinate_f130[w]),family="gaussian")
s<-summary(model)
s$coeff[grep("pca",rownames(s$coeff)),c(1,4)]

model<-glm(cad_strict[w]~scale(age[w])+as.factor(sex_f31[w])+scale(pca$PC1[w])+scale(pca$PC2[w])+scale(pca$PC3[w])+as.factor(smoking_status_f20116[w])+scale(townsend_deprivation_index_at_recruitment_f189[w])+scale(income[w])+scale(age_completed_education[w])+scale(home_location_at_assessment_north_coordinate_rounded_f20075[w])+scale(home_location_at_assessment_east_coordinate_rounded_f20074[w])+scale(place_of_birth_in_uk_north_coordinate_f129[w])+scale(place_of_birth_in_uk_east_coordinate_f130[w]),family="binomial")
s<-summary(model)
s$coeff[grep("pca",rownames(s$coeff)),c(1,4)]



model<-glm(meta_grs_norm[w]~scale(pca$PC1[w])+scale(pca$PC2[w])+scale(pca$PC3[w]),family="gaussian")
s<-summary(model)
s$coeff[grep("pca",rownames(s$coeff)),c(1,4)]

model<-glm(k2018_grs_norm[w]~scale(pca$PC1[w])+scale(pca$PC2[w])+scale(pca$PC3[w]),family="gaussian")
s<-summary(model)
s$coeff[grep("pca",rownames(s$coeff)),c(1,4)]

model<-glm(e2020_grs_norm[w]~scale(pca$PC1[w])+scale(pca$PC2[w])+scale(pca$PC3[w]),family="gaussian")
s<-summary(model)
s$coeff[grep("pca",rownames(s$coeff)),c(1,4)]




#### Looking at data by UK Biobank Assessment Center ####

w<-which(eid %in% wb)
risk_score<-"e2020_grs"
c1<-"Reading"
c2<-"Newcastle"
ks.test(get(risk_score)[intersect(w,which(uk_biobank_assessment_centre_f54==c1))],get(risk_score)[intersect(w,which(uk_biobank_assessment_centre_f54==c2))])
t.test(get(risk_score)[intersect(w,which(uk_biobank_assessment_centre_f54==c1))],get(risk_score)[intersect(w,which(uk_biobank_assessment_centre_f54==c2))])
var.test(get(risk_score)[intersect(w,which(uk_biobank_assessment_centre_f54==c1))],get(risk_score)[intersect(w,which(uk_biobank_assessment_centre_f54==c2))])


risk_score_norm_c1<-scale(get(risk_score)[intersect(w,which(uk_biobank_assessment_centre_f54==c1))])
length(which(!is.na(risk_score_norm_c1)))
mc1<-mean(get(risk_score)[intersect(w,which(uk_biobank_assessment_centre_f54==c1))],na.rm=TRUE)
sdc1<-sd(get(risk_score)[intersect(w,which(uk_biobank_assessment_centre_f54==c1))],na.rm=TRUE)
quantile(get(risk_score)[intersect(w,which(uk_biobank_assessment_centre_f54==c1))],0.9,na.rm=TRUE)
tsc1<-(test_score-mc1)/sdc1
length(which(risk_score_norm_c1>=tsc1))/length(!is.na(risk_score_norm_c1))

risk_score_norm_c2<-scale(get(risk_score)[intersect(w,which(uk_biobank_assessment_centre_f54==c2))])
mc2<-mean(get(risk_score)[intersect(w,which(uk_biobank_assessment_centre_f54==c2))],na.rm=TRUE)
sdc2<-sd(get(risk_score)[intersect(w,which(uk_biobank_assessment_centre_f54==c2))],na.rm=TRUE)
tsc2<-(test_score-mc2)/sdc2
length(which(risk_score_norm_c2>=tsc2))/length(!is.na(risk_score_norm_c2))



w<-which(eid %in% wb)
c<-"Reading"
for(risk_score in c("meta_grs","k2018_grs","e2020_grs")){
    l<-length(which(!is.na(get(risk_score)[intersect(w,which(uk_biobank_assessment_centre_f54==c))])))
    m<-mean(get(risk_score)[intersect(w,which(uk_biobank_assessment_centre_f54==c))],na.rm=TRUE)
    s<-sd(get(risk_score)[intersect(w,which(uk_biobank_assessment_centre_f54==c))],na.rm=TRUE)
    q10<-quantile(get(risk_score)[intersect(w,which(uk_biobank_assessment_centre_f54==c))],0.9,na.rm=TRUE)
    q5<-quantile(get(risk_score)[intersect(w,which(uk_biobank_assessment_centre_f54==c))],0.95,na.rm=TRUE)
    print(c(m,s,q10,q5))
}
print(l)

#Renormalize other scores for the above assessment center
for(c2 in c("Reading","Cardiff","Newcastle","Glasgow")){
    for(risk_score in c("meta_grs","k2018_grs","e2020_grs")){
        m<-mean(get(risk_score)[intersect(w,which(uk_biobank_assessment_centre_f54==c))],na.rm=TRUE)
        s<-sd(get(risk_score)[intersect(w,which(uk_biobank_assessment_centre_f54==c))],na.rm=TRUE)
        q10<-quantile(get(risk_score)[intersect(w,which(uk_biobank_assessment_centre_f54==c))],0.9,na.rm=TRUE)
        q5<-quantile(get(risk_score)[intersect(w,which(uk_biobank_assessment_centre_f54==c))],0.95,na.rm=TRUE)
        
        l2<-length(which(!is.na(get(risk_score)[intersect(w,which(uk_biobank_assessment_centre_f54==c2))])))
        m2<-(mean(get(risk_score)[intersect(w,which(uk_biobank_assessment_centre_f54==c2))],na.rm=TRUE)-m)/s
        s2<-sd(get(risk_score)[intersect(w,which(uk_biobank_assessment_centre_f54==c2))],na.rm=TRUE)/s
        # q102<-(quantile(get(risk_score)[intersect(w,which(uk_biobank_assessment_centre_f54==c2))],0.9,na.rm=TRUE)-m)/s
        q102<-length(which(get(risk_score)[intersect(w,which(uk_biobank_assessment_centre_f54==c2))]>=q10))/l2
        # q52<-(quantile(get(risk_score)[intersect(w,which(uk_biobank_assessment_centre_f54==c2))],0.95,na.rm=TRUE)-m)/s
        q52<-length(which(get(risk_score)[intersect(w,which(uk_biobank_assessment_centre_f54==c2))]>=q5))/l2
        print(c(q52))
    }
}

for(risk_score in c("meta_grs","k2018_grs","e2020_grs")){
    m<-mean(get(risk_score)[intersect(w,which(uk_biobank_assessment_centre_f54==c))],na.rm=TRUE)
    s<-sd(get(risk_score)[intersect(w,which(uk_biobank_assessment_centre_f54==c))],na.rm=TRUE)
    q10<-quantile(get(risk_score)[intersect(w,which(uk_biobank_assessment_centre_f54==c))],0.9,na.rm=TRUE)
    q5<-quantile(get(risk_score)[intersect(w,which(uk_biobank_assessment_centre_f54==c))],0.95,na.rm=TRUE)
    
    l2<-length(which(!is.na(get(risk_score)[w])))
    m2<-(mean(get(risk_score)[w],na.rm=TRUE)-m)/s
    s2<-sd(get(risk_score)[w],na.rm=TRUE)/s
    n1<-length(which(get(risk_score)[w]>=q10))/l2
    n2<-length(which(get(risk_score)[w]>=q5))/l2
    print(c(n2))
}
print(l)



n<-quantile(get(risk_score)[intersect(w,which(uk_biobank_assessment_centre_f54==c))],0.9,na.rm=TRUE)

w<-which(eid %in% wb)
for(c in c("Reading","Cardiff","Newcastle","Glasgow")){
    for(risk_score in c("meta_grs","k2018_grs","e2020_grs")){
        n1<-quantile(get(risk_score)[intersect(w,which(uk_biobank_assessment_centre_f54==c))],0.9,na.rm=TRUE)
        n2<-quantile(get(risk_score)[intersect(w,which(uk_biobank_assessment_centre_f54==c))],0.95,na.rm=TRUE)
        print(c(n1,n2))
    }
}


for(risk_score in c("meta_grs","k2018_grs","e2020_grs")){
    n1<-quantile(get(risk_score)[w],0.9,na.rm=TRUE)
    n2<-quantile(get(risk_score)[w],0.95,na.rm=TRUE)
    print(c(n1,n2))
}



#### PCA plots ####

colors<-rep(rgb(0.7,0.7,0.7),nrow(data_subs))
uk_birth_colors<-brewer.pal(10, "Paired")[c(2,3,4,6,10)]
tmp_col<-rbind(col2rgb(uk_birth_colors)/255,rep(1,5))
uk_birth_colors<-apply(tmp_col,2,function(x){return(rgb(x[1],x[2],x[3],x[4]))})

colors[which(country_of_birth_ukelsewhere_f1647=="England")]<-uk_birth_colors[4]
colors[which(country_of_birth_ukelsewhere_f1647=="Scotland")]<-uk_birth_colors[1]
colors[which(country_of_birth_ukelsewhere_f1647=="Wales")]<-uk_birth_colors[5]
colors[which(country_of_birth_ukelsewhere_f1647=="Northern Ireland")]<-uk_birth_colors[3]
colors[which(country_of_birth_ukelsewhere_f1647=="Repubic of Ireland")]<-uk_birth_colors[3]

pc_a<-1
pc_b<-3
leg.text<-c("England","Scotland","Wales","Ireland","Non-UK/Unknown")
plot_cols<-c(uk_birth_colors[c(4,1,5,3)],rgb(0.7,0.7,0.7))
png(paste0("[PATH]/UK_Biobank_white_British_PCA_PC",pc_a,"_vs_PC",pc_b,"_colored_by_country_of_birth_within_British_Isles_all.png"),height=1000,width=1000)
par(mar=c(5, 5, 2, 2)+0.1)
plot(pca[,pc_a+1],pca[,pc_b+1],col=colors,xlab=paste0("PC",pc_a," (",signif(signif(prop_variance[pc_a]*100,4)),"% of total variance)"),ylab=paste0("PC",pc_b," (",signif(signif(prop_variance[pc_b]*100,4)),"% of total variance)"),lwd=4,cex.axis=2,cex.lab=2,cex=2,xlim=c(min(pca[,pc_a+1],na.rm=TRUE),max(pca[,pc_a+1],na.rm=TRUE)),pch=1)
if(pc_b==3){
    legpos<-"topleft"
} else {
    legpos<-"bottomright"
}
legend(legpos,legend=leg.text,col=plot_cols,pch=c(rep(1,length(plot_cols))),pt.lwd=4,cex = 2)
dev.off()



#### Plotting by place of birth and coloring by PCA ####

#Color by first three PCs
rescale_min<-0.1
rescale_max<-0.7

rescaled_pc1<-(((-pca$PC1)-min(-pca$PC1,na.rm=TRUE))/((range(-pca$PC1,na.rm=TRUE)[2]-range(-pca$PC1,na.rm=TRUE)[1]))*(rescale_max-rescale_min))+rescale_min
rescaled_pc2<-((-pca$PC2-min(-pca$PC2,na.rm=TRUE))/((range(-pca$PC2,na.rm=TRUE)[2]-range(-pca$PC2,na.rm=TRUE)[1]))*(rescale_max-rescale_min))+rescale_min
rescaled_pc3<-((pca$PC3-min(pca$PC3,na.rm=TRUE))/((range(pca$PC3,na.rm=TRUE)[2]-range(pca$PC3,na.rm=TRUE)[1]))*(rescale_max-rescale_min))+rescale_min

place_of_birth_in_uk_east_coordinate_f130[which(place_of_birth_in_uk_east_coordinate_f130==(-1))]<-NA
place_of_birth_in_uk_north_coordinate_f129[which(place_of_birth_in_uk_north_coordinate_f129==(-1))]<-NA

# colors<-rep(rgb(0.8,0.8,0.8,0.2),nrow(data_subs))
colors<-rep(NA,nrow(data_subs))
colors_transparent<-rep(NA,nrow(data_subs))
colors_transparent[which(!is.na(rescaled_pc1))]<-rgb(0.8-rescaled_pc1[which(!is.na(rescaled_pc1))],rescaled_pc2[which(!is.na(rescaled_pc1))],0.8-rescaled_pc3[which(!is.na(rescaled_pc1))],0.1)
colors[which(!is.na(rescaled_pc1))]<-rgb(0.8-rescaled_pc1[which(!is.na(rescaled_pc1))],rescaled_pc2[which(!is.na(rescaled_pc1))],0.8-rescaled_pc3[which(!is.na(rescaled_pc1))])


plot_xlim<-range(place_of_birth_in_uk_east_coordinate_f130,na.rm=TRUE)+c(-100,100)
plot_ylim<-range(place_of_birth_in_uk_north_coordinate_f129,na.rm=TRUE)+c(-100,100)

png("[PATH]/Map_of_UK_Biobank_participants_place_of_birth.png",height=2034,width=1000)
plot(place_of_birth_in_uk_east_coordinate_f130,place_of_birth_in_uk_north_coordinate_f129,col=colors,lwd=4,cex=2,axes=FALSE,xlab=NA,ylab=NA,xlim=plot_xlim,ylim=plot_ylim,pch=1)
dev.off()

png("[PATH]/Map_of_UK_Biobank_participants_place_of_birth_with_density.png",height=2034,width=1000)
plot(place_of_birth_in_uk_east_coordinate_f130,place_of_birth_in_uk_north_coordinate_f129,col=colors_transparent,lwd=4,cex=2,axes=FALSE,xlab=NA,ylab=NA,xlim=plot_xlim,ylim=plot_ylim,pch=19)
dev.off()

w<-eid %in% wb
model<-glm(meta_grs_norm[w]~pca$PC1[w]+pca$PC2[w]+pca$PC3[w])
summary(model)



#### Loading missingness data and plotting ####

xlab_names<-c("MetaGRS", "K2018", "E2020")
rs<-c("metaGRS_hg19","k2018","e2020")
new_names<-c("meta_grs","k2018_","e2020")

for(risk_score in rs){
    nname<-new_names[which(c("metaGRS_hg19","K2018","E2020")==risk_score)]
    for(missingness in c(0.01,0.02,0.05,0.1)){
        for(draw in 1:10){
            tab<-read.table(paste0("[Path to missingness results/Start of missingness results file name]",nname,"_missingness_",missingness*100,"_draw_",draw,".txt"),header=TRUE,as.is=FALSE,sep="\t")
            assign(paste0(nname,"_missingness_",missingness*100,"_draw_",draw,"_tab"),tab)
        }
    }
}

for(rs in c("metaGRS","K2018","E2020")){
    nname<-new_names[which(rs==c("metaGRS","K2018","E2020"))]
    extremes<-read.table(paste0("[Path to missingness output/",rs,"_score_top_100_and_bottom_100_scores_white_British.txt"),header=FALSE,as.is=TRUE)[,1]
    assign(paste0(nname,"_ext"),extremes)
}

for(score in c("meta_grs","k2018","e2020")){
    for(missingness in c(0.01,0.02,0.05,0.1)){
        for(draw in 1:10){
            tab<-get(paste0(score,"_missingness_",missingness*100,"_draw_",draw,"_tab"))
            m<-match(get(paste0(score,"_ext")),tab[,1])
            if(draw==1){
                sum_missing_scores<-tab[m,2]
            } else {
                sum_missing_scores<-sum_missing_scores+tab[m,2]
            }
        }
        assign(paste0(score,"_ext_mean_missing_",missingness*100),sum_missing_scores/10)
    }
}

#### In the white British ####

w<-which(eid %in% wb)
wb_ukb<-eid[w]


## Normalize the missing scores ##
for(score in c("meta_grs","k2018","e2020")){
    ext<-get(paste0(score,"_ext"))
    for(missingness in c(0.01,0.02,0.05,0.1)){
        for(draw in 1:10){
            tab<-get(paste0(score,"_missingness_",missingness*100,"_draw_",draw,"_tab"))
            m<-match(wb_ukb,tab[,1])
            tab<-tab[m,]
            if(score=="meta_grs"){
                assign(paste0("meta_grs_norm_missing_",missingness*100,"_draw_",draw),(tab[,2]-mean(tab[,2],na.rm=TRUE))/sd(tab[,2],na.rm=TRUE))
                assign(paste0("meta_grs_norm_missing_",missingness*100,"_draw_",draw,"_ext"),get(paste0("meta_grs_norm_missing_",missingness*100,"_draw_",draw))[which(wb_ukb %in% ext)])
            } else {
                assign(paste0(score,"_grs_norm_missing_",missingness*100,"_draw_",draw),(tab[,2]-mean(tab[,2],na.rm=TRUE))/sd(tab[,2],na.rm=TRUE))
                assign(paste0(score,"_grs_norm_missing_",missingness*100,"_draw_",draw,"_ext"),get(paste0("meta_grs_norm_missing_",missingness*100,"_draw_",draw))[which(wb_ukb %in% ext)])
            }
        }
    }
}


for(score in c("meta_grs","k2018","e2020")){
    ext<-get(paste0(score,"_ext"))
    assign(paste0(score,"ind_ext"),which(wb_ukb %in% ext))
    for(missingness in c(0.01,0.02,0.05,0.1)){
        for(draw in 1:10){
            if(score=="meta_grs"){
                s<-get(paste0(score,"_norm_missing_",missingness*100,"_draw_",draw))[which(wb_ukb %in% ext)]
            } else {
                s<-get(paste0(score,"_grs_norm_missing_",missingness*100,"_draw_",draw))[which(wb_ukb %in% ext)]
            }
            if(draw==1){
                sum_score<-s
            } else {
                sum_score<-sum_score+s
            }
        }
        assign(paste0(score,"_norm_missing_",missingness*100,"_mean_ext"),sum_score/10)
    }
}

mat_pos<-matrix(0,ncol=3,nrow=5)
mat_neg<-matrix(0,ncol=3,nrow=5)

for(i in 1:nrow(mat_pos)){
    for(j in 1:ncol(mat_pos)){
        score<-c("meta_grs","k2018","e2020")[j]
        ext<-get(paste0(score,"_ext"))
        if(i==1){
            if(score=="meta_grs"){
                mat_pos[i,j]<-mean(get(paste0(score,"_norm"))[which(wb_ukb %in% ext & get(paste0(score,"_norm"))>0)],na.rm=TRUE)
                mat_neg[i,j]<-mean(get(paste0(score,"_norm"))[which(wb_ukb %in% ext & get(paste0(score,"_norm"))<0)],na.rm=TRUE)
            } else {
                mat_pos[i,j]<-mean(get(paste0(score,"_grs_norm"))[which(wb_ukb %in% ext & get(paste0(score,"_grs_norm"))>0)],na.rm=TRUE)
                mat_neg[i,j]<-mean(get(paste0(score,"_grs_norm"))[which(wb_ukb %in% ext & get(paste0(score,"_grs_norm"))<0)],na.rm=TRUE)
            }
        } else {
            missingness<-c(0.01,0.02,0.05,0.1)[i-1]
            mat_pos[i,j]<-mean(get(paste0(score,"_norm_missing_",missingness*100,"_mean_ext"))[which(get(paste0(score,"_norm_missing_",missingness*100,"_mean_ext"))>0)],na.rm=TRUE)
            mat_neg[i,j]<-mean(get(paste0(score,"_norm_missing_",missingness*100,"_mean_ext"))[which(get(paste0(score,"_norm_missing_",missingness*100,"_mean_ext"))<0)],na.rm=TRUE)
        }
        
    }
}


plot_cols<-brewer.pal(5,"Set1")

pdf("[PATH]/Barplot_of_mean_scores_over_missingness_for_highest_and_lowest_100_people.pdf",height=5,width=11)
par(mar=c(2,4,0.5,4.5),xpd=TRUE)
barplot(mat_pos,beside=TRUE,ylim=c(-4,4),main=NA,names.arg=c("metaGRS","K2018","E2020"),col=c(plot_cols),xpd = TRUE,axes=FALSE,ylab="Mean normalized score")
barplot(mat_neg,beside=TRUE,add=TRUE,col=c(plot_cols))
axis(2,at=seq(-4,4,by=0.5),labels=c(-4,NA,-3,NA,-2,NA,-1,NA,0,NA,1,NA,2,NA,3,NA,4))
legend(18.5,4,legend=c("Full",0.01,0.02,0.05,0.1),fill=plot_cols)
dev.off()

