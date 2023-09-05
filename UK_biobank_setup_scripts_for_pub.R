load("[UK Biobank RDATA].RData")

####Function to get minimum rows####
rowMin<-function(x){
    if(all(is.na(x))){
        return(NA)
    } else {
        return(min(x,na.rm=TRUE))
    }
}
####################################


####Functions to pull out IDs of people with give ICD10, ICD9, OPCS3, and OPCS4 code(s)####

get_icd10_patients<-function(database,icd10strings){
    subs_cols<-c(which(names(database)=="eid"),grep("diagnoses.+icd10",names(database),perl=TRUE),grep("cancer.+icd10",names(database),perl=TRUE),grep("cause_of_death.+icd10",names(ukb),perl=TRUE))
    subs<-database[,subs_cols]
    for(j in c(2:length(subs_cols))){
        for(i in 1:length(icd10strings)){
            if(j==2 && i==1){
                x<-grep(icd10strings[i],subs[,j])
            } else {
                x<-append(x,grep(icd10strings[i],subs[,j]))
                x<-unique(x)
            }
        }
    }
    x<-sort(x,decreasing=FALSE)
    return(subs[x,1])
    # rm(list=c("subs","x","subs_cols"))
}

get_icd10_dates<-function(database,icd10strings){
    dates<-rep(NA,nrow(database))
    subs_cols<-c(which(names(database)=="eid"),grep("diagnoses.+icd10",names(database),perl=TRUE),grep("cancer.+icd10",names(database),perl=TRUE),grep("cause_of_death.+icd10",names(ukb),perl=TRUE))
    subs<-database[,subs_cols]
    for(j in c(2:length(subs_cols))){
        for(i in 1:length(icd10strings)){
            x<-grep(icd10strings[i],subs[,j])
            if(length(grep("f41202",names(database)[subs_cols[j]]))>0){
                n<-unlist(strsplit(names(database)[subs_cols[j]],"_"))
                n<-n[length(n)]
                new_col<-grep(paste0("f41262_0_",n),names(database))
            } else if(length(grep("f41270",names(database)[subs_cols[j]]))>0){
                n<-unlist(strsplit(names(database)[subs_cols[j]],"_"))
                n<-n[length(n)]
                new_col<-grep(paste0("f41280_0_",n),names(database))
            } else if(length(grep("f40006",names(database)[subs_cols[j]]))>0){
                n<-unlist(strsplit(names(database)[subs_cols[j]],"_"))
                n<-n[length(n)-1]
                new_col<-grep(paste0("f40005_",n,"_0"),names(database))
            } else if(length(grep("f40001",names(database[subs_cols[j]])))>0){
                n<-unlist(strsplit(names(database)[subs_cols[j]],"_"))
                n<-n[length(n)-1]
                new_col<-grep(paste0("f40000_",n,"_0"),names(database))
            }
            if(j==2 && i==1){
                dates[x]<-format(database[x,new_col],format="%Y-%m-%d")
            } else {
                d<-cbind(format(dates[x],format="%Y-%m-%d"),format(database[x,new_col],format="%Y-%m-%d"))
                dates[x]<-apply(d,1,rowMin)
            }
        }
    }
    return(dates)
    # rm(list=c("subs","x","subs_cols"))
}

get_icd9_patients<-function(database,icd9strings){
    subs_cols<-c(which(names(database)=="eid"),grep("diagnoses.+icd9",names(database),perl=TRUE),grep("cancer.+icd9",names(database),perl=TRUE),grep("cause_of_death.+icd9",names(ukb),perl=TRUE))
    subs<-database[,subs_cols]
    for(j in c(2:length(subs_cols))){
        for(i in 1:length(icd9strings)){
            if(j==2 && i==1){
                x<-grep(icd9strings[i],subs[,j])
            } else {
                x<-append(x,grep(icd9strings[i],subs[,j]))
                x<-unique(x)
            }
        }
    }
    x<-sort(x,decreasing=FALSE)
    return(subs[x,1])
    # rm(list=c("subs","x","subs_cols"))
}

get_icd9_dates<-function(database,icd9strings){
    dates<-rep(NA,nrow(database))
    subs_cols<-c(which(names(database)=="eid"),grep("diagnoses.+icd9",names(database),perl=TRUE),grep("cancer.+icd9",names(database),perl=TRUE),grep("cause_of_death.+icd9",names(ukb),perl=TRUE))
    subs<-database[,subs_cols]
    for(j in c(2:length(subs_cols))){
        for(i in 1:length(icd9strings)){
            x<-grep(icd9strings[i],subs[,j])
            if(length(grep("f41203",names(database)[subs_cols[j]]))>0){
                n<-unlist(strsplit(names(database)[subs_cols[j]],"_"))
                n<-n[length(n)]
                new_col<-grep(paste0("f41263_0_",n),names(database))
            } else if(length(grep("f41271",names(database)[subs_cols[j]]))>0){
                n<-unlist(strsplit(names(database)[subs_cols[j]],"_"))
                n<-n[length(n)]
                new_col<-grep(paste0("f41281_0_",n),names(database))
            } else if(length(grep("f40013",names(database)[subs_cols[j]]))>0){
                n<-unlist(strsplit(names(database)[subs_cols[j]],"_"))
                n<-n[length(n)-1]
                new_col<-grep(paste0("f40005_",n,"_0"),names(database))
            }
            if(j==2 && i==1){
                dates[x]<-format(database[x,new_col],format="%Y-%m-%d")
            } else {
                d<-cbind(format(dates[x],format="%Y-%m-%d"),format(database[x,new_col],format="%Y-%m-%d"))
                dates[x]<-apply(d,1,rowMin)
            }
        }
    }
    return(dates)
    # rm(list=c("subs","x","subs_cols"))
}


get_opcs3_patients<-function(database,opcs3strings){
    subs_cols<-c(which(names(database)=="eid"),grep("^operative.+opcs3",names(database),perl=TRUE))
    subs<-database[,subs_cols]
    for(j in c(2:length(subs_cols))){
        for(i in 1:length(opcs3strings)){
            if(j==2 && i==1){
                x<-grep(opcs3strings[i],subs[,j])
            } else {
                x<-append(x,grep(opcs3strings[i],subs[,j]))
                x<-unique(x)
            }
        }
    }
    x<-sort(x,decreasing=FALSE)
    return(subs[x,1])
    # rm(list=c("subs","x","subs_cols"))
}

get_opcs3_dates<-function(database,opcs3strings){
    dates<-rep(NA,nrow(database))
    subs_cols<-c(which(names(database)=="eid"),grep("^operative.+opcs3",names(database),perl=TRUE))
    subs<-database[,subs_cols]
    for(j in c(2:length(subs_cols))){
        for(i in 1:length(opcs3strings)){
            x<-grep(opcs3strings[i],subs[,j])
            if(length(grep("f41256",names(database)[subs_cols[j]]))>0){
                n<-unlist(strsplit(names(database)[subs_cols[j]],"_"))
                n<-n[length(n)]
                new_col<-grep(paste0("f41257_0_",n),names(database))
            } else if(length(grep("f41273",names(database)[subs_cols[j]]))>0){
                n<-unlist(strsplit(names(database)[subs_cols[j]],"_"))
                n<-n[length(n)]
                new_col<-grep(paste0("f41283_0_",n),names(database))
            }
            if(j==2 && i==1){
                dates[x]<-format(database[x,new_col],format="%Y-%m-%d")
            } else {
                d<-cbind(format(dates[x],format="%Y-%m-%d"),format(database[x,new_col],format="%Y-%m-%d"))
                dates[x]<-apply(d,1,rowMin)
            }
        }
    }
    return(dates)
    # rm(list=c("subs","x","subs_cols"))
}


get_opcs4_patients<-function(database,opcs4strings){
    subs_cols<-c(which(names(database)=="eid"),grep("^operative.+opcs4",names(database),perl=TRUE))
    subs<-database[,subs_cols]
    for(j in c(2:length(subs_cols))){
        for(i in 1:length(opcs4strings)){
            if(j==2 && i==1){
                x<-grep(opcs4strings[i],subs[,j])
            } else {
                x<-append(x,grep(opcs4strings[i],subs[,j]))
                x<-unique(x)
            }
        }
    }
    x<-sort(x,decreasing=FALSE)
    return(subs[x,1])
    # rm(list=c("subs","x","subs_cols"))
}

get_opcs4_dates<-function(database,opcs4strings){
    dates<-rep(NA,nrow(database))
    subs_cols<-c(which(names(database)=="eid"),grep("^operative.+opcs4",names(database),perl=TRUE))
    subs<-database[,subs_cols]
    for(j in c(2:length(subs_cols))){
        for(i in 1:length(opcs4strings)){
            x<-grep(opcs4strings[i],subs[,j])
            if(length(grep("f41200",names(database)[subs_cols[j]]))>0){
                n<-unlist(strsplit(names(database)[subs_cols[j]],"_"))
                n<-n[length(n)]
                new_col<-grep(paste0("f41260_0_",n),names(database))
            } else if(length(grep("f41272",names(database)[subs_cols[j]]))>0){
                n<-unlist(strsplit(names(database)[subs_cols[j]],"_"))
                n<-n[length(n)]
                new_col<-grep(paste0("f41282_0_",n),names(database))
            }
            if(j==2 && i==1){
                dates[x]<-format(database[x,new_col],format="%Y-%m-%d")
            } else {
                d<-cbind(format(dates[x],format="%Y-%m-%d"),format(database[x,new_col],format="%Y-%m-%d"))
                dates[x]<-apply(d,1,rowMin)
            }
        }
    }
    return(dates)
    # rm(list=c("subs","x","subs_cols"))
}


bp_meds<-rep(NA,nrow(ukb))
diab_meds<-rep(NA,nrow(ukb))
chol_meds<-rep(NA,nrow(ukb))
for(i in 0:2){
    for(j in 0:2){
        col1<-which(names(ukb)==paste0("medication_for_cholesterol_blood_pressure_or_diabetes_f6177_",i,"_",j))
        col2<-which(names(ukb)==paste0("medication_for_cholesterol_blood_pressure_diabetes_or_take_exogenous_hormones_f6153_",i,"_",j))
        bp_meds[which(ukb[,col1]=="Blood pressure medication")]<-1
        bp_meds[which(ukb[,col2]=="Blood pressure medication")]<-1
        bp_meds[intersect(c(which(!ukb[,col1] %in% c("Do not know","Prefer not to answer")),which(!is.na(ukb[,col1]) & !bp_meds %in% c(1))),which(!bp_meds %in% c(1)))]<-0
        bp_meds[intersect(c(which(!ukb[,col2] %in% c("Do not know","Prefer not to answer")),which(!is.na(ukb[,col2]) & !bp_meds %in% c(1))),which(!bp_meds %in% c(1)))]<-0
        
        diab_meds[which(ukb[,col1]=="Insulin")]<-1
        diab_meds[which(ukb[,col2]=="Insulin")]<-1
        diab_meds[intersect(c(which(!ukb[,col1] %in% c("Do not know","Prefer not to answer")),which(!is.na(ukb[,col1]) & !diab_meds %in% c(1))),which(!diab_meds %in% c(1)))]<-0
        diab_meds[intersect(c(which(!ukb[,col2] %in% c("Do not know","Prefer not to answer")),which(!is.na(ukb[,col2]) & !diab_meds %in% c(1))),which(!diab_meds %in% c(1)))]<-0
        
        chol_meds[which(ukb[,col1]=="Cholesterol lowering medication")]<-1
        chol_meds[which(ukb[,col2]=="Cholesterol lowering medication")]<-1
        chol_meds[intersect(c(which(!ukb[,col1] %in% c("Do not know","Prefer not to answer")),which(!is.na(ukb[,col1]) & !chol_meds %in% c(1))),which(!chol_meds %in% c(1)))]<-0
        chol_meds[intersect(c(which(!ukb[,col2] %in% c("Do not know","Prefer not to answer")),which(!is.na(ukb[,col2]) & !chol_meds %in% c(1))),which(!chol_meds %in% c(1)))]<-0
    }
}


phenos<-c("age_when_attended_assessment_centre_f21003_0_0",
          "month_of_birth_f52_0_0","year_of_birth_f34_0_0",
          "sex_f31_0_0","ethnic_background_f21000_0_0",
          "genetic_sex_f22001_0_0",
          "genetic_ethnic_grouping_f22006_0_0",
          "townsend_deprivation_index_at_recruitment_f189_0_0",
          "average_total_household_income_before_tax_f738_1_0",
          "smoking_status_f20116_0_0",
          "pack_years_of_smoking_f20161_0_0",
          "hip_circumference_f49_0_0",
          "waist_circumference_f48_0_0",
          "medication_for_cholesterol_blood_pressure_diabetes_or_take_exogenous_hormones_f6153_0_0",
          "medication_for_cholesterol_blood_pressure_or_diabetes_f6177_0_0",
          "diabetes_diagnosed_by_doctor_f2443_0_0",
          "age_at_first_live_birth_f2754_0_0",
          "age_heart_attack_diagnosed_f3894_0_0",
          "age_stroke_diagnosed_f4056_0_0",
          "age_angina_diagnosed_f3627_0_0",
          "noncancer_illness_code_selfreported_f20002_0_0",
          "place_of_birth_in_uk_north_coordinate_f129_0_0",
          "place_of_birth_in_uk_east_coordinate_f130_0_0",
          "country_of_birth_ukelsewhere_f1647_0_0",
          "country_of_birth_nonuk_origin_f20115_0_0",
          "home_location_at_assessment_east_coordinate_rounded_f20074_0_0",
          "home_location_at_assessment_north_coordinate_rounded_f20075_0_0",
          "duration_of_vigorous_activity_f914_0_0",
          "number_of_daysweek_of_moderate_physical_activity_10_minutes_f884_0_0",
          "number_of_daysweek_of_vigorous_physical_activity_10_minutes_f904_0_0",
          "duration_of_vigorous_activity_f914_1_0",
          "types_of_physical_activity_in_last_4_weeks_f6164_0_4",
          "time_spent_doing_vigorous_physical_activity_f104900_2_0",
          "time_spent_doing_moderate_physical_activity_f104910_3_0",
          "time_spent_doing_light_physical_activity_f104920_1_0",
          "sleep_duration_f1160_0_0",
          "sleeplessness_insomnia_f1200_0_0",
          "daytime_dozing_sleeping_narcolepsy_f1220_0_0",
          "ldl_direct_f30780_0_0",
          "hdl_cholesterol_f30760_0_0",
          "met_minutes_per_week_for_walking_f22037_0_0",
          "met_minutes_per_week_for_moderate_activity_f22038_0_0",
          "met_minutes_per_week_for_vigorous_activity_f22039_0_0",
          "education_score_england_f26414_0_0",
          "education_score_wales_f26421_0_0",
          "education_score_scotland_f26431_0_0",
          "income_score_england_f26411_0_0",
          "income_score_scotland_f26428_0_0",
          "income_score_wales_f26418_0_0",
          "age_completed_full_time_education_f845_0_0",
          "uk_biobank_assessment_centre_f54_0_0",
          "number_of_pregnancy_terminations_f3849_0_0",
          "date_of_first_operative_procedure_main_opcs3_f41257_0_0",
          "date_lost_to_followup_f191_0_0",
          "date_of_first_operative_procedure_main_opcs4_f41260_0_0",
          "date_of_first_inpatient_diagnosis_main_icd10_f41262_0_0",
          "date_of_first_inpatient_diagnosis_main_icd9_f41263_0_0",
          "date_of_first_inpatient_diagnosis_icd10_f41280_0_0",
          "date_of_first_inpatient_diagnosis_icd9_f41281_0_0",
          "date_of_first_operative_procedure_opcs4_f41282_0_0",
          "date_of_first_operative_procedure_opcs3_f41283_0_0",
          "relative_age_of_first_facial_hair_f2375_0_0",
          "index_of_multiple_deprivation_england_f26410_0_0",
          "index_of_multiple_deprivation_wales_f26426_0_0",
          "index_of_multiple_deprivation_scotland_f26427_0_0",
          "own_or_rent_accommodation_lived_in_f680_0_0",
          "type_of_accommodation_lived_in_f670_0_0",
          "current_employment_status_f6142_0_0",
          "standing_height_f50_0_0",          
          "nitrogen_dioxide_air_pollution_2010_f24003_0_0",
          "nitrogen_oxides_air_pollution_2010_f24004_0_0",
          "particulate_matter_air_pollution_pm10_2010_f24005_0_0",
          "particulate_matter_air_pollution_pm25_2010_f24006_0_0",
          "particulate_matter_air_pollution_pm25_absorbance_2010_f24007_0_0",
          "particulate_matter_air_pollution_2510um_2010_f24008_0_0",
          "nitrogen_dioxide_air_pollution_2005_f24016_0_0",
          "nitrogen_dioxide_air_pollution_2006_f24017_0_0",
          "nitrogen_dioxide_air_pollution_2007_f24018_0_0",
          "particulate_matter_air_pollution_pm10_2007_f24019_0_0",
          "salt_added_to_food_f1478_0_0","country_of_birth_ukelsewhere_f1647_0_0",
          "vascularheart_problems_diagnosed_by_doctor_f6150_0_0",
          "alcohol_intake_frequency_f1558_0_0",
          "alcohol_usually_taken_with_meals_f1618_0_0",
          "alcohol_intake_versus_10_years_previously_f1628_0_0",
          "reason_for_reducing_amount_of_alcohol_drunk_f2664_0_0",
          "former_alcohol_drinker_f3731_0_0",
          "reason_former_drinker_stopped_drinking_alcohol_f3859_0_0",
          "alcohol_drinker_status_f20117_0_0",
          "alcohol_f100022_0_0",
          "alcohol_consumed_f100580_0_0",
          "other_alcohol_intake_f100740_0_0",
          "average_weekly_red_wine_intake_f1568_0_0",
          "average_weekly_champagne_plus_white_wine_intake_f1578_0_0",
          "average_weekly_beer_plus_cider_intake_f1588_0_0",
          "average_weekly_spirits_intake_f1598_0_0",
          "average_weekly_fortified_wine_intake_f1608_0_0",
          "average_weekly_intake_of_other_alcoholic_drinks_f5364_0_0",
          "average_monthly_red_wine_intake_f4407_0_0",
          "average_monthly_champagne_plus_white_wine_intake_f4418_0_0",
          "average_monthly_beer_plus_cider_intake_f4429_0_0",
          "average_monthly_spirits_intake_f4440_0_0",
          "average_monthly_fortified_wine_intake_f4451_0_0",
          "average_monthly_intake_of_other_alcoholic_drinks_f4462_1_0",
          "systolic_blood_pressure_manual_reading_f93_0_0",
          "systolic_blood_pressure_automated_reading_f4080_0_0",
          "diastolic_blood_pressure_manual_reading_f94_2_0",
          "diastolic_blood_pressure_automated_reading_f4079_0_1",
          "body_mass_index_bmi_f21001_0_0",
          "body_mass_index_bmi_f23104_0_0",
          "age_at_death_f40007_0_0")

cols<-c(1)
fields<-c("eid")
for(i in 1:length(phenos)){
    parts<-unlist(strsplit(phenos[i],"_"))
    field<-paste(parts[c(1:(grep("^f\\d+",parts)))],collapse="_")
    fields<-append(fields,field)
    fields<-unique(fields)
    m<-grep(paste0(field,"_"),names(ukb))
    cols<-append(cols,m)
    cols<-unique(cols)
}

cols<-setdiff(cols,grep("age_when_attended_assessment_centre_f21003_[12]_0",names(ukb)))


################################

for(f in fields){
    sf<-grep(f,names(data_subs))
    vals<-rep(NA,nrow(data_subs))
    for(j in sf){
        if(j==sf[1]){
            vals<-as.character(data_subs[,j])
        } else {
            vals[which(!is.na(data_subs[,j]) & !data_subs[,j] %in% c("None of the above","Prefer not to answer","Not sure","Do not know"))]<-as.character(data_subs[which(!is.na(data_subs[,j]) & !data_subs[,j] %in% c("None of the above","Prefer not to answer","Not sure","Do not know")),j])
        }
    }
    assign(f,vals)
    print(f)
}

ethnic_background_f21000[which(ethnic_background_f21000 %in% c("Do not know","Prefer not to answer"))]<-NA
smoking_status_f20116[which(smoking_status_f20116=="Prefer not to answer")]<-NA
smoking_status_f20116<-as.factor(smoking_status_f20116)
smoking_status_f20116<-relevel(smoking_status_f20116,ref=2)
ever_had_stillbirth_spontaneous_miscarriage_or_termination_f2774[which(ever_had_stillbirth_spontaneous_miscarriage_or_termination_f2774=="Prefer not to answer")]<-NA
ever_had_stillbirth_spontaneous_miscarriage_or_termination_f2774[which(ever_had_stillbirth_spontaneous_miscarriage_or_termination_f2774=="Do not know")]<-NA
diastolic_blood_pressure_manual_reading_f94<-as.numeric(diastolic_blood_pressure_manual_reading_f94)
systolic_blood_pressure_manual_reading_f93<-as.numeric(systolic_blood_pressure_manual_reading_f93)
diastolic_blood_pressure_automated_reading_f4079<-as.numeric(diastolic_blood_pressure_automated_reading_f4079)
systolic_blood_pressure_automated_reading_f4080<-as.numeric(systolic_blood_pressure_automated_reading_f4080)
age_when_attended_assessment_centre_f21003<-as.numeric(age_when_attended_assessment_centre_f21003)
age_heart_attack_diagnosed_f3894<-as.numeric(age_heart_attack_diagnosed_f3894)
age_heart_attack_diagnosed_f3894[which(is.na(age_heart_attack_diagnosed_f3894))]<-NA
age_stroke_diagnosed_f4056<-as.numeric(age_stroke_diagnosed_f4056)
age_stroke_diagnosed_f4056[which(age_stroke_diagnosed_f4056<0)]<-NA
age_angina_diagnosed_f3627<-as.numeric(age_angina_diagnosed_f3627)
age_angina_diagnosed_f3627[which(age_angina_diagnosed_f3627<0)]<-NA
age_at_first_live_birth_f2754<-as.numeric(age_at_first_live_birth_f2754)
age_at_first_live_birth_f2754[which(age_at_first_live_birth_f2754<0)]<-NA
age_first_birth<-age_at_first_live_birth_f2754
age_first_birth[which(is.na(age_first_birth))]<-age_of_primiparous_women_at_birth_of_child_f3872[which(is.na(age_first_birth))]
age_first_birth<-as.numeric(age_first_birth)
age_last_birth<-age_at_last_live_birth_f2764
age_last_birth[which(is.na(age_last_birth))]<-age_of_primiparous_women_at_birth_of_child_f3872[which(is.na(age_last_birth))]
age_last_birth[which(age_last_birth<0)]<-NA
age_last_birth<-as.numeric(age_last_birth)
age_first_birth<-as.numeric(age_first_birth)
age_at_death_f40007<-as.numeric(age_at_death_f40007)
dob<-as.Date(paste(month_of_birth_f52,year_of_birth_f34,"15"),format="%B %Y %d")
age<-as.numeric(as.Date("2020-08-31",format="%Y-%m-%d")-dob)/365.25
age[which(!is.na(age_at_death_f40007))]<-age_at_death_f40007[which(!is.na(age_at_death_f40007))]
age_completed_full_time_education_f845<-as.numeric(age_completed_full_time_education_f845)

body_mass_index_bmi_f21001<-as.numeric(body_mass_index_bmi_f21001)
body_mass_index_bmi_f23104<-as.numeric(body_mass_index_bmi_f23104)
diabetes_diagnosed_by_doctor_f2443[which(diabetes_diagnosed_by_doctor_f2443 %in% c("Do not know","Prefer not to answer"))]<-NA

income_score_england_f26411<-as.numeric(income_score_england_f26411)
income_score_wales_f26418<-as.numeric(income_score_wales_f26418)
income_score_scotland_f26428<-as.numeric(income_score_scotland_f26428)
income_score<-income_score_england_f26411
income_score[which(!is.na(income_score_scotland_f26428))]<-income_score_scotland_f26428[which(!is.na(income_score_scotland_f26428))]
income_score[which(!is.na(income_score_wales_f26418))]<-income_score_wales_f26418[which(!is.na(income_score_wales_f26418))]

index_of_multiple_deprivation_england_f26410<-as.numeric(index_of_multiple_deprivation_england_f26410)
index_of_multiple_deprivation_wales_f26426<-as.numeric(index_of_multiple_deprivation_wales_f26426)
index_of_multiple_deprivation_scotland_f26427<-as.numeric(index_of_multiple_deprivation_scotland_f26427)
index_of_multiple_deprivation<-index_of_multiple_deprivation_england_f26410
index_of_multiple_deprivation[which(!is.na(index_of_multiple_deprivation_wales_f26426))]<-index_of_multiple_deprivation_wales_f26426[which(!is.na(index_of_multiple_deprivation_wales_f26426))]
index_of_multiple_deprivation[which(!is.na(index_of_multiple_deprivation_scotland_f26427))]<-index_of_multiple_deprivation_scotland_f26427[which(!is.na(index_of_multiple_deprivation_scotland_f26427))]

education_score_england_f26414<-as.numeric(education_score_england_f26414)
education_score_wales_f26421<-as.numeric(education_score_wales_f26421)
education_score_scotland_f26431<-as.numeric(education_score_scotland_f26431)
education_score<-education_score_england_f26414
education_score[which(!is.na(education_score_wales_f26421))]<-education_score_wales_f26421[which(!is.na(education_score_wales_f26421))]
education_score[which(!is.na(education_score_scotland_f26431))]<-education_score_scotland_f26431[which(!is.na(education_score_scotland_f26431))]
pack_years_of_smoking_f20161<-as.numeric(pack_years_of_smoking_f20161)
pack_years_of_smoking_f20161[which(smoking_status_f20116=="Never")]<-0

salt_added_to_food_f1478[which(salt_added_to_food_f1478=="Prefer not to answer")]<-NA

whr<-as.numeric(waist_circumference_f48)/as.numeric(hip_circumference_f49)
dbp<-diastolic_blood_pressure_automated_reading_f4079
sbp<-systolic_blood_pressure_automated_reading_f4080

townsend_deprivation_index_at_recruitment_f189<-as.numeric(townsend_deprivation_index_at_recruitment_f189)
ever_smoked<-rep(0,nrow(data_subs))
ever_smoked[which(smoking_status_f20116 %in% c("Current","Previous"))]<-1
ever_smoked[which(is.na(smoking_status_f20116))]<-NA
diabetes<-rep(0,nrow(data_subs))
diabetes[which(diabetes_diagnosed_by_doctor_f2443=="Yes")]<-1
diabetes[which(is.na(diabetes_diagnosed_by_doctor_f2443))]<-NA

income<-rep(NA,length(average_total_household_income_before_tax_f738))
income[which(average_total_household_income_before_tax_f738=="Less than 18,000")]<-1
income[which(average_total_household_income_before_tax_f738=="18,000 to 30,999")]<-2
income[which(average_total_household_income_before_tax_f738=="31,000 to 51,999")]<-3
income[which(average_total_household_income_before_tax_f738=="52,000 to 100,000")]<-4
income[which(average_total_household_income_before_tax_f738=="Greater than 100,000")]<-5

qualifications<-ukb[,grep("_f6138_",names(ukb))]
for(j in 1:ncol(qualifications)){
    qualifications[,j]<-as.character(qualifications[,j])
}
nas<-apply(qualifications,1,count_nas)

qualification_none<-rep(NA,nrow(qualifications))
qualification_none[which(nas<ncol(qualifications))]<-0
qualification_none[which(qualifications[,1]=="None of the above")]<-1
qualification_none<-qualification_none[match(eid,ukb$eid)]

qualification_olevels<-rep(NA,nrow(qualifications))
qualification_olevels[which(nas<ncol(qualifications))]<-0
qualification_olevels[which(qualifications=="O levels/GCSEs or equivalent",arr.ind=TRUE)[,1]]<-1
qualification_olevels<-qualification_olevels[match(eid,ukb$eid)]

qualification_cses<-rep(NA,nrow(qualifications))
qualification_cses[which(nas<ncol(qualifications))]<-0
qualification_cses[which(qualifications=="CSEs or equivalent",arr.ind=TRUE)[,1]]<-1
qualification_cses<-qualification_cses[match(eid,ukb$eid)]

qualification_alevels<-rep(NA,nrow(qualifications))
qualification_alevels[which(nas<ncol(qualifications))]<-0
qualification_alevels[which(qualifications=="A levels/AS levels or equivalent",arr.ind=TRUE)[,1]]<-1
qualification_alevels<-qualification_alevels[match(eid,ukb$eid)]

qualification_nvq<-rep(NA,nrow(qualifications))
qualification_nvq[which(nas<ncol(qualifications))]<-0
qualification_nvq[which(qualifications=="NVQ or HND or HNC or equivalent",arr.ind=TRUE)[,1]]<-1
qualification_nvq<-qualification_nvq[match(eid,ukb$eid)]

qualification_other_professional<-rep(NA,nrow(qualifications))
qualification_other_professional[which(nas<ncol(qualifications))]<-0
qualification_other_professional[which(qualifications=="Other professional qualifications eg: nursing, teaching",arr.ind=TRUE)[,1]]<-1
qualification_other_professional<-qualification_other_professional[match(eid,ukb$eid)]

qualification_university<-rep(NA,nrow(qualifications))
qualification_university[which(nas<ncol(qualifications))]<-0
qualification_university[which(qualifications=="College or University degree",arr.ind=TRUE)[,1]]<-1
qualification_university<-qualification_university[match(eid,ukb$eid)]

qualification_school<-qualification_university+qualification_alevels+qualification_olevels
qualification_school[which(qualification_school>1)]<-1
qualification_nonschool<-qualification_nvq+qualification_cses+qualification_other_professional
qualification_nonschool[which(qualification_nonschool>1)]<-1

education_level<-qualification_university
education_level[which(qualification_none==0)]<-1
education_level[which(qualification_other_professional==1)]<-2
education_level[which(qualification_university==1)]<-2

rm(list=c("qualifications","nas"))

own_home<-rep(1,length(eid))
own_home[which(is.na(own_or_rent_accommodation_lived_in_f680))]<-NA
own_home[which(own_or_rent_accommodation_lived_in_f680 %in% c("Rent - from local authority, local council, housing association"))]<-0
own_home[which(own_or_rent_accommodation_lived_in_f680 %in% c("Own with a mortgage","Own outright (by you or someone in your household)"))]<-2

age_completed_full_time_education_f845[which(age_completed_full_time_education_f845 %in% c(-3,-1))]<-NA
age_completed_full_time_education_f845[which(age_completed_full_time_education_f845 %in% c(-3,-1))]<-NA
age_completed_full_time_education_f845[which(age_completed_full_time_education_f845==(-2))]<-0
age_completed_education<-age_completed_full_time_education_f845
age_completed_education[which(is.na(age_completed_full_time_education_f845) & qualification_university==1)]<-22

high_cholesterol<-medication_for_cholesterol_blood_pressure_or_diabetes_f6177
high_cholesterol[which(high_cholesterol %in% c("Do not know","Prefer not to answer"))]<-NA
high_cholesterol[which(high_cholesterol %in% c("None of the above", "Blood pressure medication","Insulin"))]<-0
high_cholesterol[which(medication_for_cholesterol_blood_pressure_diabetes_or_take_exogenous_hormones_f6153 %in% c("None of the above", "Blood pressure medication","Insulin","Hormone replacement therapy","Oral contraceptive pill or minipill"))]<-0
high_cholesterol[which(high_cholesterol=="Cholesterol lowering medication")]<-1
high_cholesterol[which(medication_for_cholesterol_blood_pressure_diabetes_or_take_exogenous_hormones_f6153=="Cholesterol lowering medication")]<-1
high_cholesterol<-as.numeric(high_cholesterol)
nitrogen_dioxide_air_pollution_2010_f24003<-as.numeric(nitrogen_dioxide_air_pollution_2010_f24003)
nitrogen_oxides_air_pollution_2010_f24004<-as.numeric(nitrogen_oxides_air_pollution_2010_f24004)
particulate_matter_air_pollution_pm10_2010_f24005<-as.numeric(particulate_matter_air_pollution_pm10_2010_f24005)
particulate_matter_air_pollution_pm25_2010_f24006<-as.numeric(particulate_matter_air_pollution_pm25_2010_f24006)
particulate_matter_air_pollution_pm25_absorbance_2010_f24007<-as.numeric(particulate_matter_air_pollution_pm25_absorbance_2010_f24007)
particulate_matter_air_pollution_2510um_2010_f24008<-as.numeric(particulate_matter_air_pollution_2510um_2010_f24008)

nitrogen_dioxide_air_pollution_2005_f24016<-as.numeric(nitrogen_dioxide_air_pollution_2005_f24016)
nitrogen_dioxide_air_pollution_2006_f24017<-as.numeric(nitrogen_dioxide_air_pollution_2006_f24017)
nitrogen_dioxide_air_pollution_2007_f24018<-as.numeric(nitrogen_dioxide_air_pollution_2007_f24018)
particulate_matter_air_pollution_pm10_2007_f24019<-as.numeric(particulate_matter_air_pollution_pm10_2007_f24019)

met_minutes_per_week_for_walking_f22037<-as.numeric(met_minutes_per_week_for_walking_f22037)
met_minutes_per_week_for_moderate_activity_f22038<-as.numeric(met_minutes_per_week_for_moderate_activity_f22038)
met_minutes_per_week_for_vigorous_activity_f22039<-as.numeric(met_minutes_per_week_for_vigorous_activity_f22039)

sum_met_minutes_per_week<-met_minutes_per_week_for_walking_f22037+met_minutes_per_week_for_moderate_activity_f22038+met_minutes_per_week_for_vigorous_activity_f22039

###Calculating number of alcohol units consumed per year###
average_weekly_red_wine_intake_f1568<-as.numeric(average_weekly_red_wine_intake_f1568)
average_weekly_red_wine_intake_f1568[which(average_weekly_red_wine_intake_f1568 %in% c(-1,-3))]<-NA

average_weekly_champagne_plus_white_wine_intake_f1578<-as.numeric(average_weekly_champagne_plus_white_wine_intake_f1578)
average_weekly_champagne_plus_white_wine_intake_f1578[which(average_weekly_champagne_plus_white_wine_intake_f1578 %in% c(-1,-3))]<-NA

average_weekly_beer_plus_cider_intake_f1588<-as.numeric(average_weekly_beer_plus_cider_intake_f1588)
average_weekly_beer_plus_cider_intake_f1588[which(average_weekly_beer_plus_cider_intake_f1588 %in% c(-1,-3))]<-NA

average_weekly_spirits_intake_f1598<-as.numeric(average_weekly_spirits_intake_f1598)
average_weekly_spirits_intake_f1598[which(average_weekly_spirits_intake_f1598 %in% c(-1,-3))]<-NA

average_weekly_fortified_wine_intake_f1608<-as.numeric(average_weekly_fortified_wine_intake_f1608)
average_weekly_fortified_wine_intake_f1608[which(average_weekly_fortified_wine_intake_f1608 %in% c(-1,-3))]<-NA

average_weekly_intake_of_other_alcoholic_drinks_f5364<-as.numeric(average_weekly_intake_of_other_alcoholic_drinks_f5364)
average_weekly_intake_of_other_alcoholic_drinks_f5364[which(average_weekly_intake_of_other_alcoholic_drinks_f5364 %in% c(-1,-3))]<-NA

average_monthly_red_wine_intake_f4407<-as.numeric(average_monthly_red_wine_intake_f4407)
average_monthly_red_wine_intake_f4407[which(average_monthly_red_wine_intake_f4407 %in% c(-1,-3))]<-NA

average_monthly_champagne_plus_white_wine_intake_f4418<-as.numeric(average_monthly_champagne_plus_white_wine_intake_f4418)
average_monthly_champagne_plus_white_wine_intake_f4418[which(average_monthly_champagne_plus_white_wine_intake_f4418 %in% c(-1,-3))]<-NA

average_monthly_beer_plus_cider_intake_f4429<-as.numeric(average_monthly_beer_plus_cider_intake_f4429)
average_monthly_beer_plus_cider_intake_f4429[which(average_monthly_beer_plus_cider_intake_f4429 %in% c(-1,-3))]<-NA

average_monthly_spirits_intake_f4440<-as.numeric(average_monthly_spirits_intake_f4440)
average_monthly_spirits_intake_f4440[which(average_monthly_spirits_intake_f4440 %in% c(-1,-3))]<-NA

average_monthly_fortified_wine_intake_f4451<-as.numeric(average_monthly_fortified_wine_intake_f4451)
average_monthly_fortified_wine_intake_f4451[which(average_monthly_fortified_wine_intake_f4451 %in% c(-1,-3))]<-NA

average_monthly_intake_of_other_alcoholic_drinks_f4462<-as.numeric(average_monthly_intake_of_other_alcoholic_drinks_f4462)
average_monthly_intake_of_other_alcoholic_drinks_f4462[which(average_monthly_intake_of_other_alcoholic_drinks_f4462 %in% c(-1,-3))]<-NA

alcohol_units_per_year<-rep(0,length(eid))

alcohol_units_per_year[which(!is.na(average_weekly_red_wine_intake_f1568))]<-alcohol_units_per_year[which(!is.na(average_weekly_red_wine_intake_f1568))]+((1.5*average_weekly_red_wine_intake_f1568[which(!is.na(average_weekly_red_wine_intake_f1568))])*52)

alcohol_units_per_year[which(!is.na(average_weekly_champagne_plus_white_wine_intake_f1578))]<-alcohol_units_per_year[which(!is.na(average_weekly_champagne_plus_white_wine_intake_f1578))]+((1.5*average_weekly_champagne_plus_white_wine_intake_f1578[which(!is.na(average_weekly_champagne_plus_white_wine_intake_f1578))])*52)

alcohol_units_per_year[which(!is.na(average_weekly_beer_plus_cider_intake_f1588))]<-alcohol_units_per_year[which(!is.na(average_weekly_beer_plus_cider_intake_f1588))]+((3*average_weekly_beer_plus_cider_intake_f1588[which(!is.na(average_weekly_beer_plus_cider_intake_f1588))])*52)

alcohol_units_per_year[which(!is.na(average_weekly_spirits_intake_f1598))]<-alcohol_units_per_year[which(!is.na(average_weekly_spirits_intake_f1598))]+((1.1*average_weekly_spirits_intake_f1598[which(!is.na(average_weekly_spirits_intake_f1598))])*52)

alcohol_units_per_year[which(!is.na(average_weekly_fortified_wine_intake_f1608))]<-alcohol_units_per_year[which(!is.na(average_weekly_fortified_wine_intake_f1608))]+((1.25*average_weekly_fortified_wine_intake_f1608[which(!is.na(average_weekly_fortified_wine_intake_f1608))])*52)

alcohol_units_per_year[which(!is.na(average_weekly_intake_of_other_alcoholic_drinks_f5364))]<-alcohol_units_per_year[which(!is.na(average_weekly_intake_of_other_alcoholic_drinks_f5364))]+((1.5*average_weekly_intake_of_other_alcoholic_drinks_f5364[which(!is.na(average_weekly_intake_of_other_alcoholic_drinks_f5364))])*52)

alcohol_units_per_year[which(!is.na(average_monthly_red_wine_intake_f4407))]<-alcohol_units_per_year[which(!is.na(average_monthly_red_wine_intake_f4407))]+((1.5*average_monthly_red_wine_intake_f4407[which(!is.na(average_monthly_red_wine_intake_f4407))])*12)

alcohol_units_per_year[which(!is.na(average_monthly_champagne_plus_white_wine_intake_f4418))]<-alcohol_units_per_year[which(!is.na(average_monthly_champagne_plus_white_wine_intake_f4418))]+((1.5*average_monthly_champagne_plus_white_wine_intake_f4418[which(!is.na(average_monthly_champagne_plus_white_wine_intake_f4418))])*12)

alcohol_units_per_year[which(!is.na(average_monthly_beer_plus_cider_intake_f4429))]<-alcohol_units_per_year[which(!is.na(average_monthly_beer_plus_cider_intake_f4429))]+((3*average_monthly_beer_plus_cider_intake_f4429[which(!is.na(average_monthly_beer_plus_cider_intake_f4429))])*12)

alcohol_units_per_year[which(!is.na(average_monthly_spirits_intake_f4440))]<-alcohol_units_per_year[which(!is.na(average_monthly_spirits_intake_f4440))]+((1.1*average_monthly_spirits_intake_f4440[which(!is.na(average_monthly_spirits_intake_f4440))])*12)

alcohol_units_per_year[which(!is.na(average_monthly_fortified_wine_intake_f4451))]<-alcohol_units_per_year[which(!is.na(average_monthly_fortified_wine_intake_f4451))]+((1.25*average_monthly_fortified_wine_intake_f4451[which(!is.na(average_monthly_fortified_wine_intake_f4451))])*12)

alcohol_units_per_year[which(!is.na(average_monthly_intake_of_other_alcoholic_drinks_f4462))]<-alcohol_units_per_year[which(!is.na(average_monthly_intake_of_other_alcoholic_drinks_f4462))]+((1.5*average_monthly_intake_of_other_alcoholic_drinks_f4462[which(!is.na(average_monthly_intake_of_other_alcoholic_drinks_f4462))])*12)

alcohol_units_per_year[which(alcohol_intake_frequency_f1558=="Prefer not to answer" & alcohol_units_per_year==0)]<-NA

number_pregnancies<-rowSums(cbind(number_of_live_births_f2734,number_of_spontaneous_miscarriages_f3839,number_of_stillbirths_f3829,number_of_pregnancy_terminations_f3849),na.rm=TRUE)
number_pregnancies[which(sex_f31=="Male")]<-NA

####Physical activity####
duration_of_vigorous_activity_f914<-as.numeric(duration_of_vigorous_activity_f914)
duration_of_vigorous_activity_f914[which(duration_of_vigorous_activity_f914 %in% c(-1,-3))]<-NA

number_of_daysweek_of_moderate_physical_activity_10_minutes_f884<-as.numeric(number_of_daysweek_of_moderate_physical_activity_10_minutes_f884)
number_of_daysweek_of_moderate_physical_activity_10_minutes_f884[which(number_of_daysweek_of_moderate_physical_activity_10_minutes_f884 %in% c(-3,-1))]<-NA

number_of_daysweek_of_vigorous_physical_activity_10_minutes_f904<-as.numeric(number_of_daysweek_of_vigorous_physical_activity_10_minutes_f904)
number_of_daysweek_of_vigorous_physical_activity_10_minutes_f904[which(number_of_daysweek_of_vigorous_physical_activity_10_minutes_f904 %in% c(-3,-1))]<-NA

types_of_physical_activity_in_last_4_weeks_f6164[which(types_of_physical_activity_in_last_4_weeks_f6164 %in% c("None of the above","Prefer not to answer"))]<-NA

####Distance between birth place and currently living####
home_location_at_assessment_east_coordinate_rounded_f20074<-as.numeric(home_location_at_assessment_east_coordinate_rounded_f20074)
home_location_at_assessment_north_coordinate_rounded_f20075<-as.numeric(home_location_at_assessment_north_coordinate_rounded_f20075)
place_of_birth_in_uk_east_coordinate_f130<-as.numeric(place_of_birth_in_uk_east_coordinate_f130)
place_of_birth_in_uk_north_coordinate_f129<-as.numeric(place_of_birth_in_uk_north_coordinate_f129)
home_location_at_assessment_east_coordinate_rounded_f20074[which(home_location_at_assessment_east_coordinate_rounded_f20074==(-1))]<-NA
home_location_at_assessment_north_coordinate_rounded_f20075[which(home_location_at_assessment_north_coordinate_rounded_f20075==(-1))]<-NA
place_of_birth_in_uk_north_coordinate_f129[which(place_of_birth_in_uk_north_coordinate_f129==(-1))]<-NA
place_of_birth_in_uk_east_coordinate_f130[which(place_of_birth_in_uk_east_coordinate_f130==(-1))]<-NA

birth_home_distance<-sqrt(((home_location_at_assessment_east_coordinate_rounded_f20074-place_of_birth_in_uk_east_coordinate_f130)^2)+((home_location_at_assessment_north_coordinate_rounded_f20075-place_of_birth_in_uk_north_coordinate_f129)^2))

#Diabetes from hospital records
diabetes_icd10<-rep(0,nrow(data_subs))
strings<-c("E1[0-4]","02[0-3]")
x<-get_icd10_patients(ukb,strings)
m<-match(x,data_subs$eid)
diabetes_icd10[m]<-1

#MI from hospital records
mi<-rep(0,nrow(data_subs))
strings<-c("I2[1-3]","I252")
x<-get_icd10_patients(ukb,strings)
m<-match(x,data_subs$eid)
mi[m]<-1

#Stroke from hospital records
stroke<-rep(0,nrow(data_subs))
strings<-c("I6[0-4]")
x<-get_icd10_patients(ukb,strings)
m<-match(x,data_subs$eid)
stroke[m]<-1

#Hypertension from hospital records
hypertension<-rep(0,nrow(data_subs))
strings<-c("I1[05]")
x<-get_icd10_patients(ukb,strings)
m<-match(x,data_subs$eid)
hypertension[m]<-1

#Pre-eclampsia from hospital records
pe<-rep(0,nrow(data_subs))
strings<-c("O1[456]")
x<-get_icd10_patients(ukb,strings)
m<-match(x,data_subs$eid)
pe[m]<-1

#Endometriosis from hospital records
endo<-rep(0,nrow(data_subs))
strings<-c("N80")
x<-get_icd10_patients(ukb,strings)
m<-match(x,data_subs$eid)
endo[m]<-1

#Hypercholesteremia                 
hyperchol<-rep(0,nrow(data_subs))
strings<-c("E780")
x<-get_icd10_patients(ukb,strings)
m<-match(x,data_subs$eid)
hyperchol[m]<-1

#Angina                 
angina<-rep(0,nrow(data_subs))
strings<-c("I20[0-9]")
x<-get_icd10_patients(ukb,strings)
m<-match(x,data_subs$eid)
angina[m]<-1

#CAD strict definition
cad_strict<-rep(0,nrow(data_subs))
age_cad_strict<-rep(NA,nrow(data_subs))
field_cols<-grep("_f6150_",names(ukb))
for(nh in c("Prefer not to answer",NA)){
    cad_strict[unique(which(ukb[which(ukb$eid %in% data_subs$eid),field_cols]==nh,arr.ind=TRUE)[,1])]<-NA
}
for(nh in c("None of the above","High blood pressure","Angina","Stroke")){
    cad_strict[unique(which(ukb[which(ukb$eid %in% data_subs$eid),field_cols]==nh,arr.ind=TRUE)[,1])]<-0
}
cad_strict[unique(which(ukb[which(ukb$eid %in% data_subs$eid),field_cols]=="Heart attack",arr.ind=TRUE)[,1])]<-1
field_cols<-grep("_f3894_",names(ukb))
age_cad_strict<-as.numeric(apply(ukb[which(ukb$eid %in% data_subs$eid),field_cols],1,rowMin))

field_cols<-grep("_f20002_",names(ukb))
self_report_cad_save<-which(ukb[which(ukb$eid %in% data_subs$eid),field_cols]=="1075",arr.ind=TRUE)
grp<-unique(self_report_cad_save[,1])
cad_strict[grp]<-1
field_cols<-grep("_f20009_",names(ukb))
for(i in grp){
    self_report_cols<-self_report_cad_save[which(self_report_cad_save[,1]==i),2]
    tmp_age<-setdiff(ukb[which(ukb$eid %in% data_subs$eid)[i],field_cols[min(self_report_cols)]],c(-1,-3))
    if(is.null(tmp_age)){
        tmp_age<-NA
    }
    if(all(is.na(c(tmp_age,age_cad_strict[i])))){
        next
    }
    age_cad_strict[i]<-min(c(age_cad_strict[i],tmp_age),na.rm=TRUE)
}

field_cols<-grep("_f20004_",names(ukb))
field_cols2<-grep("_f92_",names(ukb))
for(sg in c("1087","1095","1581")){
    self_report_cad_save<-which(ukb[which(ukb$eid %in% data_subs$eid),field_cols]==sg,arr.ind=TRUE)
    grp<-unique(self_report_cad_save[,1])
    cad_strict[grp]<-1
    for(i in grp){
        self_report_cols<-self_report_cad_save[which(self_report_cad_save[,1]==i),2]
        tmp_age<-setdiff(ukb[which(ukb$eid %in% data_subs$eid)[i],field_cols2[min(self_report_cols)]],c(-1,-3))
        if(is.null(tmp_age)){
            tmp_age<-NA
        }
        if(all(is.na(c(tmp_age,age_cad_strict[i])))){
            next
        }
        age_cad_strict[i]<-min(c(age_cad_strict[i],tmp_age),na.rm=TRUE)
    }
}
strings<-c("I2[1-4]","I252")
x<-get_icd10_patients(ukb,strings)
m<-match(x,data_subs$eid)
cad_strict[m]<-1
x<-get_icd10_dates(ukb,strings)
x<-as.Date(x,format="%Y-%m-%d")
m<-match(ukb$eid,data_subs$eid)
cad_dates<-format(x[m],format="%Y-%m-%d")

strings<-c("41[0-2]")
x<-get_icd9_patients(ukb,strings)
m<-match(x,data_subs$eid)
cad_strict[m]<-1
x<-get_icd9_dates(ukb,strings)
x<-as.Date(x,format="%Y-%m-%d")
m<-match(ukb$eid,data_subs$eid)
cad_dates[m]<-apply(cbind(format(as.Date(cad_dates,format="%Y-%m-%d"),format="%Y-%m-%d"),format(x[m],format="%Y-%m-%d")),1,rowMin)


strings<-c("K4[0-6]","K49","K501","K75")
x<-get_opcs4_patients(ukb,strings)
m<-match(x,data_subs$eid)
cad_strict[m]<-1
x<-get_opcs4_dates(ukb,strings)
x<-as.Date(x,format="%Y-%m-%d")
m<-match(ukb$eid,data_subs$eid)
cad_dates[m]<-apply(cbind(format(as.Date(cad_dates,format="%Y-%m-%d"),format="%Y-%m-%d"),format(x[m],format="%Y-%m-%d")),1,rowMin)

cad_ages<-(as.Date(cad_dates,format="%Y-%m-%d")-dob)/365.25

age_cad_strict[which(age_cad_strict>1000)]<-age_cad_strict[which(age_cad_strict>1000)]-as.numeric(year_of_birth_f34[which(age_cad_strict>1000)])
age_cad_strict[which(age_cad_strict<0)]<-NA

field_cols<-grep("_f40001_",names(ukb))
for(cd in c(paste0("I2",1:4),"I252")){
    self_report_cad_save<-which(ukb[which(ukb$eid %in% data_subs$eid),field_cols]==cd,arr.ind=TRUE)
    grp<-unique(self_report_cad_save[,1])
    cad_strict[grp]<-1
    for(i in grp){
        age_cad_strict[i]<-min(c(age_cad_strict[i],age_at_death_f40007[i]),na.rm=TRUE)
    }
}
age_cad_strict<-apply(cbind(age_cad_strict,cad_ages),1,rowMin)


#CAD loose definition
cad_loose<-cad_strict
field_cols<-grep("_f6150_",names(ukb))
cad_loose[unique(which(ukb[which(ukb$eid %in% data_subs$eid),field_cols]=="Angina",arr.ind=TRUE)[,1])]<-1
strings<-c("I20[0-9]","I20")
x<-get_icd10_patients(ukb,strings)
m<-match(x,data_subs$eid)
cad_loose[m]<-1
strings<-c("I25[01]","I700","I7000")
x<-get_icd10_patients(ukb,strings)
m<-match(x,data_subs$eid)
cad_loose[m]<-1

strings<-c("41[0-4]")
x<-get_icd9_patients(ukb,strings)
m<-match(x,data_subs$eid)
cad_loose[m]<-1

stillbirth_ever<-rep(NA,length(w))
stillbirth_ever[which(number_of_stillbirths_f3829==0)]<-0
stillbirth_ever[which(number_of_stillbirths_f3829>=1)]<-1
stillbirth_ever[which(ever_had_stillbirth_spontaneous_miscarriage_or_termination_f2774=="No" & number_of_live_births_f2734>=1)]<-0

miscarriage_ever<-rep(NA,length(w))
miscarriage_ever[which(number_of_spontaneous_miscarriages_f3839==0)]<-0
miscarriage_ever[which(number_of_spontaneous_miscarriages_f3839>=1)]<-1
miscarriage_ever[which(ever_had_stillbirth_spontaneous_miscarriage_or_termination_f2774=="No" & number_of_live_births_f2734>=1)]<-0

misc_lim<-2
multiple_miscarriage_ever<-rep(NA,length(w))
multiple_miscarriage_ever[which(number_of_spontaneous_miscarriages_f3839<misc_lim)]<-0
multiple_miscarriage_ever[which(number_of_spontaneous_miscarriages_f3839>=misc_lim)]<-1


#Go through stroke subtypes
for(i in 0:4){
    s<-rep(0,nrow(data_subs))
    x<-get_icd10_patients(ukb,paste0("I6",i))
    m<-match(x,data_subs$eid)
    s[m]<-1
    assign(paste0("s_i6",i),s)
}

assessment_centers<-read.table("[UK Biobank Assessmenet Center mapping tab-separated file]",header=TRUE,as.is=TRUE,sep="\t")
for(n in unique(uk_biobank_assessment_centre_f54)){
    uk_biobank_assessment_centre_f54[which(uk_biobank_assessment_centre_f54==n)]<-assessment_centers$meaning[which(assessment_centers$coding==n)]
}


#Load the results of the principal component analysis we performed on the white British subset of the UK Biobank data.
pca<-read.table("[White British PCA file]",header=TRUE,sep="\t",as.is=TRUE)
pca<-pca[,2:ncol(pca)]
names(pca)[which(names(pca)=="IID")]<-"ID"
pca<-pca[match(eid,pca$ID),]
pca$ID<-eid

prop_variance<-read.table(paste0("[White British PCA proportion of variance]"),header=FALSE,as.is=TRUE)[,1]


#Get the genotyping array (0 is Axiom, 1 is BiLEVE)
geno_array<-ukb$genotype_measurement_batch_f22000_0_0
geno_array[which(geno_array>0)]<-0
geno_array[which(geno_array<0)]<-1
geno_array<-geno_array[match(eid,ukb$eid)]

ukb_pca<-ukb[match(eid,ukb$eid),grep("f22009",names(ukb))]

wb<-read.table("[White British IDs file]")[,1]
sa<-read.table("[South Asian IDs file]")[,1]
chinese<-read.table("[Chinese IDs file]")[,1]
white<-read.table("[Other White IDs file]")[,1]
bac<-read.table("[Black, African, and Caribbean IDs file]")[,1]
irish<-read.table("[Irish IDs file]")[,1]


####Load other risk scores####

#metaGRS
for(chr in 1:22){
    tmp<-read.table(paste0("[Path to metaGRS scores/MetaGRS score per chromosome file name]",chr,".txt"),header=TRUE,as.is=TRUE)
    if(chr==1){
        grs<-tmp$SCORESUM
        names<-tmp$IID
    } else {
        grs<-grs+tmp$SCORESUM[match(tmp$IID,names)[which(!is.na(match(tmp$IID,names)))]]
    }
}
names<-as.character(names)
meta_grs<-grs[match(eid,names)]

meta_grs_norm<-rep(NA,length(eid))
meta_grs_norm[which(eid %in% wb)]<-(meta_grs[which(eid %in% wb)]-mean(meta_grs[which(eid %in% wb)],na.rm=TRUE))/sd(meta_grs[which(eid %in% wb)],na.rm=TRUE)

#K2018
for(chr in 1:22){
    tmp<-read.table(paste0("[Path to K2018 scores/K2018 score per chromosome file name]",chr,".txt"),header=TRUE,as.is=TRUE)
    if(chr==1){
        k2018_grs<-tmp$SCORESUM
        names<-tmp$IID
    } else {
        k2018_grs<-k2018_grs+tmp$SCORESUM[match(tmp$IID,names)]
    }
}
names<-as.character(names)
k2018_grs<-k2018_grs[match(eid,names)]

k2018_grs_norm<-rep(NA,length(eid))
k2018_grs_norm[which(eid %in% wb)]<-(k2018_grs[which(eid %in% wb)]-mean(k2018_grs[which(eid %in% wb)],na.rm=TRUE))/sd(k2018_grs[which(eid %in% wb)],na.rm=TRUE)

#E2020
for(chr in 1:22){
    tmp<-read.table(paste0("[Path to E2020 scores/E2020 score per chromosome file name]",chr,".txt"),header=TRUE,as.is=TRUE)
    if(chr==1){
        e2020_grs<-tmp$SCORESUM
        names<-tmp$IID
    } else {
        e2020_grs<-e2020_grs+tmp$SCORESUM[match(tmp$IID,names)]
    }
}
names<-as.character(names)
e2020_grs<-e2020_grs[match(eid,names)]*(-1)

e2020_grs_norm<-rep(NA,length(eid))
e2020_grs_norm[which(eid %in% wb)]<-(e2020_grs[which(eid %in% wb)]-mean(e2020_grs[which(eid %in% wb)],na.rm=TRUE))/sd(e2020_grs[which(eid %in% wb)],na.rm=TRUE)
