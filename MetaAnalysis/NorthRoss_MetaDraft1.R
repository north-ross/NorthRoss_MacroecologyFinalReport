library(metafor)

df <- read.csv("ExtractedData_2026-05-12.csv")
# Convert all variances to SD
table(df$Biodiversity_Variance_unit)
df$Control_SD = df$Control_var
df$Treatment_SD = df$Treatment_var
for(i in nrow(df)){
  dfi <- df[i,]
  if(dfi$Biodiversity_Variance_unit == "SE"){

    dfi$Control_SD = dfi$Control_var * sqrt(dfi$Control_n)
    dfi$Control_SD = "TEST"
    dfi$Treatment_SD = dfi$Treatment_var * sqrt(dfi$Treatment_n)
  }
  df[i,] <- dfi
}
# first set all to SD
df$Control_SD <- df$Control_var
# correct the ones that area actually SE
df[df$Biodiversity_Variance_unit == 'SE', "Control_SD"] <- df$

# Filter to just herbs
table(df$FG_Group)
herbs <- subset(df,df$FG_Group=="Herb")


datES <- escalc(measure = "SMD", # 
                m2i = Control_mean, # group 2 corresponds to the control group
                sd2i = Control_SD,
                n2i = Control_n,
                m1i = Treatment_Mean, # group 1 is the treatment group
                sd1i = Treatment_SD,
                n1i = Treatment_n,
                data = df)

grandmean1 <-rma.mv(
  yi=yi,
  V=vi, 
  random= ~1|PaperID/ID,
  method="REML",
  digits=4,
  data=datES)