library(metafor)
library(orchaRd)

df <- read.csv("ExtractedData_2026-05-12.csv")
# Convert all variances to SD
table(df$Biodiversity_Variance_unit)
df$Control_SD = df$Control_var
df$Treatment_SD = df$Treatment_var
for(i in 1:nrow(df)){
  dfi <- df[i,]
  if(dfi$Biodiversity_Variance_unit == "SE"){

    dfi$Control_SD = dfi$Control_var * sqrt(dfi$Control_n)
    dfi$Treatment_SD = dfi$Treatment_var * sqrt(dfi$Treatment_n)
  }
  df[i,] <- dfi
}
# Note -some of the variances are too high and need to be corrected

# Convert some variables to factors
df$Peatland_type <- factor(df$Peatland_type)
df$FG_Group <- factor(df$FG_Group)

# Combine peatland types
table(df$Peatland_type)
df$Peatland_type <- replace(df$Peatland_type, 
                            grepl(" fen", df$Peatland_type),
                            "Fen")

# Filter to just herbs
# Later I will compare the effect on herbs, trees and mosses
table(df$FG_Group)
herbs <- subset(df,df$FG_Group=="Herb")

# Get effect size
datES <- escalc(measure = "SMDH", # I think this is best, 
                    # since many of the studies use different units and 
                    # measurement techniques and different sample sizes
                m2i = Control_mean, 
                sd2i = Control_SD,
                n2i = Control_n,
                m1i = Treatment_mean,
                sd1i = Treatment_SD,
                n1i = Treatment_n,
                data = herbs)

# Grand mean
grandmean1 <-rma.mv(
  yi=yi,
  V=vi, 
  random= ~1|PaperID/ID,
  method="REML",
  digits=4,
  data=datES)

# Try a one-variable with Peatland type
peat1 <-rma.mv(
  yi=yi,
  V=vi, 
  mods=~Peatland_type, 
  random= ~1|PaperID/ID,
  method="REML",
  digits=4,
  data=datES)

orchard_plot(peat1,
             mod = "Peatland_type",
             group = "PaperID",
             xlab = "Standardised mean difference (with heteroscedastic population variances)")


# Let's try doing all functional group data then separating them as a variable
fg <- subset(df, df$SR_or_FG=="FG")
# remove NA
fg <- fg$FG_Group

fg$FG_Group <- replace(fg$FG_Group, 
                      grepl(" moss", fg$FG_Group),
                      "Moss")
fg <- subset(fg, fg$FG_Group != "")
fg$FG_Group <- factor(fg$FG_Group)
# Get effect size
fgES <- escalc(measure = "SMDH", # I think this is best, 
                # since many of the studies use different units and 
                # measurement techniques and different sample sizes
                m2i = Control_mean, 
                sd2i = Control_SD,
                n2i = Control_n,
                m1i = Treatment_mean,
                sd1i = Treatment_SD,
                n1i = Treatment_n,
                data = fg)
# Try a one-variable with FG
fg1 <-rma.mv(
  yi=yi,
  V=vi, 
  mods=~FG_Group, 
  random= ~1|PaperID/ID,
  method="REML",
  digits=4,
  data=fgES)

orchard_plot(fg1,
             mod = "FG_Group",
             group = "PaperID",
             xlab = "Standardised mean difference (with heteroscedastic population variances)")
# I think this looks pretty good?