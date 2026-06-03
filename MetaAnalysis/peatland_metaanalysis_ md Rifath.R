## ============================================================================
##  FINAL META-ANALYSIS SCRIPT
##  Effect of anthropogenic drainage on vascular plant functional groups
##  in boreal peatlands
##
##  Authors:  Md. Rifath Ahamed 
##  Course:   EEB 310, University of Helsinki, Spring 2026
##  Data:     ExtractedData_2026-05-15.csv
##
# install.packages(c("metafor","tidyverse"))   # uncomment if first time
library(metafor)
library(tidyverse)
## 1. LOAD DATA

dat <- read.csv("ExtractedData_2026-05-15.csv",
  stringsAsFactors = FALSE,
  fileEncoding = "latin1" #using Latin-1 encoding so Scandinavian author
  # like the names (Köster, etc.) read correctly.
)

cat("Loaded", nrow(dat), "rows from", length(unique(dat$PaperID)), "papers\n")
cat("Columns:", ncol(dat), "\n")


## ----------------------------------------------------------------------------
## 2. CLEAN COLUMNS
##    - standardise text fields
##    - coerce numeric fields to numeric (text/blank -> NA)
## ----------------------------------------------------------------------------

# Standardise the SR_or_FG flag (sometimes has white spaces / different  case)
dat$SR_or_FG <- toupper(gsub("[[:space:]]", "", dat$SR_or_FG))

# Trim whitespace on text moderators
dat$FG_Group <- trimws(dat$FG_Group)
dat$Peatland_type <- trimws(dat$Peatland_type)
dat$First_author <- trimws(dat$First_author)
dat$Exp_or_Obs <- trimws(dat$Exp_or_Obs)

# Coerce the numeric columns
numeric_cols <- c(
  "Control_mean", "Control_n", "Control_var",
  "Treatment_mean", "Treatment_n", "Treatment_var",
  "Control_WTD_mean", "Treatment_WTD_mean", "WTD_Difference"
)
numeric_cols <- intersect(numeric_cols, names(dat)) # only those that exist
# Bulk-coerce targeted character columns to numeric using list-apply
dat[numeric_cols] <- lapply(dat[numeric_cols], as.numeric)

# Quick visibility check on what's there
cat("Variance unit breakdown:\n")
print(table(dat$Biodiversity_Variance_unit, useNA = "always"))
cat("SR_or_FG breakdown:\n")
print(table(dat$SR_or_FG, useNA = "always"))
## ----------------------------------------------------------------------------
## 3. STANDARDISE SPREAD TO SD
##    escalc() needs SD. Papers report SE, SD, or variance.
##    Conversions:
##      SE  -> SD = SE * sqrt(n)
##      SD  -> SD = SD            (unchanged)
##      var -> SD = sqrt(var)

convert_to_SD <- function(value, n, unit) {
  
  if (is.na(value) | is.na(unit)) {
    return(NA_real_)
  }
  #to ensure there are no human data entry error(safety filter)
  u <- toupper(trimws(unit))
  #Now that the unit string is standardized to u, the function runs through a series of logical check
  if (u == "SE") { #If the unit is exactly "SE", R recognizes the value as a Standard Error.
    return(value * sqrt(n)) # then it converts se to sd using $SD = SE \times \sqrt{n}$.
  }
  if (u == "SD") {#
    return(value)
  }
  warning("Unknown variance unit '", unit, "' - treating as SD")
  return(value)
}
#Creates a brand-new column in dataframe and populates it with the completely standardized, uniform $SD$ values for control groups.
dat$Control_SD <- mapply(
  convert_to_SD, dat$Control_var,
  dat$Control_n, dat$Biodiversity_Variance_unit
)
#Creates a brand-new column in dataframe and populates it with the completely standardized, uniform $SD$ values for treatment groups.
dat$Treatment_SD <- mapply(
  convert_to_SD, dat$Treatment_var,
  dat$Treatment_n, dat$Biodiversity_Variance_unit
)

# datapoint_id: unique row identifier for the nested random effect
# Renamed to avoid confusion with the Exp_or_Obs column
dat$datapoint_id <- seq_len(nrow(dat))


## 4. SUBSET: FUNCTIONAL GROUPS (main) AND SPECIES RICHNESS (secondary)

# Rows usable for any effect-size calculation
#any row missing a Mean, Standard Deviation ($SD$), or Sample Size ($n$) for 
#either the Control or Treatment group is completely dropped. 
dat_complete <- dat %>%
  filter(
    !is.na(Control_mean), !is.na(Treatment_mean),
    !is.na(Control_SD), !is.na(Treatment_SD),
    !is.na(Control_n), !is.na(Treatment_n),
    #calculating pooled standard deviation or effect sizes involves dividing by
    #zero, which results in a value of infinity (Inf) or an invalid number (NaN).
    Control_SD > 0, Treatment_SD > 0
  )
# prints how many rows survived this  filter 
cat("Complete cases:", nrow(dat_complete), "of", nrow(dat), "\n")
# we got 34 data points that have mean,sd, and n for both control & treatment 

# Functional groups - excluding Laine (see header note)
dat_FG_main <- dat_complete %>%
  filter(
    SR_or_FG == "FG",
    FG_Group %in% c("Tree", "Shrub", "Herb", "Moss"),
    First_author != "Laine"
  )
cat(
  "FG main analysis k =", nrow(dat_FG_main),
  "from", length(unique(dat_FG_main$PaperID)), "papers\n"
)
print(table(dat_FG_main$FG_Group))
#FG main analysis k = 25 from 6 papers
# second subset accourding to peatland type 
print(table(dat_FG_main$Peatland_type))

# Species richness - separate analysis
dat_SR <- dat_complete %>%
  filter(SR_or_FG == "SR")
cat(
  "SR analysis k =", nrow(dat_SR),
  "from", length(unique(dat_SR$PaperID)), "papers\n"
)
# there are only  3 datapoints in 3 papers 
## ----------------------------------------------------------------------------
## 5. EFFECT SIZES
##    FG -> SMDH  (allows unequal variances between control & treatment)
##    SR -> ROM   (log response ratio: ln(treatment / control))
## ----------------------------------------------------------------------------

dat_FG_main <- escalc(
  measure = "SMDH",
  m1i = Treatment_mean, sd1i = Treatment_SD, n1i = Treatment_n,
  m2i = Control_mean, sd2i = Control_SD, n2i = Control_n,
  data = dat_FG_main, var.names = c("yi", "vi")
)
#Species Richness Analysis: ROM (Secondary Focus)
dat_SR <- escalc(
    measure = "ROM",
    m1i = Treatment_mean, sd1i = Treatment_SD, n1i = Treatment_n,
    m2i = Control_mean, sd2i = Control_SD, n2i = Control_n,
    data = dat_SR, var.names = c("yi", "vi")
)

cat("FG effect-size range:", round(range(dat_FG_main$yi, na.rm = TRUE), 2), "\n")
#FG effect-size range: -3.03 2.89


## 6. GRAND MEAN (intercept-only Multi-Level Random-Effects Meta-Analysis Model)
##    Random effects: ~ 1 | PaperID/datapoint_id
##    Why nested:  multiple effect sizes per paper are non-independent.
##    Why ~0 expected: opposing functional-group responses cancel on average.

grandmean_main <- rma.mv(yi, vi,
  random = ~ 1 | PaperID / datapoint_id,
  data = dat_FG_main, method = "REML"
)
summary(grandmean_main)


## ----------------------------------------------------------------------------
## 7. DIAGNOSTICS
##    Funnel plot: visual check for publication bias.
##      With ~6 studies, NO formal Egger-type test (Cochrane guidance).
##    QQ plot:     check normality of residuals.
## ----------------------------------------------------------------------------

# Save funnel plot
png("funnel_plot.png", width = 700, height = 500, res = 120)
funnel(grandmean_main,
  main = "Funnel plot: drainage effect on functional groups"
)
dev.off()

# Save QQ plot
png("qq_residuals.png", width = 700, height = 500, res = 120)
qqnorm(residuals(grandmean_main, type = "pearson"),
  main = "QQ plot of model residuals"
)
qqline(residuals(grandmean_main, type = "pearson"), col = "red")
dev.off()

cat("Diagnostics saved: funnel_plot.png, qq_residuals.png\n")


## ----------------------------------------------------------------------------
## 8. HETEROGENEITY (I^2 with between/within decomposition)
##    Higgins benchmarks: 25% low | 50% moderate | 75%+ high
## ----------------------------------------------------------------------------

W2 <- diag(1 / grandmean_main$vi)
X2 <- model.matrix(grandmean_main)
P2 <- W2 - W2 %*% X2 %*% solve(t(X2) %*% W2 %*% X2) %*% t(X2) %*% W2

denom <- sum(grandmean_main$sigma2) +
  (grandmean_main$k - grandmean_main$p) / sum(diag(P2))
I2_tot <- 100 * sum(grandmean_main$sigma2) / denom
I2_lvl <- 100 * grandmean_main$sigma2 / denom

cat("Total I^2:        ", round(I2_tot, 1), "%\n")
cat("Between-paper I^2:", round(I2_lvl[1], 1), "%\n")
cat("Within-paper I^2: ", round(I2_lvl[2], 1), "%\n")


## ----------------------------------------------------------------------------
## 9. MODERATOR: FUNCTIONAL GROUP
##    mod_FG          : with intercept  -> QM tests "do groups differ?"
##    mod_FG_nointcp  : without         -> each group's effect vs 0
## ----------------------------------------------------------------------------

mod_FG <- rma.mv(yi, vi,
  mods = ~FG_Group,
  random = ~ 1 | PaperID / datapoint_id,
  data = dat_FG_main, method = "REML"
)
summary(mod_FG) # read the QM(df=3) line for the omnibus test

mod_FG_nointcp <- rma.mv(yi, vi,
  mods = ~ FG_Group - 1,
  random = ~ 1 | PaperID / datapoint_id,
  data = dat_FG_main, method = "REML"
)
summary(mod_FG_nointcp) # each row = that group's effect vs zero


## ----------------------------------------------------------------------------
## 10. MODERATOR: PEATLAND TYPE
##     Caveat: some categories have small k (will print counts).
##     Reported with the n-per-cell caveat in Discussion.
## ----------------------------------------------------------------------------

cat("Sample sizes per peatland type:\n")
print(table(dat_FG_main$Peatland_type))

mod_peat <- rma.mv(yi, vi,
  mods = ~Peatland_type,
  random = ~ 1 | PaperID / datapoint_id,
  data = dat_FG_main, method = "REML"
)
summary(mod_peat) # QM here is the headline peatland-type test

mod_peat_nointcp <- rma.mv(yi, vi,
  mods = ~ Peatland_type - 1,
  random = ~ 1 | PaperID / datapoint_id,
  data = dat_FG_main, method = "REML"
)
summary(mod_peat_nointcp)


## ----------------------------------------------------------------------------
## 11. MODERATOR: WTD MAGNITUDE (continuous)
##     Tests dose-response: does larger drawdown -> larger effect?
##     Interaction with FG: do groups respond differently to magnitude?
##     Rows with missing WTD_Difference are dropped automatically by rma.mv.
## ----------------------------------------------------------------------------

dat_FG_wtd <- dat_FG_main %>% filter(!is.na(WTD_Difference))
cat(
  "WTD analysis k =", nrow(dat_FG_wtd),
  "of", nrow(dat_FG_main), "FG rows have WTD_Difference\n"
)

mod_wtd <- rma.mv(yi, vi,
  mods = ~WTD_Difference,
  random = ~ 1 | PaperID / datapoint_id,
  data = dat_FG_wtd, method = "REML"
)
summary(mod_wtd)

mod_wtd_int <- rma.mv(yi, vi,
  mods = ~ FG_Group * WTD_Difference,
  random = ~ 1 | PaperID / datapoint_id,
  data = dat_FG_wtd, method = "REML"
)
summary(mod_wtd_int)


## ----------------------------------------------------------------------------
## 12. SECONDARY: SPECIES RICHNESS (ROM)
##     Almost certainly underpowered; report descriptively.
## ----------------------------------------------------------------------------

if (nrow(dat_SR) >= 3) {
  sr_grand <- rma.mv(yi, vi,
    random = ~ 1 | PaperID / datapoint_id,
    data = dat_SR, method = "REML"
  )
  summary(sr_grand)
  cat(
    "Species richness analysis based on k =", nrow(dat_SR),
    "effect sizes - INTERPRET CAUTIOUSLY.\n"
  )
} else {
  cat("Insufficient SR data (k =", nrow(dat_SR), "); reported descriptively only.\n")
}


## ============================================================================
## 13. FIGURE 1 -- FOREST PLOT: FUNCTIONAL GROUPS
## ============================================================================

# Pull each group's estimate + CI from the no-intercept model
fg_tab <- as.data.frame(coef(summary(mod_FG_nointcp)))
fg_tab$Group <- gsub("FG_Group", "", rownames(fg_tab))
k_fg <- dat_FG_main %>% count(FG_Group, name = "k")
fg_tab <- left_join(fg_tab, k_fg, by = c("Group" = "FG_Group")) %>%
  select(Group, estimate, ci.lb, ci.ub, k)

# Overall pooled estimate as a summary row
overall_fg <- data.frame(
  Group = "Overall",
  estimate = as.numeric(grandmean_main$b),
  ci.lb = grandmean_main$ci.lb,
  ci.ub = grandmean_main$ci.ub,
  k = nrow(dat_FG_main)
)
fg_plot <- bind_rows(fg_tab, overall_fg)

fg_plot$Group <- factor(fg_plot$Group,
  levels = c("Overall", "Moss", "Herb", "Shrub", "Tree")
)
fg_plot$ypos <- as.numeric(fg_plot$Group)
fg_plot$label <- sprintf(
  "g = %+.2f [%.2f, %.2f] (k=%d)",
  fg_plot$estimate, fg_plot$ci.lb, fg_plot$ci.ub, fg_plot$k
)

# Diamond polygon for the Overall row
od <- fg_plot[fg_plot$Group == "Overall", ]
fg_diamond <- data.frame(
  x = c(od$ci.lb, od$estimate, od$ci.ub, od$estimate),
  y = c(od$ypos, od$ypos + 0.3, od$ypos, od$ypos - 0.3)
)
fg_points <- fg_plot %>% filter(Group != "Overall")

p_fg <- ggplot() +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_errorbarh(
    data = fg_points,
    aes(y = ypos, xmin = ci.lb, xmax = ci.ub), height = 0.15
  ) +
  geom_point(
    data = fg_points,
    aes(y = ypos, x = estimate, size = k, colour = estimate > 0)
  ) +
  geom_polygon(data = fg_diamond, aes(x = x, y = y), fill = "grey20") +
  geom_text(
    data = fg_plot, aes(y = ypos, label = label),
    x = 2.5, hjust = 0, size = 3.2
  ) +
  scale_colour_manual(values = c("FALSE" = "#762a83", "TRUE" = "#1b7837"), guide = "none") +
  scale_size_continuous(range = c(3, 7), name = "k (studies)") +
  scale_y_continuous(
    breaks = 1:nlevels(fg_plot$Group),
    labels = levels(fg_plot$Group)
  ) +
  coord_cartesian(xlim = c(-2.5, 4.5), clip = "off") +
  labs(
    x = "Hedges' g  (negative = drainage decreases, positive = drainage increases)",
    y = NULL, title = "Effect of drainage on vascular plant functional groups"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.y = element_blank(),
    plot.margin = margin(5, 110, 5, 5)
  )

ggsave("forest_FG.png", p_fg, width = 10, height = 4.5, dpi = 300)
cat("Saved figure: forest_FG.png\n")


## ============================================================================
## 14. FIGURE 2 -- FOREST PLOT: PEATLAND TYPE
## ============================================================================

pt_tab <- as.data.frame(coef(summary(mod_peat_nointcp)))


pt_tab$Type <- gsub("Peatland_type", "", rownames(pt_tab))
k_pt <- dat_FG_main %>% count(Peatland_type, name = "k")
pt_tab <- left_join(pt_tab, k_pt, by = c("Type" = "Peatland_type")) %>%
  select(Type, estimate, ci.lb, ci.ub, k)

# Drop any k = 1 cells (uninterpretable single observation).
# Comment out the next line if you want to show them with a caveat.
pt_tab <- pt_tab %>% filter(k > 1)

overall_pt <- data.frame(
  Type = "Overall",
  estimate = as.numeric(grandmean_main$b),
  ci.lb = grandmean_main$ci.lb,
  ci.ub = grandmean_main$ci.ub,
  k = nrow(dat_FG_main)
)
pt_plot <- bind_rows(pt_tab, overall_pt)

ord <- pt_tab$Type[order(pt_tab$estimate)]
pt_plot$Type <- factor(pt_plot$Type, levels = c("Overall", ord))
pt_plot$ypos <- as.numeric(pt_plot$Type)
pt_plot$label <- sprintf(
  "g = %+.2f [%.2f, %.2f] (k=%d)",
  pt_plot$estimate, pt_plot$ci.lb, pt_plot$ci.ub, pt_plot$k
)

od2 <- pt_plot[pt_plot$Type == "Overall", ]
pt_diamond <- data.frame(
  x = c(od2$ci.lb, od2$estimate, od2$ci.ub, od2$estimate),
  y = c(od2$ypos, od2$ypos + 0.3, od2$ypos, od2$ypos - 0.3)
)
pt_points <- pt_plot %>% filter(Type != "Overall")

p_pt <- ggplot() +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_errorbarh(
    data = pt_points,
    aes(y = ypos, xmin = ci.lb, xmax = ci.ub), height = 0.15
  ) +
  geom_point(
    data = pt_points,
    aes(y = ypos, x = estimate, size = k, colour = estimate > 0)
  ) +
  geom_polygon(data = pt_diamond, aes(x = x, y = y), fill = "grey20") +
  geom_text(
    data = pt_plot, aes(y = ypos, label = label),
    x = 2.0, hjust = 0, size = 3.2
  ) +
  scale_colour_manual(values = c("FALSE" = "#762a83", "TRUE" = "#1b7837"), guide = "none") +
  scale_size_continuous(range = c(3, 7), name = "k (studies)") +
  scale_y_continuous(breaks = 1:nlevels(pt_plot$Type),
                     labels = levels(pt_plot$Type)) +
  coord_cartesian(xlim = c(-4, 4), clip = "off") +
  labs(
    x = "Hedges' g  (negative = drainage decreases, positive = drainage increases)",
    y = NULL, title = "Effect of drainage by peatland type"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.y = element_blank(),
    plot.margin = margin(5, 110, 5, 5)
  )

ggsave("forest_peatland.png", p_pt, width = 10, height = 4.5, dpi = 300)
cat("Saved figure: forest_peatland.png\n")





## ============================================================================
##  END OF ANALYSIS
##
##  Key objects -> what they give you for the write-up:
##    grandmean_main    -> overall pooled effect (Results paragraph 2)
##    I2_tot, I2_lvl    -> heterogeneity decomposition (Results paragraph 2)
##    mod_FG            -> QM test: do FGs differ?       (Results par. 3)
##    mod_FG_nointcp    -> per-group effects vs zero      (Results par. 3, Fig 1)
##    mod_peat          -> QM test: does peat type matter?(Results par. 4, headline)
##    mod_peat_nointcp  -> per-type effects vs zero       (Results par. 4, Fig 2)
##    mod_wtd / _int    -> dose-response by WTD magnitude (Results par. 5, Fig 3)
##    sr_grand          -> SR result (descriptive only)   (Results par. 7)
##    funnel_plot.png   -> publication-bias diagnostic    (Results par. 6)
##    qq_residuals.png  -> normality diagnostic           (Results par. 6)
##
##  Run order suggestion when reading line-by-line:
##    Sections 1-5 once (setup + effect sizes)
##    Sections 6-12 individually, pausing to read each summary()
##    Sections 13-15 once each plot is wanted
## ============================================================================
