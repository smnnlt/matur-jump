# Analysis script for the manuscript
# ""
#
# written by Simon Nolte
# 
# script released under an MIT license

# load packages
library(readxl)  # data import
library(ggplot2) # visualization
library(lme4)    # multilevel modeling
library(purrr)   # functional programming
library(emmeans)
library(dplyr)

# read data
d <- read_excel("data/data.xlsx")

# rename columns
colnames(d)[c(5, 10, 22:27)] <- c("cat", "dAPHV", "FJT", "FREQ", "SUs", "SUb", "M4x60", "D4x60")
# convert from tibble to data.frame
d <- as.data.frame(d)

## PART 1.1 exploratory analysis------------------------------------------------

# check relationship between chronological and biological age
ggplot(d, aes(CA, dAPHV)) +
  geom_point(aes(color = sex)) +
  geom_smooth()
# ggsave("plots/CA-BA.png", width = 6, height = 5, dpi = 300, bg = "white")

ggplot(d, aes(CA, dAPHV)) +
  geom_point(aes(color = sex), alpha = 0.5) +
  geom_smooth(aes(color = sex), method = "lm", formula = y ~ poly(x, 2))
# maybe better to fit a second-degree polynomial instead of linear fit
# quite parallel, so we may not need an interaction here
# ggsave("plots/CA-BA_fit.png", width = 6, height = 5, dpi = 300, bg = "white")

ggplot(d, aes(CA, FOST, color = sex)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", formula = y ~ poly(x, 2))
# linear phase when 13<CA<17 but not before and after, so 2nd-degree poly might 
# be better again.
# also different trends, so use interaction effect
# ggsave("plots/CA-FOST_fit.png", width = 6, height = 5, dpi = 300, bg = "white")

# we decided to exclude athletes with small BA from the analysis 
d <- d[d$dAPHV >= 1, ]

## PART 1.2 Descriptive statistics----------------------------------------------

# get all variables
vars <- colnames(d)[12:27]
vars_jump <- vars[c(1:10, 11, 13, 14)]

# descriptive data
desc_long <- d |>
  group_by(sex) |>
  summarise(
    across(
      all_of(vars),
      list(
        n    = ~sum(!is.na(.)),
        mean = ~mean(., na.rm = TRUE),
        sd   = ~sd(., na.rm = TRUE)
      )
    ),
    .groups = "drop"
  ) |>
  tidyr::pivot_longer(
    -c(sex),
    names_to = c("variable", "stat"),
    names_sep = "_"
  ) |>
  tidyr::pivot_wider(
    names_from = stat,
    values_from = value
  )

write.csv(desc_long, "descriptive.csv")

## PART 2: modeling-------------------------------------------------------------

# calculate residuals for biological age
d$dAPHV_resid <- resid(lm(dAPHV ~ poly(CA, 2) + sex, data = d))
# alternative approach: use predicted minus mean cohort specific APHV
d$APHV_resid <- NA
d$APHV_resid[d$sex == "m" & d$cat == "J"] <- mean(d$APHV[d$sex == "m" & d$cat == "J"]) - d$APHV[d$sex == "m" & d$cat == "J"]
d$APHV_resid[d$sex == "f" & d$cat == "J"] <- mean(d$APHV[d$sex == "f" & d$cat == "J"]) - d$APHV[d$sex == "f" & d$cat == "J"]

# function to get model estimates
run_mod <- function(param, data = d, plot = FALSE) {
  # get parameter of interest
  data$out <- data[[param]]
  # compute polynoms before the model, so that emmeans work on the output
  p <- poly(data$CAc, 2)
  data$CAc_p1 <- p[,1]
  data$CAc_p2 <- p[,2]
  o <- lmerTest::lmer(out ~ CAc_p1 + CAc_p2 + sex + CAc_p1:sex + CAc_p2:sex + APHV_resid + APHV_resid:CAc + APHV_resid:sex + APHV_resid:CAc:sex + (1|ID), data = data)
  es <- summary(o)$coefficients
  ci <- confint(o)
  # calculate standardized coefficients
  bs <- effectsize::standardize_parameters(o, method = "basic")
  b <- bs[bs$Parameter == "APHV_resid",2] # Parameter for dAPHV resid
  
  # get reference value to transform to percentages
  m <- mean(data$out, na.rm = TRUE)
  
  # calculate marginal trends
  mt <- emtrends(o, specs = ~ sex | CAc, var = "APHV_resid", at = list(CAc = CAc_values))
  # transform to percentage
  mt <- as.data.frame(mt)
  mt[,c(3,4,6,7)] <- mt[,c(3,4,6,7)] / m 
  
  # calculate marginal means (sex-specific effects, averaged over CA)
  mm <- emtrends(o, specs = ~ sex, var = "APHV_resid")
  mm <- as.data.frame(mm)
  mm[,c(2,3,5,6)] <- mm[,c(2,3,5,6)] / m 
  
  df <- data.frame(
    name = param,
    es = es["APHV_resid", "Estimate"] / m,
    ci_low = ci["APHV_resid",1] / m,
    ci_up = ci["APHV_resid",2] / m,
    beta = b,
    es_int = es["APHV_resid:CAc", "Estimate"] / m,
    ci_low_int = ci["APHV_resid:CAc", 1] / m,
    ci_up_int = ci["APHV_resid:CAc", 2] / m
  )
  
  # option to create scatter plot
  if (plot) {
    ggplot(aes(APHV_resid, out_pct), data = data) +
      geom_hline(yintercept = 0, color = "grey40") +
      geom_vline(xintercept = 0, color = "grey40") +
      geom_point() +
      geom_smooth(aes(color = sex)) +
      scale_y_continuous(
        name = "Overperformance for biological age", 
        labels = scales::label_percent()
      ) +
      scale_x_continuous(
        name = "Difference from predicted biological age (y)"
      ) +
      theme_classic() +
      theme(
        panel.grid.major = element_line(color = "grey80")
      )
    ggsave(
      paste0("plots/scatter_sex_", cat, "_", param, ".png"), 
      width = 6, height = 5, dpi = 300, bg = "white"
    )
  }
  
  attr(df, "mod") <- o
  attr(df, "mt") <- mt
  attr(df, "mm") <- mm
  df
}

# center CA
d$CAc <- scale(d$CA, scale = FALSE)
# get original year values for centered CA
CAc_values <- 13:19 - mean(d$CA) 

# run all models
l_jump <- lapply(vars_jump, run_mod, data = d) |> purrr::list_rbind()

# reverse signs for scales on which lower values mean better performance
l_jump[l_jump$name %in% c("0-10m", "10-30m", "0-30m", "30-60m", "0-60m"), 2:8] <- -1 * l_jump[l_jump$name %in% c("0-10m", "10-30m", "0-30m", "30-60m", "0-60m"), 2:8]
l_jump[l_jump$name %in% c("0-10m", "10-30m", "0-30m", "30-60m", "0-60m"), c(3,4,7,8)] <- l_jump[l_jump$name %in% c("0-10m", "10-30m", "0-30m", "30-60m", "0-60m"), c(4,3,8,7)]

write.csv(l_jump, "effects_jump.csv", row.names = FALSE)

# set signs for beta thresholds
l_jump$beta_label <- ifelse(abs(l_jump$beta) > 0.3, "◊◊", ifelse(abs(l_jump$beta) > 0.1, "◊", ""))
# set levels for interaction

# plot model estimates (jump, without SU tests)
ggplot(l_jump[1:11,], aes(x = reorder(name, es, decreasing = TRUE), y = es)) +
  geom_hline(yintercept = 0, color = "grey40") +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_up), width = 0.2, linewidth = 0.7) +
  geom_text(aes(y = ci_up + 0.013, label = beta_label), size = 6) +
  scale_x_discrete(name = "Variable") +
  scale_y_continuous(
    name = "Performance difference per 1-year\n increase in maturity offset", 
    labels = scales::label_percent(),
    breaks = seq(from = -0.05, to = 0.15, by = 0.05)) +
  theme_classic(base_size = 14) +
  theme(
    panel.grid.major.y = element_line(color = "grey80"),
    axis.text.x = element_text(size = 9),
    axis.text.y = element_text(size = 9)
  )

# save plot
ggsave("plots/effects_jump_cut_three.jpg", width = 7, height = 4.6, dpi = 300, bg = "white")

# plot all model estimates (jump with SU)
ggplot(l_jump, aes(x = reorder(name, es, decreasing = TRUE), y = es)) +
  geom_hline(yintercept = 0, color = "grey40") +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_up), width = 0.2, linewidth = 0.7) +
  geom_text(aes(y = ci_up + 0.013, label = beta_label), size = 6) +
  scale_x_discrete(name = "Variable") +
  scale_y_continuous(
    name = "Performance difference per 1-year\n increase in maturity offset", 
    labels = scales::label_percent(),
    breaks = seq(from = -0.30, to = 0.15, by = 0.10)) +
  theme_classic(base_size = 14) +
  theme(
    panel.grid.major.y = element_line(color = "grey80"),
    axis.text.x = element_text(size = 9),
    axis.text.y = element_text(size = 9)
  )

# save plot
ggsave("plots/effects_jump_three.jpg", width = 7, height = 4.6, dpi = 300, bg = "white")

# analyse sex-specific effects (without interaction) ---------------------------
s_jump <- lapply(vars_jump, \(x) attr(run_mod(x, data = d), "mm")) |> purrr::list_rbind()
s_jump$name <- rep(vars_jump, each = nrow(s_jump) / length(vars_jump))

# reverse signs for scales on which lower values mean better performance
s_jump[s_jump$name %in% c("0-10m", "10-30m", "0-30m", "30-60m", "0-60m"), c(2,5,6)] <- -1 * s_jump[s_jump$name %in% c("0-10m", "10-30m", "0-30m", "30-60m", "0-60m"), c(2,5,6)]
s_jump[s_jump$name %in% c("0-10m", "10-30m", "0-30m", "30-60m", "0-60m"), c(5,6)] <- s_jump[s_jump$name %in% c("0-10m", "10-30m", "0-30m", "30-60m", "0-60m"), c(5,6)]

# set variable so that they appear in the same order as in the other plot
s_jump$name <- factor(s_jump$name, levels = l_jump$name[order(l_jump$es, decreasing = TRUE)])

pd <- position_dodge(width = 0.3)
ggplot(s_jump[!s_jump$name %in% c("SUb", "SUs"),], aes(x = name, y = APHV_resid.trend, color = sex)) +
  geom_hline(yintercept = 0, color = "grey40") +
  geom_point(size = 2, position = pd) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2, linewidth = 0.7, position = pd) +
  #geom_text(aes(y = upper.CL + 0.013, label = beta_label), size = 6, position = pd) +
  scale_x_discrete(name = "Variable") +
  scale_y_continuous(
    name = "Performance difference per 1-year\n increase in maturity offset", 
    labels = scales::label_percent(),
    breaks = seq(from = -0.30, to = 0.15, by = 0.10)) +
  theme_classic(base_size = 14) +
  theme(
    panel.grid.major.y = element_line(color = "grey80"),
    axis.text.x = element_text(size = 9),
    axis.text.y = element_text(size = 9)
  )
ggsave("plots/effects_jump_three_sex.jpg", width = 7, height = 4.6, dpi = 300, bg = "white")
write.csv(s_jump, "effects_jump_three_sex.csv")

# analyse change in BA effect -------------------------------------------------
m_jump <- lapply(vars_jump, \(x) attr(run_mod(x, data = d), "mt")) |> purrr::list_rbind()
m_jump$name <- rep(vars_jump, each = nrow(m_jump) / length(vars_jump))

# reverse signs for scales on which lower values mean better performance
m_jump[m_jump$name %in% c("0-10m", "10-30m", "0-30m", "30-60m", "0-60m"), c(3,6,7)] <- -1 * m_jump[m_jump$name %in% c("0-10m", "10-30m", "0-30m", "30-60m", "0-60m"), c(3,6,7)]
m_jump[m_jump$name %in% c("0-10m", "10-30m", "0-30m", "30-60m", "0-60m"), c(6,7)] <- m_jump[m_jump$name %in% c("0-10m", "10-30m", "0-30m", "30-60m", "0-60m"), c(6,7)]

# set variable so that they appear in the same order as in the other plot
m_jump$name <- factor(m_jump$name, levels = l_jump$name[order(l_jump$es, decreasing = TRUE)])

ggplot(m_jump[!m_jump$name %in% c("SUb", "SUs"),], aes(x = CAc, y = APHV_resid.trend, color = sex)) +
  geom_hline(yintercept = 0, color = "grey40") +
  geom_point(alpha = 0.7) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2, linewidth = 0.7, alpha = 0.7) +
  scale_x_continuous(
    name = "Chronological Age",
    breaks = CAc_values, labels = 13:19) +
  scale_y_continuous(
    name = "Performance difference per 1-year increase in maturity offset", 
    labels = scales::label_percent(),
    breaks = seq(from = -0.15, to = 0.15, by = 0.15)) +
  facet_wrap(vars(name)) +  
  theme_bw(base_size = 14) +
  theme(
    panel.grid.major.y = element_line(color = "grey80"),
    #panel.grid.minor.y = element_line()
  ) +
  coord_cartesian(ylim = c(-0.20, 0.20))
ggsave("plots/interaction_jump_three.jpg", width = 11, height = 8.7, dpi = 600, bg = "white")
write.csv(m_jump, "interaction_jump_three.csv")

#### PART 3: squad status

d$SGc <- ifelse(d$SG %in% c(1,2), "low", ifelse(d$SG %in% 3:5, "high", NA))
ggplot(d[!is.na(d$SGc),], aes(CA, dAPHV_resid, color = SGc)) +
  geom_point(alpha = 0.5) +
  theme_bw()

ggsave("plots/SGc_resid.png", width = 6, height = 5, dpi = 300, bg = "white")

lm(dAPHV_resid ~ SGc, data = d) |> summary()
mean(d$dAPHV_resid[d$SGc == "low"], na.rm = TRUE)
mean(d$dAPHV_resid[d$SGc == "high"], na.rm = TRUE)

lm(dAPHV_resid ~ SGc, data = d) |> summary()

ggplot(d, aes(CA, FOST, color = SGc)) +
  geom_point(alpha = 0.5) +
  facet_wrap(~sex)

# test for differences in performance between groups
test_SGc_diff <- function(var, data = d) {
  data$out <- data[[var]]
  p <- poly(data$CAc, 2)
  data$CAc_p1 <- p[,1]
  data$CAc_p2 <- p[,2]
  o <- lmerTest::lmer(out ~ CAc_p1 + CAc_p2 + sex + SGc + CAc_p1:sex + CAc_p2:sex + dAPHV_resid + dAPHV_resid:CAc + dAPHV_resid:sex + dAPHV_resid:CAc:sex + (1|ID), data = data)
  so <- summary(o)
  sg <- as.data.frame(t(so$coefficients["SGclow",]))
  sg$var <- var
  sg
}

# whether low performs worse than high
sgcperf_jump <- lapply(vars_jump, test_SGc_diff, data = d) |> purrr::list_rbind()
