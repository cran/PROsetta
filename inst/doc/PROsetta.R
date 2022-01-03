## ---- echo = FALSE, message = FALSE, output = 'hide'--------------------------
library(knitr)
library(kableExtra)
library(PROsetta)
library(dplyr)

## ----loadData, echo = FALSE---------------------------------------------------
d <- data_dep

## ----cfg, results = "hide", message = FALSE, eval = FALSE---------------------
#  fp <- system.file("data-raw", package = "PROsetta")
#  d <- loadData(
#    response  = "dat_DeCESD_v2.csv",
#    itemmap   = "imap_DeCESD.csv",
#    anchor    = "anchor_DeCESD.csv",
#    input_dir = fp)

## ---- eval = FALSE------------------------------------------------------------
#  file.edit(system.file("data-raw", "dat_DeCESD_v2.csv", package = "PROsetta"))

## ---- eval = FALSE------------------------------------------------------------
#  file.edit(system.file("data-raw", "imap_DeCESD.csv", package = "PROsetta"))

## ---- eval = FALSE------------------------------------------------------------
#  file.edit(system.file("data-raw", "anchor_DeCESD.csv", package = "PROsetta"))

## ----freq---------------------------------------------------------------------
freq_table <- runFrequency(d)
head(freq_table)

## ---- fig.align = "center", fig.width = 7, fig.height = 7---------------------
plot(d, scale = "combined", title = "Combined scale")

## ----eval = FALSE-------------------------------------------------------------
#  plot(d, scale = 1, title = "Scale 1") # not run
#  plot(d, scale = 2, title = "Scale 2") # not run

## ----desc---------------------------------------------------------------------
desc_table <- runDescriptive(d)
head(desc_table)

## ----alpha, cache = TRUE------------------------------------------------------
classical_table <- runClassical(d)
summary(classical_table$alpha$combined)

## ----alpha2, cache = TRUE, eval = FALSE---------------------------------------
#  classical_table <- runClassical(d, scalewise = TRUE)
#  classical_table$alpha$combined # alpha values for combined scale
#  classical_table$alpha$`1`      # alpha values for each scale, created when scalewise = TRUE
#  classical_table$alpha$`2`      # alpha values for each scale, created when scalewise = TRUE

## ----omega, cache = TRUE, eval = FALSE----------------------------------------
#  classical_table <- runClassical(d, scalewise = TRUE, omega = TRUE)
#  classical_table$omega$combined # omega values for combined scale
#  classical_table$omega$`1`      # omega values for each scale, created when scalewise = TRUE
#  classical_table$omega$`2`      # omega values for each scale, created when scalewise = TRUE

## ----eval = FALSE-------------------------------------------------------------
#  classical_table <- runClassical(d, scalewise = TRUE, omega = TRUE, nfactors = 5) # not run

## ----cfa, cache = FALSE, results = 'hide'-------------------------------------
out_cfa <- runCFA(d, scalewise = TRUE)

## ---- eval = FALSE------------------------------------------------------------
#  out_cfa <- runCFA(d, scalewise = TRUE, std.lv = TRUE) # not run

## -----------------------------------------------------------------------------
out_cfa$combined
out_cfa$`1`
out_cfa$`2`

## -----------------------------------------------------------------------------
lavaan::summary(out_cfa$combined, fit.measures = TRUE, standardized = TRUE, estimates = FALSE)

## ---- eval = FALSE------------------------------------------------------------
#  lavaan::summary(out_cfa$`1`, fit.measures = TRUE, standardized = TRUE, estimates = FALSE) # not run
#  lavaan::summary(out_cfa$`2`, fit.measures = TRUE, standardized = TRUE, estimates = FALSE) # not run

## ----calib, cache = TRUE, results = 'hide', error = TRUE, message = FALSE-----
out_calib <- runCalibration(d, technical = list(NCYCLES = 1000))

## ----calib2, cache = TRUE, results = 'hide', error = TRUE, message = FALSE----
out_calib <- runCalibration(d, technical = list(NCYCLES = 10))

## ----calib3, cache = TRUE, results = 'hide', error = TRUE, message = FALSE----
out_calib <- runCalibration(d, fixedpar = TRUE)

## ----mirt_coef, cache = TRUE--------------------------------------------------
mirt::coef(out_calib, IRTpars = TRUE, simplify = TRUE)

## ----mirt_plot, cache = TRUE, eval = FALSE, results = 'hide'------------------
#  mirt::itemfit(out_calib, empirical.plot = 1)
#  mirt::itemplot(out_calib, item = 1, type = "info")
#  mirt::itemfit(out_calib, "S_X2", na.rm = TRUE)

## ---- fig.align = 'center', fig.width = 7, fig.height = 7---------------------
plotInfo(
  out_calib, d,
  scale_label = c("PROMIS Depression", "CES-D", "Combined"),
  color = c("blue", "red", "black"),
  lty = c(1, 2, 3))

## ----fixedpar, cache = TRUE, results = 'hide', message = FALSE----------------
out_link_fixedpar <- runLinking(d, method = "FIXEDPAR")

## -----------------------------------------------------------------------------
out_link_fixedpar$ipar_linked

## ----sl, cache = TRUE, results = 'hide', error = TRUE, message = FALSE--------
out_link_sl <- runLinking(d, method = "SL", technical = list(NCYCLES = 1000))
out_link_sl

## -----------------------------------------------------------------------------
out_link_sl$ipar_linked

## -----------------------------------------------------------------------------
out_link_sl$constants

## -----------------------------------------------------------------------------
rsss_fixedpar <- runRSSS(d, out_link_fixedpar)
rsss_sl       <- runRSSS(d, out_link_sl)
round(rsss_fixedpar$`2`, 3)

## ----eqp_raw, cache = TRUE, results = 'hide', message = FALSE-----------------
out_equate <- runEquateObserved(
  d, scale_from = 2, scale_to = 1,
  eq_type = "equipercentile", smooth = "loglinear")

## -----------------------------------------------------------------------------
out_equate$concordance

## ----eqp_dir, cache = TRUE, results = 'hide', message = FALSE-----------------
out_equate_tscore <- runEquateObserved(
  d, scale_from = 2, scale_to = 1,
  type_to = "tscore", rsss = rsss_fixedpar,
  eq_type = "equipercentile", smooth = "loglinear")

## -----------------------------------------------------------------------------
out_equate_tscore$concordance

## ---- fig.align = 'center', fig.width = 7, fig.height = 7---------------------
plot(
  rsss_fixedpar$`2`$raw_2,
  rsss_fixedpar$`2`$tscore,
  xlab = "CES-D Summed Score",
  ylab = "PROMIS Depression T-score",
  type = "l", col = "blue")
lines(
  out_equate_tscore$concordance$raw_2,
  out_equate_tscore$concordance$tscore_1,
  lty = 2, col = "red")
grid()
legend(
  "topleft",
  c("Fixed-Parameter Calibration", "Equipercentile Linking"),
  lty = 1:2, col = c("blue", "red"), bg = "white"
)

## -----------------------------------------------------------------------------
scores <- getScaleSum(d, 2)
head(scores)

## ---- message = FALSE---------------------------------------------------------
eap_promis <- getTheta(d, out_link_fixedpar$ipar_linked, scale = 1)$theta
head(eap_promis)

## -----------------------------------------------------------------------------
t_promis <- data.frame(
  prosettaid = eap_promis$prosettaid,
  t_promis = round(eap_promis$theta_eap * 10 + 50, 1)
)
head(t_promis)

## -----------------------------------------------------------------------------
scores <- merge(scores, t_promis, by = "prosettaid")
head(scores)

## ---- message = FALSE---------------------------------------------------------
eap_cesd <- getTheta(d, out_link_fixedpar$ipar_linked, scale = 2)$theta
t_cesd_pattern <- data.frame(
  prosettaid = eap_cesd$prosettaid,
  t_cesd_pattern = round(eap_cesd$theta_eap * 10 + 50, 1)
)
scores <- merge(scores, t_cesd_pattern, by = "prosettaid")
head(scores)

## ---- message = FALSE---------------------------------------------------------
rsss_eap <- data.frame(
  raw_2 = rsss_fixedpar$`2`$raw_2,
  t_cesd_rsss = round(rsss_fixedpar$`2`$tscore, 1)
)
scores <- merge(scores, rsss_eap, by = "raw_2")
head(scores)

## ---- message = FALSE---------------------------------------------------------
rsss_eqp <- data.frame(
  raw_2 = out_equate_tscore$concordance$raw_2,
  t_cesd_eqp = round(out_equate_tscore$concordance$tscore_1, 1)
)
scores <- merge(scores, rsss_eqp, by = "raw_2")
head(scores)

## -----------------------------------------------------------------------------
# Reference score: IRT pattern scoring of Scale 1
c_pattern <- compareScores(
  scores$t_promis, scores$t_cesd_pattern) ## IRT response pattern EAP to T-score
c_rsss <- compareScores(
  scores$t_promis, scores$t_cesd_rsss)    ## IRT summed score EAP to T-score
c_eqp <- compareScores(
  scores$t_promis, scores$t_cesd_eqp)     ## Equipercentile summed score to T-score

stats           <- rbind(c_pattern, c_rsss, c_eqp)
rownames(stats) <- c("IRT Pattern", "IRT RSSS", "Equipercentile")
stats

