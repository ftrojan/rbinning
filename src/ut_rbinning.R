# Example usage of binning toolbox in R
# T-Mobile HR Analytics
# Deloitte 2018

source("rbinning.R")
library(tidyverse)
library(caret)

# Loading data
rdata = read_tsv("modeling_dataset_model_B.txt")

tnames = "target"

# all predictors
hnames = c("base_id","PIN","start_date","end_date","obs_date","predictor_date","target_start","target_end","tgf")
pnames = names(rdata) %>% setdiff(tnames) %>% setdiff(hnames)

# train/test SPLIT 80/20
# deterministic way which can be reproduced in Python easily
n = nrow(rdata)
ntn = round(0.8*n)
itn = seq_len(ntn)
its = setdiff(seq_len(n), itn)

#AUTOBINNING
nb = 3 # desired number of bins

b0 = autobinning(tbl = rdata, target = "target", predictors = pnames, itest = its, method = "EqualFrequency", nbins = nb,verbose = TRUE)

#MANUAL CORRECTIONS OF BINNING
b1 = b0 %>% adjust_binning(
  tbl = rdata,
  itest = its,
  ni = "month",
  target = "target",
  code = "8",
  missbin = 2
) %>% adjust_binning(
  tbl = rdata,
  itest = its,
  ni = "p_avg_fix",
  target = "target",
  code = "0.7",
  missbin = 2
) %>% adjust_binning(
  tbl = rdata,
  itest = its,
  ni = "p_bonus_fulfill_3M",
  target = "target",
  code = "0.71, 1.25",
  missbin = 3
) %>% adjust_binning(
  tbl = rdata,
  itest = its,
  ni = "worplace_churn_dif_6M",
  target = "target",
  code = "0, 0.36",
  missbin = 3
) %>% adjust_binning(
  tbl = rdata,
  itest = its,
  ni = "sal_bonus_m_dif_6M",
  target = "target",
  code = "2.3",
  missbin = 1
) %>% adjust_binning(
  tbl = rdata,
  itest = its,
  ni = "base_bon_t_3M",
  target = "target",
  code = "-0.14, 0.13",
  missbin = 3
) %>% adjust_binning(
  tbl = rdata,
  itest = its,
  ni = "base_bon_t_6M",
  target = "target",
  code = "-0.05",
  missbin = 1
) %>% adjust_binning(
  tbl = rdata,
  itest = its,
  ni = "bonus_fulfill_6M",
  target = "target",
  code = "1.06",
  missbin = 2
) %>% adjust_binning(
  tbl = rdata,
  itest = its,
  ni = "churn_rate",
  target = "target",
  code = "0.2",
  missbin = 3
) %>% adjust_binning(
  tbl = rdata,
  itest = its,
  ni = "p_sal_base_6m",
  target = "target",
  code = "1.04",
  missbin = 1
) %>% adjust_binning(
  tbl = rdata,
  itest = its,
  ni = "peer_size_dif_lag",
  target = "target",
  code = "-0.01",
  missbin = 2
) %>% adjust_binning(
  tbl = rdata,
  itest = its,
  ni = "sal_gross_6m",
  target = "target",
  code = "35000, 45000",
  missbin = 1
) %>% adjust_binning(
  tbl = rdata,
  itest = its,
  ni = "worplace_churn_v_6M",
  target = "target",
  code = "0.01",
  missbin = 1
) %>% adjust_binning(
  tbl = rdata,
  itest = its,
  ni = "base_bon_a_dif_6M",
  target = "target",
  code = "-0.14, 0.11",
  missbin = 1
) %>% adjust_binning(
  tbl = rdata,
  itest = its,
  ni = "fulfill_a_dif_6M",
  target = "target",
  code = "-0.06, 0.02",
  missbin = 1
) %>% adjust_binning(
  tbl = rdata,
  itest = its,
  ni = "p_abs_min_h",
  target = "target",
  code = "0.54, 1.09",
  missbin = 3
) %>% adjust_binning(
  tbl = rdata,
  itest = its,
  ni = "abs_h_fte_v_6M",
  target = "target",
  code = "4.07, 11.78",
  missbin = 3
) %>% adjust_binning(
  tbl = rdata,
  itest = its,
  ni = "abs_h_fte_dif_6M",
  target = "target",
  code = "0, 1.22",
  missbin = 4
)
# vname = "abs_h_fte_dif_6M"
# plot_binning_train(data.bin.final[data.bin.final$varname==vname,])
# plot_binning_test(data.bin.final[data.bin.final$varname==vname,])
# data.bin.final[data.bin.final$varname==vname,] %>% getElement("code")

# WOE TRANSFORMATION 
binned = tbl_woe(rdata, b1, keep = c(tnames,hnames), verbose = TRUE)
write_tsv(binned, "woe_dataset_R.txt")

# binning plots
binplots = b1 %>% arrange(desc(gini_train))
nbp = 10
for (i in 1:nbp){
  
  ggplot_binning(
    binplots[i,],
    "plots_R"
  )
  cat(sprintf("%d/%d\n", i, nbp))

}
