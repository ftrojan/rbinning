# Example usP001 of rbinning

source("src/rbinning.R")
library(tidyverse)
library(caret)

# 1 load data
rdata = read_tsv("data/modeling_dataset.txt")

# 2 identify which columns are predictors
tnames = "T01"
hnames = c("H01","H03","H04","H05","H06","H07","H08","H09","H02")
pnames = names(rdata) %>% setdiff(tnames) %>% setdiff(hnames)

# 3 train/test SPLIT 80/20
# deterministic way which can be reproduced in Python easily
n = nrow(rdata)
ntn = round(0.8*n)
itn = seq_len(ntn)
its = setdiff(seq_len(n), itn)

# 4 AUTOBINNING
nb = 3 # desired number of bins

b0 = autobinning(tbl = rdata, target = "T01", predictors = pnames, itest = its, method = "EqualFrequency", nbins = nb,verbose = TRUE)

# 5 MANUAL CORRECTIONS OF BINNING
b1 = b0 %>% adjust_binning(
  tbl = rdata,
  itest = its,
  ni = "P075",
  target = "T01",
  code = "8",
  missbin = 2
) %>% adjust_binning(
  tbl = rdata,
  itest = its,
  ni = "P094",
  target = "T01",
  code = "0.7",
  missbin = 2
) %>% adjust_binning(
  tbl = rdata,
  itest = its,
  ni = "P102",
  target = "T01",
  code = "0.71, 1.25",
  missbin = 3
) %>% adjust_binning(
  tbl = rdata,
  itest = its,
  ni = "P185",
  target = "T01",
  code = "0, 0.36",
  missbin = 3
) %>% adjust_binning(
  tbl = rdata,
  itest = its,
  ni = "P323",
  target = "T01",
  code = "2.3",
  missbin = 1
) %>% adjust_binning(
  tbl = rdata,
  itest = its,
  ni = "P337",
  target = "T01",
  code = "-0.14, 0.13",
  missbin = 3
) %>% adjust_binning(
  tbl = rdata,
  itest = its,
  ni = "P336",
  target = "T01",
  code = "-0.05",
  missbin = 1
) %>% adjust_binning(
  tbl = rdata,
  itest = its,
  ni = "P070",
  target = "T01",
  code = "1.06",
  missbin = 2
) %>% adjust_binning(
  tbl = rdata,
  itest = its,
  ni = "P074",
  target = "T01",
  code = "0.2",
  missbin = 3
) %>% adjust_binning(
  tbl = rdata,
  itest = its,
  ni = "P152",
  target = "T01",
  code = "1.04",
  missbin = 1
) %>% adjust_binning(
  tbl = rdata,
  itest = its,
  ni = "P202",
  target = "T01",
  code = "-0.01",
  missbin = 2
) %>% adjust_binning(
  tbl = rdata,
  itest = its,
  ni = "P049",
  target = "T01",
  code = "35000, 45000",
  missbin = 1
) %>% adjust_binning(
  tbl = rdata,
  itest = its,
  ni = "P182",
  target = "T01",
  code = "0.01",
  missbin = 1
) %>% adjust_binning(
  tbl = rdata,
  itest = its,
  ni = "P339",
  target = "T01",
  code = "-0.14, 0.11",
  missbin = 1
) %>% adjust_binning(
  tbl = rdata,
  itest = its,
  ni = "P189",
  target = "T01",
  code = "-0.06, 0.02",
  missbin = 1
) %>% adjust_binning(
  tbl = rdata,
  itest = its,
  ni = "P124",
  target = "T01",
  code = "0.54, 1.09",
  missbin = 3
) %>% adjust_binning(
  tbl = rdata,
  itest = its,
  ni = "P212",
  target = "T01",
  code = "4.07, 11.78",
  missbin = 3
) %>% adjust_binning(
  tbl = rdata,
  itest = its,
  ni = "P215",
  target = "T01",
  code = "0, 1.22",
  missbin = 4
)

# 6 WOE TRANSFORMATION 
wdata = tbl_woe(rdata, b1, keep = c(tnames,hnames), verbose = TRUE)
write_tsv(wdata, "out/woe_dataset.txt")

# 7 binning plots
binplots = b1 %>% arrange(desc(gini_train))
nbp = 10
for (i in 1:nbp){
  ggplot_binning(
    binplots[i,],
    "out"
  )
  cat(sprintf("%d/%d\n", i, nbp))
}
