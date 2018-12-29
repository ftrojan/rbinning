# rbinning
Binning as feature engineering technique for better machine learning models
You want to do four different things around binning: autobinning, 
manual adjustments, calculate WoE and plot binning graphs.
The following example shows you how to use `rbinning` in your code.

```
source('rbinning.R') # put the right path to the file
```

## Autobinning

```
binning_dataframe = autobinning(
  tbl = input_dataframe, 
  target = "T01", 
  predictors = pnames, 
  itest = its, 
  method = "EqualFrequency", 
  nbins = nb,
  verbose = TRUE
)
```

## Manual adjustments to the binning

```
adjusted_binning_dataframe = binning_dataframe %>% adjust_binning(
  tbl = input_dataframe,
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
)
```

## Score new data and calculate WoE

```
binned = tbl_woe(
  another_input_dataframe, 
  adjusted_binning_dataframe, 
  keep = c(tnames,hnames), 
  verbose = TRUE
)
```

## Plot binning graphs

```
binplots = adjusted_binning_dataframe %>% arrange(desc(gini_train))
number_of_top = 10
for (i in 1:number_of_top){
  binning_row = adjusted_binning_dataframe[i,]
  ggplot_binning(binning_row, output_directory)
  cat(sprintf("%d/%d\n", i, number_of_top))
}
```
