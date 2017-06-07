

### Limb Enhancer Genie - Source code for model tuning and combining

#### Model training / parameter tuning

* In your R-session change your working directory to the "model_tuning" folder
* Install the tunetest package

```
install.packages("./tunetest", repos=NULL, type="source")
```

* tune_models.R contains the code used to fit chromatin and sequence models
* pre-cumputed output files ('chrom.models.Rdata' & 'seq.models.Rdata') can be found in the "model_combining" folder


#### Model combining

* In your R-session change your working directory to the "model_combining" folder
* Make sure you have the tunetest package installed (see above)
* combine_models.R contains the code used to combine chromatin and sequence models on the 'leave one out' test sets
