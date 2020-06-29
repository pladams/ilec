# pa.glmtree.R
Contains miscellaneous code fragments from ILEC mortality analysis.

This code modifies the glmtree code path from the function in the R partykit package to append a weight applied to the inner GLM fitting code. This is necessary in the presence of offsets in mortality modeling. Offsets are often used as an exposure or expected claims element. In this case, using a weight for the GLM part is inappropriate. However, failing to exclude weights in the GLM confuses the parameter fluctuation tests which assume then that every row in the data is equally weighted.
