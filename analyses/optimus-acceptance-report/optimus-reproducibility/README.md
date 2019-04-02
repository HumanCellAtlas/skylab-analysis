# Reproducibility

This section tests the reproducibility of the pipeline. Using three identical runs of the t_4k data (pan T 4k dataset from 10X). We compare the output matrices and the emptyDrops output tables.

The input script is setup with the expectation that the data folder contains three directories named 'run1', 'run2', and 'run3'. These must contain two folder each, named 'MergeCountFiles' and 'RunEmptyDrops' that contain the output of the respective output folders from three different pipeline runs.

These have been removed here as they are too large. The respective run identifies on the mint cromwell deployment are:

91a04c2d-a66d-4e5a-8cf6-d6c48f5363f3
ba452651-5641-4b0e-a7d5-f828126c1452
5a2d76d8-151d-4923-8dc7-03b1bf001e85
