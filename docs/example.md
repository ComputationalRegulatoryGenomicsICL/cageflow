## Test data

The pipeline can be run on test single-end and paired-end data. The directory containing test CAGE libraries and all necessary reference genome files can be downloaded from the Sequence Read Archives. It should be placed next to this repository and uncompressed.

Single-end test reads were randomly sampled from yeast CAGE libraries "Ana" (replicate 1, ERR2495148, anaerobic conditions) and "Eth" (replicate 1, ERR2495150, ethanol limitation) generated and analyzed in ([Börlin et al., 2019](https://academic.oup.com/femsyr/article/19/2/foy128/5257840)).
The test data include two samples split into two lanes each.
One lane contains approximately 300,000 reads; therefore, the whole dataset contains around 1.2 mln reads.
See `generate_se_test_data.sh` for details on how the test dataset was generated.

To run the pipeline on these data, use the following command from outside the pipeline repository:

```
nextflow run customcageq/main.nf \
    -params-file customcageq/testdata/params_yeast_borlin_test.yaml \
    -profile singularity \
    -w work_se_test
```

adapting the suggested values for the options `-profile` (Nextflow profile name) and `-w` (Nextflow work directory) according to your system's setup.

Paired-end test reads were randomly sampled from zebrafish CAGE libraries "4-5 somites" (SRR10215487) and "prim-5" (SRR10215486) generated and analyzed in ([Nepal et al., 2020](https://doi.org/10.1038/s41467-019-13687-0)). The test data include two samples split into two lanes each. One lane contains approximately 1 mln reads; therefore, the whole dataset contains around 4 mln reads. See `generate_pe_test_data.sh` for details on how the test dataset was generated.

To run the pipeline on these data, use the following command from outside the pipeline repository:

```
nextflow run customcageq/main.nf \
    -params-file customcageq/testdata/params_danio_nepal_test.yaml \
    -profile singularity \
    -w work_pe_test
```

adapting the suggested values for the options `-profile` (Nextflow profile name) and `-w` (Nextflow work directory) according to your system's setup.
