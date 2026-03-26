# Supplemental Figure S3: comparison with differential-correlation methods

Count simulation comparing MOSAIC, MOSAIC without projection, DGCA and memento.

## Requirements

- R packages: Seurat (>= 5), RSpectra, Matrix, DGCA, PRROC, foreach, doParallel, digest, jsonlite, ggplot2.
- Python with numpy, pandas, scipy, anndata and memento from commit 224c5c9
  (`pip install "git+https://github.com/yelabucsf/scrna-parameter-estimation.git@224c5c982f952cb1fd0f0de4623852dd435dbc1d"`);
  set `S3_PYTHON` to its interpreter.

## Run

Run from this folder (paths are relative to the working directory). `s3_frozen.rds` is the configuration used for the figure; Stage 3 overwrites it.

```bash
export S3_PYTHON=/path/to/memento-env/bin/python
Rscript s3_stage3_calibrate.R --workers 70
Rscript s3_stage3c_predict.R --workers 70
Rscript s3_stage4_run.R --mode controls --workers 70
Rscript s3_stage4_run.R --mode formal --workers 70
Rscript s3_plot.R
Rscript s3_stage4b_memento_coef.R --workers 70   # optional
```

Outputs: `stage4/figS3.{png,pdf}` (Supplemental Figure S3), `stage4/formal_auprc.csv` (AUPRC per cohort and method), `stage3/` (calibration logs).
