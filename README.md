# LIV MSM simulation workflow

This guide runs the complete simulation workflow from the command line, with every long-running task kept in the background. Run every command from the project root:

```bash
cd /Users/main/Coding/LIV_MSM-main
mkdir -p logs
```

The parameter files in `Script/par/` define the DGPs and are loaded automatically by the worker scripts. Do not run them as standalone simulation jobs.

## 1. Prepare the semi-synthetic JTPA DGP

The JTPA continuous-outcome DGP requires two fixed quantile maps. Generate the `L1` map first, then the terminal-earnings map. The following job also runs the DGP diagnostic after both maps have been created.

```bash
nohup bash -lc '
  Rscript Script/par/calibrate_jtpa_l1_quantile_map.R &&
  Rscript Script/par/calibrate_jtpa_earnings_quantile_map.R &&
  Rscript Script/par/check_jtpa_dgps.R
' > logs/01_prepare_jtpa.log 2>&1 &
```

Monitor the job with:

```bash
tail -f logs/01_prepare_jtpa.log
```

Wait until it finishes successfully before starting the JTPA simulations. The two calibration scripts use 250,000 simulated observations per DGP by default, so this step can take time. Running `source("...")` in R only loads the calibration functions; it does not generate the CSV files. Use the `Rscript` commands above, or call the corresponding `calibrate_*()` functions explicitly from R.

## 2. Run the ordinary simulation study with XGBoost

Each worker performs 500 replications at sample sizes 2,000 and 5,000, writing results to `data/xgb2000/` and `data/xgb5000/`.

```bash
for estimand in Complier Nudge Overall; do
  for policy in 00 01 10 11; do
    job="xgb_${estimand}_${policy}"
    nohup Rscript "Script/xgb/${estimand}${policy}.R" \
      > "logs/${job}.log" 2>&1 &
  done
done
```

## 3. Run the semi-synthetic JTPA study with XGBoost

Each worker performs 500 replications at sample size 10,000. Results are written to `data/jtpa10000/`.

```bash
for estimand in Complier Overall; do
  for policy in 00 01 0p 0n 10 11 1p 1n; do
    job="jtpa_xgb_${estimand}_${policy}"
    nohup Rscript "Script/xgb/${estimand}_jtpa${policy}.R" \
      > "logs/${job}.log" 2>&1 &
  done
done
```

## 4. Run the ordinary simulation study with MGCV

Each worker performs 500 replications at sample sizes 2,000 and 5,000. Results are written to `data/mgcv2000/` and `data/mgcv5000/`.

```bash
for estimand in Complier Nudge Overall; do
  for policy in 00 01 10 11; do
    job="mgcv_${estimand}_${policy}"
    nohup Rscript "Script/mgcv/${estimand}${policy}.R" \
      > "logs/${job}.log" 2>&1 &
  done
done
```

The three blocks above launch 40 long-running R jobs. On a machine with limited memory or CPU, run one block at a time, or edit the loops to submit only a subset of policies before launching the next subset.

## 5. Monitor the background jobs

```bash
jobs -l
rg -n "Error|Execution halted" logs
tail -f logs/jtpa_xgb_Overall_00.log
```

Do not start the summary scripts until all simulation workers have completed and the error scan is empty.

## 6. Create the tables and figures

The ordinary-simulation summary creates four PDF figures and four LaTeX tables. The JTPA summary creates one PDF figure and one LaTeX table.

```bash
nohup env SUMMARY_TRUTH_ACCURACY=1000000 \
  Rscript Script/summary_ggplot.R \
  > logs/summary_ordinary.log 2>&1 &

nohup env JTPA_TRUTH_ACCURACY=500000 \
  Rscript Script/summary_ggplot_jtpa.R \
  > logs/summary_jtpa.log 2>&1 &
```

The ordinary-simulation outputs are written to `outputs/simulation/`:

- `plot_xgb_complier_overall_ggplot.pdf` and
  `summary_xgb_complier_overall_table_integrated.tex`;
- `plot_xgb_nudge_ggplot.pdf` and `summary_xgb_nudge_table.tex`;
- `plot_mgcv_complier_overall_ggplot.pdf` and
  `summary_mgcv_complier_overall_table_integrated.tex`;
- `plot_mgcv_nudge_ggplot.pdf` and `summary_mgcv_nudge_table.tex`;

The JTPA outputs are written to `outputs/jtpa/`:

- `jtpa_12_dgp_ggplot.pdf` and `jtpa_12_dgp_table.tex`.

For a faster diagnostic summary only, lower `SUMMARY_TRUTH_ACCURACY` or `JTPA_TRUTH_ACCURACY`. Do not use those smaller values for the final reported results.









