#!/usr/bin/env bash

echo "# ---------> NORMALIZATION: Start"
Rscript 1_normalization.R
echo "# <--------- NORMALIZATION: End"

echo "# ---------> DIFFERENTIAL EXPRESSION: Start"
Rscript 2_0_diff_exp.R
python 2_A_diff_exp_plots.py
echo "# <--------- DIFFERENTIAL EXPRESSION: End"

echo "# ---------> GSEA: Start"
Rscript 3_A_0_gsea_piano.R
python 3_A_1_gsea_piano_plots.py
Rscript 3_B_0_gsea_viper.R
python 3_B_1_gsea_viper_plots.py
echo "# <--------- GSEA: End"

