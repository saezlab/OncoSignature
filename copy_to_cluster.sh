#!/usr/bin/env bash

mkdir out


# 2+2_AML_patients
mkdir out/2+2_AML_patients
mkdir out/2+2_AML_patients/data

echo "Entering 2+2_AML_patients directory"
cd 2+2_AML_patients
cp -r data/msigdb_path_to_gsymbol.csv ../out/2+2_AML_patients/data/msigdb_path_to_gsymbol.csv
cp -r data/raw_data.txt ../out/2+2_AML_patients/data/raw_data.txt
cp -r 1_normalization.R ../out/2+2_AML_patients/1_normalization.R
cp -r 2_0_diff_exp.R ../out/2+2_AML_patients/2_0_diff_exp.R
cp -r 2_A_diff_exp_plots.py ../out/2+2_AML_patients/2_A_diff_exp_plots.py
cp -r 3_A_0_gsea_piano.R ../out/2+2_AML_patients/3_A_0_gsea_piano.R
cp -r 3_A_1_gsea_piano_plots.py ../out/2+2_AML_patients/3_A_1_gsea_piano_plots.py
cp -r 3_B_0_gsea_viper.R ../out/2+2_AML_patients/3_B_0_gsea_viper.R
cp -r 3_B_1_gsea_viper_plots.py ../out/2+2_AML_patients/3_B_1_gsea_viper_plots.py
cp -r run_all.sh ../out/2+2_AML_patients/run_all.sh
echo "Leaving 2+2_AML_patients directory"
cd ..

# 4x_AML_cell_lines
mkdir out/4x_AML_cell_lines
mkdir out/4x_AML_cell_lines/data

echo "Entering 4x_AML_cell_lines directory"
cd 4x_AML_cell_lines
cp -r data/msigdb_path_to_gsymbol.csv ../out/4x_AML_cell_lines/data/msigdb_path_to_gsymbol.csv
cp -r data/raw_data.txt ../out/4x_AML_cell_lines/data/raw_data.txt
cp -r 1_normalization.R ../out/4x_AML_cell_lines/1_normalization.R
cp -r 2_0_diff_exp.R ../out/4x_AML_cell_lines/2_0_diff_exp.R
cp -r 2_A_diff_exp_plots.py ../out/4x_AML_cell_lines/2_A_diff_exp_plots.py
cp -r 3_A_0_gsea_piano.R ../out/4x_AML_cell_lines/3_A_0_gsea_piano.R
cp -r 3_A_1_gsea_piano_plots.py ../out/4x_AML_cell_lines/3_A_1_gsea_piano_plots.py
cp -r 3_B_0_gsea_viper.R ../out/4x_AML_cell_lines/3_B_0_gsea_viper.R
cp -r 3_B_1_gsea_viper_plots.py ../out/4x_AML_cell_lines/3_B_1_gsea_viper_plots.py
cp -r run_all.sh ../out/4x_AML_cell_lines/run_all.sh
echo "Leaving 4x_AML_cell_lines directory"
cd ..

# 44_AML_ex_vivo
mkdir out/44_AML_ex_vivo
mkdir out/44_AML_ex_vivo/data

echo "Entering 44_AML_ex_vivo directory"
cd 44_AML_ex_vivo
cp -r data/msigdb_path_to_gsymbol.csv ../out/44_AML_ex_vivo/data/msigdb_path_to_gsymbol.csv
cp -r "data/Phospho (STY)Sites-for-normalization.txt" "../out/44_AML_ex_vivo/data/Phospho (STY)Sites-for-normalization.txt"
cp -r 1_normalization.R ../out/44_AML_ex_vivo/1_normalization.R
cp -r 2_0_diff_exp.R ../out/44_AML_ex_vivo/2_0_diff_exp.R
cp -r 2_A_diff_exp_plots.py ../out/44_AML_ex_vivo/2_A_diff_exp_plots.py
cp -r 3_A_0_gsea_piano.R ../out/44_AML_ex_vivo/3_A_0_gsea_piano.R
cp -r 3_A_1_gsea_piano_plots.py ../out/44_AML_ex_vivo/3_A_1_gsea_piano_plots.py
cp -r 3_B_0_gsea_viper.R ../out/44_AML_ex_vivo/3_B_0_gsea_viper.R
cp -r 3_B_1_gsea_viper_plots.py ../out/44_AML_ex_vivo/3_B_1_gsea_viper_plots.py
cp -r run_all.sh ../out/44_AML_ex_vivo/run_all.sh
echo "Leaving 44_AML_ex_vivo directory"
cd ..

scp -r out/ np176757@cluster.rz.rwth-aachen.de:OncoSignature/
rm -r out/

echo "===== DONE! ====="
