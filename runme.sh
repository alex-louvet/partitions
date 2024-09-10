g++ -fopenmp -Ofast main.cpp

mkdir run_artifacts

for i in {1..10}
do
	echo "====== RUN " $i " ======"
	./a.out -a '1 1 4 1 10 1' -n 8192 -d 2 -t 16 -f grid -k 150 -s -e
	./a.out -a '1 1 4 1 10 1' -n 8192 -d 2 -t 32 -f grid -k 150 -s -e
	./a.out -a '1 1 4 1 10 1' -n 8192 -d 2 -t 64 -f grid -k 150 -s -e
	./a.out -a '1 1 4 1 10 1' -n 8192 -d 2 -t 128 -f grid -k 150 -s -e
	./a.out -a '1 1 4 1 10 1' -n 8192 -d 2 -t 256 -f grid -k 150 -s -e
	./a.out -a '1 1 4 1 10 1' -n 8192 -d 2 -t 512 -f grid -k 150 -s -e
	./a.out -a '1 1 4 1 10 1' -n 8192 -d 2 -t 1024 -f grid -k 150 -s -e
	./a.out -a '1 1 4 1 10 1' -n 8192 -d 2 -t 2048 -f grid -k 150 -s -e
	./a.out -a '1 1 4 1 10 1' -n 8192 -d 3 -t 16 -f grid -k 150 -s -e
	./a.out -a '1 1 4 1 10 1' -n 8192 -d 3 -t 32 -f grid -k 150 -s -e
	./a.out -a '1 1 4 1 10 1' -n 8192 -d 3 -t 64 -f grid -k 150 -s -e
	./a.out -a '1 1 4 1 10 1' -n 8192 -d 3 -t 128 -f grid -k 150 -s -e
	./a.out -a '1 1 4 1 10 1' -n 8192 -d 3 -t 256 -f grid -k 150 -s -e
	./a.out -a '1 1 4 1 10 1' -n 8192 -d 3 -t 512 -f grid -k 150 -s -e
	./a.out -a '1 1 4 1 10 1' -n 8192 -d 4 -t 16 -f grid -k 150 -s -e
	./a.out -a '1 1 4 1 10 1' -n 8192 -d 4 -t 32 -f grid -k 150 -s -e
	./a.out -a '1 1 4 1 10 1' -n 8192 -d 4 -t 64 -f grid -k 150 -s -e
	./a.out -a '1 1 4 1 10 1' -n 8192 -d 4 -t 128 -f grid -k 150 -s -e
	./a.out -a '1 1 4 1 10 1' -n 8192 -d 4 -t 256 -f grid -k 150 -s -e
	./a.out -a '1 1 4 1 10 1' -n 8192 -d 4 -t 512 -f grid -k 150 -s -e
	./a.out -a '1 1 4 1 10 1' -n 8192 -d 5 -t 16 -f grid -k 150 -s -e
	./a.out -a '1 1 4 1 10 1' -n 8192 -d 5 -t 32 -f grid -k 150 -s -e
	./a.out -a '1 1 4 1 10 1' -n 8192 -d 5 -t 64 -f grid -k 150 -s -e
	./a.out -a '1 1 4 1 10 1' -n 8192 -d 5 -t 128 -f grid -k 150 -s -e
	./a.out -a '1 1 4 1 10 1' -n 8192 -d 5 -t 256 -f grid -k 150 -s -e
	./a.out -a '1 1 4 1 10 1' -n 8192 -d 5 -t 512 -f grid -k 150 -s -e
	./a.out -a '1 1 4 1 10 1' -n 8192 -d 10 -t 16 -f grid -k 150 -s -e
	./a.out -a '1 1 4 1 10 1' -n 8192 -d 10 -t 32 -f grid -k 150 -s -e
	./a.out -a '1 1 4 1 10 1' -n 8192 -d 10 -t 64 -f grid -k 150 -s -e
	./a.out -a '1 1 4 1 10 1' -n 8192 -d 10 -t 128 -f grid -k 150 -s -e
	./a.out -a '1 1 4 1 10 1' -n 8192 -d 10 -t 256 -f grid -k 150 -s -e
	./a.out -a '1 1 4 1 10 1' -n 8192 -d 10 -t 512 -f grid -k 150 -s -e
	./a.out -a '1 1 4 1 10 1' -n 2048 -d 2 -t 512 -f grid -k 150 -s -e
	./a.out -a '1 1 4 1 10 1' -n 4096 -d 2 -t 512 -f grid -k 150 -s -e
	./a.out -a '1 1 4 1 10 1' -n 16384 -d 2 -t 512 -f grid -k 150 -s -e
	./a.out -a '1 1 10 1' -n 4000 -d 1 -p 2 -t 32 -f power_law -k 150 -s -e
	./a.out -a '1 1 10 1' -n 4000 -d 1 -p 2.5 -t 32 -f power_law -k 150 -s -e
	./a.out -a '1 1 10 1' -n 4000 -d 1 -p 3 -t 32 -f power_law -k 150 -s -e
	./a.out -a '1 1 10 1' -n 4000 -d 1 -p 2 -t 128 -f power_law -k 150 -s -e
	./a.out -a '1 1 10 1' -n 4000 -d 1 -p 2.5 -t 128 -f power_law -k 150 -s -e
	./a.out -a '1 1 10 1' -n 4000 -d 1 -p 3 -t 128 -f power_law -k 150 -s -e
	./a.out -a '1 1 10 1' -n 4000 -d 1 -p 2 -t 512 -f power_law -k 150 -s -e
	./a.out -a '1 1 10 1' -n 4000 -d 1 -p 2.5 -t 512 -f power_law -k 150 -s -e
	./a.out -a '1 1 10 1' -n 4000 -d 1 -p 3 -t 512 -f power_law -k 150 -s -e
	./a.out -a '1 1 10 1' -n 10000 -d 1 -p 2 -t 32 -f power_law -k 150 -s -e
	./a.out -a '1 1 10 1' -n 10000 -d 1 -p 2.5 -t 32 -f power_law -k 150 -s -e
	./a.out -a '1 1 10 1' -n 10000 -d 1 -p 3 -t 32 -f power_law -k 150 -s -e
	./a.out -a '1 1 10 1' -n 10000 -d 1 -p 2 -t 128 -f power_law -k 150 -s -e
	./a.out -a '1 1 10 1' -n 10000 -d 1 -p 2.5 -t 128 -f power_law -k 150 -s -e
	./a.out -a '1 1 10 1' -n 10000 -d 1 -p 3 -t 128 -f power_law -k 150 -s -e
	./a.out -a '1 1 10 1' -n 10000 -d 1 -p 2 -t 512 -f power_law -k 150 -s -e
	./a.out -a '1 1 10 1' -n 10000 -d 1 -p 2.5 -t 512 -f power_law -k 150 -s -e
	./a.out -a '1 1 10 1' -n 10000 -d 1 -p 3 -t 512 -f power_law -k 150 -s -e
done

mkdir latex_report
mkdir latex_report/img
touch latex_report/report.tex
echo "\documentclass[10pt]{article}
\usepackage{graphicx}
\usepackage{amsmath}
\usepackage{amsfonts}
\usepackage{multirow}
\usepackage{booktabs}
\usepackage[section]{placeins}


\\newcommand{\\F}{\\mathcal{F}}

\begin{document}
\title{A Greedy Algorithm for Low-Crossing Partitions for General Set Systems}
\date{\today}

\maketitle" >> latex_report/report.tex

python3 visualization/draw_activity.py 8192 2 16 latex_report/img/activity_16.png
python3 visualization/draw_activity.py 8192 2 32 latex_report/img/activity_32.png
python3 visualization/draw_activity.py 8192 2 64 latex_report/img/activity_64.png

python3 visualization/draw_comparison.py n grid latex_report/img/crossing_n.png
python3 visualization/draw_comparison.py t grid latex_report/img/crossing_t.png
python3 visualization/draw_comparison.py d grid latex_report/img/crossing_d.png

python3 visualization/draw_runtime.py n grid latex_report/img/runtime_n.png
python3 visualization/draw_runtime.py t grid latex_report/img/runtime_t.png
python3 visualization/draw_runtime.py d grid latex_report/img/runtime_d.png

echo "\\section{Potential function bound violation}

\\begin{figure}[!htb]
\\centering
\\includegraphics[width=0.8\\textwidth]{img/activity_16.png}\\\\
\\includegraphics[width=0.8\\textwidth]{img/activity_32.png}\\\\
\\includegraphics[width=0.8\\textwidth]{img/activity_64.png}\\\\
\\caption{Number of potential function violations.}
\\end{figure}
" >> latex_report/report.tex

echo "

\\begin{table}[!htb]
  \\centering
      \\begin{tabular}{@{}ccc@{}}


          \\multicolumn{1}{c|}{Input}
          & \\multicolumn{1}{c|}{GreedyPotential}
          & \\multicolumn{1}{c}{MinWeight}
          \\\\
          \\multicolumn{1}{c|}{\$n,d,t\$}
          & \\multicolumn{1}{c|}{\# violations}
          & \\multicolumn{1}{c}{\# violations}
          \\\\
          \\midrule" >> latex_report/report.tex

python3 visualization/draw_table.py 8192 2 128 grid rate >> latex_report/report.tex
python3 visualization/draw_table.py 8192 2 256 grid rate >> latex_report/report.tex
python3 visualization/draw_table.py 8192 2 512 grid rate >> latex_report/report.tex
python3 visualization/draw_table.py 8192 2 1024 grid rate >> latex_report/report.tex
python3 visualization/draw_table.py 8192 2 2048 grid rate >> latex_report/report.tex
echo "\\midrule" >> latex_report/report.tex
python3 visualization/draw_table.py 8192 2 512 grid rate >> latex_report/report.tex
python3 visualization/draw_table.py 8192 3 512 grid rate >> latex_report/report.tex
python3 visualization/draw_table.py 8192 4 512 grid rate >> latex_report/report.tex
python3 visualization/draw_table.py 8192 5 512 grid rate >> latex_report/report.tex
python3 visualization/draw_table.py 8192 10 512 grid rate >> latex_report/report.tex
echo "\\midrule" >> latex_report/report.tex
python3 visualization/draw_table.py 2048 2 512 grid rate >> latex_report/report.tex
python3 visualization/draw_table.py 4096 2 512 grid rate >> latex_report/report.tex
python3 visualization/draw_table.py 8192 2 512 grid rate >> latex_report/report.tex
python3 visualization/draw_table.py 16384 2 512 grid rate >> latex_report/report.tex
echo "\\end{tabular}
  \\caption{Potential function violations}
  \\vspace{-15pt}
\\end{table}" >> latex_report/report.tex

echo "\\newpage
\\section{Crossing number on the grid set system}

\\begin{figure}[!htb]
\\centering
\\includegraphics[width=0.48\\textwidth]{img/crossing_d.png}
\\includegraphics[width=0.48\\textwidth]{img/runtime_d.png}
\\includegraphics[width=0.48\\textwidth]{img/crossing_n.png}
\\includegraphics[width=0.48\\textwidth]{img/runtime_n.png}
\\includegraphics[width=0.48\\textwidth]{img/crossing_t.png}
\\includegraphics[width=0.48\\textwidth]{img/runtime_t.png}
\\caption{Average crossing numbers and runtimes of the 3 variations of the algorithm depending on the parameters \$n,d,t\$ of the set system.}
\\end{figure}
" >> latex_report/report.tex

echo "
\\begin{table}[!htb]
    \\centering
        \\begin{tabular}{@{}ccccccccccc@{}}


            \\multicolumn{1}{c|}{Input}
            & \\multicolumn{1}{c|}{}
            & \\multicolumn{3}{c|}{MinWeight}
            & \\multicolumn{3}{c|}{GreedyPotential}
            & \\multicolumn{3}{c}{PartitionSetsAtOnce}


            \\\\


            \\multicolumn{1}{c|}{\$n,d,t\$}
            & \\multicolumn{1}{c|}{\$t^{1-1/d}\$}
            & $\\kappa_{\\mathcal{F}}$
            & $\\bar{\\kappa_F}$
            & \\multicolumn{1}{c|}{Runtime}
            & $\\kappa_{\\mathcal{F}}$
            & $\\bar{\\kappa_F}$
            & \\multicolumn{1}{c|}{Runtime}
            & $\\kappa_{\\mathcal{F}}$
            & $\\bar{\\kappa_F}$
            & Runtime


            \\\\
            \\midrule" >> latex_report/report.tex

python3 visualization/draw_table.py 8192 2 128 grid grid >> latex_report/report.tex
python3 visualization/draw_table.py 8192 2 256 grid grid >> latex_report/report.tex
python3 visualization/draw_table.py 8192 2 512 grid grid >> latex_report/report.tex
python3 visualization/draw_table.py 8192 2 1024 grid grid >> latex_report/report.tex
python3 visualization/draw_table.py 8192 2 2048 grid grid >> latex_report/report.tex
echo "\\midrule" >> latex_report/report.tex
python3 visualization/draw_table.py 8192 2 512 grid grid >> latex_report/report.tex
python3 visualization/draw_table.py 8192 3 512 grid grid >> latex_report/report.tex
python3 visualization/draw_table.py 8192 4 512 grid grid >> latex_report/report.tex
python3 visualization/draw_table.py 8192 5 512 grid grid >> latex_report/report.tex
python3 visualization/draw_table.py 8192 10 512 grid grid >> latex_report/report.tex
echo "\\midrule" >> latex_report/report.tex
python3 visualization/draw_table.py 2048 2 512 grid grid >> latex_report/report.tex
python3 visualization/draw_table.py 4096 2 512 grid grid >> latex_report/report.tex
python3 visualization/draw_table.py 8192 2 512 grid grid >> latex_report/report.tex
python3 visualization/draw_table.py 16384 2 512 grid grid >> latex_report/report.tex

echo "\\end{tabular}
    
    \\caption{ : \$\\kappa_\\F, \\bar{\\kappa_\\F}\$ and runtime of our algorithms on the grid set system.}
\\end{table}" >> latex_report/report.tex

echo "
\\newpage
\\section{Crossing number on power-law graph induced set system}
\\begin{table}[!htb]
    \\centering

        \\begin{tabular}{@{}ccccccc@{}}


            \\multicolumn{1}{c|}{Input}
            & \\multicolumn{1}{c|}{\$\\mathbb{E}\$[VC-dim]}
            & \\multicolumn{1}{c|}{}
            & \\multicolumn{2}{c|}{MinWeight}
            & \\multicolumn{2}{c}{PartitionSetsAtOnce}
            
            \\\\
            
            \\multicolumn{1}{c|}{\$n,\beta,t\$}
            & \\multicolumn{1}{c|}{}
            & \\multicolumn{1}{c|}{\$t^{1-1/d}\$}
            & \$\\kappa_\\F\$
            & \\multicolumn{1}{c|}{runtime (s)}
            & \$\\kappa_\\F\$
            & \\multicolumn{1}{c}{runtime (s)}
            \\\\
            \\midrule" >> latex_report/report.tex

python3 visualization/draw_table.py 4000 2 32 power_law pl >> latex_report/report.tex
python3 visualization/draw_table.py 4000 2.5 32 power_law pl >> latex_report/report.tex
python3 visualization/draw_table.py 4000 3 32 power_law pl >> latex_report/report.tex
echo "\\midrule" >> latex_report/report.tex
python3 visualization/draw_table.py 4000 2 128 power_law pl >> latex_report/report.tex
python3 visualization/draw_table.py 4000 2.5 128 power_law pl >> latex_report/report.tex
python3 visualization/draw_table.py 4000 3 128 power_law pl >> latex_report/report.tex
echo "\\midrule" >> latex_report/report.tex
python3 visualization/draw_table.py 4000 2 512 power_law pl >> latex_report/report.tex
python3 visualization/draw_table.py 4000 2.5 512 power_law pl >> latex_report/report.tex
python3 visualization/draw_table.py 4000 3 512 power_law pl >> latex_report/report.tex
echo "\\midrule" >> latex_report/report.tex
echo "\\midrule" >> latex_report/report.tex
python3 visualization/draw_table.py 10000 2 32 power_law pl >> latex_report/report.tex
python3 visualization/draw_table.py 10000 2.5 32 power_law pl >> latex_report/report.tex
python3 visualization/draw_table.py 10000 3 32 power_law pl >> latex_report/report.tex
echo "\\midrule" >> latex_report/report.tex
python3 visualization/draw_table.py 10000 2 128 power_law pl >> latex_report/report.tex
python3 visualization/draw_table.py 10000 2.5 128 power_law pl >> latex_report/report.tex
python3 visualization/draw_table.py 10000 3 128 power_law pl >> latex_report/report.tex
echo "\\midrule" >> latex_report/report.tex
python3 visualization/draw_table.py 10000 2 512 power_law pl >> latex_report/report.tex
python3 visualization/draw_table.py 10000 2.5 512 power_law pl >> latex_report/report.tex
python3 visualization/draw_table.py 10000 3 512 power_law pl >> latex_report/report.tex

echo "\\end{tabular}
    \\caption{ : Crossing number and runtime of our algorithms on the power-law graph neighborhood set system.}
\\end{table}" >> latex_report/report.tex


echo "
\\newpage
\\section{Error factor comparison with random samples}
\\begin{table}[!htb]
    \\centering

        \\begin{tabular}{@{}ccccc@{}}


            \\multicolumn{1}{c|}{Input}
            & \\multicolumn{1}{c|}{Random sample}
            & \\multicolumn{1}{c|}{MinWeight}
            & \\multicolumn{1}{c|}{GreedyPotential}
            & \\multicolumn{1}{c}{PartitionSetsAtOnce}


            \\\\


            \\multicolumn{1}{c|}{\$n,d,t\$}
            & \\multicolumn{1}{c|}{Error factor}
            & \\multicolumn{1}{c|}{Error factor}
            & \\multicolumn{1}{c|}{Error factor}
            & \\multicolumn{1}{c}{Error factor}


            \\\\
            \\midrule" >> latex_report/report.tex

python3 visualization/draw_table.py 8192 2 16 grid approx >> latex_report/report.tex
python3 visualization/draw_table.py 8192 2 32 grid approx >> latex_report/report.tex
python3 visualization/draw_table.py 8192 2 64 grid approx >> latex_report/report.tex
python3 visualization/draw_table.py 8192 2 128 grid approx >> latex_report/report.tex
python3 visualization/draw_table.py 8192 2 256 grid approx >> latex_report/report.tex
echo "\\midrule" >> latex_report/report.tex
python3 visualization/draw_table.py 8192 3 16 grid approx >> latex_report/report.tex
python3 visualization/draw_table.py 8192 3 32 grid approx >> latex_report/report.tex
python3 visualization/draw_table.py 8192 3 64 grid approx >> latex_report/report.tex
python3 visualization/draw_table.py 8192 3 128 grid approx >> latex_report/report.tex
python3 visualization/draw_table.py 8192 3 256 grid approx >> latex_report/report.tex
echo "\\midrule" >> latex_report/report.tex
python3 visualization/draw_table.py 8192 4 16 grid approx >> latex_report/report.tex
python3 visualization/draw_table.py 8192 4 32 grid approx >> latex_report/report.tex
python3 visualization/draw_table.py 8192 4 64 grid approx >> latex_report/report.tex
python3 visualization/draw_table.py 8192 4 128 grid approx >> latex_report/report.tex
python3 visualization/draw_table.py 8192 4 256 grid approx >> latex_report/report.tex
echo "\\midrule" >> latex_report/report.tex
python3 visualization/draw_table.py 8192 5 16 grid approx >> latex_report/report.tex
python3 visualization/draw_table.py 8192 5 32 grid approx >> latex_report/report.tex
python3 visualization/draw_table.py 8192 5 64 grid approx >> latex_report/report.tex
python3 visualization/draw_table.py 8192 5 128 grid approx >> latex_report/report.tex
python3 visualization/draw_table.py 8192 5 256 grid approx >> latex_report/report.tex
echo "\\midrule" >> latex_report/report.tex
python3 visualization/draw_table.py 8192 10 16 grid approx >> latex_report/report.tex
python3 visualization/draw_table.py 8192 10 32 grid approx >> latex_report/report.tex
python3 visualization/draw_table.py 8192 10 64 grid approx >> latex_report/report.tex
python3 visualization/draw_table.py 8192 10 128 grid approx >> latex_report/report.tex
python3 visualization/draw_table.py 8192 10 256 grid approx >> latex_report/report.tex
echo "\\midrule" >> latex_report/report.tex

echo "\\end{tabular}
    \\caption{ : \$\\max\\limits_{F \\in \\F} \\left\\lvert\\frac{|F|}{|X|} - \\frac{|F\\cap A|}{|A|}\\right\\rvert\$ for our algorithms on the grid set system averaged over 10 runs.}
\\end{table}" >> latex_report/report.tex

echo "\\end{document}" >> latex_report/report.tex

cd latex_report
latexmk -pdf -silent report.tex
latexmk -c -silent report.tex
