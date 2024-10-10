g++ -fopenmp -Ofast main.cpp
g++ -fopenmp -Ofast main_fb.cpp -o fb

if [ ! -d "run_artifacts" ]; then
  mkdir run_artifacts
fi
	
./a.out -a '1 1 4 1 10 1' -n 8192 -d 2 -t 16 -f grid -k 150 -s -e


if [ ! -d "latex_report" ]; then
    mkdir latex_report
fi

if [ ! -d "latex_report/img" ]; then
    mkdir latex_report/img
fi

if [ ! -d "latex_report/report.tex" ]; then
    touch latex_report/report.tex
fi

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

echo "\\section{Potential function bound violation}

\\begin{figure}[!htb]
\\centering
\\includegraphics[width=0.8\\textwidth]{img/activity_16.png}\\\\
\\caption{Number of potential function violations.}
\\end{figure}
" >> latex_report/report.tex

echo "\\end{document}" >> latex_report/report.tex

cd latex_report
latexmk -pdf -silent report.tex
latexmk -c -silent report.tex
