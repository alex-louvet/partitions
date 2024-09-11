if [ -d "a.out" ]; then
    rm -r a.out
fi
if [ -d "latex_report" ]; then
    rm -r latex_report
fi
if [ -d "run_artifacts" ]; then
    rm -r run_artifacts
fi
if [ -d "potential.csv" ]; then
    rm -r potential.csv
fi
touch potential.csv
echo "algo;n;t;d;table" > potential.csv
rm potential2.csv
if [ -d "potential2.csv" ]; then
    rm -r potential2.csv
fi
echo "algo;n;t;d;table" > potential2.csv
if [ -d "results.csv" ]; then
    rm -r results.csv
fi
touch results.csv
echo "algo;n;t;ss_type;m;d;p;max_crossing;avg_crossing;min_crossing;rate_violation;time;max_approx_partition;max_random_sample" > results.csv

