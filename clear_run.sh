rm a.out
rm -r run_artifacts
rm -r latex_report
rm potential.csv
touch potential.csv
echo "algo;n;t;d;table" > potential.csv
rm potential2.csv
touch potential2.csv
echo "algo;n;t;d;table" > potential2.csv
rm results.csv
touch results.csv
echo "algo;n;t;ss_type;m;d;p;max_crossing;avg_crossing;min_crossing;rate_violation;time;max_approx_partition;max_random_sample" > results.csv

