#include <cstdlib>
#include <ios>
#include <iostream>
#include <fstream>
#include <pthread.h>
#include <string>

#include "partition.cpp"

using namespace std;

void writeCSVFile(SetSystem res, string fileName){
    ofstream MyFile(fileName);
    for (Point& p : res.points) {
        for (const float& f : p.coordinates) {
            MyFile << f << ",";
        }
        MyFile << "\n";
    }
    MyFile << "sets\n";
    for (Set& s : res.sets) {
        for (const bool& b : s.points) {
            MyFile << b << ",";
        }
        MyFile << "\n";
    }  
}

/* a.out argv 
 1 -> algorithm (1: min, 2: min + rate, 3: uniform rate)
 2 -> n (number of points)
 3 -> t (number of partitions)
 4 -> type of set system (1: grid, 2: random hs, 3: grid graph)
 5 -> d (dimension)
 6 -> save (1 to save)
*/
int main(int argc, char** argv){
    srand(time(NULL));
    int n = 1024;
    int d = 2;
    int m = d*sqrt(n);
    int t = 16;
    SetSystem test;
    string ss_type = "grid";
    if (argc >= 3){
        n = stoi(argv[2]);

    }
    if (argc >= 4){
        t = stoi(argv[3]);

    }
    if (argc >= 5){
        if (stoi(argv[4]) == 1){
            if (argc >= 6){
                d = stoi(argv[5]);
                test = Grid(n,stoi(argv[5]));
            }
            else{
                test = Grid(n,2);
            }
            
        }
        if (stoi(argv[4]) == 2){
            ss_type = "random_hs";
            if (argc >= 6){
                d = stoi(argv[5]);
                test = RandomHyperplanes(n,stoi(argv[5]),n*log(n));
            }
            else{
                test = RandomHyperplanes(n,2, n*log(n));
            }
            
        }
        if (stoi(argv[4]) == 3){
            ss_type = "grid_graph";
            if (argc >= 6){
                d = stoi(argv[5]);
                test = GridGraph(sqrt(n),stoi(argv[5]));
            }
            else{
                test = GridGraph(sqrt(n),1);
            }        
        }
    } else{
        test = Grid(n,d);
    }
    cout << "algo "  << argv[1] << ", n = "<< n << ", d = " << d << ", t = " << t << " " << ss_type << endl;
    clock_t tStart = clock();
    Result res;
    if (argc >= 2){
        if (stoi(argv[1]) == 1){
            res = partition_min_violation(test,t);
            if (argc >= 7 && stoi(argv[6]) == 1){
                writeCSVFile(res, to_string(time(NULL)) + "_" + ss_type + "_algo_min.csv");
            }
        }

        if (stoi(argv[1]) == 2){
            res = partition_min_rate_violation(test,t);
            if (argc >= 7 && stoi(argv[6]) == 1){
                writeCSVFile(res, to_string(time(NULL)) + "_" + ss_type + "_algo_min_rate.csv");
            }
        }

        if (stoi(argv[1]) == 3){
            res = partition_min_rate_violation(test,t);
            if (argc >= 7 && stoi(argv[6]) == 1){
                writeCSVFile(res, to_string(time(NULL)) + "_" + ss_type + "_algo_rate.csv");
            }
        }

        if (stoi(argv[1]) == 4){
            res = partition_sampling_violation(test,t,1.0/2);
            if (argc >= 7 && stoi(argv[6]) == 1){
                writeCSVFile(res, to_string(time(NULL)) + "_" + ss_type + "_sampling1.csv");
            }
        }

        if (stoi(argv[1]) == 5){
            res = partition_sampling_fixed_violation(test,t,m/2);
            if (argc >= 7 && stoi(argv[6]) == 1){
                writeCSVFile(res, to_string(time(NULL)) + "_" + ss_type + "_sampling2.csv");
            }
        }

        if (stoi(argv[1]) == 6){
            res = partition_sampling_fixed_violation2(test,t,m);
            if (argc >= 7 && stoi(argv[6]) == 1){
                writeCSVFile(res, to_string(time(NULL)) + "_" + ss_type + "_sampling2.csv");
            }
        }

    } else {
        res = partition_min(test,t);
        if (argc >= 7 && stoi(argv[6]) == 1){
            writeCSVFile(res, to_string(time(NULL)) + "_" + ss_type + "_algo_min.csv");
        }
    }
    printf("Time taken: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
    double time = (double)(clock() - tStart)/CLOCKS_PER_SEC;

    int maxcrossing = res.intersections.at(0);
    int mincrossing = res.intersections.at(0);
    float avgcrossing = 0.;

    for (int j = 0; j < res.intersections.size(); j++){
        if (res.intersections.at(j) > maxcrossing){
            maxcrossing = res.intersections.at(j);
        }
        if (res.intersections.at(j) < mincrossing){
            mincrossing = res.intersections.at(j);
        }
        avgcrossing += static_cast<float>(res.intersections.at(j));
    }
    avgcrossing /= static_cast<float>(res.intersections.size());

    int rate_violation = 0;
    for (int j = 0; j < res.weights.size(); j++){
        if (res.weights.at(j) > 2*pow(static_cast<float>(j%(n/t-1)+1),1.0/d)){
            rate_violation++;
        };
    }

    ofstream MyFile("results.csv",std::ios_base::app);

    MyFile << argv[1] << ";" << n << ";" << t << ";" << argv[4] << ";" << m << ";" << d << ";" << maxcrossing << ";" << avgcrossing << ";" << mincrossing << ";" << rate_violation << ";" << time << endl;

    return 0;
}