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
    if (argc >= 2){
        if (stoi(argv[1]) == 1){
            SetSystem res = partition_min(test,t);
            if (argc >= 7 && stoi(argv[6]) == 1){
                writeCSVFile(res, to_string(time(NULL)) + "_" + ss_type + "_algo_min.csv");
            }
        }

        if (stoi(argv[1]) == 2){
            SetSystem res = partition_min_rate(test,t);
            if (argc >= 7 && stoi(argv[6]) == 1){
                writeCSVFile(res, to_string(time(NULL)) + "_" + ss_type + "_algo_min_rate.csv");
            }
        }

        if (stoi(argv[1]) == 3){
            SetSystem res = partition_rate(test,t);
            if (argc >= 7 && stoi(argv[6]) == 1){
                writeCSVFile(res, to_string(time(NULL)) + "_" + ss_type + "_algo_rate.csv");
            }
        }

    } else {
        SetSystem res = partition_min(test,t);
        if (argc >= 7 && stoi(argv[6]) == 1){
            writeCSVFile(res, to_string(time(NULL)) + "_" + ss_type + "_algo_min.csv");
        }
    }
    printf("Time taken: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);

    return 0;
}