#include <cmath>
#include <cstdlib>
#include <ios>
#include <iostream>
#include <fstream>
#include <pthread.h>
#include <string>
#include <vector>
#include <chrono>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

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

vector<int> simple_tokenizer(string s)
{
    vector<int> res;
    stringstream ss(s);
    string word;
    while (ss >> word) {
        res.push_back(stoi(word));
    }
    return res;
}

bool isPrime(int n){
    for (int i = 2; i <= sqrt(n); i++){
        if (n%i == 0){
            return false;
        }
    }
    return true;
}

int main(int argc, char** argv){
    srand(time(NULL));
    int n = 1024;
    int d = 2;
    int m = d*sqrt(n);
    int t = 16;
    float p = .1;
    SetSystem test;
    string ss_type = "grid";
    bool save = false;
    int c;
    vector<int> algoList;
    bool hask = false;
    int warmup = floor(log(m*n));
    while ((c = getopt (argc, argv, "a:n:t:d:f:r:k:s")) != -1){
        switch (c)
        {
        case 'a':
            algoList = simple_tokenizer(optarg);
            break;
        case 'n':
            n = stoi(optarg);
            break;
        case 't':
            t = stoi(optarg);
            break;
        case 'd':
            d = stoi(optarg);
            break;
        case 'f':
            ss_type = optarg;
            break;
        case 's':
            save = true;
            break;
        case 'r':
            srand(stoi(optarg));
            break;
        case 'k':
            warmup = stoi(optarg);
            hask = true;
            break;
        case '?':
            if (optopt == 'a' || optopt == 'n' || optopt == 't' || optopt == 'd' || optopt == 'f' || optopt == 'r' || optopt == 'k')
            fprintf (stderr, "Option -%c requires an argument.\n", optopt);
            else if (isprint (optopt))
            fprintf (stderr, "Unknown option `-%c'.\n", optopt);
            else
            fprintf (stderr,
                    "Unknown option character `\\x%x'.\n",
                    optopt);
            return 1;
        default:
            abort ();
        }
    }
    if (!hask){
        warmup = floor(log(m*n));
    }
    if (ss_type == "grid"){
        test = Grid(n,d);
    }
    if (ss_type == "random_hs"){
        m = n*log(n);
        test = RandomHyperplanes(n,d,m);
    }
    if (ss_type == "grid_graph"){
        test = GridGraph(sqrt(n),d);      
    }
    if (ss_type == "linear_grid"){
        test = LinearGrid(n,d);
    }
    if (ss_type == "exponential_grid"){
        test = ExponentialGrid(n,d);
    }
    if (ss_type == "directional_grid"){
        test = DirectionalGrid(n,d);
    }
    if (ss_type == "projective_plane"){
        if (!isPrime(n)){
            fprintf (stderr, "Projective plane needs n to be prime, given: '%d'.\n", n);
            return 1;
        }
        test = ProjectivePlane(n);
    }
    if (ss_type == "ERGGraph"){
        test = ERGGraph(n,d,p);
    }
    m = test.sets.size();
    n = test.points.size();
    test.buildAdjacency(false);
    for (Set& s : test.sets){
        cout << s.points_indices.size() << endl;
    }
    
    for (int k = 0; k < algoList.size()/2;k++){
        for (int ite = 0; ite < algoList.at(2*k+1); ite ++){
            cout << "algo "  << algoList.at(2*k) << ", n = "<< n << ", d = " << d << ", t = " << t << " " << ss_type << endl;
            test.resetWeight();
            auto start_time = chrono::high_resolution_clock::now();
            Result res;
            if (argc >= 2){
                if (algoList.at(2*k) == 1){
                    res = partition_min_stats(test,t);
                    if (save){
                        writeCSVFile(res, to_string(time(NULL)) + "_" + ss_type + "_algo_min.csv");
                    }
                }

                if (algoList.at(2*k) == 2){
                    res = partition_rate_stats(test,t);
                    if (save){
                        writeCSVFile(res, to_string(time(NULL)) + "_" + ss_type + "_algo_rate.csv");
                    }
                }

                if (algoList.at(2*k) == 3){
                    res = partition_sampling(test,t,pow(static_cast<float>(m),2/3));
                    if (save){
                        writeCSVFile(res, to_string(time(NULL)) + "_" + ss_type + "_sampling.csv");
                    }
                }

                if (algoList.at(2*k) == 4){
                    res = no_weight_update_deque_insert_middle(test,t);
                    if (save){
                        writeCSVFile(res, to_string(time(NULL)) + "_" + ss_type + "_no_update_deque_middle.csv");
                    }
                }

                if (algoList.at(2*k) == 5){
                    res = partition_distance_set_weight_par(test,t,l1,0);
                    if (save){
                        writeCSVFile(res, to_string(time(NULL)) + "_" + ss_type + "_par_l1.csv");
                    }
                }

                if (algoList.at(2*k) == 6){
                    res = partition_distance_set_weight_par(test,t,l2,warmup);
                    if (save){
                        writeCSVFile(res, to_string(time(NULL)) + "_" + ss_type + "_par_l2.csv");
                    }
                }

                if (algoList.at(2*k) == 7){
                    res = partition_distance_set_weight_par(test,t,sw,warmup);
                    if (save){
                        writeCSVFile(res, to_string(time(NULL)) + "_" + ss_type + "_par_sw.csv");
                    }
                }

                if (algoList.at(2*k) == 8){
                    res = partition_distance_set_weight_par(test,t,dw,warmup);
                    if (save){
                        writeCSVFile(res, to_string(time(NULL)) + "_" + ss_type + "_par_dw.csv");
                    }
                }
                
                if (algoList.at(2*k) == 9){
                    res = partition_distance_set_weight_par(test,t,sw_weighted,warmup);
                    if (save){
                        writeCSVFile(res, to_string(time(NULL)) + "_" + ss_type + "_par_sw_weighted.csv");
                    }
                }

                if (algoList.at(2*k) == 10){
                    res = partition_distance_set_weight_par(test,t,sw_weighted_w_sample, warmup);
                    if (save){
                        writeCSVFile(res, to_string(time(NULL)) + "_" + ss_type + "_par_sw_sample.csv");
                    }
                }

            } else {
                fprintf (stderr, "Algorithm -%c not defined.\n", algoList.at(2*k));
                return 1;
            }
            auto end_time = chrono::high_resolution_clock::now();
            std::chrono::duration<double> duration = end_time - start_time;
            printf("Time taken: %.2fs\n", duration.count());
            auto time = duration.count();

            cout << "Computing intersection number" << endl;
            if (res.intersections.size() == 0){
                for (int j = 0; j < m; j++){
                    res.intersections.push_back(0);
                }
                #pragma omp parallel for
                for (int j = 0; j < m; j ++){
                    for (int i = 0; i < t; i++){
                        int start = -1;
                        for (int k = 0; k < n; k++){
                            if (res.sets.at(i).points.at(k)){
                                if (start == -1){
                                    start = k;
                                } else {
                                    if (intersects(Edge(start,k),test.sets.at(j))){
                                        res.intersections.at(j) += 1;
                                        break;
                                    }
                                }
                            }
                        }
                    }
                }
            }

            int maxcrossing = 0;
            int mincrossing = 0;
            float avgcrossing = 0.;
            
            if (res.intersections.size() > 0){
                maxcrossing = res.intersections.at(0);
                mincrossing = res.intersections.at(0);

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
            }

            int rate_stats = 0;
            for (int j = 0; j < res.weights.size(); j++){
                if (res.weights.at(j) > 2*pow(static_cast<float>(j%(n/t-1)+1),1.0/d)){
                    rate_stats++;
                };
            }

            ofstream MyFile("results.csv",std::ios_base::app);

            MyFile << algoList.at(2*k) << ";" << n << ";" << t << ";" << ss_type << ";" << m << ";" << d << ";" << maxcrossing << ";" << avgcrossing << ";" << mincrossing << ";" << rate_stats << ";" << time << endl;
        }
    }
    
    return 0;
}
