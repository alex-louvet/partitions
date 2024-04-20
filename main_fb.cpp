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
    cout << "Writing " << fileName << endl;
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

SetSystem replacePoints(SetSystem res){
    int side = floor(sqrt(res.sets.size()));
    float cl = 0;
    float cside = 0;
    for (Set& s : res.sets){
        for (int i = 0; i < res.points.size(); i++){
            if (s.points.at(i)){
                res.points.at(i).coordinates.at(0) += cl;
                res.points.at(i).coordinates.at(1) += cside;
            }
        }
        cl += 1.2;
        if (cl/1.2 > side){
            cl = 0;
            cside += 1.2;
        }
    }
    return res;
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
    string line;
    string file = "./data/facebook/facebook_combined.txt";
    int max = 0;
    vector<vector<int>> edges;
    vector<Point> pt;
    vector<Set> s;
    int d = 1;
    int c;
    vector<int> algoList;
    int t = 32;
    float p = .1;
    string ss_type = "dakota";
    bool save = true;
    bool hask = false;
    int warmup;
    int r;

    while ((c = getopt (argc, argv, "a:n:t:d:r:e:k:s")) != -1){
        switch (c)
        {
        case 'a':
            algoList = simple_tokenizer(optarg);
            break;
        case 'n':
            ss_type = optarg;
            break;
        case 't':
            t = stoi(optarg);
            break;
        case 'd':
            d = stoi(optarg);
            break;
        case 's':
            save = true;
            break;
        case 'r':
            r = stoi(optarg);
            srand(r);
            break;
        case 'e':
            file = optarg;
            break;
        case 'k':
            warmup = stoi(optarg);
            hask = true;
            break;
        case '?':
            if (optopt == 'a' || optopt == 'n' || optopt == 't' || optopt == 'd' || optopt == 'r' || optopt == 'e' || optopt == 'k')
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

    //building set system from network files
    ifstream road_net (file);

    if (road_net.is_open())
    {
        while ( getline (road_net,line) )
        {
            int a = stoi(line.substr(0,line.find(" ")));
            line.erase(0,line.find(" ") + 1);
            int b = stoi(line);
            if (a > max){
                max = a;
                while (edges.size() < max){
                    vector<int> temp;
                    temp.clear();
                    edges.push_back(temp);
                }
            }
            if (b > max){
                max = b;
                while (edges.size() < max){
                    vector<int> temp;
                    temp.clear();
                    edges.push_back(temp);
                }
            }
            edges.at(a-1).push_back(b-1);
            edges.at(b-1).push_back(a-1);
        }
        road_net.close();

    }

    else{
        cout << "Unable to open file";
        return 1;
    }

    for (int i = 0; i < max; i++){
        pt.push_back(Point(2));
    }

    int n = max;
    #pragma omp for
    for (int i = 0; i < max; i++){
        vector<int> steps(max,-1);
        vector<bool> visited(n,false);
        deque<int> file = {i};
        steps.at(i) = 0;
        visited.at(i) = true;
        while(file.size() > 0){
            int a = file.at(0);
            file.pop_front();
            for (int& x : edges.at(a)){
                if (!visited.at(x)){
                    visited.at(x) = true;
                    if (steps.at(x) == -1){
                        steps.at(x) = steps.at(a) + 1;
                    } else {
                        steps.at(x) = min(steps.at(x), steps.at(a) + 1);
                    }
                    if (steps.at(x) < d){
                        file.push_back(x);
                    }
                }
            }
        }
        vector<bool> temp;
        for (int i = 0; i < n; i++){
            if (steps.at(i) > -1 && steps.at(i) <= d){
                temp.push_back(1);
            } else {
                temp.push_back(0);
            }
        }
        s.push_back(temp);
    }

    SetSystem test = SetSystem(pt,s);

    int m = test.sets.size();
    n = test.points.size();
    if (!hask){
        warmup = floor(log(m*n));
    }
    test.buildAdjacency(false);

    for (int k = 0; k < algoList.size()/2;k++){
        for (int ite = 0; ite < algoList.at(2*k+1); ite ++){
            cout << "algo "  << algoList.at(2*k) << ", n = "<< n << ", d = " << d << ", t = " << t << " " << ss_type << endl;
            test.resetWeight();
            auto start_time = chrono::high_resolution_clock::now();
            Result res;
            if (algoList.at(2*k) == 1){
                res = partition_min_stats(test,t);
                if (save){
                    res = replacePoints(res);
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
                    res = replacePoints(res);
                    writeCSVFile(res, to_string(time(NULL)) + "_" + ss_type + "_par_sw_sample.csv");
                }
            }
            auto end_time = chrono::high_resolution_clock::now();
            std::chrono::duration<double> duration = end_time - start_time;
            printf("Runtime: %.2fs\n", duration.count());
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

            MyFile << algoList.at(2*k) << ";" << n << ";" << t << ";" << ss_type << ";" << m << ";" << d << ";" << p << ";" << maxcrossing << ";" << avgcrossing << ";" << mincrossing << ";" << rate_stats << ";" << time << endl;
        }
    }
    
    return 0;
}
