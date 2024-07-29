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
#include "part_size.cpp"

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
    int r = time(NULL);
    srand(r);
    int n = 1024;
    int d = 2;
    int m = d*sqrt(n);
    int t = 16;
    float p = .1;
    SetSystem test;
    string ss_type = "grid";
    bool save = false;
    int c;
    int centers = 3;
    vector<int> algoList;
    bool hask = false;
    int warmup = floor(log(m*n));
    string filename = "";
    float constant = 2.0;
    bool approx = false;
    while ((c = getopt (argc, argv, "a:n:t:d:f:r:p:m:i:c:k:es")) != -1){
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
        case 'm':
            m = stoi(optarg);
            break;
        case 'f':
            ss_type = optarg;
            break;
        case 's':
            save = true;
            break;
        case 'r':
            r = stoi(optarg);
            srand(r);
            break;
        case 'p':
            p = stof(optarg);
            break;
        case 'i':
            filename = optarg;
            break;
        case 'k':
            warmup = stoi(optarg);
            hask = true;
            break;
        case 'g':
            centers = stoi(optarg);
            break;
        case 'c':
            constant = stof(optarg);
            break;
        case 'e':
            approx = true;
            break;
        case '?':
            if (optopt == 'a' || optopt == 'n' || optopt == 't' || optopt == 'd' || optopt == 'f' || optopt == 'r' || optopt == 'p' || optopt == 'm' || optopt == 'i' || optopt == 'c' || optopt == 'k' || optopt == 'g')
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

    if (ss_type == "grid"){
        test = Grid(n,d);
    }
    else if (ss_type == "grid_centers"){
        test = GridWithCenters(n,d,centers);
    }
    else if (ss_type == "random_hs"){
        m = n*log(n);
        test = RandomHyperplanes(n,d,m);
    }
    else if (ss_type == "grid_graph"){
        test = GridGraph(sqrt(n),d);      
    }
    else if (ss_type == "linear_grid"){
        test = LinearGrid(n,d);
    }
    else if (ss_type == "exponential_grid"){
        test = ExponentialGrid(n,d);
    }
    else if (ss_type == "directed_grid"){
        test = DirectionalGrid(n,d);
    }
    else if (ss_type == "random"){
        test = Random(n,d,m,p);
    }
    else if (ss_type == "projective_plane"){
        if (!isPrime(n)){
            fprintf (stderr, "Projective plane needs n to be prime, given: '%d'.\n", n);
            return 1;
        }
        test = ProjectivePlane(n);
    }
    else if (ss_type == "ERGGraph"){
        test = ERGGraph(n,d,p);
    }
    else if (ss_type == "power_law"){
        test = PowerLaw(n,d,p,r);
    }
    else if (ss_type == "concentric_circles"){
        test = ConcentricCircles(n,m);
    }
    else if (ss_type == "file"){
        test = SetSystem(filename);
    }
    else {
        fprintf (stderr, "Unknown set system, given: '%s'.\n", ss_type);
        return 1;
    }
    m = test.sets.size();
    n = test.points.size();
    vector<int> list = equi_distro(n,t);
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
            if (argc >= 2){
                if (algoList.at(2*k) == 1){
                    res = partition_min_stats(test,t, equi_distro(n,t));
                    if (save){
                        writeCSVFile(res, to_string(time(NULL)) + "_" + ss_type + "_algo_min.csv");
                    }
                }

                if (algoList.at(2*k) == 2){
                    res = partition_rate_stats(test,t, constant, equi_distro(n,t));
                    if (save){
                        writeCSVFile(res, to_string(time(NULL)) + "_" + ss_type + "_algo_rate.csv");
                    }
                }

                if (algoList.at(2*k) == 3){
                    res = partition_sampling(test,t,pow(static_cast<float>(m),2/3), constant, equi_distro(n,t));
                    if (save){
                        writeCSVFile(res, to_string(time(NULL)) + "_" + ss_type + "_sampling.csv");
                    }
                }

                if (algoList.at(2*k) == 4){
                    res = no_weight_update_deque_insert_middle(test,t, constant, equi_distro(n,t));
                    if (save){
                        writeCSVFile(res, to_string(time(NULL)) + "_" + ss_type + "_no_update_deque_middle.csv");
                    }
                }

                if (algoList.at(2*k) == 5){
                    res = partition_distance_set_weight_par(test,t,l1,0, equi_distro(n,t));
                    if (save){
                        writeCSVFile(res, to_string(time(NULL)) + "_" + ss_type + "_par_l1.csv");
                    }
                }

                if (algoList.at(2*k) == 6){
                    res = partition_distance_set_weight_par(test,t,l2,warmup, equi_distro(n,t));
                    if (save){
                        writeCSVFile(res, to_string(time(NULL)) + "_" + ss_type + "_par_l2.csv");
                    }
                }

                if (algoList.at(2*k) == 7){
                    res = partition_distance_set_weight_par(test,t,sw,warmup, equi_distro(n,t));
                    if (save){
                        writeCSVFile(res, to_string(time(NULL)) + "_" + ss_type + "_par_sw.csv");
                    }
                }

                if (algoList.at(2*k) == 8){
                    res = partition_distance_set_weight_par(test,t,dw,warmup, equi_distro(n,t));
                    if (save){
                        writeCSVFile(res, to_string(time(NULL)) + "_" + ss_type + "_par_dw.csv");
                    }
                }
                
                if (algoList.at(2*k) == 9){
                    res = partition_distance_set_weight_par(test,t,sw_weighted,warmup, equi_distro(n,t));
                    if (save){
                        writeCSVFile(res, to_string(time(NULL)) + "_" + ss_type + "_par_sw_weighted.csv");
                    }
                }

                if (algoList.at(2*k) == 10){
                    res = partition_distance_set_weight_par(test,t,sw_weighted_w_sample, warmup, equi_distro(n,t));
                    if (save){
                        writeCSVFile(res, to_string(time(NULL)) + "_" + ss_type + "_par_sw_sample_equi.csv");
                    }
                }

                if (algoList.at(2*k) == 11){
                    list = interval_size(n,t,n/t/2,4*n/t);
                    res = partition_distance_set_weight_par(test,t,sw_weighted_w_sample, warmup, list);
                    if (save){
                        writeCSVFile(res, to_string(time(NULL)) + "_" + ss_type + "_par_sw_sample_interval.csv");
                    }
                }

                if (algoList.at(2*k) == 12){
                    list = same_number_different_size_3(n,t);
                    res = partition_distance_set_weight_par(test,t,sw_weighted_w_sample, warmup, list);
                    if (save){
                        writeCSVFile(res, to_string(time(NULL)) + "_" + ss_type + "_par_sw_sample_3.csv");
                    }
                }
                if (algoList.at(2*k) == 13){
                    list = same_number_different_size_linear(n,t, 0.9);
                    res = partition_distance_set_weight_par(test,t,sw_weighted_w_sample, warmup, list);
                    if (save){
                        writeCSVFile(res, to_string(time(NULL)) + "_" + ss_type + "_par_sw_sample_linear.csv");
                    }
                }
                if (algoList.at(2*k) == 14){
                    list = same_number_different_size_3(n,t);
                    res = partition_min_stats(test, t, list);
                    if (save){
                        writeCSVFile(res, to_string(time(NULL)) + "_" + ss_type + "_min_3.csv");
                    }
                }
                if (algoList.at(2*k) == 15){
                    list = same_number_different_size_linear(n,t,0.9);
                    res = partition_min_stats(test, t, list);
                    if (save){
                        writeCSVFile(res, to_string(time(NULL)) + "_" + ss_type + "_min_linear.csv");
                    }
                }
                if (algoList.at(2*k) == 16){
                    list = ninty_percent(n,t,(n/(8*t)));
                    res = partition_min_stats(test, t, list);
                    if (save){
                        writeCSVFile(res, to_string(time(NULL)) + "_" + ss_type + "_min_90percent.csv");
                    }
                }


            } else {
                fprintf (stderr, "Algorithm -%c not defined.\n", algoList.at(2*k));
                return 1;
            }
            auto end_time = chrono::high_resolution_clock::now();
            std::chrono::duration<double> duration = end_time - start_time;
            printf("Runtime: %.2fs\n", duration.count());
            auto time = duration.count();

            float max_random_sample = 0.;
            float max_approx_partition = 0.;

            if (approx){
                vector<int> random_sample;
                while (random_sample.size() < t){
                    int a = rand()%n;
                    bool test = false;
                    for (int& i : random_sample){
                        if (i == a){
                            test = true;
                            break;
                        }
                    }
                    if (!test){
                        random_sample.push_back(a);
                    }
                }

                vector<int> approx_partition;
                for (Set& s : res.sets){
                    vector<int> temp;
                    for (int i = 0; i < n; i++){
                        if (s.points.at(i)){
                            temp.push_back(i);
                        }
                    }
                    approx_partition.push_back(temp.at(rand()%temp.size()));
                }

                vector<float> approx_partition_count(m,0);
                vector<float> random_sample_count(m,0);
                #pragma omp parallel for
                for (int j = 0; j < m; j ++){
                    int set_size = 0;
                    int approx_partition_inter = 0;
                    int random_sample_inter = 0;
                    for (int i = 0; i < n; i++){
                        if (test.sets.at(j).points.at(i)){
                            set_size ++;
                        }
                    }
                    for (int& i : approx_partition){
                        if (test.sets.at(j).points.at(i)){
                            approx_partition_inter ++;
                        }
                    }                
                    for (int& i : random_sample){
                        if (test.sets.at(j).points.at(i)){
                            random_sample_inter ++;
                        }
                    }
                    approx_partition_count.at(j) = abs(static_cast<float>(set_size)/n - static_cast<float>(approx_partition_inter)/t);
                    random_sample_count.at(j) = abs(static_cast<float>(set_size)/n - static_cast<float>(random_sample_inter)/t);
                }

                for (int j = 0; j < m; j ++){
                    if (max_approx_partition < approx_partition_count.at(j)){
                        max_approx_partition = approx_partition_count.at(j);
                    }
                    if (max_random_sample < random_sample_count.at(j)){
                        max_random_sample = random_sample_count.at(j);
                    }
                }
                cout << "Max epsilon: "<< max_random_sample << " (random) - " << max_approx_partition << " (with simplicial partition)" << endl;

            }

            cout << "Computing intersection number" << endl;
            if (res.intersections.size() == 0){
                for (int j = 0; j < m; j++){
                    res.intersections.push_back(0);
                }
                #pragma omp parallel for
                for (int j = 0; j < m; j ++){
                    for (int i = 0; i < res.sets.size(); i++){
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
            /*if (res.weights.size() > 0){
                ofstream pot("potential.csv",std::ios_base::app);
                pot << algoList.at(2*k) << ";" << n << ";" << t << ";" << d << ";";
                for (int j = 0; j < list.size(); j++){
                    int part_stat = 0;
                    for (int jp = 0; jp < list.at(j) - 1; jp++){
                        if (res.weights.at(partialSum(list,j) - j+jp) > constant*pow(static_cast<float>(jp+1),1.0/d)){
                            part_stat++;
                        };
                    }
                    pot << part_stat << ",";
                }
                pot << endl;

                int rate_stats2 = 0;
                ofstream pot2("potential2.csv",std::ios_base::app);
                pot2 << algoList.at(2*k) << ";" << n << ";" << t << ";" << d << ";";
                for (int j = 0; j < list.size(); j++){
                    int part_stat = 0;
                    for (int jp = 0; jp < list.at(j) - 1; jp++){
                        if (res.weights.at(partialSum(list,j) - j+jp) > t/(t-j)*constant*pow(static_cast<float>(jp+1),1.0/d)){
                            part_stat++;
                        };
                    }
                    pot2 << part_stat << ",";
                }
                pot2 << endl;
            }*/
            
            for (int j = 0; j < list.size(); j++){
                for (int jp = 0; jp < list.at(j) - 1; jp++){
                    if (res.weights.at(partialSum(list,j-1) + jp - j +1) > constant*pow(jp + 1,1.0/d)){
                        rate_stats++;
                    };
                }
                //cout << res.weights.at(j) << " ";
            }

            ofstream MyFile("results.csv",std::ios_base::app);

            MyFile << algoList.at(2*k) << ";" << n << ";" << t << ";" << ss_type << ";" << m << ";" << d << ";" << p << ";" << maxcrossing << ";" << avgcrossing << ";" << mincrossing << ";" << rate_stats << ";" << time << ";" << max_approx_partition << ";" << max_random_sample << endl;
        }
    }
    
    return 0;
}
