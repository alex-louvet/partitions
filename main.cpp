#include <cstdlib>
#include <ios>
#include <iostream>
#include <fstream>
#include <pthread.h>
#include <string>
#include <vector>

#include "partition.cpp"

using namespace std;

void writeCSVFile(SetSystem res, string fileName){
    ofstream MyFile("results/" + fileName);
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

/* a.out argv 
 1 -> algorithm 
        1: min,
        2: min + rate, 
        3: uniform rate,
        4: sampling with proba for each set to be sampled,
        5: sampling w/fixed number of sets sampled every rounds and partition weight added over time,
        6: sampling w/fixed number of sets sampled every rounds and partition weight recomputed with sample)
        7: Rate + as soon as 1 pt gets below rate no more weight update and updates high weight sets first
    Can be combined in the first argument by giving algo number and # of iterations, ex: '1 2 2 1' runs min twice and then min+rate once on the same set system
 2 -> n (number of points)
 3 -> t (number of partitions)
 4 -> type of set system (1: grid, 2: random hs, 3: grid graph)
 5 -> d (dimension)
 6 -> save (1 to save)*/

int main(int argc, char** argv){
    srand(time(NULL));
    int n = 1024;
    int d = 2;
    int m = d*sqrt(n);
    int t = 16;
    SetSystem test;
    string ss_type = "grid";
    vector<int> algoList = simple_tokenizer(argv[1]);
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
            m = n*log(n);
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

        if (stoi(argv[4]) == 4){
            ss_type = "linear_grid";
            if (argc >= 6){
                d = stoi(argv[5]);
                test = LinearGrid(n,stoi(argv[5]));
            }
            else{
                test = LinearGrid(n,2);
            }
        }

        if (stoi(argv[4]) == 5){
            ss_type = "exponential_grid";
            if (argc >= 6){
                d = stoi(argv[5]);
                test = ExponentialGrid(n,stoi(argv[5]));
            }
            else{
                test = ExponentialGrid(n,2);
            }
        }

        if (stoi(argv[4]) == 6){
            ss_type = "directional_grid";
            if (argc >= 6){
                d = stoi(argv[5]);
                test = DirectionalGrid(n,stoi(argv[5]));
            }
            else{
                test = DirectionalGrid(n,2);
            }
        }
        
    } else{
        test = Grid(n,d);
    }
    m = test.sets.size();
    test.buildAdjacency();
    for (int k = 0; k < algoList.size()/2;k++){
        for (int ite = 0; ite < algoList.at(2*k+1); ite ++){
            cout << "algo "  << algoList.at(2*k) << ", n = "<< n << ", d = " << d << ", t = " << t << " " << ss_type << endl;
            clock_t tStart = clock();
            Result res;
            if (argc >= 2){
                if (algoList.at(2*k) == 1){
                    res = partition_min_stats(test,t);
                    if (argc >= 7 && stoi(argv[6]) == 1){
                        writeCSVFile(res, to_string(time(NULL)) + "_" + ss_type + "_algo_min.csv");
                    }
                }

                if (algoList.at(2*k) == 2){
                    res = partition_min_rate_stats(test,t);
                    if (argc >= 7 && stoi(argv[6]) == 1){
                        writeCSVFile(res, to_string(time(NULL)) + "_" + ss_type + "_algo_min_rate.csv");
                    }
                }

                if (algoList.at(2*k) == 3){
                    res = partition_min_rate_stats(test,t);
                    if (argc >= 7 && stoi(argv[6]) == 1){
                        writeCSVFile(res, to_string(time(NULL)) + "_" + ss_type + "_algo_rate.csv");
                    }
                }

                if (algoList.at(2*k) == 4){
                    res = partition_sampling_stats(test,t,1.0/pow(static_cast<float>(m),2/3));
                    if (argc >= 7 && stoi(argv[6]) == 1){
                        writeCSVFile(res, to_string(time(NULL)) + "_" + ss_type + "_sampling1.csv");
                    }
                }

                if (algoList.at(2*k) == 5){
                    res = partition_sampling_fixed_stats(test,t,pow(static_cast<float>(m),2/3));
                    if (argc >= 7 && stoi(argv[6]) == 1){
                        writeCSVFile(res, to_string(time(NULL)) + "_" + ss_type + "_sampling2.csv");
                    }
                }

                if (algoList.at(2*k) == 6){
                    res = partition_sampling_fixed_stats2(test,t,pow(static_cast<float>(m),2/3));
                    if (argc >= 7 && stoi(argv[6]) == 1){
                        writeCSVFile(res, to_string(time(NULL)) + "_" + ss_type + "_sampling2.csv");
                    }
                }

                if (algoList.at(2*k) == 7){
                    res = no_weight_update_deque_insert_end(test,t);
                    if (argc >= 7 && stoi(argv[6]) == 1){
                        writeCSVFile(res, to_string(time(NULL)) + "_" + ss_type + "_no_update_deque_end.csv");
                    }
                }

                if (algoList.at(2*k) == 8){
                    res = no_weight_update_deque_insert_middle(test,t);
                    if (argc >= 7 && stoi(argv[6]) == 1){
                        writeCSVFile(res, to_string(time(NULL)) + "_" + ss_type + "_no_update_deque_middle.csv");
                    }
                }

                if (algoList.at(2*k) == 9){
                    res = no_weight_update_insert_sorted(test,t);
                    if (argc >= 7 && stoi(argv[6]) == 1){
                        writeCSVFile(res, to_string(time(NULL)) + "_" + ss_type + "_no_update_sorted.csv");
                    }
                }

                if (algoList.at(2*k) == 10){
                    res = no_weight_update_point_sort_set_sort(test,t);
                    if (argc >= 7 && stoi(argv[6]) == 1){
                        writeCSVFile(res, to_string(time(NULL)) + "_" + ss_type + "_no_update_sorted_pt.csv");
                    }
                }

                if (algoList.at(2*k) == 11){
                    res = no_weight_update_point_sort_set_deque(test,t);
                    if (argc >= 7 && stoi(argv[6]) == 1){
                        writeCSVFile(res, to_string(time(NULL)) + "_" + ss_type + "_no_update_deque_middle_pt.csv");
                    }
                }

                if (algoList.at(2*k) == 12){
                    res = no_weight_update_insert_sorted2(test,t);
                    if (argc >= 7 && stoi(argv[6]) == 1){
                        writeCSVFile(res, to_string(time(NULL)) + "_" + ss_type + "_no_update_sorted.csv");
                    }
                }

                if (algoList.at(2*k) == 13){
                    res = partition_no_set(test,t,dw);
                    if (argc >= 7 && stoi(argv[6]) == 1){
                        writeCSVFile(res, to_string(time(NULL)) + "_" + ss_type + "_no_set.csv");
                    }
                }

            } else {
                res = partition_min_stats(test,t);
                if (argc >= 7 && stoi(argv[6]) == 1){
                    writeCSVFile(res, to_string(time(NULL)) + "_" + ss_type + "_algo_min.csv");
                }
            }
            printf("Time taken: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
            double time = (double)(clock() - tStart)/CLOCKS_PER_SEC;

            if (res.intersections.size() == 0){
                for (int j = 0; j < m; j++){
                    res.intersections.push_back(0);
                }
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

            MyFile << algoList.at(2*k) << ";" << n << ";" << t << ";" << argv[4] << ";" << m << ";" << d << ";" << maxcrossing << ";" << avgcrossing << ";" << mincrossing << ";" << rate_stats << ";" << time << endl;
        }
    }
    
    return 0;
}

/*int main(){
    vector<Set> sets;
    vector<int> L;
    for (int i = 0; i < 1; i++){
        Set s(0);
        s.weight = 0*i;
        sets.push_back(s);
    }
    Set s(0);
    s.weight = 2;
    cout << "S ";
    for (Set& s : sets){
        cout << s.weight << " ";
    }
    cout << endl;
    cout << "Weight s " << s.weight << endl;
    cout << "Res " << insertAt(sets,s) << endl;
    sets.insert(next(sets.begin(),insertAt(sets,s)),s);
    cout << "S ";
    for (Set& s : sets){
        cout << s.weight << " ";
    }
    return 0;
}*/
