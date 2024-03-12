#include <cstdlib>
#include <ios>
#include <iostream>
#include <fstream>
#include <pthread.h>
#include <string>
#include <vector>
#include <chrono>

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

int main(int argc, char** argv){
    SetSystem test;
    
    ifstream fin;
    fin.open("book_tags.csv", ios::in);
    string line;
    string delimiter = ",";
    vector<tuple<int,int,int>> ss;
    int c = 0;
    while(!fin.eof() && c < stoi(argv[2])){
        fin>>line;
        size_t pos = 0;
        string token;
        vector<int> temp;
        while ((pos = line.find(delimiter)) != string::npos) {
            token = line.substr(0, pos);
            temp.push_back(stoi(token));
            line.erase(0, pos + delimiter.length());
        }
        //line.pop_back();
        temp.push_back(stoi(line));
        ss.push_back(make_tuple(temp.at(0),temp.at(1),temp.at(2)));
        c ++;
    }
    int n = 0;
    int m = 0;
    int t = stoi(argv[1]);


    for (tuple<int,int,int>& t : ss){
        if (get<0>(t) > n){
            n = get<0>(t);
        }
        if (get<1>(t) > m){
            m = get<1>(t);
        }
    }

    vector<Point> pts;
    for (int i = 0; i < n+1; i++){
        pts.push_back(Point(1));
    }
    vector<vector<bool>> sets;
    for (int i = 0; i < m+1; i++){
        vector<bool> temp(n+1,0);
        sets.push_back(temp);
    }
    for (tuple<int,int,int>& t : ss){
        sets.at(get<1>(t)).at(get<0>(t)) = 1;
    }

    vector<Set> set;
    for (vector<bool> s : sets){
        set.push_back(Set(s));
    }

    test = SetSystem(pts,set);

    m = test.sets.size();
    test.buildAdjacency(false);
    auto start_time = chrono::high_resolution_clock::now();
    Result res;
    res = partition_distance_set_weight_par(test,t,sw_weighted_w_sample);
    writeCSVFile(res, to_string(time(NULL)) + "_parallel" + "_set_weight_sample_batch.csv");
    auto end_time = chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end_time - start_time;
    printf("Time taken: %.2fs\n", duration.count());
    auto time = duration.count();

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

    ofstream MyFile("results.csv",std::ios_base::app);

    MyFile << "17;" << "book" << ";" << maxcrossing << ";" << avgcrossing << ";" << mincrossing << ";" << time << endl;

    
    return 0;
}