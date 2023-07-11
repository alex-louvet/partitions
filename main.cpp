#include <iostream>
#include <fstream>
#include <string>

#include "partition.cpp"

using namespace std;

void writePythonFile(Result res, string fileName){
    ofstream MyFile(fileName);
    for (Point& p : res.ss.points) {
        for (const float& f : p.coordinates) {
            MyFile << f << ",";
        }
        MyFile << "\n";
    }
    MyFile << "sets\n";
    for (Set& s : res.ss.sets) {
        for (const bool& b : s.points) {
            MyFile << b << ",";
        }
        MyFile << "\n";
    }
    MyFile << "histogram\n";
    for (const int& x : res.histogram) {
        MyFile << x << ",";
    }
    MyFile << "\n";
    MyFile << "weight\n";
    for (const float& x : res.weight) {
        MyFile << x << ",";
    }
    MyFile << "\n";
    
}

int main(int argc, char** argv){
    int n = 1024;
    int d = 2;
    int m = d*sqrt(n);
    int t = 16;
    if (argc >= 3){
        n = stoi(argv[2]);

    }
    if (argc >= 4){
        t = stoi(argv[3]);

    }
    if (argc >= 5){
        d = stoi(argv[4]);

    }
    if (argc >= 6){
        m = stoi(argv[5]);

    }
    srand(time(NULL));
    SetSystem test = Grid(n,d);
    if (argc >= 2){
        if (stoi(argv[1]) == 1){
            Result res = mwu_min(test,t);
            writePythonFile(res,to_string(time(NULL)) + " min");
        }
        if (stoi(argv[1]) == 2){
            Result res = mwu_random_uniform(test,t);
            writePythonFile(res,to_string(time(NULL)) + " uniform");
        }
        if (stoi(argv[1]) == 3){
            Result res = mwu_random_weight(test,t);
            writePythonFile(res,to_string(time(NULL)) + " weight");
        }

        if (stoi(argv[1]) == 4){
            Result res = mwu_bfs(test,t);
            writePythonFile(res,to_string(time(NULL)) + " bfs");
        }

        if (stoi(argv[1]) == 5){
            Result res = mwu_bfs_min(test,t);
            writePythonFile(res,to_string(time(NULL)) + " bfs_min");
        }

        if (stoi(argv[1]) == 6){
            Result res = mwu_min_bfs(test,t);
            writePythonFile(res,to_string(time(NULL)) + " min_bfs");
        }

    } else {
        Result res = mwu_min(test,t);
        writePythonFile(res,to_string(time(NULL)) + " min");
    } 
    return 0;
}