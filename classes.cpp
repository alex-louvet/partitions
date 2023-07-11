#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

class Point {
    public:
        vector<float> coordinates;

        Point(int dimension){
            for (int i = 0; i < dimension; i++){
                coordinates.push_back(static_cast <float> (rand()) / static_cast <float> (RAND_MAX));
            }
        }
};

class Edge {
    public:
        int points[2];
        int weight;
        int next_ite_weight;

        Edge(int index1, int index2){
            points[0] = index1;
            points[1] = index2;
            weight = 0;
            next_ite_weight = 0;
        }

        void increase(){
            weight ++;
        }
        void decrease(){
            weight --;
        }
        void resetWeight(){
            weight = 0;
        }
};

class Set {
    public:
        vector<bool> points;
        int weight;


        Set(vector<bool> indices){
            points = indices;
            weight = 0;
        }
        Set(int size){
            for (int i = 0; i < size; i++){
                points.push_back(0);
            }
            weight = 0;
        }

        void increase(){
            weight ++;
        }
        void decrease(){
            weight --;
        }
        void resetWeight(){
            weight = 0;
        }
};

class SetSystem {
    public:
        vector<Point> points;
        vector<Set> sets;

        SetSystem(int n,int d,int m){
            for (int i = 0; i < n; i++){
                points.push_back(Point(d));
            }
            for (int i = 0; i < n; i++){
                sets.push_back(Set(n));
            }
        }

        SetSystem(vector<Point> p,vector<Set> s){
            points = p;
            sets = s;
        }

        SetSystem(){
            vector<Point> p;
            vector<Set> s;
            points = p;
            sets = s;
        }
};

Set RandomSet(int size){
    vector<bool> points;
    for (int i = 0; i < size; i++){
        points.push_back(rand() > 0.5*RAND_MAX);
    }
    return Set(points);
}

SetSystem Grid(int n, int d){
    vector<Point> p;
    vector<Set> s;

    for (int i = 0; i < n; i++){
        p.push_back(Point(d));
    }

    vector<bool> temp;
    for (int dim = 0; dim < d; dim++){
        for (int i = 0; i < sqrt(n); i++){
            temp.clear();
            for (const Point& pt : p) {
                if (pt.coordinates[dim] > i/sqrt(n)){
                    temp.push_back(1);
                } else {
                    temp.push_back(0);
                }
            }
            s.push_back(Set(temp));
        }
    }

    return SetSystem(p,s);
}

//TODO maybe
SetSystem RandomHyperplanes(int n, int d, int m){
    vector<Point> p;
    vector<Set> s;
    for (int i = 0; i < n; i++){
        p.push_back(Point(d));
    }

    return SetSystem(p,s);     
}

class Result {
    
    public:
        SetSystem ss;
        vector<int> histogram;
        vector<float> weight;

        Result(SetSystem setsys, vector<int> histo, vector<float> weights){
            ss = setsys;
            histogram = histo;
            weight = weights;
        }
        Result(){}
};
