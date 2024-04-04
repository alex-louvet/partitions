#include <cstddef>
#include <cstdlib>
#include <tuple>
#include <vector>
#include <cmath>
#include <math.h>
#include <random>

#include "Eigen/Dense"

using namespace std;

class Point {
    public:
        vector<float> coordinates;
        vector<int> belongs_to;
        vector<int> not_in;

        Point(int dimension){
            for (int i = 0; i < dimension; i++){
                coordinates.push_back(static_cast <float> (rand()) / static_cast <float> (RAND_MAX));
            }
        }

        Point(vector<float> coords){
            coordinates = coords;
        }
};

struct
    {
        bool operator()(tuple<int,int> a, tuple<int,int> b) const { return get<1>(a) < get<1>(b); }
    }
    indexWeightOrder;

struct
{
    bool operator()(tuple<int,float> a, tuple<int,float> b) const { return get<1>(a) < get<1>(b); }
}
floatWeightOrder;

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
        vector<int> points_indices;
        vector<int> complement_indices;


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
        void buildAdjacency(){
            points_indices.clear();
            complement_indices.clear();
            for (int i = 0; i < points.size(); i++){
                if (points.at(i)){
                    points_indices.push_back(i);
                } else {
                    complement_indices.push_back(i);
                }
            }
        }
};

struct
    {
        bool operator()(Set a, Set b) const { return a.weight < b.weight; }
    }
    setWeightOrder;

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
        void buildAdjacency(bool pts){
            for (Point& p : points){
                p.belongs_to.clear();
            }
            #pragma omp parallel for
            for (int j = 0; j < sets.size(); j++){
                sets.at(j).buildAdjacency();
                if (pts){
                    for (int i = 0 ; i < sets.at(j).points.size(); i++){
                        if (sets.at(j).points.at(i)){
                            points.at(i).belongs_to.push_back(j);
                        } else{
                            points.at(i).not_in.push_back(j);
                        }
                    }
                }
            }
        }
        
        void resetWeight(){
            for (Set& s : sets){
                s.resetWeight();
            }
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
        for (int i = 0; i < pow(n,1.0/d); i++){
            temp.clear();
            for (const Point& pt : p) {
                if (pt.coordinates[dim] > i/(pow(n,1.0/d))){
                    temp.push_back(1);
                } else {
                    temp.push_back(0);
                }
            }
            s.push_back(Set(temp));
        }
    }

		random_device rd;
		mt19937 g(rd());
		for (int i = 0; i < 10*s.size(); i++){
			shuffle(s.begin(), s.end(), g);
		}

    return SetSystem(p,s);
}

float dot(vector<float> p, vector<float> q)	//dot product of two vectors
{
	assert(p.size() == q.size());
	float sum = 0.0f;
	for (size_t i = 0; i < p.size(); i++) sum += p[i] * q[i];
	return sum;
}

SetSystem RandomHyperplanes(int n, int d, int m){
    vector<Point> p;
    vector<Set> s;
    vector<int> sample;
    for (int i = 0; i < n; i++){
        p.push_back(Point(d));
    }
    for (int z = 0; z < m ; z++){
        sample.clear();
        while (sample.size() < d){
            bool push = true;
            int temp = rand()%n;
            for (const int& c : sample) {
                if (c == temp){
                    push = false;
                }
            }
            if (push){
                sample.push_back(rand()%n);
            }
        }
        Eigen::MatrixXd A(d, d);
		for (int i = 0; i < d; i++)
			for (int j = 0; j < d; j++)
				A(i, j) = p.at(sample.at(i)).coordinates.at(j);
		Eigen::VectorXd b(d);
		for (size_t i = 0; i < d; i++) b(i) = 1.0f;
		//VectorXd ans = A.colPivHouseholderQr().solve(b);
		Eigen::VectorXd ans = A.partialPivLu().solve(b);
        vector<float> N;
        N.reserve(d);
		for (size_t i = 0; i < d; i++) N.push_back(ans(i));
        Set set = Set(n);
        for (int i = 0; i < n ; i++){
            set.points.at(i) = (dot(p.at(i).coordinates,N) > 1.0) ? 1 : 0;
        }
        for (const int& index : sample) {
            set.points.at(index) = 1;
        }
        s.push_back(set);
    }

    return SetSystem(p,s);
}

SetSystem GridGraph(int n, int d){
    vector<Point> pts;
    for (int i = 0; i < n ; i++){
        for (int j = 0; j < n ; j++){
            pts.push_back(Point(vector<float> {static_cast<float>(i)/n,static_cast<float>(j)/n}));
        }
    }

    const int pts_len = pts.size();  

    vector<Set> sets;
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            Set temp = Set(vector<bool>(pts_len,0));
            temp.points.at(i*n+j) = 1;
            for (int d1 = 0; d1 <= d ; d1++){
                for (int d2 = 0; d2 <= d1; d2++){
                    int d3 = d1-d2;
                    if (j >= d2){
                        if (i >= d3){
                            temp.points.at(i*n+j-d2-n*d3) = 1;
                        }
                        if (i < n-d3){
                            temp.points.at(i*n+j-d2+n*d3) = 1;
                        }
                            
                    }
                    if (j < n-d2){
                        if (i >= d3){
                            temp.points.at(i*n+j+d2-n*d3) = 1;
                        }
                        if (i < n-d3){
                            temp.points.at(i*n+j+d2+n*d3) = 1;
                        }
                            
                    }
                } 
            }
            sets.push_back(temp);
        }
    }
    
    return SetSystem(pts,sets);
}

SetSystem LinearGrid(int n, int d){
    vector<Point> p;
    vector<Set> s;

    for (int i = 0; i < n; i++){
        p.push_back(Point(d));
    }

    vector<bool> temp;
    for (int dim = 0; dim < d; dim++){
        for (int i = 0; i < pow(n,1.0/d); i++){
            temp.clear();
            for (const Point& pt : p) {
                if (pt.coordinates[dim] > i/(pow(n,1.0/d))){
                    temp.push_back(1);
                } else {
                    temp.push_back(0);
                }
            }
            for (int j = 0; j < i; j++){
                s.push_back(Set(temp));
            }
        }
    }

		random_device rd;
		mt19937 g(rd());
		for (int i = 0; i < 10*s.size(); i++){
			shuffle(s.begin(), s.end(), g);
		}

    return SetSystem(p,s);
}

SetSystem ExponentialGrid(int n, int d){
    vector<Point> p;
    vector<Set> s;

    for (int i = 0; i < n; i++){
        p.push_back(Point(d));
    }

    vector<bool> temp;
    for (int dim = 0; dim < d; dim++){
        for (int i = 0; i < pow(n,1.0/d); i++){
            temp.clear();
            for (const Point& pt : p) {
                if (pt.coordinates[dim] > i/(pow(n,1.0/d))){
                    temp.push_back(1);
                } else {
                    temp.push_back(0);
                }
            }
            for (int j = 0; j < pow(2,i); j++){
                s.push_back(Set(temp));
            }
        }
    }

		random_device rd;
		mt19937 g(rd());
		for (int i = 0; i < 10*s.size(); i++){
			shuffle(s.begin(), s.end(), g);
		}

    return SetSystem(p,s);
}

SetSystem DirectionalGrid(int n, int d){
    vector<Point> p;
    vector<Set> s;

    for (int i = 0; i < n; i++){
        p.push_back(Point(d));
    }

    vector<bool> temp;
    for (int dim = 0; dim < d; dim++){
        for (int i = 0; i < pow(n,1.0/d); i++){
            temp.clear();
            for (const Point& pt : p) {
                if (pt.coordinates[dim] > i/(pow(n,1.0/d))){
                    temp.push_back(1);
                } else {
                    temp.push_back(0);
                }
            }
            if (dim == 0){
                for (int j = 0; j < pow(n,1.0/d); j++){
                    s.push_back(Set(temp));
                }
            } else {
                s.push_back(Set(temp));
            }
        }
    }

		random_device rd;
		mt19937 g(rd());
		for (int i = 0; i < 10*s.size(); i++){
			shuffle(s.begin(), s.end(), g);
		}

    return SetSystem(p,s);
}

/*
Credits to Salmelu for its code that I just adapted to my data model:
https://github.com/Salmelu/ProjectivePlane
*/
SetSystem ProjectivePlane(int o){
    vector<Point> p;
    vector<Set> s;

    p.push_back(Point({1.,0.,0.}));

    for (float i = 0; i < o; i++){
        p.push_back(Point({i,1.,0.}));
    }

    for (float i = 0; i < o; i++){
        for (float j = 0; j < o; j++){
            p.push_back(Point({i,j,1.}));
        }
    }

    vector<bool> unique(o*o*o, true);
    unique.at(0) = false;

    for (int a = 0; a < o; a++) {
		for (int b = 0; b < o; b++) {
			for (int c = 0; c < o; c++) {
                if (!unique[(a * o * o) + (b * o) + c]) continue;
                vector<bool> temp;
                temp.clear();
                for (int i = 0; i < p.size(); i++){
                    if ((static_cast<int>(a*p.at(i).coordinates.at(0) + b*p.at(i).coordinates.at(1) + c*p.at(i).coordinates.at(2)) % o) == 0) {
                        temp.push_back(1);
                    } else{
                        temp.push_back(0);
                    }
                }
                for (int i = 2; i < o; i++) {			
					unique.at(((i * a) % o) * o * o + ((i * b) % o) * o + (i * c) % o) = false;
				}
                s.push_back(temp);
            }
        }
    }

    return SetSystem(p,s);
}

class Result: public SetSystem {
    public:
        vector<float> weights;
        vector<int> intersections;
        Result(int n,int d,int m){
            for (int i = 0; i < n; i++){
                points.push_back(Point(d));
            }
            for (int i = 0; i < n; i++){
                sets.push_back(Set(n));
            }
        }

        Result(vector<Point> p,vector<Set> s){
            points = p;
            sets = s;
        }

        Result(){
            vector<Point> p;
            vector<Set> s;
            points = p;
            sets = s;
        }

        Result(SetSystem ss){
            points = ss.points;
            sets = ss.sets;
        }
};
