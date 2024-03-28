#include <cmath>
#include <cstdlib>
#include <iterator>
#include <vector>
#include <iostream>

#include "classes.cpp"

bool intersects(int a, int b, Set s){
    return s.points.at(a) != s.points.at(b);
}

float distance(Point a, Point b){
    float res = 0;
    for (int d = 0; d < a.coordinates.size(); d++){
        res += pow((a.coordinates.at(d) - b.coordinates.at(d)),2);
    }
    return sqrt(res);
}

int pick_according_to_weights(vector<float> weights){
    float total = 0;
    for (int i = 0; i < weights.size(); i++){
        total += weights.at(i);
    }
    float stop_at = static_cast<float>(rand())/RAND_MAX*total;
    float partial_sum = 0;
    int i = 0;
    while(i < weights.size() && partial_sum + weights.at(i) < stop_at){
        partial_sum += weights.at(i);
        i++;
    }
    return i < weights.size() ? i : weights.size() - 1;
}

int pick_according_to_weights_exponentially(vector<float> weights, float coeff){
    float total = 0;
    for (int i = 0; i < weights.size(); i++){
        total += (weights.at(i) == 0 ? 0 : pow(coeff,weights.at(i)));
    }
    if (total == INFINITY){
        return pick_according_to_weights_exponentially(weights, (coeff+1)/2);
    }
    float stop_at = static_cast<float>(rand())/RAND_MAX*total;
    float partial_sum = 0.;
    int i = 0;
    while(i < weights.size() && partial_sum + (weights.at(i) == 0 ? 0 : pow(coeff,weights.at(i))) < stop_at){
        partial_sum += (weights.at(i) == 0 ? 0 : pow(coeff,weights.at(i)));
        i++;
    }
    return i < weights.size() ? i : weights.size() - 1;
}


vector<float> l2(vector<Point> pts, vector<bool> available, int start, vector<Set> sets, int k){
    int n = pts.size();
    int d = pts.at(0).coordinates.size();
    vector<float> distances(n,0);
    for (int i = 0; i < n; i++){
        float temp = 0;
        for (int k = 0; k < d; k++){
            temp += abs(pts.at(i).coordinates.at(k) - pts.at(start).coordinates.at(k));
        }
        distances.at(i) = temp;
    }
    return distances;
}

vector<float> l1(vector<Point> pts, vector<bool> available, int start, vector<Set> sets, int k){
    int n = pts.size();
    int d = pts.at(0).coordinates.size();
    vector<float> distances(n,0);
    for (int i = 0; i < n; i++){
        float temp = 0;
        for (int k = 0; k < d; k++){
            temp += sqrt(pow(pts.at(i).coordinates.at(k) - pts.at(start).coordinates.at(k), 2.0));
        }
        distances.at(i) = temp;
    }
    return distances;
}

vector<float> dw(vector<Point> pts, vector<bool> available, int start, vector<Set> sets, int k){
    int n = pts.size();
    int d = pts.at(0).coordinates.size();
    int m = sets.size();
    vector<float> distances(n,1);
    vector<float> set_weight(m,log(n*m)+1);
    for (int i = 0; i < n; i++){
        if (!available.at(i)){
            distances.at(i) = 0;
        }
    }

    for (int t = 0; t < k; t++){
        int s = pick_according_to_weights_exponentially(set_weight, 2);
        int p = pick_according_to_weights_exponentially(distances, 2);

        if (sets.at(s).points.at(start)){
            for (int& i : sets.at(s).complement_indices){
                if (available.at(i)){
                    distances.at(i) += 1;
                }
            }
        } else {
            for (int& i : sets.at(s).points_indices){
                if (available.at(i)){
                    distances.at(i) += 1;
                }
            }
        }

        for (int j = 0 ; j < m; j++){
            if (intersects(start, p, sets.at(j))){
                set_weight.at(j) -= 1;
            }
        }
    }

    return distances;
}

vector<float> sw(vector<Point> pts, vector<bool> available, int start, vector<Set> sets, int k){
    int n = pts.size();
    int d = pts.at(0).coordinates.size();
    int m = sets.size();
    vector<float> distances(n,1);
    for (int i = 0; i < n; i++){
        if (!available.at(i)){
            distances.at(i) = 0;
        }
    }
    for (int t = 0; t < k; t++){
        int s = floor(m*static_cast<float>(rand())/RAND_MAX);

        if (sets.at(s).points.at(start)){
            for (int& i : sets.at(s).complement_indices){
                if (available.at(i)){
                    distances.at(i) += 1;
                }
            }
        } else {
            for (int& i : sets.at(s).points_indices){
                if (available.at(i)){
                    distances.at(i) += 1;
                }
            }
        }
    }
    
    return distances;
}

vector<float> sw_weighted(vector<Point> pts, vector<bool> available, int start, vector<Set> sets, int k){
    int n = pts.size();
    int d = pts.at(0).coordinates.size();
    int m = sets.size();
    vector<float> distances(n,1);
    for (int i = 0; i < n; i++){
        if (!available.at(i)){
            distances.at(i) = 0;
        }
    }
    for (int t = 0; t < k; t++){
        int s = floor(m*static_cast<float>(rand())/RAND_MAX);
        if (sets.at(s).points.at(start)){
            for (int& i : sets.at(s).complement_indices){
                if (available.at(i)){
                    distances.at(i) += 1 << sets.at(s).weight;
                }
            }
        } else {
            for (int& i : sets.at(s).points_indices){
                if (available.at(i)){
                    distances.at(i) += 1 << sets.at(s).weight;
                }
            }
        }
    }
    
    return distances;
}

vector<float> sw_weighted_w_sample(vector<Point> pts, vector<bool> available, int start, vector<Set> sets, int k){
    int n = pts.size();
    int d = pts.at(0).coordinates.size();
    int m = sets.size();
    vector<float> distances(n,1);
    for (int i = 0; i < n; i++){
        if (!available.at(i)){
            distances.at(i) = 0;
        }
    }
    vector<float> set_weight;
    for (int j = 0; j < m; j++){
        set_weight.push_back(sets.at(j).weight);
    }
    for (int t = 0; t < k; t++){
        int s = pick_according_to_weights_exponentially(set_weight, 2);
        if (sets.at(s).points.at(start)){
            for (int& i : sets.at(s).complement_indices){
                if (available.at(i)){
                    distances.at(i) += 1 << sets.at(s).weight;
                }
            }
        } else {
            for (int& i : sets.at(s).points_indices){
                if (available.at(i)){
                    distances.at(i) += 1 << sets.at(s).weight;
                }
            }
        }
    }
    
    return distances;
}
