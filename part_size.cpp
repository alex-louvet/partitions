#include <vector>
#include <cmath>

vector<int> equi_distro(int n, int t){
    vector<int> res;
    int c = n;
    while (c >= n/t){
        res.push_back(n/t);
        c -= n/t;
    }
    return res;
}

vector<int> reduce_size(int n, int t, int min){
    vector<int> res;
    int c = n;
    int floor_current = n/2;
    int current = n/t;
    while (current >= min){
        while (c > floor_current){
            res.push_back(current);
            c -= current;
        }
        floor_current /= 2;
        current /= 2;
    }
    while (c >= current){
        res.push_back(current);
        c -= current;
    }

    return res;
}

vector<int> interval_size(int n, int t, int min, int max){
    vector<int> res;
    int c = n;
    int floor_current = n/2;
    int current = max;
    while (current >= min){
        while (c > floor_current){
            res.push_back(current);
            c -= current;
        }
        floor_current /= 2;
        current /= 2;
    }
    while (c >= current){
        res.push_back(current);
        c -= current;
    }

    return res;
}

vector<int> same_number_different_size_3(int n, int t){
    vector<int> res;
    for (int i = 0; i < t/3 ; i++){
        res.push_back(3*n/(2*t));
    }
    for (int i = 0; i < t/3 + t%3 ; i++){
        res.push_back(n/t);
    }
    for (int i = 0; i < t/3 ; i++){
        res.push_back(n/(2*t));
    }
    return res;
}

vector<int> same_number_different_size_linear(int n, int t, float mineps){
    vector<int> res;
    float incr = 2*mineps/t;
    float current = mineps;
    for (int i = 0 ; i < t/2; i++){
        res.push_back((1+current)*n/t);
        current -= incr;
    }
    if (t%2 == 1){
        res.push_back(n/t);
    }
    for (int i = t/2 - 1 ; i >= 0; i--){
        res.push_back(2*n/t - res.at(i));
    }
    return res;
}

vector<int> ninty_percent(int n, int t, int min){
    vector<int> res;
    int c = n;
    int current = n/(2*t);
    while (c >= 0.1*n){
        res.push_back(n/t);
        c -= n/t;
    }
    while (current > min){
        int atStart = c;
        while (c >= atStart/2){
            res.push_back(current);
            c -= current;
        }
        current /= 2;
    }
    while (c >= current){
        res.push_back(current);
        c -= current;
    }
    return res;
}