#include <algorithm>
//#include <bits/fs_fwd.h>
//#include <bits/types/clock_t.h>
#include <cstdlib>
#include <math.h>
#include <tuple>
#include <vector>
#include <cmath>
#include <random>
#include <iostream>
#include <deque>
#include <list>

#include "Eigen/Core"
#include "Eigen/src/Core/util/Constants.h"
#include "distances.cpp"

using namespace std;

bool intersects(Edge e, Set s){
    return s.points.at(e.points[0]) != s.points.at(e.points[1]);
}

int sumSet(Set s){
    int res = 0;
    for (const bool& value : s.points) {
        res += value;
    }
    return res;
}

int posNum(vector<unsigned long> weight){
    int res = 0;
    for (const float& value : weight) {
        if (value >= 0){
            res ++;
        }
    }
    return res;
}

Result partition_min_stats(SetSystem ss, int t){
    
    const int n = ss.points.size();
    const int d = ss.points.at(0).coordinates.size();
    const int m = ss.sets.size();

    //Store partition weight at each iteration
    vector<unsigned long> weights;

    //Store the number of partitions intersected
    vector<int> intersections(m,0);

    //Initialize weight vector
    vector<unsigned long> weight(n,0);
    vector<bool> available(n,true);

    //reset sets weight
    for (Set& s : ss.sets) {
        s.resetWeight();
    }

    vector<Set> s;
    Result res = Result(ss.points, s);
    vector<bool> intersect_partition(m,0);

    for (int i = 0 ; i < t; i++){
        //cout << "\nPartition " << i+1 << "\n";

        Set partition = Set(n);

        unsigned long setsWeight = 0;
        for (Set& s : ss.sets) {
            setsWeight += 1 << s.weight;
        }
        unsigned long partitionWeight = 0;


        vector<bool> intersect_partition(m,0);

        vector<int> admissible_start;
        for (int k = 0; k < n; k++){
            if (available.at(k)){
                admissible_start.push_back(k);
                weight.at(k) = 0;
            }
        }
        int start = admissible_start.at(rand()%admissible_start.size());
        partition.points.at(start) = 1;
        available.at(start) = false;
        for (Set& s : ss.sets) {
            if (s.points.at(start)){
                for (int& k : s.complement_indices){
                    weight.at(k) += 1 << s.weight;
                }
            } else {
                for (int& k : s.points_indices){
                    weight.at(k) += 1 << s.weight;
                }
            }
        }

        for (int k = 1 ; k < n / t; k++){
            int min = -1;
            for (auto j = 0; j < n; j++) {
                // Replace previous minimum if the weight is strictly smaller in all valid cases (no restriction if the partition has no edge yet, connectiveness otherwise)
                if (available.at(j) && (min == -1 || weight.at(j) < weight.at(min))){
                    min = j;
                } else {
                    // If the algorithm finds a valid edge with the same weight it swaps it with probability 1/2 to avoid obteining the same result for every algo because of order
                    if (available.at(j) && (min == -1 || (weight.at(j) <= weight.at(min) && rand() > 0.5*RAND_MAX))){
                        min = j;
                    }
                }
            }

            //Add selected edge to the partition
            partition.points.at(min) = 1;

            //Update partition weight
            partitionWeight += weight.at(min);
            available.at(min) = false;
            weights.push_back(partitionWeight*pow(static_cast<float>(n)-(static_cast<float>(n)/static_cast<float>(t))*static_cast<float>(i),1.0/d)/setsWeight);

            // Remove from the weight of edges the weight of sets that intersect the selected edge
            for (int j = 0; j < ss.sets.size(); j++) {
                if (!intersect_partition.at(j) && intersects(Edge(start,min),ss.sets.at(j))){
                    intersect_partition.at(j) = 1;
                    if (ss.sets.at(j).points.at(start)){
                        for (int& pt : ss.sets.at(j).complement_indices) {
                            if (available.at(pt) && intersects(Edge(start,pt),ss.sets.at(j))){
                                weight.at(pt) -= (1 << ss.sets.at(j).weight);
                            }
                        }
                    } else {
                        for (int& pt : ss.sets.at(j).points_indices) {
                            if (available.at(pt) && intersects(Edge(start,pt),ss.sets.at(j))){
                                weight.at(pt) -= (1 << ss.sets.at(j).weight);
                            }
                        }
                    }
                    // Double set weight for next iteration
                    ss.sets.at(j).increase();
                }
            }
        }

        // Store results
        res.sets.push_back(partition);
        for (int j = 0; j < m; j++){
            intersections.at(j) += intersect_partition.at(j);
        }
    }

    //fill the last partition with the remaining points
    for (int i = 0; i < n ; i++){
        if (available.at(i)){
            res.sets.at(res.sets.size()-1).points.at(i) = 1;
        }
    }
    res.intersections = intersections;
    for (unsigned long& w : weights){
        res.weights.push_back(static_cast<float>(w));
    }
    return res;
}

Result partition_min_rate_stats(SetSystem ss, int t){
    
    const int n = ss.points.size();
    const int d = ss.points.at(0).coordinates.size();
    const int m = ss.sets.size();

    vector<unsigned long> weights;
    vector<int> intersections(m,0);

    //Initialize weight vector
    vector<unsigned long> weight(n,0);
    vector<bool> available(n,true);

    //reset sets weight
    for (Set& s : ss.sets) {
        s.resetWeight();
    }

    vector<Set> s;
    Result res = Result(ss.points, s);
    vector<bool> intersect_partition(m,0);

    for (int i = 0 ; i < t; i++){
        //cout << "\nPartition " << i+1 << "\n";

        Set partition = Set(n);

        unsigned long setsWeight = 0;
        for (Set& s : ss.sets) {
            setsWeight += 1 << s.weight;
        }
        unsigned long partitionWeight = 0;


        vector<bool> intersect_partition(m,0);

        vector<int> admissible_start;
        for (int k = 0; k < n; k++){
            if (weight.at(k) >= 0){
                admissible_start.push_back(k);
                weight.at(k) = 0;
            }
        }
        int start = admissible_start.at(rand()%admissible_start.size());
        available.at(start) = false;
        partition.points.at(start) = 1;
        for (Set& s : ss.sets) {
            if (s.points.at(start)){
                for (int& k : s.complement_indices){
                    weight.at(k) += 1 << s.weight;
                }
            } else {
                for (int& k : s.points_indices){
                    weight.at(k) += 1 << s.weight;
                }
            }
        }

        for (int k = 1 ; k < n / t; k++){
            int min = -1;
            bool test = false;
            int j = 0;
            while (!test && j < n) {
                if (available.at(j) && (partitionWeight + weight.at(j))*pow(static_cast<float>(n)-(static_cast<float>(n)/static_cast<float>(t))*static_cast<float>(i),1.0/d)/setsWeight <= 1*pow(static_cast<float>(k),1/d)){
                    min = j;
                    test = true;
                } else {
                    // Replace previous minimum if the weight is strictly smaller in all valid cases (no restriction if the partition has no edge yet, connectiveness otherwise)
                    if (available.at(j) && (min == -1 || weight.at(j) < weight.at(min))){
                        min = j;
                    } else {
                        // If the algorithm finds a valid edge with the same weight it swaps it with probability 1/2 to avoid obteining the same result for every algo because of order
                        if (available.at(j) && (min == -1 || weight.at(j) <= weight.at(min) && rand() > 0.5*RAND_MAX)){
                            min = j;
                        }
                    }
                }
                j++;
            }
            //Add selected edge to the partition
            partition.points.at(min) = 1;

            //Update partition weight
            partitionWeight += weight.at(min);
            available.at(min) = false;
            weights.push_back(partitionWeight*pow(n-(n/t)*i,1/d)/setsWeight);

            // Remove from the weight of edges the weight of sets that intersect the selected edge
            for (int j = 0; j < ss.sets.size(); j++) {
                if (!intersect_partition.at(j) && intersects(Edge(start,min),ss.sets.at(j))){
                    intersect_partition.at(j) = 1;
                    if (ss.sets.at(j).points.at(start)){
                        for (int& pt : ss.sets.at(j).complement_indices) {
                            if (available.at(pt) && intersects(Edge(start,pt),ss.sets.at(j))){
                                weight.at(pt) -= (1 << ss.sets.at(j).weight);
                            }
                        }
                    } else {
                        for (int& pt : ss.sets.at(j).points_indices) {
                            if (available.at(pt) && intersects(Edge(start,pt),ss.sets.at(j))){
                                weight.at(pt) -= (1 << ss.sets.at(j).weight);
                            }
                        }
                    }
                    // Double set weight for next iteration
                    ss.sets.at(j).increase();
                }
            }
        }

        // Store results
        res.sets.push_back(partition);
        for (int j = 0; j < m; j++){
            intersections.at(j) += intersect_partition.at(j);
        }
    }

    for (int i = 0; i < n ; i++){
        if (available.at(i)){
            res.sets.at(res.sets.size()-1).points.at(i) = 1;
        }
    }
    res.intersections = intersections;
    for (unsigned long& w : weights){
        res.weights.push_back(static_cast<float>(w));
    }    
    return res;
}

Result partition_rate_stats(SetSystem ss, int t){
    
    const int n = ss.points.size();
    const int d = ss.points.at(0).coordinates.size();
    const int m = ss.sets.size();

    vector<unsigned long> weights;
    vector<int> intersections(m,0);

    //Initialize weight vector
    vector<unsigned long> weight(n,0);
    vector<bool> available(n,true);

    //reset sets weight
    for (Set& s : ss.sets) {
        s.resetWeight();
    }

    vector<Set> s;
    Result res = Result(ss.points, s);
    vector<bool> intersect_partition(m,0);

    for (int i = 0 ; i < t; i++){
        //cout << "\nPartition " << i+1 << "\n";

        Set partition = Set(n);

        unsigned long setsWeight = 0;
        for (Set& s : ss.sets) {
            setsWeight += 1 << s.weight;
        }
        unsigned long partitionWeight = 0.;


        vector<bool> intersect_partition(m,0);

        vector<int> admissible_start;
        for (int k = 0; k < n; k++){
            if (available.at(k)){
                admissible_start.push_back(k);
                weight.at(k) = 0;
            }
        }
        int start = admissible_start.at(rand()%admissible_start.size());
        available.at(start) = false;
        partition.points.at(start) = 1;
        for (Set& s : ss.sets) {
            if (s.points.at(start)){
                for (int& k : s.complement_indices){
                    weight.at(k) += 1 << s.weight;
                }
            } else {
                for (int& k : s.points_indices){
                    weight.at(k) += 1 << s.weight;
                }
            }
        }

        for (int k = 1 ; k < n / t; k++){
            int min = -1;
            int j = 0;
            vector<int> candidates;
            candidates.clear();
            while (j < n) {
                if (available.at(j) && (partitionWeight + weight.at(j))*pow(static_cast<float>(n)-(static_cast<float>(n)/static_cast<float>(t))*static_cast<float>(i),1.0/d)/setsWeight <= 2*pow(static_cast<float>(k),1.0/d)){
                    candidates.push_back(j);
                    break;
                } else {
                    // Replace previous minimum if the weight is strictly smaller in all valid cases (no restriction if the partition has no edge yet, connectiveness otherwise)
                    if (available.at(j) && (min == -1 || weight.at(j) < weight.at(min))){
                        min = j;
                    } else {
                        // If the algorithm finds a valid edge with the same weight it swaps it with probability 1/2 to avoid obteining the same result for every algo because of order
                        if (available.at(j) && (min == -1 || weight.at(j) <= weight.at(min) && rand() > 0.5*RAND_MAX)){
                            min = j;
                        }
                    }
                }
                j++;
            }
            if (candidates.size() > 0){
                min = candidates.at(rand()%candidates.size());
            }

            //Add selected edge to the partition
            partition.points.at(min) = 1;

            //Update partition weight
            partitionWeight += weight.at(min);
            available.at(min) = false;
            weights.push_back(partitionWeight*pow(static_cast<float>(n)-(static_cast<float>(n)/static_cast<float>(t))*static_cast<float>(i),1.0/d)/setsWeight);

            // Remove from the weight of edges the weight of sets that intersect the selected edge
            for (int j = 0; j < ss.sets.size(); j++) {
                if (!intersect_partition.at(j) && intersects(Edge(start,min),ss.sets.at(j))){
                    intersect_partition.at(j) = 1;
                    if (ss.sets.at(j).points.at(start)){
                        for (int& pt : ss.sets.at(j).complement_indices) {
                            if (available.at(pt) && intersects(Edge(start,pt),ss.sets.at(j))){
                                weight.at(pt) -= (1 << ss.sets.at(j).weight);
                            }
                        }
                    } else {
                        for (int& pt : ss.sets.at(j).points_indices) {
                            if (available.at(pt) && intersects(Edge(start,pt),ss.sets.at(j))){
                                weight.at(pt) -= (1 << ss.sets.at(j).weight);
                            }
                        }
                    }
                    
                    // Double set weight for next iteration
                    ss.sets.at(j).increase();
                }
            }
        }

        // Store results
        res.sets.push_back(partition);
        for (int j = 0; j < m; j++){
            intersections.at(j) += intersect_partition.at(j);
        }
                
    }

    for (int i = 0; i < n ; i++){
        if (available.at(i)){
            res.sets.at(res.sets.size()-1).points.at(i) = 1;
        }
    }

    res.intersections = intersections;
    for (unsigned long& w : weights){
        res.weights.push_back(static_cast<float>(w));
    }
    return res;
}


Result partition_sampling_stats(SetSystem ss, int t, float p){
    
    const int n = ss.points.size();
    const int d = ss.points.at(0).coordinates.size();
    const int m = ss.sets.size();

    vector<unsigned long> weights;
    vector<int> intersections(m,0);

    //Initialize weight vector
    vector<unsigned long> weight(n,0);
    vector<bool> available(n,true);

    //reset sets weight
    for (Set& s : ss.sets) {
        s.resetWeight();
    }

    vector<Set> s;
    Result res = Result(ss.points, s);
    vector<bool> intersect_partition(m,0);

    for (int i = 0 ; i < t; i++){
        //cout << "\nPartition " << i+1 << "\n";

        Set partition = Set(n);

        unsigned long setsWeight = 0;
        for (Set& s : ss.sets) {
            setsWeight += 1 << s.weight;
        }
        unsigned long partitionWeight = 0.;


        vector<bool> intersect_partition(m,0);

        vector<int> admissible_start;
        for (int k = 0; k < n; k++){
            if (available.at(k)){
                admissible_start.push_back(k);
                weight.at(k) = 0;
            }
        }
        int start = admissible_start.at(rand()%admissible_start.size());
        available.at(start) = false;
        partition.points.at(start) = 1;
        /*for (Set& s : ss.sets) {
            for (int k = 0; k < n; k++){
                if (available.at(k) && intersects(Edge(start,k), s)){
                    weight.at(k) += 1 << s.weight;
                }
            }
        }*/

        for (int k = 1 ; k < n / t; k++){
            int min = -1;
            vector<int> sampled;
            for (int j = 0; j < m; j++){
                if (rand() < p*RAND_MAX){
                    sampled.push_back(rand()%m);
                    
                }
            }

            for (int j = 0; j < n; j++){
                if (available.at(j)){
                    for (int& samp : sampled) {
                        if (intersects(Edge(start,j), ss.sets.at(samp))){
                            if (!intersect_partition.at(samp)){
                                weight.at(j) += 1 << ss.sets.at(samp).weight;
                            }
                        }
                    }

                    if (1/p*(partitionWeight + weight.at(j))/static_cast<float>(setsWeight) <= 0*pow(static_cast<float>(k),1.0/d)/pow(static_cast<float>(n-n/static_cast<float>(t)*i),1.0/d)){
                        min = j;
                        break;
                    } else {
                        // Replace previous minimum if the weight is strictly smaller in all valid cases (no restriction if the partition has no edge yet, connectiveness otherwise)
                        if (min == -1 || weight.at(j) < weight.at(min)){
                            min = j;
                        } else {
                            // If the algorithm finds a valid edge with the same weight it swaps it with probability 1/2 to avoid obteining the same result for every algo because of order
                            if (min == -1 || weight.at(j) <= weight.at(min) && rand() > 0.5*RAND_MAX){
                                min = j;
                            }
                        }
                    }
                }
            }
            //Add selected edge to the partition
            partition.points.at(min) = 1;

            //Update partition weight
            partitionWeight += weight.at(min);
            available.at(min) = false;
            weights.push_back(partitionWeight*pow(static_cast<float>(n)-(static_cast<float>(n)/static_cast<float>(t))*static_cast<float>(i),1.0/d)/setsWeight);

            // Remove from the weight of edges the weight of sets that intersect the selected edge
            for (int j = 0; j < ss.sets.size(); j++) {
                if (!intersect_partition.at(j) && intersects(Edge(start,min),ss.sets.at(j))){
                    intersect_partition.at(j) = 1;
                    /*for (int pt = 0; pt < n; pt++) {
                        if (available.at(pt) && intersects(Edge(start,pt),ss.sets.at(j))){
                            weight.at(pt) -= (1 << ss.sets.at(j).weight);
                        }
                    }*/
                    // Double set weight for next iteration
                    ss.sets.at(j).increase();
                }
            }

        }
        

        // Store results
        res.sets.push_back(partition);
        for (int j = 0; j < m; j++){
            intersections.at(j) += intersect_partition.at(j);
        }
                
    }

    for (int i = 0; i < n ; i++){
        if (available.at(i)){
            res.sets.at(res.sets.size()-1).points.at(i) = 1;
        }
    }

    res.intersections = intersections;
        for (unsigned long& w : weights){
        res.weights.push_back(static_cast<float>(w));
    }
    return res; 
}

Result partition_sampling_fixed_stats(SetSystem ss, int t, int sample_size){
    
    const int n = ss.points.size();
    const int d = ss.points.at(0).coordinates.size();
    const int m = ss.sets.size();

    vector<unsigned long> weights;
    vector<int> intersections(m,0);

    //Initialize weight vector
    vector<unsigned long> weight(n,0);
    vector<bool> available(n,true);

    //reset sets weight
    for (Set& s : ss.sets) {
        s.resetWeight();
    }

    vector<Set> s;
    Result res = Result(ss.points, s);
    vector<bool> intersect_partition(m,0);

    for (int i = 0 ; i < t; i++){
        //cout << "\nPartition " << i+1 << "\n";

        Set partition = Set(n);

        unsigned long setsWeight = 0;
        for (Set& s : ss.sets) {
            setsWeight += 1 << s.weight;
        }
        unsigned long partitionWeight = 0.;


        vector<bool> intersect_partition(m,0);

        vector<int> admissible_start;
        for (int k = 0; k < n; k++){
            if (available.at(k)){
                admissible_start.push_back(k);
                weight.at(k) = 0;
            }
        }

        int start = admissible_start.at(rand()%admissible_start.size());
        available.at(start) = false;
        partition.points.at(start) = 1;
        /*for (Set& s : ss.sets) {
            for (int k = 0; k < n; k++){
                if (available.at(k) && intersects(Edge(start,k), s)){
                    weight.at(k) += 1 << s.weight;
                }
            }
        }*/

        for (int k = 1 ; k < n / t; k++){
            int min = -1;

            //get sampled_size distinct random number between 0 and m-1
            vector<int> urn;
            for (int j = 0; j < m; j++){
                urn.push_back(j);
            }
            random_device rd;
            mt19937 g(rd());
            shuffle(urn.begin(), urn.end(),g);
            vector<int>::const_iterator first = urn.begin();
            vector<int>::const_iterator last = urn.begin() + sample_size;
            vector<int> sampled(first,last);
            setsWeight=0;
            for (int& samp : sampled) {
                setsWeight += ss.sets.at(samp).weight;
            }
            for (int j = 0; j < n; j++){
                if (available.at(j)){
                    for (int& samp : sampled) {
                        if (intersects(Edge(start,j), ss.sets.at(samp))){
                            if (!intersect_partition.at(samp)){
                                weight.at(j) += 1 << ss.sets.at(samp).weight;
                            }
                        }
                    }

                    if (m/static_cast<float>(sample_size)*(partitionWeight + weight.at(j))/static_cast<float>(setsWeight) <= 1*pow(static_cast<float>(k),1.0/d)/pow(static_cast<float>(n-n/static_cast<float>(t)*i),1.0/d)){
                        min = j;
                        break;
                    } else {
                        // Replace previous minimum if the weight is strictly smaller in all valid cases (no restriction if the partition has no edge yet, connectiveness otherwise)
                        if (min == -1 || weight.at(j) < weight.at(min)){
                            min = j;
                        } else {
                            // If the algorithm finds a valid edge with the same weight it swaps it with probability 1/2 to avoid obteining the same result for every algo because of order
                            if (min == -1 || weight.at(j) <= weight.at(min) && rand() > 0.5*RAND_MAX){
                                min = j;
                            }
                        }
                    }
                }
            }
            //Add selected edge to the partition
            partition.points.at(min) = 1;

            //Update partition weight
            partitionWeight += weight.at(min);
            available.at(min) = false;
            weights.push_back(partitionWeight*pow(static_cast<float>(n)-(static_cast<float>(n)/static_cast<float>(t))*static_cast<float>(i),1.0/d)/setsWeight);

            // Remove from the weight of edges the weight of sets that intersect the selected edge
            for (int j = 0; j < ss.sets.size(); j++) {
                if (!intersect_partition.at(j) && intersects(Edge(start,min),ss.sets.at(j))){
                    intersect_partition.at(j) = 1;
                    /*for (int pt = 0; pt < n; pt++) {
                        if (available.at(pt) && intersects(Edge(start,pt),ss.sets.at(j))){
                            weight.at(pt) -= (1 << ss.sets.at(j).weight);
                        }
                    }*/
                    // Double set weight for next iteration
                    ss.sets.at(j).increase();
                }
            }

        }
        

        // Store results
        res.sets.push_back(partition);
        for (int j = 0; j < m; j++){
            intersections.at(j) += intersect_partition.at(j);
        }
                
    }

    for (int i = 0; i < n ; i++){
        if (available.at(i)){
            res.sets.at(res.sets.size()-1).points.at(i) = 1;
        }
    }

    res.intersections = intersections;
    for (unsigned long& w : weights){
        res.weights.push_back(static_cast<float>(w));
    }
    return res;
}

Result partition_sampling_fixed_stats2(SetSystem ss, int t, int sample_size){
    
    const int n = ss.points.size();
    const int d = ss.points.at(0).coordinates.size();
    const int m = ss.sets.size();

    vector<unsigned long> weights;
    vector<int> intersections(m,0);

    //Initialize weight vector
    vector<unsigned long> weight(n,0);
    vector<bool> available(n,true);

    //reset sets weight
    for (Set& s : ss.sets) {
        s.resetWeight();
    }

    vector<Set> s;
    Result res = Result(ss.points, s);
    vector<bool> intersect_partition(m,0);

    for (int i = 0 ; i < t; i++){
        //cout << "\nPartition " << i+1 << "\n";

        Set partition = Set(n);

        unsigned long setsWeight = 0;
        for (Set& s : ss.sets) {
            setsWeight += 1 << s.weight;
        }
        unsigned long partitionWeight = 0.;


        vector<bool> intersect_partition(m,0);

        vector<int> admissible_start;
        for (int k = 0; k < n; k++){
            if (available.at(k)){
                admissible_start.push_back(k);
                weight.at(k) = 0;
            }
        }

        int start = admissible_start.at(rand()%admissible_start.size());
        available.at(start) = false;
        partition.points.at(start) = 1;
        /*for (Set& s : ss.sets) {
            for (int k = 0; k < n; k++){
                if (available.at(k) && intersects(Edge(start,k), s)){
                    weight.at(k) += 1 << s.weight;
                }
            }
        }*/

        for (int k = 1 ; k < n / t; k++){
            int min = -1;

            //get sampled_size distinct random number between 0 and m-1
            vector<int> urn;
            for (int j = 0; j < m; j++){
                urn.push_back(j);
            }
            random_device rd;
            mt19937 g(rd());
            shuffle(urn.begin(), urn.end(),g);
            vector<int>::const_iterator first = urn.begin();
            vector<int>::const_iterator last = urn.begin() + sample_size;
            vector<int> sampled(first,last);
            setsWeight=0;
            for (int& samp : sampled) {
                if (intersect_partition.at(samp)){
                    setsWeight += ss.sets.at(samp).weight;
                }
            }
            for (int j = 0; j < n; j++){
                if (available.at(j)){
                    for (int& samp : sampled) {
                        if (intersects(Edge(start,j), ss.sets.at(samp))){
                            if (!intersect_partition.at(samp)){
                                weight.at(j) += 1 << ss.sets.at(samp).weight;
                            }
                        }
                    }

                    if (m/static_cast<float>(sample_size)*(partitionWeight + weight.at(j))/static_cast<float>(setsWeight) <= 1*pow(static_cast<float>(k),1.0/d)/pow(static_cast<float>(n-n/static_cast<float>(t)*i),1.0/d)){
                        min = j;
                        break;
                    } else {
                        // Replace previous minimum if the weight is strictly smaller in all valid cases (no restriction if the partition has no edge yet, connectiveness otherwise)
                        if (min == -1 || weight.at(j) < weight.at(min)){
                            min = j;
                        } else {
                            // If the algorithm finds a valid edge with the same weight it swaps it with probability 1/2 to avoid obteining the same result for every algo because of order
                            if (min == -1 || weight.at(j) <= weight.at(min) && rand() > 0.5*RAND_MAX){
                                min = j;
                            }
                        }
                    }
                }
            }
            //Add selected edge to the partition
            partition.points.at(min) = 1;

            //Update partition weight
            partitionWeight += weight.at(min);
            available.at(min) = false;
            weights.push_back(partitionWeight*pow(static_cast<float>(n)-(static_cast<float>(n)/static_cast<float>(t))*static_cast<float>(i),1.0/d)/setsWeight);

            // Remove from the weight of edges the weight of sets that intersect the selected edge
            for (int j = 0; j < ss.sets.size(); j++) {
                if (!intersect_partition.at(j) && intersects(Edge(start,min),ss.sets.at(j))){
                    intersect_partition.at(j) = 1;
                    /*for (int pt = 0; pt < n; pt++) {
                        if (available.at(pt) && intersects(Edge(start,pt),ss.sets.at(j))){
                            weight.at(pt) -= (1 << ss.sets.at(j).weight);
                        }
                    }*/
                    // Double set weight for next iteration
                    ss.sets.at(j).increase();
                }
            }

        }
        

        // Store results
        res.sets.push_back(partition);
        for (int j = 0; j < m; j++){
            intersections.at(j) += intersect_partition.at(j);
        }
                
    }

    for (int i = 0; i < n ; i++){
        if (available.at(i)){
            res.sets.at(res.sets.size()-1).points.at(i) = 1;
        }
    }

    res.intersections = intersections;
    for (unsigned long& w : weights){
        res.weights.push_back(static_cast<float>(w));
    }
    return res;
}

Result no_weight_update_deque_insert_end(SetSystem ss, int t){
    const int n = ss.points.size();
    const int d = ss.points.at(0).coordinates.size();
    const int m = ss.sets.size();

    vector<unsigned long> weights;
    vector<int> intersections(m,0);

    //Initialize weight vector
    vector<unsigned long> weight(n,0);
    vector<bool> available(n,true);

    //reset sets weight
    for (Set& s : ss.sets) {
        s.resetWeight();
    }

    vector<Set> s;
    Result res = Result(ss.points, s);
    vector<bool> intersect_partition(m,0);

    for (int i = 0 ; i < t; i++){
        //cout << "\nPartition " << i+1 << "\n";

        Set partition = Set(n);

        unsigned long setsWeight = 0;
        for (Set& s : ss.sets) {
            setsWeight += 1 << s.weight;
        }
        unsigned long partitionWeight = 0.;


        vector<bool> intersect_partition(m,0);

        vector<int> admissible_start;
        for (int k = 0; k < n; k++){
            if (available.at(k)){
                admissible_start.push_back(k);
                weight.at(k) = 0;
            }
        }
        int start = admissible_start.at(rand()%admissible_start.size());
        available.at(start) = false;
        partition.points.at(start) = 1;
        for (Set& s : ss.sets) {
            if (s.points.at(start)){
                for (int& k : s.complement_indices){
                    weight.at(k) += 1 << s.weight;
                }
            } else {
                for (int& k : s.points_indices){
                    weight.at(k) += 1 << s.weight;
                }
            }
        }

        deque<int> update_weight;
        update_weight.clear();

        for (int k = 1 ; k < n / t; k++){
            int min = -1;
            int j = 0;
            while (j < n) {
                if (available.at(j) && (partitionWeight + weight.at(j))*pow(static_cast<float>(n)-(static_cast<float>(n)/static_cast<float>(t))*static_cast<float>(i),1.0/d)/setsWeight <= 2*pow(static_cast<float>(k),1.0/d)){
                    min=j;
                    break;
                } else {
                    // Replace previous minimum if the weight is strictly smaller in all valid cases (no restriction if the partition has no edge yet, connectiveness otherwise)
                    if (available.at(j) && (min == -1 || weight.at(j) < weight.at(min))){
                        min = j;
                    } else {
                        // If the algorithm finds a valid edge with the same weight it swaps it with probability 1/2 to avoid obteining the same result for every algo because of order
                        if (available.at(j) && (min == -1 || weight.at(j) <= weight.at(min) && rand() > 0.5*RAND_MAX)){
                            min = j;
                        }
                    }
                }
                j++;
            }
            //Add selected edge to the partition
            partition.points.at(min) = 1;

            //Update partition weight
            partitionWeight += weight.at(min);
            available.at(min) = false;
            weights.push_back(partitionWeight*pow(static_cast<float>(n)-(static_cast<float>(n)/static_cast<float>(t))*static_cast<float>(i),1.0/d)/setsWeight);

            // Remove from the weight of edges the weight of sets that intersect the selected edge
            for (int j = 0; j < ss.sets.size(); j++) {
                if (!intersect_partition.at(j) && intersects(Edge(start,min),ss.sets.at(j))){
                    intersect_partition.at(j) = 1;
                    if (update_weight.size() > 0 && ss.sets.at(j).weight >= ss.sets.at(update_weight.at(update_weight.size()-1)).weight){
                        update_weight.push_front(j);
                    } else {
                        update_weight.push_back(j);
                    }
                    // Double set weight for next iteration
                    ss.sets.at(j).increase();
                }
            }
            bool test = false;
            while (!test and update_weight.size() > 0){
                int set = update_weight.at(0);
                if (ss.sets.at(set).points.at(start)){
                    for (int& pt : ss.sets.at(set).complement_indices) {
                        if (available.at(pt) && intersects(Edge(start,pt),ss.sets.at(set))){
                            weight.at(pt) -= (1 << (ss.sets.at(set).weight - 1));
                            if ((partitionWeight + weight.at(pt))*pow(static_cast<float>(n)-(static_cast<float>(n)/static_cast<float>(t))*static_cast<float>(i),1.0/d)/setsWeight <= 2*pow(static_cast<float>(k+1),1.0/d)){
                                test = true;
                            }
                        }
                    }
                } else {
                    for (int& pt : ss.sets.at(set).points_indices) {
                        if (available.at(pt) && intersects(Edge(start,pt),ss.sets.at(set))){
                            weight.at(pt) -= (1 << (ss.sets.at(set).weight - 1));
                            if ((partitionWeight + weight.at(pt))*pow(static_cast<float>(n)-(static_cast<float>(n)/static_cast<float>(t))*static_cast<float>(i),1.0/d)/setsWeight <= 2*pow(static_cast<float>(k+1),1.0/d)){
                                test = true;
                            }
                        }
                    }
                }
                update_weight.pop_front();
            }
        }

        // Store results
        res.sets.push_back(partition);
        for (int j = 0; j < m; j++){
            intersections.at(j) += intersect_partition.at(j);
        }
                
    }

    for (int i = 0; i < n ; i++){
        if (available.at(i)){
            res.sets.at(res.sets.size()-1).points.at(i) = 1;
        }
    }

    res.intersections = intersections;
    for (unsigned long& w : weights){
        res.weights.push_back(static_cast<float>(w));
    }
    return res;
}

Result no_weight_update_deque_insert_middle(SetSystem ss, int t){
    const int n = ss.points.size();
    const int d = ss.points.at(0).coordinates.size();
    const int m = ss.sets.size();

    vector<unsigned long> weights;
    vector<int> intersections(m,0);

    //Initialize weight vector
    vector<unsigned long> weight(n,0);
    vector<bool> available(n,true);

    //reset sets weight
    for (Set& s : ss.sets) {
        s.resetWeight();
    }

    vector<Set> s;
    Result res = Result(ss.points, s);
    vector<bool> intersect_partition(m,0);
    deque<int> update_weight;

    for (int i = 0 ; i < t; i++){
        //cout << "\nPartition " << i+1 << "\n";

        Set partition = Set(n);

        unsigned long setsWeight = 0;
        for (Set& s : ss.sets) {
            setsWeight += 1 << s.weight;
        }
        unsigned long partitionWeight = 0.;


        vector<bool> intersect_partition(m,0);

        vector<int> admissible_start;
        for (int k = 0; k < n; k++){
            if (available.at(k)){
                admissible_start.push_back(k);
                weight.at(k) = 0;
            }
        }
        int start = admissible_start.at(rand()%admissible_start.size());
        available.at(start) = false;
        partition.points.at(start) = 1;
        for (Set& s : ss.sets) {
            if (s.points.at(start)){
                for (int& k : s.complement_indices){
                    weight.at(k) += 1 << s.weight;
                }
            } else {
                for (int& k : s.points_indices){
                    weight.at(k) += 1 << s.weight;
                }
            }
        }

        update_weight.clear();

        for (int k = 1 ; k < n / t; k++){
            int min = -1;
            int j = 0;
            while (j < n) {
                if (available.at(j) && (partitionWeight + weight.at(j))*pow(static_cast<float>(n)-(static_cast<float>(n)/static_cast<float>(t))*static_cast<float>(i),1.0/d)/setsWeight <= 2*pow(static_cast<float>(k),1.0/d)){
                    min=j;
                    break;
                } else {
                    // Replace previous minimum if the weight is strictly smaller in all valid cases (no restriction if the partition has no edge yet, connectiveness otherwise)
                    if (available.at(j) && (min == -1 || weight.at(j) < weight.at(min))){
                        min = j;
                    } else {
                        // If the algorithm finds a valid edge with the same weight it swaps it with probability 1/2 to avoid obteining the same result for every algo because of order
                        if (available.at(j) && (min == -1 || weight.at(j) <= weight.at(min) && rand() > 0.5*RAND_MAX)){
                            min = j;
                        }
                    }
                }
                j++;
            }

            //Add selected edge to the partition
            partition.points.at(min) = 1;

            //Update partition weight
            partitionWeight += weight.at(min);
            available.at(min) = false;
            weights.push_back(partitionWeight*pow(static_cast<float>(n)-(static_cast<float>(n)/static_cast<float>(t))*static_cast<float>(i),1.0/d)/setsWeight);

            // Remove from the weight of edges the weight of sets that intersect the selected edge
            for (int j = 0; j < ss.sets.size(); j++) {
                if (!intersect_partition.at(j) && intersects(Edge(start,min),ss.sets.at(j))){
                    intersect_partition.at(j) = 1;
                    if (update_weight.size() > 0 && ss.sets.at(j).weight >= ss.sets.at(update_weight.at((update_weight.size()-1)/2)).weight){
                        update_weight.push_front(j);
                    } else {
                        update_weight.push_back(j);
                    }
                    // Double set weight for next iteration
                    ss.sets.at(j).increase();
                }
            }
            bool test = false;
            while (!test and update_weight.size() > 0){
                int set = update_weight.at(0);
                if (ss.sets.at(set).points.at(start)){
                    for (int& pt : ss.sets.at(set).complement_indices) {
                        if (available.at(pt) && intersects(Edge(start,pt),ss.sets.at(set))){
                            weight.at(pt) -= (1 << (ss.sets.at(set).weight - 1));
                            if ((partitionWeight + weight.at(pt))*pow(static_cast<float>(n)-(static_cast<float>(n)/static_cast<float>(t))*static_cast<float>(i),1.0/d)/setsWeight <= 2*pow(static_cast<float>(k+1),1.0/d)){
                                test = true;
                            }
                        }
                    }
                } else {
                    for (int& pt : ss.sets.at(set).points_indices) {
                        if (available.at(pt) && intersects(Edge(start,pt),ss.sets.at(set))){
                            weight.at(pt) -= (1 << (ss.sets.at(set).weight - 1));
                            if ((partitionWeight + weight.at(pt))*pow(static_cast<float>(n)-(static_cast<float>(n)/static_cast<float>(t))*static_cast<float>(i),1.0/d)/setsWeight <= 2*pow(static_cast<float>(k+1),1.0/d)){
                                test = true;
                            }
                        }
                    }
                }
                update_weight.pop_front();
            }
        }

        // Store results
        res.sets.push_back(partition);
        for (int j = 0; j < m; j++){
            intersections.at(j) += intersect_partition.at(j);
        }
                
    }

    for (int i = 0; i < n ; i++){
        if (available.at(i)){
            res.sets.at(res.sets.size()-1).points.at(i) = 1;
        }
    }

    res.intersections = intersections;
    for (unsigned long& w : weights){
        res.weights.push_back(static_cast<float>(w));
    }
    return res;
}

Result no_weight_update_insert_sorted(SetSystem ss, int t){
const int n = ss.points.size();
    const int d = ss.points.at(0).coordinates.size();
    const int m = ss.sets.size();

    vector<unsigned long> weights;
    vector<int> intersections(m,0);

    //Initialize weight vector
    vector<unsigned long> weight(n,0);
    vector<bool> available(n,true);

    //reset sets weight
    for (Set& s : ss.sets) {
        s.resetWeight();
    }

    vector<Set> s;
    Result res = Result(ss.points, s);
    vector<bool> intersect_partition(m,0);
    vector<Set> update_weight;

    for (int i = 0 ; i < t; i++){
        //cout << "\nPartition " << i+1 << "\n";

        Set partition = Set(n);

        unsigned long setsWeight = 0;
        for (Set& s : ss.sets) {
            setsWeight += 1 << s.weight;
        }
        unsigned long partitionWeight = 0.;


        vector<bool> intersect_partition(m,0);

        vector<int> admissible_start;
        for (int k = 0; k < n; k++){
            if (available.at(k)){
                admissible_start.push_back(k);
                weight.at(k) = 0;
            }
        }
        int start = admissible_start.at(rand()%admissible_start.size());
        available.at(start) = false;
        partition.points.at(start) = 1;
        for (Set& s : ss.sets) {
            if (s.points.at(start)){
                for (int& k : s.complement_indices){
                    weight.at(k) += 1 << s.weight;
                }
            } else {
                for (int& k : s.points_indices){
                    weight.at(k) += 1 << s.weight;
                }
            }
        }

        update_weight.clear();

        for (int k = 1 ; k < n / t; k++){
            int min = -1;
            int j = 0;
            while (j < n) {
                if (available.at(j) && (partitionWeight + weight.at(j))*pow(static_cast<float>(n)-(static_cast<float>(n)/static_cast<float>(t))*static_cast<float>(i),1.0/d)/setsWeight <= 1*pow(static_cast<float>(k),1.0/d)){
                    min=j;
                    break;
                } else {
                    // Replace previous minimum if the weight is strictly smaller in all valid cases (no restriction if the partition has no edge yet, connectiveness otherwise)
                    if (available.at(j) && (min == -1 || weight.at(j) < weight.at(min))){
                        min = j;
                    } else {
                        // If the algorithm finds a valid edge with the same weight it swaps it with probability 1/2 to avoid obteining the same result for every algo because of order
                        if (available.at(j) && (min == -1 || weight.at(j) <= weight.at(min) && rand() > 0.5*RAND_MAX)){
                            min = j;
                        }
                    }
                }
                j++;
            }

            //Add selected edge to the partition
            partition.points.at(min) = 1;

            //Update partition weight
            partitionWeight += weight.at(min);
            available.at(min) = false;
            weights.push_back(partitionWeight*pow(static_cast<float>(n)-(static_cast<float>(n)/static_cast<float>(t))*static_cast<float>(i),1.0/d)/setsWeight);

            // Remove from the weight of edges the weight of sets that intersect the selected edge
            for (int j = 0; j < ss.sets.size(); j++) {
                if (!intersect_partition.at(j) && intersects(Edge(start,min),ss.sets.at(j))){
                    intersect_partition.at(j) = 1;
                    update_weight.push_back(ss.sets.at(j));
                    // Double set weight for next iteration
                    ss.sets.at(j).increase();
                }
            }

            sort(update_weight.begin(), update_weight.end(),setWeightOrder);
            bool test = false;
            while (!test and update_weight.size() > 0){
                Set set = update_weight.at(update_weight.size()-1);
                if (set.points.at(start)){
                    for (int& pt : set.complement_indices) {
                        if (available.at(pt) && intersects(Edge(start,pt),set)){
                            weight.at(pt) -= (1 << (set.weight - 1));
                            if ((partitionWeight + weight.at(pt))*pow(static_cast<float>(n)-(static_cast<float>(n)/static_cast<float>(t))*static_cast<float>(i),1.0/d)/setsWeight <= 1*pow(static_cast<float>(k+1),1.0/d)){
                                test = true;
                            }
                        }
                    }
                } else {
                    for (int& pt : set.points_indices) {
                        if (available.at(pt) && intersects(Edge(start,pt),set)){
                            weight.at(pt) -= (1 << (set.weight - 1));
                            if ((partitionWeight + weight.at(pt))*pow(static_cast<float>(n)-(static_cast<float>(n)/static_cast<float>(t))*static_cast<float>(i),1.0/d)/setsWeight <= 1*pow(static_cast<float>(k+1),1.0/d)){
                                test = true;
                            }
                        }
                    }
                }
                update_weight.pop_back();
            }
        }

        // Store results
        res.sets.push_back(partition);
        for (int j = 0; j < m; j++){
            intersections.at(j) += intersect_partition.at(j);
        }
                
    }

    for (int i = 0; i < n ; i++){
        if (available.at(i)){
            res.sets.at(res.sets.size()-1).points.at(i) = 1;
        }
    }

    res.intersections = intersections;
    for (unsigned long& w : weights){
        res.weights.push_back(static_cast<float>(w));
    }
    return res;
}

Result no_weight_update_point_sort_set_sort(SetSystem ss, int t){
    const int n = ss.points.size();
    const int d = ss.points.at(0).coordinates.size();
    const int m = ss.sets.size();

    vector<float> weights;
    vector<int> intersections(m,0);

    //Initialize weight vector
    vector<unsigned long> weight(n,0);
    vector<bool> available(n,true);

    //reset sets weight
    for (Set& s : ss.sets) {
        s.resetWeight();
    }

    vector<Set> s;
    Result res = Result(ss.points, s);
    vector<bool> intersect_partition(m,0);

    for (int i = 0 ; i < t; i++){
        //cout << "\nPartition " << i+1 << "\n";

        Set partition = Set(n);

        unsigned long setsWeight = 0;
        for (Set& s : ss.sets) {
            setsWeight += 1 << s.weight;
        }
        unsigned long partitionWeight = 0.;


        vector<bool> intersect_partition(m,0);

        vector<int> admissible_start;
        for (int k = 0; k < n; k++){
            if (available.at(k)){
                admissible_start.push_back(k);
                weight.at(k) = 0;
            }
        }
        int start = admissible_start.at(rand()%admissible_start.size());
        available.at(start) = false;
        partition.points.at(start) = 1;
        for (Set& s : ss.sets) {
            if (s.points.at(start)){
                for (int& k : s.complement_indices){
                    weight.at(k) += 1 << s.weight;
                }
            } else {
                for (int& k : s.points_indices){
                    weight.at(k) += 1 << s.weight;
                }
            }
        }

        vector<tuple<int,int>> point_order;
        for (int pt = 0; pt < n; pt ++){
            if (available.at(pt)){
                point_order.push_back(make_tuple(pt,weight.at(pt)));
            }
        }
        sort(point_order.begin(),point_order.end(),indexWeightOrder);

        vector<Set> update_weight;
        update_weight.clear();

        for (int k = 1 ; k < n / t; k++){
            int minm = -1;
            int j = 0;
            while (j < point_order.size()) {
                if (available.at(get<0>(point_order.at(j))) && (partitionWeight + get<1>(point_order.at(j)))*pow(static_cast<float>(n)-(static_cast<float>(n)/static_cast<float>(t))*static_cast<float>(i),1.0/d)/setsWeight <= 1*pow(static_cast<float>(k),1.0/d)){
                    minm=j;
                    break;
                } else {
                    // Replace previous minimum if the weight is strictly smaller in all valid cases (no restriction if the partition has no edge yet, connectiveness otherwise)
                    if (available.at(get<0>(point_order.at(j))) && (minm == -1 || get<0>(point_order.at(j)) < get<1>(point_order.at(minm)))){
                        minm = j;
                    } else {
                        // If the algorithm finds a valid edge with the same weight it swaps it with probability 1/2 to avoid obteining the same result for every algo because of order
                        if (available.at(get<0>(point_order.at(j))) && (minm == -1 || get<1>(point_order.at(j)) <= get<1>(point_order.at(minm)) && rand() > 0.5*RAND_MAX)){
                            minm = j;
                        }
                    }
                }
                j++;
            }      
            int min = get<0>(point_order.at(minm));
            point_order.erase(point_order.begin()+minm);
            //Add selected edge to the partition
            partition.points.at(min) = 1;

            //Update partition weight
            partitionWeight += weight.at(min);
            available.at(min) = false;
            weights.push_back(partitionWeight*pow(static_cast<float>(n)-(static_cast<float>(n)/static_cast<float>(t))*static_cast<float>(i),1.0/d)/setsWeight);
            // Remove from the weight of edges the weight of sets that intersect the selected edge
            for (int j = 0; j < ss.sets.size(); j++) {
                if (!intersect_partition.at(j) && intersects(Edge(start,min),ss.sets.at(j))){
                    intersect_partition.at(j) = 1;
                    update_weight.push_back(ss.sets.at(j));

                    // Double set weight for next iteration
                    ss.sets.at(j).increase();
                }
            }
            sort(update_weight.begin(),update_weight.end(),setWeightOrder);
            bool test = false;
            
            while (!test and update_weight.size() > 0){
                Set set = update_weight.at(update_weight.size()-1);
                for (tuple<int,int> pt : point_order) {
                    if (intersects(Edge(start,get<0>(pt)),set)){
                        int newWeight = get<1>(pt) - (1 << (set.weight - 1));
                        weight.at(get<0>(pt)) = newWeight;
                        pt = make_tuple(get<0>(pt), newWeight);
                        if ((partitionWeight + newWeight)*pow(static_cast<float>(n)-(static_cast<float>(n)/static_cast<float>(t))*static_cast<float>(i),1.0/d)/setsWeight <= 1*pow(static_cast<float>(k+1),1.0/d)){
                            test = true;
                        }
                    }
                }
                update_weight.pop_back();
            }
        }

        // Store results
        res.sets.push_back(partition);
        for (int j = 0; j < m; j++){
            intersections.at(j) += intersect_partition.at(j);
        }
                
    }

    for (int i = 0; i < n ; i++){
        if (available.at(i)){
            res.sets.at(res.sets.size()-1).points.at(i) = 1;
        }
    }

    res.intersections = intersections;
    res.weights=weights;
    return res;
}

Result no_weight_update_point_sort_set_deque(SetSystem ss, int t){
    const int n = ss.points.size();
    const int d = ss.points.at(0).coordinates.size();
    const int m = ss.sets.size();

    vector<unsigned long> weights;
    vector<int> intersections(m,0);

    //Initialize weight vector
    vector<unsigned long> weight(n,0);
    vector<bool> available(n,true);

    //reset sets weight
    for (Set& s : ss.sets) {
        s.resetWeight();
    }

    vector<Set> s;
    Result res = Result(ss.points, s);
    vector<bool> intersect_partition(m,0);

    for (int i = 0 ; i < t; i++){
        //cout << "\nPartition " << i+1 << "\n";

        Set partition = Set(n);

        unsigned long setsWeight = 0;
        for (Set& s : ss.sets) {
            setsWeight += 1 << s.weight;
        }
        unsigned long partitionWeight = 0.;


        vector<bool> intersect_partition(m,0);

        vector<int> admissible_start;
        for (int k = 0; k < n; k++){
            if (available.at(k)){
                admissible_start.push_back(k);
                weight.at(k) = 0;
            }
        }
        int start = admissible_start.at(rand()%admissible_start.size());
        available.at(start) = false;
        partition.points.at(start) = 1;
        for (Set& s : ss.sets) {
            if (s.points.at(start)){
                for (int& k : s.complement_indices){
                    weight.at(k) += 1 << s.weight;
                }
            } else {
                for (int& k : s.points_indices){
                    weight.at(k) += 1 << s.weight;
                }
            }
        }

        vector<tuple<int,int>> point_order;
        for (int pt = 0; pt < n; pt ++){
            if (available.at(pt)){
                point_order.push_back(make_tuple(pt,weight.at(pt)));
            }
        }
        sort(point_order.begin(),point_order.end(),indexWeightOrder);
        
				deque<int> update_weight;
        update_weight.clear();

        for (int k = 1 ; k < n / t; k++){
            int minm = -1;
            int j = 0;
            while (j < point_order.size()) {
                if (available.at(get<0>(point_order.at(j))) && (partitionWeight + get<1>(point_order.at(j)))*pow(static_cast<float>(n)-(static_cast<float>(n)/static_cast<float>(t))*static_cast<float>(i),1.0/d)/setsWeight <= 1*pow(static_cast<float>(k),1.0/d)){
                    minm=j;
                    break;
                } else {
                    // Replace previous minimum if the weight is strictly smaller in all valid cases (no restriction if the partition has no edge yet, connectiveness otherwise)
                    if (available.at(get<0>(point_order.at(j))) && (minm == -1 || get<0>(point_order.at(j)) < get<1>(point_order.at(minm)))){
                        minm = j;
                    } else {
                        // If the algorithm finds a valid edge with the same weight it swaps it with probability 1/2 to avoid obteining the same result for every algo because of order
                        if (available.at(get<0>(point_order.at(j))) && (minm == -1 || get<1>(point_order.at(j)) <= get<1>(point_order.at(minm)) && rand() > 0.5*RAND_MAX)){
                            minm = j;
                        }
                    }
                }
                j++;
            }            
            int min = get<0>(point_order.at(minm));
            point_order.erase(point_order.begin()+minm);
            //Add selected edge to the partition
            partition.points.at(min) = 1;

            //Update partition weight
            partitionWeight += weight.at(min);
            available.at(min) = false;
            weights.push_back(partitionWeight*pow(static_cast<float>(n)-(static_cast<float>(n)/static_cast<float>(t))*static_cast<float>(i),1.0/d)/setsWeight);

            // Remove from the weight of edges the weight of sets that intersect the selected edge
            for (int j = 0; j < ss.sets.size(); j++) {
                if (!intersect_partition.at(j) && intersects(Edge(start,min),ss.sets.at(j))){
                    intersect_partition.at(j) = 1;
                    if (update_weight.size() > 0 && ss.sets.at(j).weight >= ss.sets.at(update_weight.at((update_weight.size()-1)/2)).weight){
                        update_weight.push_front(j);
                    } else {
                        update_weight.push_back(j);
                    }
                    // Double set weight for next iteration
                    ss.sets.at(j).increase();
                }
            }

            bool test = false;
            while (!test and update_weight.size() > 0){
                Set set = ss.sets.at(update_weight.at(update_weight.size()-1));
                for (tuple<int,int> pt : point_order) {
                    if (intersects(Edge(start,get<0>(pt)),set)){
                        int newWeight = get<1>(pt) - (1 << (set.weight - 1));
                        weight.at(get<0>(pt)) = newWeight;
                        pt = make_tuple(get<0>(pt), newWeight);
                        if ((partitionWeight + newWeight)*pow(static_cast<float>(n)-(static_cast<float>(n)/static_cast<float>(t))*static_cast<float>(i),1.0/d)/setsWeight <= 1*pow(static_cast<float>(k+1),1.0/d)){
                            test = true;
                        }
                    }
                }
                update_weight.pop_back();
            }
        }

        // Store results
        res.sets.push_back(partition);
        for (int j = 0; j < m; j++){
            intersections.at(j) += intersect_partition.at(j);
        }
                
    }

    for (int i = 0; i < n ; i++){
        if (available.at(i)){
            res.sets.at(res.sets.size()-1).points.at(i) = 1;
        }
    }

    res.intersections = intersections;
    for (unsigned long& w : weights){
        res.weights.push_back(static_cast<float>(w));
    }
    return res;
}

Result no_weight_update_point_sort_invert(SetSystem ss, int t){
    const int n = ss.points.size();
    const int d = ss.points.at(0).coordinates.size();
    const int m = ss.sets.size();

    vector<float> weights;
    vector<int> intersections(m,0);

    //Initialize weight vector
    vector<unsigned long> weight(n,0);
    vector<bool> available(n,true);

    //reset sets weight
    for (Set& s : ss.sets) {
        s.resetWeight();
    }

    vector<Set> s;
    Result res = Result(ss.points, s);
    vector<bool> intersect_partition(m,0);

    for (int i = 0 ; i < t; i++){
        //cout << "\nPartition " << i+1 << "\n";

        Set partition = Set(n);

        unsigned long setsWeight = 0;
        for (Set& s : ss.sets) {
            setsWeight += 1 << s.weight;
        }
        unsigned long partitionWeight = 0.;


        vector<bool> intersect_partition(m,0);

        vector<int> admissible_start;
        for (int k = 0; k < n; k++){
            if (available.at(k)){
                admissible_start.push_back(k);
                weight.at(k) = 0;
            }
        }
        int start = admissible_start.at(rand()%admissible_start.size());
        available.at(start) = false;
        partition.points.at(start) = 1;
        for (Set& s : ss.sets) {
            if (s.points.at(start)){
                for (int& k : s.complement_indices){
                    weight.at(k) += 1 << s.weight;
                }
            } else {
                for (int& k : s.points_indices){
                    weight.at(k) += 1 << s.weight;
                }
            }
        }

        vector<tuple<int,unsigned long>> point_order;
        for (int pt = 0; pt < n; pt ++){
            if (available.at(pt)){
                point_order.push_back(make_tuple(pt,weight.at(pt)));
            }
        }
        sort(point_order.begin(),point_order.end(),indexWeightOrder);

        vector<vector<Set>> update_weight;
        update_weight.clear();
        vector<int> last_updated(n,0);

        for (int k = 1 ; k < n / t; k++){
            int minm = -1;
            int j = 0;
            while (j < point_order.size()) {
                if (available.at(get<0>(point_order.at(j))) && (partitionWeight + get<1>(point_order.at(j)))*pow(static_cast<float>(n)-(static_cast<float>(n)/static_cast<float>(t))*static_cast<float>(i),1.0/d)/setsWeight <= 1*pow(static_cast<float>(k),1.0/d)){
                    minm=j;
                    break;
                } else {
                    // Replace previous minimum if the weight is strictly smaller in all valid cases (no restriction if the partition has no edge yet, connectiveness otherwise)
                    if (available.at(get<0>(point_order.at(j))) && (minm == -1 || get<0>(point_order.at(j)) < get<1>(point_order.at(minm)))){
                        minm = j;
                    } else {
                        // If the algorithm finds a valid edge with the same weight it swaps it with probability 1/2 to avoid obteining the same result for every algo because of order
                        if (available.at(get<0>(point_order.at(j))) && (minm == -1 || get<1>(point_order.at(j)) <= get<1>(point_order.at(minm)) && rand() > 0.5*RAND_MAX)){
                            minm = j;
                        }
                    }
                }
                j++;
            }      
            int min = get<0>(point_order.at(minm));
            point_order.erase(point_order.begin()+minm);
            //Add selected edge to the partition
            partition.points.at(min) = 1;

            //Update partition weight
            partitionWeight += weight.at(min);
            available.at(min) = false;
            weights.push_back(partitionWeight*pow(static_cast<float>(n)-(static_cast<float>(n)/static_cast<float>(t))*static_cast<float>(i),1.0/d)/setsWeight);
            // Remove from the weight of edges the weight of sets that intersect the selected edge
            vector<Set> temp;
            for (int j = 0; j < ss.sets.size(); j++) {
                if (!intersect_partition.at(j) && intersects(Edge(start,min),ss.sets.at(j))){
                    intersect_partition.at(j) = 1;
                    temp.push_back(ss.sets.at(j));

                    // Double set weight for next iteration
                    ss.sets.at(j).increase();
                }
            }
            update_weight.push_back(temp);
            bool test = false;
            
            int p = 0;
            while (!test and p < point_order.size()){
                tuple<int,unsigned long> pt = point_order.at(p);
                for (Set set : update_weight.at(last_updated.at(get<0>(pt)))) {
                    if (intersects(Edge(start,get<0>(pt)),set)){
                        weight.at(get<0>(pt)) -= (1 << (set.weight - 1));
                        if ((partitionWeight + weight.at(get<0>(pt)))*pow(static_cast<float>(n)-(static_cast<float>(n)/static_cast<float>(t))*static_cast<float>(i),1.0/d)/setsWeight < 1*pow(static_cast<float>(k+1),1.0/d)){
                            test = true;
                        }
                    }
                }
                if (last_updated.at(get<0>(pt)) == update_weight.size() - 1){
                    p++;
                } else {
                    last_updated.at(get<0>(pt)) ++;
                }
            }
        }

        // Store results
        res.sets.push_back(partition);
        for (int j = 0; j < m; j++){
            intersections.at(j) += intersect_partition.at(j);
        }
                
    }

    for (int i = 0; i < n ; i++){
        if (available.at(i)){
            res.sets.at(res.sets.size()-1).points.at(i) = 1;
        }
    }

    res.intersections = intersections;
    res.weights=weights;
    return res;
}

Result partition_no_set(SetSystem ss, int t, vector<float> (*lf)(vector<Point>, vector<bool>, int, vector<Set>)){
    
    const int n = ss.points.size();
    const int d = ss.points.at(0).coordinates.size();
    const int m = ss.sets.size();

    //Initialize weight vector
    vector<bool> available(n,true);

    vector<Set> s;
    Result res = Result(ss.points, s);

    for (int i = 0 ; i < t; i++){
        //cout << "\nPartition " << i+1 << "\n";

        Set partition = Set(n);

        vector<int> admissible_start;
        for (int k = 0; k < n; k++){
            if (available.at(k)){
                admissible_start.push_back(k);
            }
        }
        int start = admissible_start.at(rand()%admissible_start.size());
        partition.points.at(start) = 1;
        available.at(start) = false;

        vector<float> distances = lf(ss.points, available, start, ss.sets);
        vector<tuple<int,float>> tosort;
        for (int i = 0; i < n ; i++){
            if (available.at(i)){
                tosort.push_back(make_tuple(i,distances.at(i)));
            }
        }

        sort(tosort.begin(),tosort.end(),floatWeightOrder);

        for (int i = 0; i < n/t-1; i++){
            //Add selected edge to the partition
            partition.points.at(get<0>(tosort.at(i))) = 1;

            //Update partition weight
            available.at(get<0>(tosort.at(i))) = false;
        
        }
        // Store results
        res.sets.push_back(partition);
    }

    //fill the last partition with the remaining points
    for (int i = 0; i < n ; i++){
        if (available.at(i)){
            res.sets.at(res.sets.size()-1).points.at(i) = 1;
        }
    }

    return res;
}
