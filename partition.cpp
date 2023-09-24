#include <bits/types/clock_t.h>
#include <cstdlib>
#include <iostream>
#include <ostream>
#include <type_traits>
#include <vector>
#include <cmath>
#include <ctime>

#include "classes.cpp"

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

SetSystem partition_min(SetSystem ss, int t){
    
    const int n = ss.points.size();
    const int d = ss.points.at(0).coordinates.size();
    const int m = ss.sets.size();

    //Store partition weight at each iteration
    vector<float> weights;

    //Initialize weight vector
    vector<unsigned long> weight(n,0);
    vector<bool> available(n,true);

    //reset sets weight
    for (Set& s : ss.sets) {
        s.resetWeight();
    }

    vector<Set> s;
    SetSystem res = SetSystem(ss.points, s);
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
            for (int k = 0; k < n; k++){
                if (available.at(k) && intersects(Edge(start,k), s)){
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
            weights.push_back(partitionWeight*sqrt(n-(n/t)*i)/setsWeight);

            // Remove from the weight of edges the weight of sets that intersect the selected edge
            for (int j = 0; j < ss.sets.size(); j++) {
                if (!intersect_partition.at(j) && intersects(Edge(start,min),ss.sets.at(j))){
                    intersect_partition.at(j) = 1;
                    for (int pt = 0; pt < n; pt++) {
                        if (available.at(pt) && intersects(Edge(start,pt),ss.sets.at(j))){
                            weight.at(pt) -= (1 << ss.sets.at(j).weight);
                        }
                    }
                    // Double set weight for next iteration
                    ss.sets.at(j).increase();
                }
            }
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

SetSystem partition_min_rate(SetSystem ss, int t){
    
    const int n = ss.points.size();
    const int d = ss.points.at(0).coordinates.size();
    const int m = ss.sets.size();

    vector<float> weights;

    //Initialize weight vector
    vector<unsigned long> weight(n,0);
    vector<bool> available(n,true);

    //reset sets weight
    for (Set& s : ss.sets) {
        s.resetWeight();
    }

    vector<Set> s;
    SetSystem res = SetSystem(ss.points, s);
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
            for (int k = 0; k < n; k++){
                if (available.at(k) && intersects(Edge(start,k), s)){
                    weight.at(k) += 1 << s.weight;
                }
            }
        }

        for (int k = 1 ; k < n / t; k++){
            int min = -1;
            bool test = false;
            int j = 0;
            while (!test && j < n) {
                if (available.at(j) && (partitionWeight + weight.at(j))*sqrt(n-(n/t)*i)/setsWeight <= 2*pow(k,1/d)){
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
            weights.push_back(partitionWeight*sqrt(n-(n/t)*i)/setsWeight);

            // Remove from the weight of edges the weight of sets that intersect the selected edge
            for (int j = 0; j < ss.sets.size(); j++) {
                if (!intersect_partition.at(j) && intersects(Edge(start,min),ss.sets.at(j))){
                    intersect_partition.at(j) = 1;
                    for (int pt = 0; pt < n; pt++) {
                        if (available.at(pt) && intersects(Edge(start,pt),ss.sets.at(j))){
                            weight.at(pt) -= (1 << ss.sets.at(j).weight);
                        }
                    }
                    // Double set weight for next iteration
                    ss.sets.at(j).increase();
                }
            }
        }

        // Store results
        res.sets.push_back(partition);
                
    }

    for (int i = 0; i < n ; i++){
        if (available.at(i)){
            res.sets.at(res.sets.size()-1).points.at(i) = 1;
        }
    }

    return res;
}

SetSystem partition_rate(SetSystem ss, int t){
    
    const int n = ss.points.size();
    const int d = ss.points.at(0).coordinates.size();
    const int m = ss.sets.size();

    vector<float> weights;

    //Initialize weight vector
    vector<unsigned long> weight(n,0);
    vector<bool> available(n,true);

    //reset sets weight
    for (Set& s : ss.sets) {
        s.resetWeight();
    }

    vector<Set> s;
    SetSystem res = SetSystem(ss.points, s);
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
            for (int k = 0; k < n; k++){
                if (available.at(k) && intersects(Edge(start,k), s)){
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
                if (available.at(j) && (partitionWeight + weight.at(j))*sqrt(n-(n/t)*i)/setsWeight <= 2*pow(k,1/d)){
                    candidates.push_back(j);
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
            weights.push_back(partitionWeight*sqrt(n-(n/t)*i)/setsWeight);

            // Remove from the weight of edges the weight of sets that intersect the selected edge
            for (int j = 0; j < ss.sets.size(); j++) {
                if (!intersect_partition.at(j) && intersects(Edge(start,min),ss.sets.at(j))){
                    intersect_partition.at(j) = 1;
                    for (int pt = 0; pt < n; pt++) {
                        if (available.at(pt) && intersects(Edge(start,pt),ss.sets.at(j))){
                            weight.at(pt) -= (1 << ss.sets.at(j).weight);
                        }
                    }
                    // Double set weight for next iteration
                    ss.sets.at(j).increase();
                }
            }
        }

        // Store results
        res.sets.push_back(partition);
                
    }

    for (int i = 0; i < n ; i++){
        if (available.at(i)){
            res.sets.at(res.sets.size()-1).points.at(i) = 1;
        }
    }

    return res;
}