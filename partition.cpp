#include <bits/types/clock_t.h>
#include <cstdlib>
#include <iostream>
#include <ostream>
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

Result mwu_min(SetSystem ss, int t){
    
    const int n = ss.points.size();
    const int d = ss.points.at(0).coordinates.size();
    const int m = ss.sets.size();

    vector<int> histogram(m,0);
    vector<float> weights;
    time_t timer;
    time(&timer);

    //Initialize edges
    vector<Edge> edges;
    for (int i = 0; i < n - 1; i++){
        for (int j = i + 1; j < n; j++){
            edges.push_back(Edge(i, j));
        }
    }

    //reset sets weight
    for (Set& s : ss.sets) {
        s.resetWeight();
    }

    for (Edge& e : edges) {
    e.weight = 0;
    e.next_ite_weight = 0;
        for (Set& s : ss.sets){
            if (intersects(e,s)){
                e.weight += (1 << s.weight);
                e.next_ite_weight += (1 << s.weight);
            }
        }
    }

    vector<Set> s;
    SetSystem res = SetSystem(ss.points, s);
    vector<bool> intersect_partition(m,0);

    for (int i = 0 ; i < t; i++){
        cout << "\nPartition " << i+1 << "\n";
        Set partition = Set(n);

        int setsWeight = 0;
        for (Set& s : ss.sets) {
            setsWeight += 1 << s.weight;
        }
        float partitionWeight = 0.;


        vector<bool> intersect_partition(m,0);

        // Find edge with minimum weight among all
        for (int k = 1 ; k < n / t; k++){
            time_t temp;
            time(&temp);
            cout << k << " " << (temp-timer) <<endl;
            time(&timer);
            int min = 0;
            for (auto j = 0; j < edges.size(); j++) {
                // Replace previous minimum if the weight is strictly smaller in all valid cases (no restriction if the partition has no edge yet, connectiveness otherwise)
                if (edges.at(j).weight < edges.at(min).weight && (k == 1 || partition.points.at(edges.at(j).points[0]) || partition.points.at(edges.at(j).points[1]))){
                    min = j;
                } else {
                    // If the algorithm finds a valid edge with the same weight it swaps it with probability 1/2 to avoid obteining the same result for every algo because of order
                    if (edges.at(j).weight == edges.at(min).weight && (k == 1 || partition.points.at(edges.at(j).points[0]) || partition.points.at(edges.at(j).points[1])) && rand() > 0.5*RAND_MAX){
                        min = j;
                    }
                }
            }

            //Add selected edge to the partition
            partition.points.at(edges.at(min).points[0]) = 1;
            partition.points.at(edges.at(min).points[1]) = 1;

            //Update partition weight
            partitionWeight += edges.at(min).weight;
            weights.push_back(partitionWeight*sqrt(n-(n/t)*i)/setsWeight);

            // Remove from the weight of edges the weight of sets that intersect the selected edge
            for (int j = 0; j < ss.sets.size(); j++) {
                if (intersects(edges.at(min),ss.sets.at(j)) and !intersect_partition.at(j)){
                    intersect_partition.at(j) = 1;
                    for (Edge& e : edges) {
                        if (intersects(e,ss.sets.at(j))){
                            e.weight -= (1 << ss.sets.at(j).weight);
                            e.next_ite_weight += (1 << ss.sets.at(j).weight);
                        }
                    }
                    // Double set weight for next iteration
                    ss.sets.at(j).increase();
                }
            }

            //Remove edges that connect two points that have been added to the partition by different partitions
            edges.erase(edges.begin() + min);
            for (auto j = edges.size() - 1; j > 0; j--) {
                if (partition.points.at(edges.at(j).points[0]) && partition.points.at(edges.at(j).points[1])){
                    edges.erase(edges.begin() + j);
                }
            }
        }

        // Remove all edges linking to points added to the previous partition
        for (auto j = edges.size() - 1; j > 0; j--) {
                if (partition.points.at(edges.at(j).points[0]) || partition.points.at(edges.at(j).points[1])){
                    edges.erase(edges.begin() + j);
                }
            }
        
        // Store results
        res.sets.push_back(partition);
        for (int j = 0; j < intersect_partition.size();j++){
            histogram.at(j) += intersect_partition.at(j);
        }
                
        //Recompute edge weights
        for (Edge& e : edges) {
            e.weight = e.next_ite_weight;
        }
    }

    return Result(res,histogram,weights);
}

Result mwu_random_uniform(SetSystem ss, int t){
    
    const int n = ss.points.size();
    const int d = ss.points.at(0).coordinates.size();
    const int m = ss.sets.size();

    vector<int> histogram(m,0);
    vector<float> weights;

    //Initialize edges
    vector<Edge> edges;
    for (int i = 0; i < n - 1; i++){
        for (int j = i + 1; j < n; j++){
            edges.push_back(Edge(i, j));
        }
    }

    //reset sets weight
    for (Set& s : ss.sets) {
        s.resetWeight();
    }

    for (Edge& e : edges) {
        e.weight = 0;
        e.next_ite_weight = 0;
            for (Set& s : ss.sets){
                if (intersects(e,s)){
                    e.weight += (1 << s.weight);
                    e.next_ite_weight += (1 << s.weight);
                }
            }
    }

    vector<Set> s;
    SetSystem res = SetSystem(ss.points, s);

    vector<int> candidates;

    for (int i = 0 ; i < t; i++){
        //cout << "\nPartition " << i+1 << "\n";
        Set partition = Set(n);

        int setsWeight = 0;
        for (Set& s : ss.sets) {
            setsWeight += 1 << s.weight;
        }
        float partitionWeight = 0.;

        // Find a random edge from all edges enforcing rate
        vector<bool> intersect_partition(m,0);
        for (int k = 1 ; k < n / t; k++){
            candidates.clear();
            int min = 0;
            for (auto j = 0; j < edges.size(); j++) {
                // Replace previous minimum if the weight is strictly smaller in all valid cases (no restriction if the partition has no edge yet, connectiveness otherwise)
                if (edges.at(j).weight < edges.at(min).weight && (k == 1 || partition.points.at(edges.at(j).points[0]) || partition.points.at(edges.at(j).points[1]))){
                    min = j;
                } else {
                    // If the algorithm finds a valid edge with the same weight it swaps it with probability 1/2 to avoid obteining the same result for every algo because of order
                    if (edges.at(j).weight == edges.at(min).weight && (k == 1 || partition.points.at(edges.at(j).points[0]) || partition.points.at(edges.at(j).points[1])) && rand() > 0.5*RAND_MAX){
                        min = j;
                    }
                }
                if (partitionWeight + edges.at(j).weight <= setsWeight*sqrt(k+1)/sqrt(n-(n/t)*i) && (k == 1 || partition.points.at(edges.at(j).points[0]) || partition.points.at(edges.at(j).points[1]))){
                    candidates.push_back(j);
                }
            }
            if (candidates.size() > 0){
                min = candidates.at(candidates.size()*rand()/RAND_MAX);
            }

            //Add selected edge to the partition
            partition.points.at(edges.at(min).points[0]) = 1;
            partition.points.at(edges.at(min).points[1]) = 1;

            //Update partition weight
            partitionWeight += edges.at(min).weight;
            weights.push_back(partitionWeight*sqrt(n-(n/t)*i)/setsWeight);

            // Remove from the weight of edges the weight of sets that intersect the selected edge
            for (int j = 0; j < ss.sets.size(); j++) {
                if (intersects(edges.at(min),ss.sets.at(j)) and !intersect_partition.at(j)){
                    intersect_partition.at(j) = 1;
                    for (Edge& e : edges) {
                        if (intersects(e,ss.sets.at(j))){
                            e.weight -= 1 << ss.sets.at(j).weight;
                            e.next_ite_weight += (1 << ss.sets.at(j).weight);
                        }
                    }
                    // Double set weight for next iteration
                    ss.sets.at(j).increase();
                }
            }


            //Remove edges that connect two points that have been added to the partition by different partitions
            edges.erase(edges.begin() + min);
            for (auto j = edges.size() - 1; j > 0; j--) {
                if (partition.points.at(edges.at(j).points[0]) && partition.points.at(edges.at(j).points[1])){
                    edges.erase(edges.begin() + j);

                }
            }
        }

        // Remove all edges linking to points added to the previous partition
        for (auto j = edges.size() - 1; j > 0; j--) {
                if (partition.points.at(edges.at(j).points[0]) || partition.points.at(edges.at(j).points[1])){
                    edges.erase(edges.begin() + j);
                }
            }

        // Store results
        res.sets.push_back(partition);
        for (int j = 0; j < intersect_partition.size();j++){
            histogram.at(j) += intersect_partition.at(j);
        }

        //Recompute edge weights
        for (Edge& e : edges) {
            e.weight = e.next_ite_weight;
        }
    }

    return Result(res,histogram,weights);
}

Result mwu_random_weight(SetSystem ss, int t){
    
    const int n = ss.points.size();
    const int d = ss.points.at(0).coordinates.size();
    const int m = ss.sets.size();

    vector<int> histogram(m,0);
    vector<float> weights;

    //Initialize edges
    vector<Edge> edges;
    for (int i = 0; i < n - 1; i++){
        for (int j = i + 1; j < n; j++){
            edges.push_back(Edge(i, j));
        }
    }

    //reset sets weight
    for (Set& s : ss.sets) {
        s.resetWeight();
    }

    for (Edge& e : edges) {
        e.weight = 0;
        e.next_ite_weight = 0;
            for (Set& s : ss.sets){
                if (intersects(e,s)){
                    e.weight += (1 << s.weight);
                    e.next_ite_weight += (1 << s.weight);
                }
            }
    }

    vector<Set> s;
    SetSystem res = SetSystem(ss.points, s);

    vector<int> candidates;

    for (int i = 0 ; i < t; i++){
        //cout << "\nPartition " << i+1 << "\n";
        Set partition = Set(n);

        int setsWeight = 0;
        for (Set& s : ss.sets) {
            setsWeight += 1 << s.weight;
        }
        float partitionWeight = 0.;

        // Find a random edge from all edges enforcing rate with respect to exp(-w)
        vector<bool> intersect_partition(m,0);
        for (int k = 1 ; k < n / t; k++){
            candidates.clear();
            int min = 0;
            for (auto j = 0; j < edges.size(); j++) {
                // Replace previous minimum if the weight is strictly smaller in all valid cases (no restriction if the partition has no edge yet, connectiveness otherwise)
                if (edges.at(j).weight < edges.at(min).weight && (k == 1 || partition.points.at(edges.at(j).points[0]) || partition.points.at(edges.at(j).points[1]))){
                    min = j;
                } else {
                    // If the algorithm finds a valid edge with the same weight it swaps it with probability 1/2 to avoid obteining the same result for every algo because of order
                    if (edges.at(j).weight == edges.at(min).weight && (k == 1 || partition.points.at(edges.at(j).points[0]) || partition.points.at(edges.at(j).points[1])) && rand() > 0.5*RAND_MAX){
                        min = j;
                    }
                }
                if (partitionWeight + edges.at(j).weight <= setsWeight*sqrt(k+1)/sqrt(n-(n/t)*i) && (k == 1 || partition.points.at(edges.at(j).points[0]) || partition.points.at(edges.at(j).points[1]))){
                    candidates.push_back(j);
                }
            }
            if (candidates.size() > 0){
                float sum = 0;
                float partial = 0;
                for (auto j = 0; j < candidates.size(); j++) {
                    sum += exp(-1*edges.at(candidates.at(j)).weight);
                }
                int j=0;
                float coeff = static_cast<float>(rand())/static_cast<float>(RAND_MAX);
                while(partial < sum*coeff && j < candidates.size() - 1){
                    partial += exp(-1*edges.at(candidates.at(j)).weight);
                    j++;
                }
                min = candidates.at(j);
            }

            //Add selected edge to the partition
            partition.points.at(edges.at(min).points[0]) = 1;
            partition.points.at(edges.at(min).points[1]) = 1;

            //Update partition weight
            partitionWeight += edges.at(min).weight;
            weights.push_back(partitionWeight*sqrt(n-(n/t)*i)/setsWeight);

            // Remove from the weight of edges the weight of sets that intersect the selected edge
            for (int j = 0; j < ss.sets.size(); j++) {
                if (intersects(edges.at(min),ss.sets.at(j)) and !intersect_partition.at(j)){
                    intersect_partition.at(j) = 1;
                    for (Edge& e : edges) {
                        if (intersects(e,ss.sets.at(j))){
                            e.weight -= 1 << ss.sets.at(j).weight;
                            e.next_ite_weight += (1 << ss.sets.at(j).weight);
                        }
                    }
                    // Double set weight for next iteration
                    ss.sets.at(j).increase();
                }
            }


            //Remove edges that connect two points that have been added to the partition by different partitions
            edges.erase(edges.begin() + min);
            for (auto j = edges.size() - 1; j > 0; j--) {
                if (partition.points.at(edges.at(j).points[0]) && partition.points.at(edges.at(j).points[1])){
                    edges.erase(edges.begin() + j);

                }
            }
        }

        // Remove all edges linking to points added to the previous partition
        for (auto j = edges.size() - 1; j > 0; j--) {
                if (partition.points.at(edges.at(j).points[0]) || partition.points.at(edges.at(j).points[1])){
                    edges.erase(edges.begin() + j);
                }
            }

        // Store results
        res.sets.push_back(partition);
        for (int j = 0; j < intersect_partition.size();j++){
            histogram.at(j) += intersect_partition.at(j);
        }

        //Recompute edge weights
        for (Edge& e : edges) {
            e.weight = e.next_ite_weight;
        }
    }

    return Result(res,histogram,weights);
}

Result mwu_bfs(SetSystem ss, int t){
    
    const int n = ss.points.size();
    const int d = ss.points.at(0).coordinates.size();
    const int m = ss.sets.size();

    vector<int> histogram(m,0);
    vector<float> weights;


    //Initialize edges
    vector<Edge> edges;
    for (int i = 0; i < n - 1; i++){
        for (int j = i + 1; j < n; j++){
            edges.push_back(Edge(i, j));
        }
    }

    vector<int> bfs_table(n,t+1);

    //reset sets weight
    for (Set& s : ss.sets) {
        s.resetWeight();
    }

    //Recompute edge weights
    for (Edge& e : edges) {
        e.weight = 0;
        e.next_ite_weight = 0;
            for (Set& s : ss.sets){
                if (intersects(e,s)){
                    e.weight += (1 << s.weight);
                    e.next_ite_weight += (1 << s.weight);
                }
            }
    }

    vector<Set> s;
    SetSystem res = SetSystem(ss.points, s);

    vector<int> candidates;

    for (int i = 0 ; i < t; i++){
        //cout << "\nPartition " << i+1 << "\n";
        Set partition = Set(n);

        int setsWeight = 0;
        for (Set& s : ss.sets) {
            setsWeight += 1 << s.weight;
        }
        float partitionWeight = 0.;

        // Find a random edge from all edges enforcing rate with respect to bfs distance
        vector<bool> intersect_partition(m,0);
        for (int k = 1 ; k < n / t; k++){
            candidates.clear();
            int min = 0;
            for (auto j = 0; j < edges.size(); j++) {
                // Replace previous minimum if the weight is strictly smaller in all valid cases (no restriction if the partition has no edge yet, connectiveness otherwise)
                if (edges.at(j).weight < edges.at(min).weight && (k == 1 || partition.points.at(edges.at(j).points[0]) || partition.points.at(edges.at(j).points[1]))){
                    min = j;
                } else {
                    // If the algorithm finds a valid edge with the same weight it swaps it with probability 1/2 to avoid obteining the same result for every algo because of order
                    if (edges.at(j).weight == edges.at(min).weight && (k == 1 || partition.points.at(edges.at(j).points[0]) || partition.points.at(edges.at(j).points[1])) && rand() > 0.5*RAND_MAX){
                        min = j;
                    }
                }
                if (partitionWeight + edges.at(j).weight <= setsWeight*sqrt(k+1)/sqrt(n-(n/t)*i) && (k == 1 || partition.points.at(edges.at(j).points[0]) || partition.points.at(edges.at(j).points[1]))){
                    candidates.push_back(j);
                }
            }

            int min_index = partition.points.at(edges.at(min).points[0]) == 1 ? 0 : 1;
            if (candidates.size() > 0 && k > 1){
                min = candidates.at(0);
                min_index = partition.points.at(edges.at(candidates.at(0)).points[0]) == 1 ? 0 : 1;
                for (int& j : candidates) {
                        int indexlookup = partition.points.at(edges.at(j).points[0]) == 1 ? 0 : 1;
                        if (bfs_table.at(edges.at(j).points[indexlookup]) < bfs_table.at(edges.at(min).points[min_index])){
                            min = j;
                            min_index = indexlookup;
                        } else {
                            if (bfs_table.at(edges.at(j).points[indexlookup]) == bfs_table.at(edges.at(min).points[min_index]) && rand() > 0.5*RAND_MAX){
                                min = j;
                                min_index = indexlookup;
                            }
                        }
                    }
            }

            if (k == 1){
                bfs_table.at(edges.at(min).points[0]) = 0;
                bfs_table.at(edges.at(min).points[1]) = 0;
            } else {
                bfs_table.at(edges.at(min).points[1 - min_index]) = bfs_table.at(edges.at(min).points[min_index]) + 1;
            }
            
            //Add selected edge to the partition
            partition.points.at(edges.at(min).points[0]) = 1;
            partition.points.at(edges.at(min).points[1]) = 1;

            //Update partition weight
            partitionWeight += edges.at(min).weight;
            weights.push_back(partitionWeight*sqrt(n-(n/t)*i)/setsWeight);

            // Remove from the weight of edges the weight of sets that intersect the selected edge
            for (int j = 0; j < ss.sets.size(); j++) {
                if (intersects(edges.at(min),ss.sets.at(j)) and !intersect_partition.at(j)){
                    intersect_partition.at(j) = 1;
                    for (Edge& e : edges) {
                        if (intersects(e,ss.sets.at(j))){
                            e.weight -= 1 << ss.sets.at(j).weight;
                            e.next_ite_weight += (1 << ss.sets.at(j).weight);
                        }
                    }
                    // Double set weight for next iteration
                    ss.sets.at(j).increase();
                }
            }


            //Remove edges that connect two points that have been added to the partition by different partitions
            edges.erase(edges.begin() + min);
            for (auto j = edges.size() - 1; j > 0; j--) {
                if (partition.points.at(edges.at(j).points[0]) && partition.points.at(edges.at(j).points[1])){
                    edges.erase(edges.begin() + j);

                }
            }
        }

        // Remove all edges linking to points added to the previous partition
        for (auto j = edges.size() - 1; j > 0; j--) {
                if (partition.points.at(edges.at(j).points[0]) || partition.points.at(edges.at(j).points[1])){
                    edges.erase(edges.begin() + j);
                }
            }

        // Store results
        res.sets.push_back(partition);
        for (int j = 0; j < intersect_partition.size();j++){
            histogram.at(j) += intersect_partition.at(j);
        }

        //Recompute edge weights
        for (Edge& e : edges) {
            e.weight = e.next_ite_weight;
        }
    }

    return Result(res,histogram,weights);
}

Result mwu_bfs_min(SetSystem ss, int t){
    
    const int n = ss.points.size();
    const int d = ss.points.at(0).coordinates.size();
    const int m = ss.sets.size();

    vector<int> histogram(m,0);
    vector<float> weights;


    //Initialize edges
    vector<Edge> edges;
    for (int i = 0; i < n - 1; i++){
        for (int j = i + 1; j < n; j++){
            edges.push_back(Edge(i, j));
        }
    }

    vector<int> bfs_table(n,t+1);

    //reset sets weight
    for (Set& s : ss.sets) {
        s.resetWeight();
    }

    //Recompute edge weights
    for (Edge& e : edges) {
        e.weight = 0;
        e.next_ite_weight = 0;
            for (Set& s : ss.sets){
                if (intersects(e,s)){
                    e.weight += (1 << s.weight);
                    e.next_ite_weight += (1 << s.weight);
                }
            }
    }

    vector<Set> s;
    SetSystem res = SetSystem(ss.points, s);

    vector<int> candidates;

    for (int i = 0 ; i < t; i++){
        //cout << "\nPartition " << i+1 << "\n";
        Set partition = Set(n);

        int setsWeight = 0;
        for (Set& s : ss.sets) {
            setsWeight += 1 << s.weight;
        }
        float partitionWeight = 0.;

        // Find a random edge from all edges enforcing rate with respect to bfs distance
        vector<bool> intersect_partition(m,0);
        for (int k = 1 ; k < n / t; k++){
            candidates.clear();
            int min = 0;
            for (auto j = 0; j < edges.size(); j++) {
                // Replace previous minimum if the weight is strictly smaller in all valid cases (no restriction if the partition has no edge yet, connectiveness otherwise)
                if (edges.at(j).weight < edges.at(min).weight && (k == 1 || partition.points.at(edges.at(j).points[0]) || partition.points.at(edges.at(j).points[1]))){
                    min = j;
                } else {
                    // If the algorithm finds a valid edge with the same weight it swaps it with probability 1/2 to avoid obteining the same result for every algo because of order
                    if (edges.at(j).weight == edges.at(min).weight && (k == 1 || partition.points.at(edges.at(j).points[0]) || partition.points.at(edges.at(j).points[1])) && rand() > 0.5*RAND_MAX){
                        min = j;
                    }
                }
                if (partitionWeight + edges.at(j).weight <= setsWeight*sqrt(k+1)/sqrt(n-(n/t)*i) && (k == 1 || partition.points.at(edges.at(j).points[0]) || partition.points.at(edges.at(j).points[1]))){
                    candidates.push_back(j);
                }
            }

            int min_index = partition.points.at(edges.at(min).points[0]) == 1 ? 0 : 1;
            if (candidates.size() > 0 && k > 1){
                min = candidates.at(0);
                min_index = partition.points.at(edges.at(candidates.at(0)).points[0]) == 1 ? 0 : 1;
                for (int& j : candidates) {
                        int indexlookup = partition.points.at(edges.at(j).points[0]) == 1 ? 0 : 1;
                        if (bfs_table.at(edges.at(j).points[indexlookup]) < bfs_table.at(edges.at(min).points[min_index])){
                            min = j;
                            min_index = indexlookup;
                        } else {
                            if (bfs_table.at(edges.at(j).points[indexlookup]) == bfs_table.at(edges.at(min).points[min_index]) && edges.at(j).weight < edges.at(min).weight){
                                min = j;
                                min_index = indexlookup;
                            } else {
                                if (bfs_table.at(edges.at(j).points[indexlookup]) == bfs_table.at(edges.at(min).points[min_index]) && edges.at(j).weight == edges.at(min).weight && rand() < 0.5*RAND_MAX){
                                    min = j;
                                    min_index = indexlookup;
                                }
                        }
                    }
                }
            }

            if (k == 1){
                bfs_table.at(edges.at(min).points[0]) = 0;
                bfs_table.at(edges.at(min).points[1]) = 0;
            } else {
                bfs_table.at(edges.at(min).points[1 - min_index]) = bfs_table.at(edges.at(min).points[min_index]) + 1;
            }
            
            //Add selected edge to the partition
            partition.points.at(edges.at(min).points[0]) = 1;
            partition.points.at(edges.at(min).points[1]) = 1;

            //Update partition weight
            partitionWeight += edges.at(min).weight;
            weights.push_back(partitionWeight*sqrt(n-(n/t)*i)/setsWeight);

            // Remove from the weight of edges the weight of sets that intersect the selected edge
            for (int j = 0; j < ss.sets.size(); j++) {
                if (intersects(edges.at(min),ss.sets.at(j)) and !intersect_partition.at(j)){
                    intersect_partition.at(j) = 1;
                    for (Edge& e : edges) {
                        if (intersects(e,ss.sets.at(j))){
                            e.weight -= 1 << ss.sets.at(j).weight;
                            e.next_ite_weight += (1 << ss.sets.at(j).weight);
                        }
                    }
                    // Double set weight for next iteration
                    ss.sets.at(j).increase();
                }
            }


            //Remove edges that connect two points that have been added to the partition by different partitions
            edges.erase(edges.begin() + min);
            for (auto j = edges.size() - 1; j > 0; j--) {
                if (partition.points.at(edges.at(j).points[0]) && partition.points.at(edges.at(j).points[1])){
                    edges.erase(edges.begin() + j);

                }
            }
        }

        // Remove all edges linking to points added to the previous partition
        for (auto j = edges.size() - 1; j > 0; j--) {
                if (partition.points.at(edges.at(j).points[0]) || partition.points.at(edges.at(j).points[1])){
                    edges.erase(edges.begin() + j);
                }
            }

        // Store results
        res.sets.push_back(partition);
        for (int j = 0; j < intersect_partition.size();j++){
            histogram.at(j) += intersect_partition.at(j);
        }

        //Recompute edge weights
        for (Edge& e : edges) {
            e.weight = e.next_ite_weight;
        }
    }

    return Result(res,histogram,weights);
}

Result mwu_min_bfs(SetSystem ss, int t){
    
    const int n = ss.points.size();
    const int d = ss.points.at(0).coordinates.size();
    const int m = ss.sets.size();

    vector<int> histogram(m,0);
    vector<float> weights;


    //Initialize edges
    vector<Edge> edges;
    for (int i = 0; i < n - 1; i++){
        for (int j = i + 1; j < n; j++){
            edges.push_back(Edge(i, j));
        }
    }

    vector<int> bfs_table(n,t+1);

    //reset sets weight
    for (Set& s : ss.sets) {
        s.resetWeight();
    }

    //Recompute edge weights
    for (Edge& e : edges) {
        e.weight = 0;
        e.next_ite_weight = 0;
            for (Set& s : ss.sets){
                if (intersects(e,s)){
                    e.weight += (1 << s.weight);
                    e.next_ite_weight += (1 << s.weight);
                }
            }
    }


    vector<Set> s;
    SetSystem res = SetSystem(ss.points, s);

    vector<int> candidates;

    for (int i = 0 ; i < t; i++){
        //cout << "\nPartition " << i+1 << "\n";
        Set partition = Set(n);

        int setsWeight = 0;
        for (Set& s : ss.sets) {
            setsWeight += 1 << s.weight;
        }
        float partitionWeight = 0.;

        // Find a random edge from all edges enforcing rate with respect to bfs distance
        vector<bool> intersect_partition(m,0);
        for (int k = 1 ; k < n / t; k++){
            candidates.clear();
            int min = 0;
            for (auto j = 0; j < edges.size(); j++) {
                // Replace previous minimum if the weight is strictly smaller in all valid cases (no restriction if the partition has no edge yet, connectiveness otherwise)
                if (edges.at(j).weight < edges.at(min).weight && (k == 1 || partition.points.at(edges.at(j).points[0]) || partition.points.at(edges.at(j).points[1]))){
                    min = j;
                } else {
                    // If the algorithm finds a valid edge with the same weight it swaps it with probability 1/2 to avoid obteining the same result for every algo because of order
                    if (edges.at(j).weight == edges.at(min).weight && (k == 1 || partition.points.at(edges.at(j).points[0]) || partition.points.at(edges.at(j).points[1])) && rand() > 0.5*RAND_MAX){
                        min = j;
                    }
                }
                if (partitionWeight + edges.at(j).weight <= setsWeight*sqrt(k+1)/sqrt(n-(n/t)*i) && (k == 1 || partition.points.at(edges.at(j).points[0]) || partition.points.at(edges.at(j).points[1]))){
                    candidates.push_back(j);
                }
            }

            int min_index = partition.points.at(edges.at(min).points[0]) == 1 ? 0 : 1;
            if (candidates.size() > 0 && k > 1){
                min = candidates.at(0);
                min_index = partition.points.at(edges.at(candidates.at(0)).points[0]) == 1 ? 0 : 1;
                for (int& j : candidates) {
                        int indexlookup = partition.points.at(edges.at(j).points[0]) == 1 ? 0 : 1;
                        if (edges.at(j).weight < edges.at(min).weight){
                            min = j;
                            min_index = indexlookup;
                        } else {
                            if (edges.at(j).weight == edges.at(min).weight && bfs_table.at(edges.at(j).points[indexlookup]) < bfs_table.at(edges.at(min).points[min_index])){
                                min = j;
                                min_index = indexlookup;
                            } else {
                                if (bfs_table.at(edges.at(j).points[indexlookup]) == bfs_table.at(edges.at(min).points[min_index]) && edges.at(j).weight == edges.at(min).weight && rand() < 0.5*RAND_MAX){
                                    min = j;
                                    min_index = indexlookup;
                                }
                        }
                    }
                }
            }

            if (k == 1){
                bfs_table.at(edges.at(min).points[0]) = 0;
                bfs_table.at(edges.at(min).points[1]) = 0;
            } else {
                bfs_table.at(edges.at(min).points[1 - min_index]) = bfs_table.at(edges.at(min).points[min_index]) + 1;
            }
            
            //Add selected edge to the partition
            partition.points.at(edges.at(min).points[0]) = 1;
            partition.points.at(edges.at(min).points[1]) = 1;

            //Update partition weight
            partitionWeight += edges.at(min).weight;
            weights.push_back(partitionWeight*sqrt(n-(n/t)*i)/setsWeight);

            // Remove from the weight of edges the weight of sets that intersect the selected edge
            for (int j = 0; j < ss.sets.size(); j++) {
                if (intersects(edges.at(min),ss.sets.at(j)) and !intersect_partition.at(j)){
                    intersect_partition.at(j) = 1;
                    for (Edge& e : edges) {
                        if (intersects(e,ss.sets.at(j))){
                            e.weight -= 1 << ss.sets.at(j).weight;
                            e.next_ite_weight += (1 << ss.sets.at(j).weight);
                        }
                    }
                    // Double set weight for next iteration
                    ss.sets.at(j).increase();
                }
            }


            //Remove edges that connect two points that have been added to the partition by different partitions
            edges.erase(edges.begin() + min);
            for (auto j = edges.size() - 1; j > 0; j--) {
                if (partition.points.at(edges.at(j).points[0]) && partition.points.at(edges.at(j).points[1])){
                    edges.erase(edges.begin() + j);

                }
            }
        }

        // Remove all edges linking to points added to the previous partition
        for (auto j = edges.size() - 1; j > 0; j--) {
                if (partition.points.at(edges.at(j).points[0]) || partition.points.at(edges.at(j).points[1])){
                    edges.erase(edges.begin() + j);
                }
            }

        // Store results
        res.sets.push_back(partition);
        for (int j = 0; j < intersect_partition.size();j++){
            histogram.at(j) += intersect_partition.at(j);
        }

        //Recompute edge weights
        for (Edge& e : edges) {
            e.weight = e.next_ite_weight;
        }
    }

    return Result(res,histogram,weights);
}
