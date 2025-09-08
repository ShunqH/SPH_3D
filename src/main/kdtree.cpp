#include <vector>  
#include <iostream>
#include <algorithm>   
#include <memory>
#include "particles.h"
#include "kdtree.h"

using namespace std;

// Constructor 
Tree::Tree(Particles& pts):pts(pts){
    int n = pts.endid;          // n of particles

    // indices of particles
    indices.resize(n);
    for (int i=0; i<n; i++){
        indices[i] = i; 
    }
    
    // build up tree
    root = build(0, n, 0);
}

// build function
unique_ptr<Tree::Node> Tree::build(int left, int right, int depth){
    // end of tree
    if (left>=right){
        return nullptr; 
    }

    // devide by axis
    int ax = depth%3; 
    int mid = partition(left, right, ax);

    // setup this node 
    // auto node = unique_ptr<Node>(new Node(indices[mid], ax))
    auto node = make_unique<Node>(indices[mid], ax);
    node->left = build(left, mid, depth+1); 
    node->right = build(mid+1, right, depth+1); 
    return node; 
}

vector<int> Tree::search(double x0, double y0, double z0, double r){
    vector<int> results; 
    double r2 = r*r; 
    search_rec(root.get(), x0, y0, z0, r2, results); 
    return results; 
}

void Tree::search_rec(Node* node, const double& x0, const double& y0, const double& z0, const double& r2, std::vector<int>& results){
    if (node == nullptr) return;        // end of tree
    int id = node->index;               // id of this node
    // distance from node to target
    double dist2 = (pts.x1[id]-x0)*(pts.x1[id]-x0) 
                 + (pts.x2[id]-y0)*(pts.x2[id]-y0) 
                 + (pts.x3[id]-z0)*(pts.x3[id]-z0);     
    if (dist2<=r2){                     // this node is inside of range
        results.push_back(id); 
    }
    int ax = node->axis; 
    double p, q; 
    // |p-q|: distance from target to the divided ax
    switch (ax) {
        case 0: 
            q = x0; 
            p = pts.x1[id];
            break; 
        case 1: 
            q = y0; 
            p = pts.x2[id];
            break; 
        case 2: 
            q = z0; 
            p = pts.x3[id]; 
            break; 
    }
    if ((p-q)*(p-q)<=r2){       // both side 
        search_rec(node->left.get(), x0, y0, z0, r2, results); 
        search_rec(node->right.get(), x0, y0, z0, r2, results); 
    }else if (q<p){             // left side
        search_rec(node->left.get(), x0, y0, z0, r2, results); 
    }else{                      // right side 
        search_rec(node->right.get(), x0, y0, z0, r2, results); 
    }
    return; 
}

int Tree::partition(int left, int right, int ax){
    auto cmp = [&](int a, int b) {
        if (ax == 0) return pts.x1[a] < pts.x1[b];
        if (ax == 1) return pts.x2[a] < pts.x2[b];
        return pts.x3[a] < pts.x3[b];
    };
    int mid = left + (right - left) / 2;
    std::nth_element(indices.begin() + left, indices.begin() + mid, 
                    indices.begin() + right, cmp);
    return mid;
}