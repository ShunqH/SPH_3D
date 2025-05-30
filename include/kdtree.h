#ifndef KDTREE_H
#define KDTREE_H
#include <vector>
#include <memory>
#include "particles.h"

class Tree{
public:
    // Constructor 
    Tree(Particles& pts); 
    
    //search function
    std::vector<int> search(double x0, double y0, double z0, double r); 

private:
    std::vector<int> indices; 
    Particles& pts;
    // node structure
    struct Node{
        int index;      // index of the particle
        int axis;       // based on which ax
        std::unique_ptr<Node> left;
        std::unique_ptr<Node> right;
    
        Node(int id, int ax) : index(id), axis(ax), left(nullptr), right(nullptr) {}
    }; 
    std::unique_ptr<Node> root; 
    std::unique_ptr<Node> build(int left, int right, int depth); 

    int partition(int left, int right, int axis); 

    void search_rec(Node* node, const double& x0, const double& y0, const double& z0, const double& r2, std::vector<int>& results);


}; 

#endif 