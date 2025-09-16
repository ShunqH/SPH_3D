#include <vector>
#include <fstream>
#include <sstream>
#include <iomanip> 
#include <iostream>
#include <string>
#include "utils.h" 
#include "particles.h"

using namespace std;

void WriteParticles(Particles& pts, int step){
    int n = pts.ghostid; 

    // string Path = "./output/"; 
    ostringstream filename_stream;
    filename_stream << "output_" << setw(5) << setfill('0') << step;  
    string filename = filename_stream.str(); 

    ofstream outFile(filename, ios::binary);

    // write rows and cols
    outFile.write(reinterpret_cast<char*>(&n), sizeof(n));
    
    // write data
    outFile.write(reinterpret_cast<char*>(pts.x1.data()), n*sizeof(double));
    outFile.write(reinterpret_cast<char*>(pts.x2.data()), n*sizeof(double));
    outFile.write(reinterpret_cast<char*>(pts.x3.data()), n*sizeof(double));
    outFile.write(reinterpret_cast<char*>(pts.vel1.data()), n*sizeof(double));
    outFile.write(reinterpret_cast<char*>(pts.vel2.data()), n*sizeof(double));
    outFile.write(reinterpret_cast<char*>(pts.vel3.data()), n*sizeof(double));
    outFile.write(reinterpret_cast<char*>(pts.acc1.data()), n*sizeof(double));
    outFile.write(reinterpret_cast<char*>(pts.acc2.data()), n*sizeof(double));
    outFile.write(reinterpret_cast<char*>(pts.acc3.data()), n*sizeof(double));
    outFile.write(reinterpret_cast<char*>(pts.accu.data()), n*sizeof(double));
    outFile.write(reinterpret_cast<char*>(pts.den.data()), n*sizeof(double));
    outFile.write(reinterpret_cast<char*>(pts.mas.data()), n*sizeof(double));
    outFile.write(reinterpret_cast<char*>(pts.len.data()), n*sizeof(double));
    outFile.write(reinterpret_cast<char*>(pts.ene.data()), n*sizeof(double));
    outFile.write(reinterpret_cast<char*>(pts.status.data()), n*sizeof(int));

    outFile.close();
    return; 
}

double dot(const vector<double>& va, const vector<double>& vb){
    if (va.size()!=3 || vb.size()!=3){
        cerr << "vector size error in dot(va, vb), va.size = " << va.size()
             << ", vb.size = " << vb.size() << endl;
        return 0;
    }
    return va[0]*vb[0] + va[1]*vb[1] + va[2]*vb[2]; 
}

double power(double base, int order) {
    double result = 1.0;
    bool negative = false;

    if (order < 0) {
        negative = true;
        order = -order;
    }
    while (order > 0) {
        if (order % 2 == 1) {
            result *= base;
        }
        base *= base;
        order /= 2;
    }
    if (negative){
        return 1.0/result;
    }else{
        return result; 
    }
}
