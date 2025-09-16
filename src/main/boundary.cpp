#include <stdexcept>  
#include <iostream> 
#include <cmath>
#include "boundary.h" 
#include "particles.h"

using namespace std;


// Free boundaries for 6 directions
BoundaryX1L_Free::BoundaryX1L_Free(): Boundary() {}
void BoundaryX1L_Free::SetBoundary(Particles& pts){ return; }

BoundaryX1R_Free::BoundaryX1R_Free(): Boundary() {}
void BoundaryX1R_Free::SetBoundary(Particles& pts){ return; }

BoundaryX2L_Free::BoundaryX2L_Free(): Boundary() {}
void BoundaryX2L_Free::SetBoundary(Particles& pts){ return; }

BoundaryX2R_Free::BoundaryX2R_Free(): Boundary() {}
void BoundaryX2R_Free::SetBoundary(Particles& pts){ return; }

BoundaryX3L_Free::BoundaryX3L_Free(): Boundary() {}
void BoundaryX3L_Free::SetBoundary(Particles& pts){ return; }

BoundaryX3R_Free::BoundaryX3R_Free(): Boundary() {}
void BoundaryX3R_Free::SetBoundary(Particles& pts){ return; }

// Periodic boundaries for 6 directions
// Copying a layer from the left activate zone to the right side as ghost particles. The layer thickness is dx. 
BoundaryX1L_Periodic::BoundaryX1L_Periodic(double xl, double xr, double dx): Boundary(xl, xr, dx) {}
void BoundaryX1L_Periodic::SetBoundary(Particles& pts){
    int n = pts.ghostid; 
    for (int i=0; i<n; i++){
        if (pts.x1[i] < xl - dx){
            pts.x1[i] += xr - xl;
        }else if (pts.x1[i] <= xl){
            Particle pt = pts.extract(i); 
            pt.status = -1; 
            pts.add_particle(pt);
            pts.x1[i] += xr - xl;
        }else if ((pts.x1[i] >= xr - dx)&&(pts.x1[i] <= xr)){
            Particle pt = pts.extract(i); 
            pt.x1 = pt.x1 - (xr - xl); 
            pt.status = -1; 
            pts.add_particle(pt);  
        } 
    }
    return ;
}

BoundaryX1R_Periodic::BoundaryX1R_Periodic(double xl, double xr, double dx): Boundary(xl, xr, dx) {}
void BoundaryX1R_Periodic::SetBoundary(Particles& pts){
    int n = pts.ghostid; 
    for (int i=0; i<n; i++){
        if (pts.x1[i] > xr + dx){
            pts.x1[i] -= xr - xl;
        }else if (pts.x1[i] >= xr){
            Particle pt = pts.extract(i);
            pt.status = -1; 
            pts.add_particle(pt); 
            pts.x1[i] -= xr - xl;
        }else if ((pts.x1[i] <= xl + dx)&&(pts.x1[i] >= xl)){
            Particle pt = pts.extract(i); 
            pt.x1 = pt.x1 + (xr - xl); 
            pt.status = -1; 
            pts.add_particle(pt);  
        } 
    }
    return ;
}

BoundaryX2L_Periodic::BoundaryX2L_Periodic(double xl, double xr, double dx): Boundary(xl, xr, dx) {}
void BoundaryX2L_Periodic::SetBoundary(Particles& pts){
    int n = pts.ghostid; 
    for (int i=0; i<n; i++){
        if (pts.x2[i] < xl - dx){
            pts.x2[i] += xr - xl;
        }else if (pts.x2[i] <= xl){
            Particle pt = pts.extract(i); 
            pt.status = -1; 
            pts.add_particle(pt);
            pts.x2[i] += xr - xl;
        }else if ((pts.x2[i] >= xr - dx)&&(pts.x2[i] <= xr)){
            Particle pt = pts.extract(i); 
            pt.x2 = pt.x2 - (xr - xl); 
            pt.status = -1; 
            pts.add_particle(pt);  
        } 
    }
    return ;
}

BoundaryX2R_Periodic::BoundaryX2R_Periodic(double xl, double xr, double dx): Boundary(xl, xr, dx) {}
void BoundaryX2R_Periodic::SetBoundary(Particles& pts){
    // cout<<"set bd"<<endl; 
    int n = pts.ghostid; 
    for (int i=0; i<n; i++){
        if (pts.x2[i] > xr + dx){
            pts.x2[i] -= xr - xl;
        }else if (pts.x2[i] >= xr){
            Particle pt = pts.extract(i);
            pt.status = -1; 
            pts.add_particle(pt); 
            pts.x2[i] -= xr - xl;
        }else if ((pts.x2[i] <= xl + dx)&&(pts.x2[i] >= xl)){
            Particle pt = pts.extract(i); 
            pt.x2 = pt.x2 + (xr - xl); 
            pt.status = -1; 
            pts.add_particle(pt);  
        } 
    }
    return ;
}

BoundaryX3L_Periodic::BoundaryX3L_Periodic(double xl, double xr, double dx): Boundary(xl, xr, dx) {}
void BoundaryX3L_Periodic::SetBoundary(Particles& pts){
    int n = pts.ghostid; 
    for (int i=0; i<n; i++){
        if (pts.x3[i] < xl - dx){
            pts.x3[i] += xr - xl;
        }else if (pts.x3[i] <= xl){
            Particle pt = pts.extract(i); 
            pt.status = -1; 
            pts.add_particle(pt);
            pts.x3[i] += xr - xl;
        }else if ((pts.x3[i] >= xr - dx)&&(pts.x3[i] <= xr)){
            Particle pt = pts.extract(i); 
            pt.x3 = pt.x3 - (xr - xl); 
            pt.status = -1; 
            pts.add_particle(pt);  
        } 
    }
    return ;
}

BoundaryX3R_Periodic::BoundaryX3R_Periodic(double xl, double xr, double dx): Boundary(xl, xr, dx) {}
void BoundaryX3R_Periodic::SetBoundary(Particles& pts){
    int n = pts.ghostid; 
    for (int i=0; i<n; i++){
        if (pts.x3[i] > xr + dx){
            pts.x3[i] -= xr - xl;
        }else if (pts.x3[i] >= xr){
            Particle pt = pts.extract(i);
            pt.status = -1; 
            pts.add_particle(pt); 
            pts.x3[i] -= xr - xl;
        }else if ((pts.x3[i] <= xl + dx)&&(pts.x3[i] >= xl)){
            Particle pt = pts.extract(i); 
            pt.x3 = pt.x3 + (xr - xl); 
            pt.status = -1; 
            pts.add_particle(pt);  
        } 
    }
    return ;
}

// Static boundaries for 6 directions
// Mirror a layer from the left activate zone to the left side as ghost particles. The layer thickness is dx. 
BoundaryX1L_Static::BoundaryX1L_Static(double xl, double xr, double dx): Boundary(xl, xr, dx) {}
void BoundaryX1L_Static::SetBoundary(Particles& pts){
    int n = pts.ghostid; 
    for (int i=0; i<n; i++){
        if (pts.x1[i] < xl - dx){
            cerr << "particle cross through static boundary at x, left."
                 << endl;
        }else if (pts.x1[i] < xl){
            Particle pt = pts.extract(i); 
            pt.status = -1; 
            pts.add_particle(pt);
            pts.x1[i] = 2*xl - pts.x1[i];
            pts.vel1[i] = -pts.vel1[i]; 
        }else if (pts.x1[i] <= xl + dx){
            Particle pt = pts.extract(i); 
            pt.x1 = 2*xl - pt.x1; 
            pt.vel1 = -pt.vel1;
            pt.status = -1; 
            pts.add_particle(pt);  
        } 
    }
    return ;
}

BoundaryX1R_Static::BoundaryX1R_Static(double xl, double xr, double dx): Boundary(xl, xr, dx) {}
void BoundaryX1R_Static::SetBoundary(Particles& pts){
    int n = pts.ghostid; 
    for (int i=0; i<n; i++){
        if (pts.x1[i] > xr + dx){
            cerr << "particle cross through static boundary at x, right."
                 << endl;
        }else if (pts.x1[i] > xr){
            Particle pt = pts.extract(i); 
            pt.status = -1; 
            pts.add_particle(pt);
            pts.x1[i] = 2*xr - pts.x1[i];
            pts.vel1[i] = - pts.vel1[i]; 
        }else if (pts.x1[i] >= xr - dx){
            Particle pt = pts.extract(i); 
            pt.x1 = 2*xr - pt.x1; 
            pt.vel1 = -pt.vel1;
            pt.status = -1; 
            pts.add_particle(pt);  
        } 
    }
    return ;
}

BoundaryX2L_Static::BoundaryX2L_Static(double xl, double xr, double dx): Boundary(xl, xr, dx) {}
void BoundaryX2L_Static::SetBoundary(Particles& pts){
    int n = pts.ghostid; 
    for (int i=0; i<n; i++){
        if (pts.x2[i] < xl - dx){
            cerr << "particle cross through static boundary at y, left."
                 << endl;
        }else if (pts.x2[i] < xl){
            Particle pt = pts.extract(i); 
            pt.status = -1; 
            pts.add_particle(pt);
            pts.x2[i] = 2*xl - pts.x2[i];
            pts.vel2[i] = -pts.vel2[i]; 
        }else if (pts.x2[i] <= xl + dx){
            Particle pt = pts.extract(i); 
            pt.x2 = 2*xl - pt.x2; 
            pt.vel2 = -pt.vel2;
            pt.status = -1; 
            pts.add_particle(pt);  
        } 
    }
    return ;
}

BoundaryX2R_Static::BoundaryX2R_Static(double xl, double xr, double dx): Boundary(xl, xr, dx) {}
void BoundaryX2R_Static::SetBoundary(Particles& pts){
    int n = pts.ghostid; 
    for (int i=0; i<n; i++){
        if (pts.x2[i] > xr + dx){
            cerr << "particle cross through static boundary at y, right."
                 << endl;
        }else if (pts.x2[i] > xr){
            Particle pt = pts.extract(i); 
            pt.status = -1; 
            pts.add_particle(pt);
            pts.x2[i] = 2*xr - pts.x2[i];
            pts.vel2[i] = - pts.vel2[i]; 
        }else if (pts.x2[i] >= xr - dx){
            Particle pt = pts.extract(i); 
            pt.x2 = 2*xr - pt.x2; 
            pt.vel2 = -pt.vel2;
            pt.status = -1; 
            pts.add_particle(pt);  
        } 
    }
    return ;
}

BoundaryX3L_Static::BoundaryX3L_Static(double xl, double xr, double dx): Boundary(xl, xr, dx) {}
void BoundaryX3L_Static::SetBoundary(Particles& pts){
    int n = pts.ghostid; 
    for (int i=0; i<n; i++){
        if (pts.x3[i] < xl - dx){
            cerr << "particle cross through static boundary at z, left."
                 << endl;
        }else if (pts.x3[i] < xl){
            Particle pt = pts.extract(i); 
            pt.status = -1; 
            pts.add_particle(pt);
            pts.x3[i] = 2*xl - pts.x3[i];
            pts.vel3[i] = -pts.vel3[i]; 
        }else if (pts.x3[i] <= xl + dx){
            Particle pt = pts.extract(i); 
            pt.x3 = 2*xl - pt.x3; 
            pt.vel3 = -pt.vel3;
            pt.status = -1; 
            pts.add_particle(pt);  
        } 
    }
    return ;
}

BoundaryX3R_Static::BoundaryX3R_Static(double xl, double xr, double dx): Boundary(xl, xr, dx) {}
void BoundaryX3R_Static::SetBoundary(Particles& pts){
    int n = pts.ghostid; 
    for (int i=0; i<n; i++){
        if (pts.x3[i] > xr + dx){
            cerr << "particle cross through static boundary at z, right."
                 << endl;
        }else if (pts.x3[i] > xr){
            Particle pt = pts.extract(i); 
            pt.status = -1; 
            pts.add_particle(pt);
            pts.x3[i] = 2*xr - pts.x3[i];
            pts.vel3[i] = - pts.vel3[i]; 
        }else if (pts.x3[i] >= xr - dx){
            Particle pt = pts.extract(i); 
            pt.x3 = 2*xr - pt.x3; 
            pt.vel3 = -pt.vel3;
            pt.status = -1; 
            pts.add_particle(pt);  
        } 
    }
    return ;
}

unique_ptr<Boundary> Set_BDX1L(int type, double xl, double xr, double dx){
    if (type == 0){
        return make_unique<BoundaryX1L_Free> (); 
    }else if (type == 1){
        return make_unique<BoundaryX1L_Periodic> (xl, xr, dx); 
    }else if (type == 2){
        return make_unique<BoundaryX1L_Static> (xl, xr, dx); 
    }else{
        throw invalid_argument("unknow boundary condition for x1 left. ");
    }
}

unique_ptr<Boundary> Set_BDX1R(int type, double xl, double xr, double dx){
    if (type == 0){
        return make_unique<BoundaryX1R_Free> (); 
    }else if (type == 1){
        return make_unique<BoundaryX1R_Periodic> (xl, xr, dx); 
    }else if (type == 2){
        return make_unique<BoundaryX1R_Static> (xl, xr, dx); 
    }else{
        throw invalid_argument("unknow boundary condition for x1 right. ");
    }
}

unique_ptr<Boundary> Set_BDX2L(int type, double xl, double xr, double dx){
    if (type == 0){
        return make_unique<BoundaryX2L_Free> (); 
    }else if (type == 1){
        return make_unique<BoundaryX2L_Periodic> (xl, xr, dx); 
    }else if (type == 2){
        return make_unique<BoundaryX2L_Static> (xl, xr, dx); 
    }else{
        throw invalid_argument("unknow boundary condition for x2 left. ");
    }
}

unique_ptr<Boundary> Set_BDX2R(int type, double xl, double xr, double dx){
    if (type == 0){
        return make_unique<BoundaryX2R_Free> (); 
    }else if (type == 1){
        return make_unique<BoundaryX2R_Periodic> (xl, xr, dx); 
    }else if (type == 2){
        return make_unique<BoundaryX2R_Static> (xl, xr, dx); 
    }else{
        throw invalid_argument("unknow boundary condition for x2 right. ");
    }
}

unique_ptr<Boundary> Set_BDX3L(int type, double xl, double xr, double dx){
    if (type == 0){
        return make_unique<BoundaryX3L_Free> (); 
    }else if (type == 1){
        return make_unique<BoundaryX3L_Periodic> (xl, xr, dx); 
    }else if (type == 2){
        return make_unique<BoundaryX3L_Static> (xl, xr, dx); 
    }else{
        throw invalid_argument("unknow boundary condition for x3 left. ");
    }
}

unique_ptr<Boundary> Set_BDX3R(int type, double xl, double xr, double dx){
    if (type == 0){
        return make_unique<BoundaryX3R_Free> (); 
    }else if (type == 1){
        return make_unique<BoundaryX3R_Periodic> (xl, xr, dx); 
    }else if (type == 2){
        return make_unique<BoundaryX3R_Static> (xl, xr, dx); 
    }else{
        throw invalid_argument("unknow boundary condition for x3 right. ");
    }
}