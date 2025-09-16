#ifndef BOUNDARY_H
#define BOUNDARY_H

#include <memory>

class Particles;

class Boundary{
protected:
    double xl, xr, dx;
public: 
    Boundary(double xl=0, double xr=0, double dx=0): 
        xl(xl), xr(xr), dx(dx) {}
    virtual void SetBoundary(Particles& pts) = 0; 
    virtual ~Boundary() {}; 
}; 

class BoundaryX1L_Free: public Boundary{
public: 
    BoundaryX1L_Free(); 
    void SetBoundary(Particles& pts) override; 
}; 

class BoundaryX1R_Free: public Boundary{
public: 
    BoundaryX1R_Free(); 
    void SetBoundary(Particles& pts) override; 
}; 

class BoundaryX2L_Free: public Boundary{
public: 
    BoundaryX2L_Free(); 
    void SetBoundary(Particles& pts) override; 
}; 

class BoundaryX2R_Free: public Boundary{
public: 
    BoundaryX2R_Free(); 
    void SetBoundary(Particles& pts) override; 
}; 

class BoundaryX3L_Free: public Boundary{
public: 
    BoundaryX3L_Free(); 
    void SetBoundary(Particles& pts) override; 
}; 

class BoundaryX3R_Free: public Boundary{
public: 
    BoundaryX3R_Free(); 
    void SetBoundary(Particles& pts) override; 
}; 


class BoundaryX1L_Periodic: public Boundary{
public: 
    BoundaryX1L_Periodic(double xl, double xr, double dx); 
    void SetBoundary(Particles& pts) override; 
}; 

class BoundaryX1R_Periodic: public Boundary{
public: 
    BoundaryX1R_Periodic(double xl, double xr, double dx); 
    void SetBoundary(Particles& pts) override; 
}; 

class BoundaryX2L_Periodic: public Boundary{
public: 
    BoundaryX2L_Periodic(double xl, double xr, double dx); 
    void SetBoundary(Particles& pts) override; 
}; 

class BoundaryX2R_Periodic: public Boundary{
public: 
    BoundaryX2R_Periodic(double xl, double xr, double dx); 
    void SetBoundary(Particles& pts) override; 
}; 

class BoundaryX3L_Periodic: public Boundary{
public: 
    BoundaryX3L_Periodic(double xl, double xr, double dx); 
    void SetBoundary(Particles& pts) override; 
}; 

class BoundaryX3R_Periodic: public Boundary{
public: 
    BoundaryX3R_Periodic(double xl, double xr, double dx); 
    void SetBoundary(Particles& pts) override; 
}; 


class BoundaryX1L_Static: public Boundary{
public: 
    BoundaryX1L_Static(double xl, double xr, double dx); 
    void SetBoundary(Particles& pts) override; 
}; 

class BoundaryX1R_Static: public Boundary{
public: 
    BoundaryX1R_Static(double xl, double xr, double dx); 
    void SetBoundary(Particles& pts) override; 
}; 

class BoundaryX2L_Static: public Boundary{
public: 
    BoundaryX2L_Static(double xl, double xr, double dx); 
    void SetBoundary(Particles& pts) override; 
}; 

class BoundaryX2R_Static: public Boundary{
public: 
    BoundaryX2R_Static(double xl, double xr, double dx); 
    void SetBoundary(Particles& pts) override; 
}; 

class BoundaryX3L_Static: public Boundary{
public: 
    BoundaryX3L_Static(double xl, double xr, double dx); 
    void SetBoundary(Particles& pts) override; 
}; 

class BoundaryX3R_Static: public Boundary{
public: 
    BoundaryX3R_Static(double xl, double xr, double dx); 
    void SetBoundary(Particles& pts) override; 
}; 

std::unique_ptr<Boundary> Set_BDX1L(int type, double xl=0, double xr=0, double dx=0); 
std::unique_ptr<Boundary> Set_BDX1R(int type, double xl=0, double xr=0, double dx=0); 
std::unique_ptr<Boundary> Set_BDX2L(int type, double xl=0, double xr=0, double dx=0); 
std::unique_ptr<Boundary> Set_BDX2R(int type, double xl=0, double xr=0, double dx=0); 
std::unique_ptr<Boundary> Set_BDX3L(int type, double xl=0, double xr=0, double dx=0); 
std::unique_ptr<Boundary> Set_BDX3R(int type, double xl=0, double xr=0, double dx=0); 


#endif