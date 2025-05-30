#ifndef EOS_H
#define EOS_H

#include <memory>

class EoS{
public: 
    virtual double pressure(double den, double ene) const = 0; 
    virtual double cs(double den, double ene) const = 0; 
    virtual ~EoS() {}; 
}; 

class EoS_adi: public EoS{
private:
    double gamma; 
public: 
    EoS_adi(double gamma); 
    double pressure(double den, double ene) const override;
    double cs(double den, double ene) const override;
}; 

class EoS_iso: public EoS{
private:
    double sound_speed; 
public: 
    EoS_iso(double sound_speed); 
    double pressure(double den, double ene) const override;
    double cs(double den, double ene) const override;
}; 

std::unique_ptr<EoS> set_eos(int type, double parm); 

#endif