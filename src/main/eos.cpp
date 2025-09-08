#include <cmath>
#include <stdexcept>  
#include "eos.h"

using namespace std;

EoS_adi::EoS_adi(double gamma): gamma(gamma) {}; 

double EoS_adi::pressure(double den, double ene) const {
    return (gamma-1.)*den*ene; 
}

double EoS_adi::cs(double den, double ene) const {
    return sqrt(gamma*(gamma-1.)*ene); 
}

EoS_iso::EoS_iso(double sound_speed): sound_speed(sound_speed) {};

double EoS_iso::pressure(double den, double ene) const {
    return sound_speed*sound_speed*den; 
}

double EoS_iso::cs(double den, double ene) const {
    return sound_speed; 
}

unique_ptr<EoS> set_eos(int type, double parm){
    if (type == 0){
        return make_unique<EoS_adi> (parm); 
    }else if (type == 1){
        return make_unique<EoS_iso> (parm); 
    }else{
        throw invalid_argument("unknow equation of state. ");
    }
}