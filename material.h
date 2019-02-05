#ifndef MATERIAL_H
#define MATERIAL_H

struct isotropic_material{
    double m_C11;
    double m_C12;
    double m_C33;
    double m_rho;
    //TODO: revisar la el valor de C33 en notacion de voight
    isotropic_material(double lamda,double mu,double rho):m_C11(lamda+2*mu),m_C12(lamda),m_C33(0.5*mu),m_rho(rho){}
};

#endif // MATERIAL_H
