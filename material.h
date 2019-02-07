#ifndef MATERIAL_H
#define MATERIAL_H

struct isotropic_material{
    double m_C11;
    double m_C12;
    double m_C33;
    double m_rho;
    //TODO: revisar la el valor de C33 en notacion de voight
    isotropic_material(double E,double nu,double rho):m_C11(E*(1-nu)/(1+nu)/(1-2*nu)),m_C12(m_C11*nu/(1-nu)),m_C33(0.5*E/(nu+1)),m_rho(rho){}
};
#endif // MATERIAL_H
