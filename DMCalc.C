#include <iostream>
#include <map>
#include <string>
#include "TMath.h"
#include "DMCalc.h"

DMCalc::DMCalc(const char* Material="Ar", bool IsSI=kTRUE, int RateMode=4) {
    // initialization
    const int mattypes= 2;
    std::string mat[mattypes]= {"Ar", "Xe"};
    double tmp_A[mattypes]= {40, 132};
    for (int i= 0; i< mattypes; ++i) {
        matdb[mat[i]]= tmp_A[i];
    }
    if ( matdb[Material]==0 ) {
        std::cout << "WARNING: Undefined material: " << Material << std::endl;
        std::cout << "Set to default material 'Ar40'." << std::endl;
        A= 40;
    } else {
        std::cout << "Target material: " << Material << std::endl;
        A= matdb[Material];
    }
    IsSpinIndependent= IsSI;

    c       = 2.99792458e10;// speed of light
    c2      = c*c;      
    Avogadro= 6.0221413e23; // Avogadro number
    rho0    = 0.3e9/c2;     // dark matter particle density: eV/c2/cm3
    v0      = 2.2e7;        // cm/s
    sigma_n = 1e-42;        // WIMP-nucleon cross section
    m_N     = A*9.315e8/c2; // eV/c2
    m_n     = 9.315e8/c2;   // eV/c2
    v_esc   = 6e7;          // escape velocity: cm/s
    v_e     = 2.44e7;       // velocity of earth: cm/s

    mode= RateMode;
    switch(mode) {
        case 1:
            std::cout << "mode 1: V_earth=0, V_escape=inf." << std::endl;
            break;
        case 2:
            std::cout << "mode 2: V_earth=0, V_escape=" <<  v_esc << " cm/s." << std::endl;
            break;
        case 3:
            std::cout << "mode 3: V_earth=" << v_e << " cm/s, V_escape=inf." << std::endl;
            break;
        case 4:
            std::cout << "mode 4: V_earth=" << v_e << " cm/s, V_escape=" << v_esc << " cm/s." << std::endl;
            break;
        default:
            std::cout << "WARNING: Undefined mode " << mode << std::endl;
            std::cout << "Set to mode 4: " << "V_earth=" << v_e << ", V_escape=" << v_esc << std::endl;
    }
}

DMCalc::~DMCalc() {
}

double DMCalc::sigma0 (double m_w) {
    double tmp_sigma0;
    if (IsSpinIndependent) {
        double mu_N= m_w*m_N/(m_w + m_N);
        double mu_n= m_w*m_n/(m_w + m_n);
        tmp_sigma0= sigma_n*mu_N*mu_N*A*A/mu_n/mu_n;
    } else {
        tmp_sigma0= 0;
    }

    return tmp_sigma0;
}

double DMCalc::FormFactor(double Er) {
    double tmp_F;
    // The SI form factor is essentially a Fourier transform of the mass distribution of the nucleus. A reasonably accurate approximation is the Helm form factor
    if (IsSpinIndependent) {
        // q(MeVc^-1)=[2*0.932(GeVc^-2)AEr(keV)]^0.5, hbarc=197.3MeVfm
        // qr(dimensionless)= 6.92e-3*sqrt(A*Er)(a*A^1/3+b) approxi. equal to 6.92e-3*sqrt(A*Er)(1.14*A^1/3)
        // Er in keV, a/b in fm
        double qr =6.92e-3*TMath::Sqrt(A*Er/1000.)*(1.14*TMath::Power(A,1.0/3.0));
        double qs =6.92e-3*TMath::Sqrt(A*Er/1000.)*0.9;
        tmp_F= 3*(TMath::Sin(qr)/(qr*qr*qr)-TMath::Cos(qr)/(qr*qr))*TMath::Exp(-qs*qs/2);
    } else {
        tmp_F= 0;
    }

    return tmp_F;
}

double DMCalc::dRdEr(double m_w, double Er) {
    double F= FormFactor(Er);
    switch(mode) {
        case 1: 
            return F*F*dRdEr0inf(m_w, Er);
        case 2:
            return F*F*dRdEr0esc(m_w, Er);
        case 3:
            return F*F*dRdEreinf(m_w, Er);
        case 4:
            return F*F*dRdEreesc(m_w, Er);
        default:
            return F*F*dRdEreesc(m_w, Er);
    }
}

double DMCalc::dRdEr0inf(double m_w, double Er) {
    double R0= 2*Avogadro*rho0*sigma0(m_w)*v0/TMath::Sqrt(TMath::Pi())/A/m_w;
    double r= 4*m_w*m_N/(m_w+m_N)/(m_w+m_N);
    double E0 = (0.5*(m_w*c2*1.79e-36)*(v0/100)*(v0/100))/(1.609e-19);
//    double E0= m_w*v0*v0*c2/2;
    double tmp_dRdEr0inf= R0*TMath::Exp(-Er/E0/r)/E0/r;

    return tmp_dRdEr0inf;
}

double DMCalc::dRdEr0esc(double m_w, double Er) {
    double k1_k0= TMath::Erf(v_esc/v0)-2*v_esc*TMath::Exp(-v_esc*v_esc/v0/v0)/TMath::Sqrt(TMath::Pi())/v0;
    double R0= 2*Avogadro*rho0*sigma0(m_w)*v0/TMath::Sqrt(TMath::Pi())/A/m_w;
    double r= 4*m_w*m_N/(m_w+m_N)/(m_w+m_N);
    double E0 = (0.5*(m_w*c2*1.79e-36)*(v0/100)*(v0/100))/(1.609e-19);
//    double E0= m_w*v0*v0*c2/2;
    double tmp_dRdEr0esc= (dRdEr0inf(m_w, Er) - R0*TMath::Exp(-v_esc*v_esc/v0/v0)/E0/r)/k1_k0;

    return tmp_dRdEr0esc;
}

double DMCalc::dRdEreinf(double m_w, double Er) {
    double R0= 2*Avogadro*rho0*sigma0(m_w)*v0/TMath::Sqrt(TMath::Pi())/A/m_w;
    double r= 4*m_w*m_N/(m_w+m_N)/(m_w+m_N);
    double E0 = (0.5*(m_w*c2*1.79e-36)*(v0/100)*(v0/100))/(1.609e-19);
//    double E0= m_w*v0*v0*c2/2;
    double vmin= TMath::Sqrt(Er/E0/r)*v0;
    double Part_I= R0*TMath::Sqrt(TMath::Pi())*v0*1000/E0/r/4/v_e;
    double Part_II= TMath::Erf((vmin+v_e)/v0) - TMath::Erf((vmin-v_e)/v0);
    double tmp_dRdEreinf= Part_I*Part_II;

    return tmp_dRdEreinf;
}

double DMCalc::dRdEreesc(double m_w, double Er) {
    double k1_k0= TMath::Erf(v_esc/v0)-2*v_esc*TMath::Exp(-v_esc*v_esc/v0/v0)/TMath::Sqrt(TMath::Pi())/v0;
    double R0= 2*Avogadro*rho0*sigma0(m_w)*v0/TMath::Sqrt(TMath::Pi())/A/m_w;
    double r= 4*m_w*m_N/(m_w+m_N)/(m_w+m_N);
    double E0 = (0.5*(m_w*c2*1.79e-36)*(v0/100)*(v0/100))/(1.609e-19);
//    double E0= m_w*v0*v0*c2/2;
    double tmp_dRdEreesc= (dRdEreinf(m_w, Er) - R0*TMath::Exp(-v_esc*v_esc/v0/v0)/E0/r)/k1_k0;

    return tmp_dRdEreesc;
}

double DMCalc::DetResponse(double Evis, double Er) {
    double tsigma= TMath::Sqrt(0.965 + 0.787*Evis + 5.33e-3*Evis*Evis);
    double tmp_DetResponse= TMath::Exp(-(Er-Evis)*(Er-Evis)/2/tsigma/tsigma)/tsigma/TMath::Sqrt(2*TMath::Pi());

    return tmp_DetResponse;
}
