#ifndef DMCALC_H
#define DMCALC_H
#include "TF1.h"

class DMCalc {
    public:
        DMCalc();
        ~DMCalc();

        double sigma0(double m_w);
        double dRdEr0inf(double m_w, double Er);
        double dRdEr0esc(double m_w, double Er);
        double dRdEreinf(double m_w, double Er);
        double dRdEreesc(double m_w, double Er);
        double FormFactor(double Er);
        double DetResponse(double Evis, double Er);

    private:
        double c;
        double c2;
        double Avogadro;
        double rho0;
        double sigma_n;
        double v0;
        double A;
        double m_w;
        double m_N;
        double m_n;
        double v_esc;
        double v_e;
        double quench;
};

#endif
