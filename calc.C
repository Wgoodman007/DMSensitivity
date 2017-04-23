#include <iostream>
#include "TMath.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TApplication.h"
#include "TGraph.h"
#include "DMCalc.h"

using namespace std;

int main(int argc, char** argv) {
    TApplication* tap= new TApplication("DMCalc", &argc, argv);

    double c= 2.99792458e10;
    double c2= c*c;
    double sigma_n = 1e-42; // WIMP-Nucleon cross-section (cmˆ2)
    double t= 86400*365.25*3; // exposure time in seconds
    double mDEAP= 1000000; // target massing
    double quench  = 0.25;         // recoil energy quenching factor
    double Ethresh = 15 ; // energy threshold in keVee
    double Ewin = 25 ; // energy window in keVee

    DMCalc* deap= new DMCalc("Xe");

    double wimpmass[1000] ;
    double wimpcsc[1000] ;
    int n = 0 ;
    // measured energy distribution (keVee)
    TH1D* teff = new TH1D("teff" ,"teff" ,500 ,0 ,500) ;
    // WIMP recoil energy distribution (keVr)
    TH1D* DiffSpec= new TH1D("DiffSpec", "Differential spectrum", 500, 0, 500); 
    for (double i= 0; i<= 4.2; i+=0.1) {
        wimpmass[n]= TMath::Power(10, i);
        double m_w= TMath::Power(10, i)*1e9/c2;
        teff->Reset(); // reset for each mass loop
        for (int j= 1; j<= 500; ++j) {
            // recoil energy in eV
            double Er= (j-0.5)*1000; 
            double dRdEr= TMath::Power(deap->FormFactor(Er),2)*deap->dRdEreinf(m_w, Er)*t*mDEAP;
            DiffSpec->SetBinContent(j, dRdEr);
            double Evis= Er*quench/1000.;
            for (int k= 1; k<= 500; ++k) {
                double tmp_Er= k - 0.5;
                teff->AddBinContent(k, dRdEr*deap->DetResponse(Evis, tmp_Er));
            }
        }
        // - - - - - - - - rate as integral over energy plot
        int binlo= teff->GetXaxis()->FindBin(Ethresh);
        int binhi= teff->GetXaxis()->FindBin(Ethresh+Ewin);
        double rate= teff->Integral(binlo , binhi);
        // 10% probability to see 0 events for mean=2.3
        wimpcsc[n]= 2.3*sigma_n/rate;
        n++;
    }

    TCanvas* cc= new TCanvas("cc", "cc", 800, 600);
    cc->SetLogy();
    cc->SetLogx();
    TGraph* graph = new TGraph(n, wimpmass, wimpcsc);
    graph->Draw("AL");
    graph->SetMinimum(1e-48);
    graph->SetMaximum(1e-40);
    graph->GetXaxis()->SetLimits(8,10000);
    graph->SetLineWidth(1);

    tap->Run();

    return 0;
}
