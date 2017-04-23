#include <iostream>
#include "TMath.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TApplication.h"
#include "TGraph.h"
#include "TROOT.h"
#include "DMCalc.h"

using namespace std;

int main(int argc, char** argv) {
    gROOT->ForceStyle();
    TApplication* tap= new TApplication("DMCalc", &argc, argv);

    double c= 2.99792458e10;
    double c2= c*c;
    double sigma_n = 1e-42;   // WIMP-Nucleon cross-section (cmˆ2)
    double t= 86400*365.25*3; // exposure time in seconds
    double mDEAP= 1000000;    // target massing
    double quench  = 0.25;    // recoil energy quenching factor
    double Ethresh = 15 ;     // energy threshold in keVee
    double Ewin = 25 ;        // energy window in keVee

    // Arguments: Target, IsSpinIndependent, mode of event rate
    DMCalc* deap= new DMCalc("Ar", "", 3);
//    DMCalc* deap2= new DMCalc("Xe", "", 4);

    double wimpmass[1000] ;
    double wimpcsc[1000] ;
//    double wimpcsc2[1000] ;
    int n= 0 ;
    // measured energy distribution (keVee)
    TH1D* teff = new TH1D("teff" ,"teff" ,500 ,0 ,500) ;
    TH1D* teff2 = new TH1D("teff2" ,"teff2" ,500 ,0 ,500) ;
    for (double i= 0; i<= 4.2; i+=0.1) {
        wimpmass[n]= TMath::Power(10, i);
        double m_w= TMath::Power(10, i)*1e9/c2;
        teff->Reset();                  // reset for each mass loop
//        teff2->Reset();                  // reset for each mass loop
        for (int j= 1; j<= 500; ++j) {
            double Er= (j-0.5)*1000;    // recoil energy in eV
            double DRdRdEr= deap->dRdEr(m_w, Er)*t*mDEAP;  // differential event rate
//            double DRdRdEr2= deap2->dRdEr(m_w, Er)*t*mDEAP;  // differential event rate
            double Evis= Er*quench/1000.;
            for (int k= 1; k<= 500; ++k) {
                double tmp_Er= k - 0.5;
                teff->AddBinContent(k, DRdRdEr*deap->DetResponse(Evis, tmp_Er));
//                teff2->AddBinContent(k, DRdRdEr2*deap2->DetResponse(Evis, tmp_Er));
            }
        }
        // - - - - - - - - rate as integral over energy plot
        int binlo= teff->GetXaxis()->FindBin(Ethresh);
        int binhi= teff->GetXaxis()->FindBin(Ethresh+Ewin);
        double rate= teff->Integral(binlo , binhi);
//        double rate2= teff2->Integral(binlo , binhi);
        // 10% probability to see 0 events for mean=2.3
        wimpcsc[n]= 2.3*sigma_n/rate; // don't understand!!!
//        wimpcsc2[n]= 2.3*sigma_n/rate2; // don't understand!!!
        n++;
    }

    TCanvas* cc= new TCanvas("cc", "cc", 800, 600);
    cc->SetGrid();
    cc->SetLogy();
    cc->SetLogx();
    TGraph* graph = new TGraph(n, wimpmass, wimpcsc);
    graph->SetMinimum(1e-48);
    graph->SetMaximum(1e-40);
    graph->GetXaxis()->SetLimits(8,10000);
    graph->SetLineWidth(1);
    graph->GetXaxis()->SetTitle("m_{#chi} [GeV]");
    graph->GetXaxis()->SetLabelSize(0.03);
    graph->GetXaxis()->SetLabelOffset(0.012);
    graph->GetXaxis()->SetTitleSize(0.035);
    graph->GetXaxis()->SetTitleOffset(1.30);
    graph->GetXaxis()->SetTickLength(0.02);
    graph->GetYaxis()->SetTitle("#sigma [cm^{2}]");
    graph->GetYaxis()->SetLabelSize(0.03);
    graph->GetYaxis()->SetLabelOffset(0.012);
    graph->GetYaxis()->SetTitleSize(0.035);
    graph->GetYaxis()->SetTitleOffset(1.25);
    graph->GetYaxis()->SetTickLength(0.01);
    graph->SetLineColor(kPink-9);
    graph->Draw("AL");
//    TGraph* graph2 = new TGraph(n, wimpmass, wimpcsc2);
//    graph2->SetLineColor(kBlue);
//    graph2->Draw("SAME");

    tap->Run();

    return 0;
}
