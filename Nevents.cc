#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "TCanvas.h"

#include "RooRealVar.h"
#include "RooHistPdf.h"
#include "RooDataSet.h"
#include "RooPlot.h"
#include "TH1D.h"
#include "TFile.h"
#include "TClass.h"
#include "TList.h"
#include <iostream>

using namespace RooFit;

void Nevents() {

    TFile *file = TFile::Open("Run2data.root");
    TH1F *h_irr = (TH1F*)file->Get("h_irr");
    TH1F *h_red = (TH1F*)file->Get("h_red");

    TFile *fileRatio = TFile::Open("h_m_yy_ratio_jose.root");
    TH1F *hRatio = (TH1F*)fileRatio->Get("h_m_yy_ratio");


    std::cout << "Il y a " << h_red->GetNbinsX() << "bins dans h_red "<< std::endl;
    std::cout << "Il y a " << h_irr->GetNbinsX() << "bins dans h_irr "<< std::endl;
    std::cout << "Il y a " << hRatio->GetNbinsX() << "bins dans hRatio "<< std::endl;

    TH1F *h_irr_run3 = new TH1F("h_irr_run3", "h_irr_run3", 30, 5, 35);
    TH1F *h_red_run3 = new TH1F("h_red_run3", "h_red_run3", 30, 5, 35);

    for (int i = 1; i <= h_red->GetNbinsX(); ++i) {

        double red = h_red->GetBinContent(i);
        double err_red = h_red->GetBinError(i);
        double irr = h_irr->GetBinContent(i);
        double err_irr = h_irr->GetBinError(i);

        double mass = h_red->GetBinCenter(i);

        // Chercher quel bin de hRatio correspond au bin central de h_red
        int binIndex_ratio = hRatio->FindBin(mass);

        double ratio = hRatio->GetBinContent(binIndex_ratio);
        double err_ratio = hRatio->GetBinError(binIndex_ratio);

        double redRun3 = ratio*red;
        double irrRun3 = ratio*irr;

        double err_irrRun3 = irrRun3 * sqrt(
            (err_ratio / ratio) * (err_ratio / ratio) +
            (err_irr / irr) * (err_irr / irr)
        );  

        double err_redRun3 = redRun3 * sqrt(
            (err_ratio / ratio) * (err_ratio / ratio) +
            (err_red / red) * (err_red / red)
        );  

        // Ajouter le bruit de fond réductible et irréductible au bin correspondant
        h_irr_run3->SetBinContent(i, irrRun3);
        h_irr_run3->SetBinError(i, err_irrRun3);
        h_red_run3->SetBinContent(i, redRun3);
        h_red_run3->SetBinError(i, err_redRun3);

    }

    TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);

    h_red->SetLineColor(kBlack);  
    h_red->GetXaxis()->SetTitle("mass (GeV)");
    h_red->GetYaxis()->SetTitle("Event / (1GeV)");
    h_red->SetStats(0);
    h_red->Draw("E"); 
    h_red_run3->SetLineColor(kBlue);
    h_red_run3->Draw("Same E");
    c1->SaveAs("bkg_red_run3.png");

    TCanvas *c2 = new TCanvas("c2", "c2", 800, 600);


    // Dessiner le bruit de fond réductible (BFR) avec une couleur différente
    h_irr->SetLineColor(kBlack);  
    h_irr->GetXaxis()->SetTitle("mass (GeV)");
    h_irr->GetYaxis()->SetTitle("Event / (1GeV)");
    //h_irr->SetStats(0);
    h_irr->Draw("E"); 
    h_irr_run3->SetLineColor(kBlue);
    h_irr_run3->Draw("Same E");
    c2->SaveAs("bkg_irr_run3.png");


    double eventRed = 0;
    for (int i = 1; i <= h_red_run3->GetNbinsX(); ++i) {
        eventRed += h_red_run3->GetBinContent(i);
    }
    std::cout << "Nombre d'évènements réductibles : " << eventRed << std::endl;

    double eventIrr = 0;
    for (int i = 1; i <= h_irr_run3->GetNbinsX(); ++i) {
        eventIrr += h_irr_run3->GetBinContent(i);
    }
    std::cout << "Nombre d'évènements irréductibles : " << eventIrr << std::endl;

    std::cout << "Nombre d'évènements au total : " << eventIrr + eventRed<< std::endl;







    //int nEventsRed = h_red_run3->GetEntries();
    //std::cout << "Nombre d'événements réductibles : " << nEventsRed << std::endl;
    //int nEventsIrr = h_irr_run3->GetEntries();
    //std::cout << "Nombre d'événements irréductibles : " << nEventsIrr << std::endl;
    //std::cout << "Nombre d'événements total : " << nEventsIrr + nEventsRed << std::endl;



    //auto *outputFile = TFile::Open("bkg_red_run3.root", "RECREATE");
    //h_red_run3->Write();
    //outputFile->Close();
//
    //auto *outputFile2 = TFile::Open("bkg_irr_run3.root", "RECREATE");
    //h_irr_run3->Write();
    //outputFile2->Close();

    file->Close();
    fileRatio->Close();

}