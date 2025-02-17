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

void bkg_tot() {

    //TDatime startTime;
    //Int_t theSeed = 0;
    //theSeed = startTime.GetTime();
    //RooRandom::randomGenerator()->SetSeed(theSeed);
    //cout << "Will use " << theSeed << " as random seed ... " << endl;


    // Load histogram
    TFile *file = TFile::Open("h_m_yy_JZ/h_m_yy_JZ_barrel_20.root");
    TH1D *h_red = (TH1D*)file->Get("h_m_yy_tot");
    h_red->GetXaxis()->SetRangeUser(5, 35); // je m'intéresse qu'aux faibles masses


    TFile *file_run2 = TFile::Open("run2data.root"); 
    TH1F *hRatio = (TH1F*)file_run2->Get("hRatio");

    std::cout << "Il y a " << h_red->GetNbinsX() << "bins dans h_red "<< std::endl;
    std::cout << "Il y a " << hRatio->GetNbinsX() << "bins dans hRatio "<< std::endl;


    TH1F *h_tot = new TH1F("h_tot", "Bkg tot", 30, 5, 35);

    //for (int i = 1; i <= h_red->GetNbinsX(); ++i) {
    //std::cout << "Bin " << i << " : "
    //          << "mass value = " << h_red->GetBinCenter(i) 
    //          << ", low edge = " << h_red->GetBinLowEdge(i)
    //          << ", high edge = " << h_red->GetBinLowEdge(i+1)
    //          << ", content = " << h_red->GetBinContent(i)
    //          << ", error = " << h_red->GetBinError(i) << std::endl;
    //}
//
    //for (int i = 1; i <= hRatio->GetNbinsX(); ++i) {
    //std::cout << "Bin Ratio " << i << " : "
    //          << " mass value ratio = " << hRatio->GetBinCenter(i) 
    //          << ", low edge ratio = " << hRatio->GetBinLowEdge(i)
    //          << ", high edge ratio = " << hRatio->GetBinLowEdge(i+1)
    //          << ", content ratio = " << hRatio->GetBinContent(i)
    //          << ", error ratio = " << hRatio->GetBinError(i) << std::endl;
    //}




    for (int i = 1; i <= hRatio->GetNbinsX(); ++i) {

        double ratio = hRatio->GetBinContent(i);
        double err_ratio = hRatio->GetBinError(i);

        // Trouver la position du bin central de h_red (par exemple pour l'index i)
        double mass = hRatio->GetBinCenter(i);

        // Chercher quel bin de hRatio correspond au bin central de h_red
        int binIndex_red = h_red->FindBin(mass);

        double red = h_red->GetBinContent(binIndex_red);
        double err_red = h_red->GetBinError(binIndex_red);

        double irr = ratio * red;

        double err_irr = irr * sqrt(
            (err_ratio / ratio) * (err_ratio / ratio) +
            (err_red / red) * (err_red / red)
        );  

        double tot = irr + red;
        double err_tot = err_red + err_irr;

        std::cout << "pour mass = " << mass << " red = " << red << " +- "<< err_red <<", irr = " << irr <<" +- " << err_irr << ", et total bkg =  "<< tot<< " +- "<< err_tot<< std::endl;

        // Ajouter le bruit de fond réductible et irréductible au bin correspondant
        h_tot->SetBinContent(i, tot);
        h_tot->SetBinError(i, err_tot);
    }

    // Créer un canvas pour afficher les résultats
    TCanvas *c1 = new TCanvas("c1", "Superposition du BFR et du BDTot", 800, 600);
    c1->SetLogy();

    // Dessiner le bruit de fond réductible (BFR) avec une couleur différente
    h_red->SetLineColor(kRed);  
    h_red->GetXaxis()->SetTitle("mass (GeV)");
    h_red->GetYaxis()->SetTitle("Event / (1GeV)");
    h_red->SetStats(0);  // Mettre BFR en rouge
    h_red->Draw("P");                // Dessiner l'histogramme réductible
//
    // Dessiner le bruit de fond total (BFR + BDI) sur le même canvas avec l'option "SAME"
    h_tot->SetLineColor(kBlue); 
    h_tot->SetStats(0);
    h_tot->Draw("Same P");        // Dessiner l'histogramme total

    // Ajouter des légendes pour identifier chaque histogramme
    TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9); // Position de la légende
    legend->AddEntry(h_red, "Bkg reductible", "l");  // Ajouter BFR à la légende
    legend->AddEntry(h_tot, "Bkg total", "l"); // Ajouter BDTot à la légende
    legend->Draw();

    // Sauvegarder l'image
    c1->SaveAs("bkgtot.png");

    auto *outputFile = TFile::Open("bkg_tot.root", "RECREATE");
    h_tot->Write();
    outputFile->Close();

}