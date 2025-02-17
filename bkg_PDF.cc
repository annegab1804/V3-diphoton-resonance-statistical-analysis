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

void bkg_PDF() {


    // Load histogram
    TFile *file = TFile::Open("h_m_yy_JZ/h_m_yy_JZ_barrel_20.root");
    TH1D *h = (TH1D*)file->Get("h_m_yy_tot");
    h->GetXaxis()->SetRangeUser(5, 35); // je m'intéresse qu'aux faibles masses


    TF1 * myTwoExp = new TF1( "myTwoExp" , "[0]*TMath::Exp(-x/[1])+[2]*TMath::Exp(-x/[3])" , 5 , 35 )  ;
    myTwoExp->SetParameters(1.0e8,2.0,1.0e5,10.0) ;
    h->Fit(myTwoExp,"","",5,35);


    TCanvas *c1 = new TCanvas("c1", "Fit bkg", 800, 600);
    h->GetXaxis()->SetTitle("mass (GeV)");
    h->GetYaxis()->SetTitle("Event / (1GeV)");
    h->SetMarkerStyle(20); 
    h->Draw("E");

    // Dessiner la fonction ajustée
    myTwoExp->SetLineColor(kRed);
    myTwoExp->SetLineWidth(2);
    myTwoExp->Draw("SAME");

    c1->SaveAs("H_yy_JZ_fit_bkg.pdf");

    // Extraire les paramètres du fit
    double A1 = myTwoExp->GetParameter(0);  
    double bkg_x1_init = myTwoExp->GetParameter(1);  
    double A2 = myTwoExp->GetParameter(2);  
    double bkg_x2_init = myTwoExp->GetParameter(3);  

    // Calcul de la fraction f
    double integral1 = bkg_x1_init*(TMath::Exp(-5/bkg_x1_init)-TMath::Exp(-35/bkg_x1_init));
    double integral2 = bkg_x2_init*(TMath::Exp(-5/bkg_x2_init)-TMath::Exp(-35/bkg_x2_init));
    double f_init = A1*integral1 / (A1*integral1 + A2*integral2);
    std::cout << "f=" << f_init << std::endl;

    //// on crée la PDF
    RooRealVar mass("mass", "mass", 5, 35, "GeV");
    RooRealVar bkg_x1("bkg_x1", "bkg_x1", bkg_x1_init, 0 , 100 , "GeV" );
    RooFormulaVar bkg_a1("bkg_a1", "bkg_a1", "-1/@0", bkg_x1 );
    RooRealVar bkg_x2("bkg_x2", "bkg_x2", bkg_x2_init, 0 , 100 , "GeV" );
    RooFormulaVar bkg_a2("bkg_a2", "bkg_a2", "-1/@0", bkg_x2 );
    RooRealVar f("f", "f", f_init, 0, 1 );

	RooExponential exp1_PDF( "exp1_PDF" , "exp1_PDF" , mass , bkg_a1 ) ;
    RooExponential exp2_PDF( "exp2_PDF" , "exp2_PDF" , mass , bkg_a2 ) ;

    RooAddPdf bkg_PDF("bkg_PDF", "bkg_PDF", exp1_PDF, exp2_PDF, f);

    RooDataHist dataHist("dataHist", "dataHist", mass, Import(*h));

    RooPlot* frame = mass.frame();
    dataHist.plotOn(frame);
    bkg_PDF.plotOn(frame);

    TCanvas * Canvas = new TCanvas( "Canvas" , "Canvas" ) ;
    Canvas->SetLogy();


    frame->Draw() ;

    Canvas -> SaveAs("H_yy_JZ_PDF_bkg_log.pdf");

}