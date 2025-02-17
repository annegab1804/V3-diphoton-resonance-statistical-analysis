#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "TCanvas.h"

#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooExponential.h"
#include "RooGaussian.h"
#include "RooAddPdf.h"
#include "RooFitResult.h"
#include "RooDataSet.h"
#include "RooPlot.h"
#include "TH1D.h"
#include "RooDataHist.h"
#include "RooFormulaVar.h"
#include "RooCrystalBall.h"


using namespace RooFit ;

void DSCB() {

    RooRealVar x("x", "x", -10, 10);

    RooRealVar mean ( "mean" , "mean" , 0, -1, 1) ;
    RooRealVar sigmaL ("sigmaL", "sigmaL", 1, 0.1, 2);
    RooRealVar sigmaR ("sigmaR","sigmaR", 1.5, 0.1, 2);
    RooRealVar alphaL ("alphaL", "alphaL", 1.5, 0.1, 5);
    RooRealVar alphaR ("alphaR", "alphaR", 1.2, 0.1, 5);
    RooRealVar nL("nL", "nL", 2 , 0.1, 15);
    RooRealVar nR("nR", "nR", 3 , 0.1, 15);

	RooCrystalBall PDF( "PDF" , "PDF" , x , mean , sigmaL, sigmaR, alphaL, nL, alphaR, nR ) ;

    TCanvas *c = new TCanvas("c", "Double-Sided Crystal Ball with RooFit", 800, 600);
    RooPlot *frame = x.frame();
    c->SetTitle("Double-Sided Crystal Ball ");
    
    // Tracer la fonction DSCB
    PDF.plotOn(frame, LineColor(kBlue));


    x.setRange("full", -10, 10);

    x.setRange("droite", 1.2*1.5, 10); 
    x.setRange("gauche", -10, -1.5); 
    PDF.plotOn(frame, LineColor(kRed), Range("droite"), NormRange("full"));
    PDF.plotOn(frame, LineColor(kRed), Range("gauche"), NormRange("full"));

    // Affichage
    frame->Draw();



    // Ajouter les annotations pour sigma, alpha et n

    TLine *lineM = new TLine(0, 0, 0, 0.053);
    lineM->SetLineColor(kBlack);
    lineM->SetLineStyle(2);
    lineM->Draw();

    TLine *lineSr = new TLine(0, 0.032, 1.5, 0.032);
    lineSr->SetLineColor(kRed);
    lineSr->SetLineStyle(2);
    lineSr->Draw();

    TLine *lineSl = new TLine(-1, 0.032, 0, 0.032);
    lineSl->SetLineColor(kRed);
    lineSl->SetLineStyle(2);
    lineSl->Draw();


    TLine *lineAl = new TLine(0, 0.0175, -1.5, 0.0175);
    lineAl->SetLineColor(kBlue);
    lineAl->SetLineStyle(2);
    lineAl->Draw();

    TLine *lineAr = new TLine(0, 0.026, 1.2*1.5, 0.026);
    lineAr->SetLineColor(kBlue);
    lineAr->SetLineStyle(2);
    lineAr->Draw();

    TLatex latex;
    latex.SetTextSize(0.03);
    latex.SetTextAlign(22); // Alignement centré

    // Ajouter les annotations pour sigma, alpha et n

    latex.DrawLatex(0.7, 0.034, "#bf{#sigma_{R}}"); // Pour sigmaR
    latex.DrawLatex(-0.5, 0.034, "#bf{#sigma_{L}}"); // Pour sigmaL
    latex.DrawLatex(-0.8, 0.02, "#bf{#alpha_{L}#sigma_{L}}"); // Pour alphaL
    latex.DrawLatex(1, 0.028, "#bf{#alpha_{R}#sigma_{R}}"); // Pour alphaR
    c->SetTitle("Double-Sided Crystal Ball ");


    c->SaveAs("DSCB_RooFit.png");

}