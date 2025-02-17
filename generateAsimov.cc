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

void generateAsimov() {

    // Définir les valeurs initiales 
	const Float_t m =  0.0622905   ;
    const Float_t nL_init = 4.47038  ;
    const Float_t nR_init = 5.07865 ;
    const Float_t aSigmaL = 0.012622 ;
    const Float_t bSigmaL = 0.0759829   ;
    const Float_t aSigmaR =  0.0108157   ;
    const Float_t bSigmaR = 0.157893 ;
    const Float_t aAlphaL =  0.0193975;
    const Float_t bAlphaL =  1.32431 ;
    const Float_t aAlphaR = 0.0399246  ;
    const Float_t bAlphaR = 1.07553  ;

    const Float_t m0 =  0.308234   ;
    const Float_t m1 =  -0.0439373   ;
    const Float_t m2 =  0.0022759    ;
    const Float_t m3 =  -3.52928e-05    ;

    const Float_t Nb_init = 100000;
    const Float_t Ns_init = 0;

    const Float_t bkg_x1_init = 2.60979  ;
    const Float_t bkg_x2_init = 10.5988  ;
    const Float_t f_init = 0.681738;

    // Créer Bkg pdf

    RooRealVar mass("mass", "mass", 5, 35, "GeV");
    RooRealVar bkg_x1("bkg_x1", "bkg_x1", bkg_x1_init, 0 , 100 , "GeV" );
    RooFormulaVar bkg_a1("bkg_a1", "bkg_a1", "-1/@0", bkg_x1 );
    RooRealVar bkg_x2("bkg_x2", "bkg_x2", bkg_x2_init, 0 , 100 , "GeV" );
    RooFormulaVar bkg_a2("bkg_a2", "bkg_a2", "-1/@0", bkg_x2 );
    RooRealVar f("f", "f", f_init, 0, 1 );

	RooExponential exp1_PDF( "exp1_PDF" , "exp1_PDF" , mass , bkg_a1 ) ;
    RooExponential exp2_PDF( "exp2_PDF" , "exp2_PDF" , mass , bkg_a2 ) ;

    RooAddPdf bkg_PDF("bkg_PDF", "bkg_PDF", exp1_PDF, exp2_PDF, f);


    // Créer un Asimov dataset (données binned sans fluctuations)
    int nbins = 300; // Nombre de bins pour l'histogramme
    mass.setBins(nbins);
    int n_events = 13153;  // Nombre total d'événements attendus
    RooDataHist* asimovData = bkg_PDF.generateBinned(mass, RooFit::ExpectedData(), RooFit::NumEvents(n_events));


    // Optionnel : Tracer les données Asimov
    TCanvas* c = new TCanvas("c", "Asimov Dataset", 800, 600);
    RooPlot* frame = mass.frame();
    asimovData->plotOn(frame); // Afficher les points de l'Asimov dataset
    bkg_PDF.plotOn(frame, RooFit::LineColor(kRed)); // Superposer la PDF
    frame->Draw();
    c->SaveAs("asimov_dataset.png");

    // Define the sig PDF

    RooRealVar mean ( "mean" , "mean" , 20, 0, 40, "GeV" ) ;
    RooRealVar sigmaL ("sigmaL", "sigmaL", aSigmaL*20 + bSigmaL, 1e-3, 2);
    RooRealVar sigmaR ("sigmaR","sigmaR", aSigmaR*20 + bSigmaR, 1e-3, 2);
    RooRealVar alphaL ("alphaL", "alphaL", aAlphaL*20 + bAlphaL, 1e-3, 10);
    RooRealVar alphaR ("alphaR", "alphaR", aAlphaR*20 + bAlphaR, 1e-3, 10);
    RooRealVar nL("nL", "nL", 6.882 , 1e-3, 15);
    RooRealVar nR("nR", "nR", 6.36303 , 1e-3, 15);

	RooCrystalBall sig_PDF( "sig_PDF" , "sig_PDF" , mass , mean , sigmaL, sigmaR, alphaL, nL, alphaR, nR ) ;

    // Define model PDF
    RooRealVar Nb("Nb", "Number of Background events", 13153 );
    RooRealVar Ns("Ns", "Number of Signal events", 0  );

    RooAddPdf model_PDF("model_PDF", "model_PDF", RooArgList(sig_PDF, bkg_PDF), RooArgList(Ns, Nb));

    TGraph graph_sigmas;
    graph_sigmas.SetTitle("Ns uncertainty");
    graph_sigmas.GetXaxis()->SetTitle("Mx (GeV)");  
    graph_sigmas.GetYaxis()->SetTitle("Sigma Ns"); 
    graph_sigmas.SetMarkerStyle(20);       
    graph_sigmas.SetMarkerSize(0.25); 

    TGraph graph_Ns;
    graph_Ns.SetTitle("Number of signal events");
    graph_Ns.GetXaxis()->SetTitle("Mx (GeV)");  
    graph_Ns.GetYaxis()->SetTitle("Ns");  
    graph_Ns.SetMarkerStyle(20);       
    graph_Ns.SetMarkerSize(0.25);

    // Scan 

    for (int i = 6000; i <= 34000 ; i+=10) {

    // Reset parameter values

    mean.setVal(m3*TMath::Power(i*1e-3, 3)+ m2*TMath::Power(i*1e-3, 2)+ m1*(i*1e-3) + m0 + i*1e-3);
    sigmaL.setVal(aSigmaL*i*1e-3+bSigmaL);
    sigmaR.setVal(aSigmaR*i*1e-3+bSigmaR);
    alphaL.setVal(aAlphaL*i*1e-3+bAlphaL);
    alphaR.setVal(aAlphaR*i*1e-3+bAlphaR);
    nL.setVal(nL_init);
    nR.setVal(nR_init);
    Ns.setVal(Ns_init);
    Nb.setVal(Nb_init);

    bkg_x1.setVal(bkg_x1_init);
    bkg_x2.setVal(bkg_x2_init);
    f.setVal(f_init);

        // Fit

        sigmaL.setConstant( kTRUE ) ;
        sigmaR.setConstant( kTRUE ) ;
        alphaL.setConstant( kTRUE ) ;
        alphaR.setConstant( kTRUE ) ;
        mean.setConstant(kTRUE);
        nL.setConstant(kTRUE);
        nR.setConstant(kTRUE);

        bkg_x1.setConstant(kFALSE);
        bkg_x2.setConstant(kFALSE);
        f.setConstant(kFALSE);

        Nb.setConstant( kFALSE);
        Ns.setConstant( kFALSE);

        RooFitResult * myFit = model_PDF.fitTo( *asimovData , Save()) ;

        myFit->Print( "v" ) ;

        double NsValue = Ns.getVal();
        double NsError = Ns.getError();

        int pointIndex = (i - 6000) / 10;

        graph_sigmas.SetPoint(pointIndex, i * 1e-3, NsError);
        graph_Ns.SetPoint(pointIndex, i * 1e-3, NsValue);


        // Libérer la mémoire du dataset
        delete myFit;
        

    }

    TCanvas *Canvas = new TCanvas("Canvas", "Canvas", 800, 400);
    Canvas->Divide(2, 1);

    Canvas->cd(1);
 
    graph_Ns.Draw("AP");

    Canvas->cd(2);

    graph_sigmas.Draw("AP");

    Canvas->SaveAs("Ns_asimov.pdf");

    graph_sigmas.SetName("sigmaNs");

    auto *outputFile = TFile::Open("sigmaNs.root", "RECREATE");
    graph_sigmas.Write();
    outputFile->Close();




}
