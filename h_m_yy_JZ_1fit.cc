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

void h_m_yy_JZ_1fit() {

    //TDatime startTime;
    //Int_t theSeed = 0;
    //theSeed = startTime.GetTime();
    //RooRandom::randomGenerator()->SetSeed(theSeed);
    //cout << "Will use " << theSeed << " as random seed ... " << endl;

    // Définir les valeurs initiales 
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

    const Float_t bkg_x1_init = 2.35988 ;
    const Float_t bkg_x2_init = 9.72257 ;
    const Float_t f_init = 0.68181;


    // Define the bkg PDF

    RooRealVar mass("mass", "mass", 6, 35, "GeV");
    RooRealVar bkg_x1("bkg_x1", "bkg_x1", bkg_x1_init, 0 , 100 , "GeV" );
    RooFormulaVar bkg_a1("bkg_a1", "bkg_a1", "-1/@0", bkg_x1 );
    RooRealVar bkg_x2("bkg_x2", "bkg_x2", bkg_x2_init, 0 , 100 , "GeV" );
    RooFormulaVar bkg_a2("bkg_a2", "bkg_a2", "-1/@0", bkg_x2 );
    RooRealVar f("f", "f", f_init, 0, 1 );

	RooExponential exp1_PDF( "exp1_PDF" , "exp1_PDF" , mass , bkg_a1 ) ;
    RooExponential exp2_PDF( "exp2_PDF" , "exp2_PDF" , mass , bkg_a2 ) ;

    RooAddPdf bkg_PDF("bkg_PDF", "bkg_PDF", exp1_PDF, exp2_PDF, f);


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
    RooRealVar Nb("Nb", "Number of Background events", 100000 );
    RooRealVar Ns("Ns", "Number of Signal events", 0  );

    RooAddPdf model_PDF("model_PDF", "model_PDF", RooArgList(sig_PDF, bkg_PDF), RooArgList(Ns, Nb));

    // create a toy dataset

    RooRealVar nEvt( "nEvt" , "nEvt" ,  100000 ) ;
	RooDataSet * data = bkg_PDF.generate( mass , nEvt.getVal() ) ;

    // Fit Bkg Only

    sigmaL.setConstant( kTRUE ) ;
    sigmaR.setConstant( kTRUE ) ;
    alphaL.setConstant( kTRUE ) ;
    alphaR.setConstant( kTRUE ) ;
    mean.setConstant(kTRUE);
    nL.setConstant(kTRUE);
    nR.setConstant(kTRUE);

    Ns.setConstant(kTRUE);
    Nb.setConstant(kFALSE);

    bkg_x1.setConstant( kFALSE ) ;
    bkg_x2.setConstant( kFALSE ) ;
    f.setConstant( kFALSE ) ;

    RooFitResult* myFit_bkg = model_PDF.fitTo(*data, Save());

	myFit_bkg -> Print( "v" ) ;

    // log likelihood 

    double nllValue_bkg = myFit_bkg->minNll();


    // Reset parameter values

    double m = 5;

    mean.setVal(m3*TMath::Power(m, 3)+ m2*TMath::Power(m, 2)+ m1*(m) + m0 + m);
    sigmaL.setVal(aSigmaL*m+bSigmaL);
    sigmaR.setVal(aSigmaR*m+bSigmaR);
    alphaL.setVal(aAlphaL*m+bAlphaL);
    alphaR.setVal(aAlphaR*m+bAlphaR);
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

        RooFitResult * myFit = model_PDF.fitTo( *data , Save()) ;

        myFit->Print( "v" ) ;

        
        double nllValue = myFit->minNll();
        double NsValue = Ns.getVal();
        double NsError = Ns.getError();
        double q = -2*(nllValue-nllValue_bkg);

        // Print nll
        std::cout << "Mass: " << 15<< " GeV -> NLL: " << nllValue << std::endl;
        std::cout << "Mass: " << 15<< " GeV -> q: " << q<< std::endl;
        

    RooPlot* frame = mass.frame();
    data->plotOn(frame, Name("data")); // Tracer l'histogramme
    model_PDF.plotOn(frame, Name("model"));
    model_PDF.plotOn( frame , Components( bkg_PDF ) , LineStyle( kDashed ), Name("background-only")) ;
    model_PDF.plotOn( frame , Components( sig_PDF ) , LineStyle( kDashed) , LineColor(kRed), Name("signal-only")) ;
    frame->GetXaxis()->SetRangeUser(m-2, m+4);
    frame->GetYaxis()->SetRangeUser(0, 5000);

    // Créer une légende
    TLegend* legend = new TLegend(0.7, 0.7, 0.9, 0.9); // Position (x1, y1, x2, y2)
    legend->SetBorderSize(0);                          // Pas de bordure
    legend->SetFillStyle(0);                           // Pas de fond
    legend->AddEntry(frame->findObject("data"), "Toy data", "lep");   // L'histogramme des données
    legend->AddEntry(frame->findObject("model"), "Model", "l"); // Modèle total
    legend->AddEntry(frame->findObject("background-only"), "Background", "l");    // Fond
    legend->AddEntry(frame->findObject("signal-only"), "Signal", "l");        // Signal



    TCanvas * Canvas = new TCanvas( "Canvas" , "Canvas" ) ;

    frame->Draw() ;
    legend->Draw();

    Canvas -> SaveAs("H_yy_JZ_fit_5.pdf");


}



