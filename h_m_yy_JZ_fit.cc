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

void h_m_yy_JZ_fit() {

    //TDatime startTime;
    //Int_t theSeed = 0;
    //theSeed = startTime.GetTime();
    //RooRandom::randomGenerator()->SetSeed(theSeed);
    //cout << "Will use " << theSeed << " as random seed ... " << endl;

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

    RooPlot* frame = mass.frame();
    data->plotOn(frame); // Tracer l'histogramme
    model_PDF.plotOn(frame);

    TCanvas * Canvas = new TCanvas( "Canvas" , "Canvas" ) ;

    frame->Draw() ;

    Canvas -> SaveAs("H_yy_JZ_fit_bkg_only.pdf");




    // Create Graphics

    TGraph graph_p;
    graph_p.SetTitle("p_value");
    graph_p.GetXaxis()->SetTitle("M_{X} [GeV]");  
    graph_p.GetYaxis()->SetTitle("p_value");
    graph_p.SetMarkerStyle(20);       
    graph_p.SetMarkerSize(0.25);  

    TGraph graph_NLL2;
    graph_NLL2.SetTitle("Fluctuations");
    graph_NLL2.GetXaxis()->SetTitle("M_{X} [GeV]");  
    graph_NLL2.GetYaxis()->SetTitle("q");
    graph_NLL2.SetMarkerStyle(20);       
    graph_NLL2.SetMarkerSize(0.25);  

    TGraph graph_NLL;
    graph_NLL.SetTitle("Delta Log-likelihood");
    graph_NLL.GetXaxis()->SetTitle("Mx (GeV)");  
    graph_NLL.GetYaxis()->SetTitle("Delta Log-likelihood");
    graph_NLL.SetMarkerStyle(20);       
    graph_NLL.SetMarkerSize(0.25);  

    TGraph graph_Ns;
    graph_Ns.SetTitle("Number of signal events");
    graph_Ns.GetXaxis()->SetTitle("Mx (GeV)");  
    graph_Ns.GetYaxis()->SetTitle("Ns");  
    graph_Ns.SetMarkerStyle(20);       
    graph_Ns.SetMarkerSize(0.25);

    TGraph graph_sigmas;
    graph_sigmas.SetTitle("Ns uncertainty");
    graph_sigmas.GetXaxis()->SetTitle("Mx (GeV)");  
    graph_sigmas.GetYaxis()->SetTitle("Sigma Ns"); 
    graph_sigmas.SetMarkerStyle(20);       
    graph_sigmas.SetMarkerSize(0.25); 

    TGraph graph_Ns_sigmas;
    graph_Ns_sigmas.SetTitle("Ns/SigmaS");
    graph_Ns_sigmas.GetXaxis()->SetTitle("Mx (GeV)");  
    graph_Ns_sigmas.GetYaxis()->SetTitle("Ns/SigmaS"); 
    graph_Ns_sigmas.SetMarkerStyle(20);       
    graph_Ns_sigmas.SetMarkerSize(0.25);


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

        RooFitResult * myFit = model_PDF.fitTo( *data , Save()) ;

        myFit->Print( "v" ) ;

        
        double nllValue = myFit->minNll();
        double NsValue = Ns.getVal();
        double NsError = Ns.getError();

        // Print nll
        std::cout << "Mass: " << i*1e-3 << " GeV -> NLL: " << nllValue << std::endl;

        double q;

        if(NsValue<0){
            q = 0;
        } else{
            q = -2*(nllValue- nllValue_bkg);
        }

        double p_value = TMath::Prob(q, 1);





        // Taking away outliers

        if (NsError < 1000) {  

            int pointIndex = (i - 6000) / 10;


            graph_p.SetPoint(pointIndex, i * 1e-3, p_value);
            graph_NLL2.SetPoint(pointIndex, i * 1e-3, q);
            graph_NLL.SetPoint(pointIndex, i * 1e-3, nllValue - nllValue_bkg);
            graph_Ns.SetPoint(pointIndex, i * 1e-3, NsValue);
            graph_sigmas.SetPoint(pointIndex, i * 1e-3, NsError);
            graph_Ns_sigmas.SetPoint(pointIndex, i * 1e-3, NsValue / NsError);

        }


        // Libérer la mémoire du dataset
        delete myFit;
        

    }


    //Tracer le graphe
    TCanvas *Canvas2 = new TCanvas("Canvas2", "Canvas2", 400, 400);
    graph_NLL2.Draw("AP");

    Canvas2 -> SaveAs("H_yy_JZ_fluctuations.pdf");

    TCanvas *Canvas4 = new TCanvas("Canvas4", "Canvas4", 400, 400);
    Canvas4->SetLogy();
    graph_p.GetYaxis()->SetRangeUser(1e-7, 1);
    graph_p.Draw("AP");

    double sigmaLevels[5] = {0.3173, 0.0455, 0.0027, 0.000063, 2.87e-7};
    int colors[5] = {kBlack, kRed, kGreen+2, kMagenta, kCyan+2};

    const char* sigmaLabels[5] = {"1#sigma", "2#sigma", "3#sigma", "4#sigma", "5#sigma"};

    TLine *lines[5];
    TLatex *labels[5];

    for (int i = 0; i < 5; i++) {
        // Tracer la ligne horizontale
        lines[i] = new TLine(4, sigmaLevels[i], 36, sigmaLevels[i]);
        lines[i]->SetLineColor(colors[i]);
        lines[i]->SetLineStyle(2); // Pointillé
        lines[i]->Draw("Same");

        // Ajouter un label juste à côté de l'axe des ordonnées
        labels[i] = new TLatex(3 - 0.2, sigmaLevels[i], sigmaLabels[i]); 
        labels[i]->SetTextSize(0.035);  // Taille du texte
        labels[i]->SetTextColor(colors[i]); // Couleur du texte
        labels[i]->SetTextAlign(32);  // Alignement à droite
        labels[i]->Draw();
    }



    //Canvas4 -> SaveAs("H_yy_JZ_p_value.pdf");

    TCanvas *Canvas3 = new TCanvas("Canvas3", "Canvas3", 1600, 400);
    Canvas3->Divide(4, 1);

    Canvas3->cd(1);
 
    graph_NLL.Draw("AP");

    Canvas3->cd(2);

    graph_Ns.Draw("AP");

    Canvas3->cd(3);

    graph_sigmas.Draw("AP");

    Canvas3->cd(4);

    graph_Ns_sigmas.Draw("AP");

    // save graphics

    //Canvas3 -> SaveAs("H_yy_JZ_fit_Ns.pdf");

    
        


}



