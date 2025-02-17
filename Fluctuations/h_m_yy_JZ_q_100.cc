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
#include "TFile.h"
#include "TH1F.h"
#include "RooDataHist.h"
#include "RooFormulaVar.h"
#include "RooCrystalBall.h"
#include "TLatex.h"
#include "vector"

using namespace RooFit ;

void h_m_yy_JZ_q_100() {

    TFile *fileq = new TFile("fluctuations_100.root", "RECREATE");
    TFile *filep = new TFile("p-value_100.root", "RECREATE");

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

    const Float_t Nb_init = 100000;
    const Float_t Ns_init = 0;

    const Float_t bkg_x1_init = 2.35988 ;
    const Float_t bkg_x2_init = 9.72257 ;
    const Float_t f_init = 0.68181;

    const Float_t m0 =  0.308234   ;
    const Float_t m1 =  -0.0439373   ;
    const Float_t m2 =  0.0022759    ;
    const Float_t m3 =  -3.52928e-05    ;


    for (int j = 0; j < 100 ; j++) {

        TGraph graphp;
        graphp.SetTitle("p_value");
        graphp.GetXaxis()->SetTitle("M_{X} [GeV]");  
        graphp.GetYaxis()->SetTitle("p_value");
        graphp.SetMarkerStyle(20);       
        graphp.SetMarkerSize(0.25);  

        TGraph graphq;
        graphq.SetTitle("Fluctuations");
        graphq.GetXaxis()->SetTitle("M_{X} [GeV]");  
        graphq.GetYaxis()->SetTitle("q");
        graphq.SetMarkerStyle(20);       
        graphq.SetMarkerSize(0.25);  

        TDatime startTime;
        Int_t theSeed = 0;
        theSeed = startTime.GetTime();
        RooRandom::randomGenerator()->SetSeed(theSeed);
        cout << "Will use " << theSeed << " as random seed ... " << endl;

        //// on crée la PDF
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

        RooRealVar mean( "mean" , "mean" , 20, 0, 40, "GeV" ) ;
        RooRealVar sigmaL ("sigmaL", "sigmaL", aSigmaL*20 + bSigmaL, 1e-3, 2);
        RooRealVar sigmaR ("sigmaR","sigmaR", aSigmaR*20 + bSigmaR, 1e-3, 2);
        RooRealVar alphaL ("alphaL", "alphaL", aAlphaL*20 + bAlphaL, 1e-3, 10);
        RooRealVar alphaR ("alphaR", "alphaR", aAlphaR*20 + bAlphaR, 1e-3, 10);
        RooRealVar nL("nL", "nL", 2.30541 , 1e-3, 10);
        RooRealVar nR("nR", "nR", 8.29111 , 1e-3, 10);

	    RooCrystalBall sig_PDF( "sig_PDF" , "sig_PDF" , mass , mean , sigmaL, sigmaR, alphaL, nL, alphaR, nR ) ;

        // Define model PDF
        RooRealVar Nb("Nb", "Number of Background events", 100000 );
        RooRealVar Ns("Ns", "Number of Signal events", 0);
        RooAddPdf model_PDF("model_PDF", "model_PDF", RooArgList(sig_PDF, bkg_PDF), RooArgList(Ns, Nb));

        // create a toy dataset

        RooRealVar nEvt( "nEvt" , "nEvt" ,  100000 ) ;
	    RooDataSet * data = bkg_PDF.generate( mass , nEvt.getVal() ) ;

        // Ajuster la PDF sur les données

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

        // On extrait la log vraissemblance
        double nllValue_bkg = myFit_bkg->minNll();


        // on ajoute la pdf du signal pour différentes masses

        for (int i = 6000; i <= 34000 ; i+=10) {

        // Réinitialiser les valeurs des paramètres

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

            // On extrait la log vraissemblance
            double nllValue = myFit->minNll();
            double NsValue = Ns.getVal();
            double NsError = Ns.getError();

            double q;

            if(NsValue<0){
                q = 0;
            } else{
                q = -2*(nllValue- nllValue_bkg);
            }

            double p_value = TMath::Prob(q, 1);



            // Filtrer les valeurs de sigma
            if (NsError<1000){

                int pointIndex = (i - 6000) / 10;

                graphp.SetPoint(pointIndex, i * 1e-3, p_value);
                graphq.SetPoint(pointIndex, i * 1e-3, q);
            
            }


            // Libérer la mémoire du dataset
            delete myFit;



    }

    // Enregistrer chaque TGraph dans le fichier ROOT avec un nom unique
    TString graphNameq = Form("graphq_%d", j+1);  // Nom du graphique: graph_1, graph_2, ..., graph_20
    fileq->cd();
    graphq.Write(graphNameq); 

    TString graphNamep = Form("graphp_%d", j+1);  // Nom du graphique: graph_1, graph_2, ..., graph_20
    filep->cd();
    graphp.Write(graphNamep); 

    }


    // Fermeture du fichier ROOT après l'écriture
    fileq->Close();
    filep->Close();
    
    // Nettoyage de la mémoire
    delete fileq;
    delete filep;



}