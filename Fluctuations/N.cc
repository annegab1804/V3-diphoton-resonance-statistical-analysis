#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "TCanvas.h"
#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooExponential.h"
#include "RooGaussian.h"
#include "RooCrystalBall.h"
#include "RooPolynomial.h"
#include "RooAddPdf.h"
#include "RooFitResult.h"
#include "RooDataSet.h"
#include "RooPlot.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TBox.h"
#include "TAxis.h"


void N() {

    TF1* linearFitSigmaL = new TF1("linearFitSigmaL",  "pol1", 5, 35);
    TF1* linearFitSigmaR = new TF1("linearFitSigmaR",  "pol1", 5, 35);

    linearFitSigmaL->SetParameters(0.0759829 ,0.012622) ;
    linearFitSigmaR->SetParameters(0.157893 ,0.0108157) ;

    double m = 6;
    int n = 0;

    std::vector<double> centres;   
    std::vector<double> erreursL;  
    std::vector<double> erreursR;


    while(m<=34) {

        n++; 
        std::cout << "centre de la région statistique n°" << n<< ":" << m << std::endl;

        // Calcul des sigma gauche et droit pour m
        double sigmaL = linearFitSigmaL->Eval(m);
        double sigmaR = linearFitSigmaR->Eval(m);

        // Stocker les données pour le graphe
        centres.push_back(m);

        double x = m ;
        
        while(x-m < (1.5*sigmaR + 1.5*sigmaL)){
            x += 0.001;
            sigmaL = linearFitSigmaL->Eval(x);
        }

        m = x;  
 

    }

    std::cout << "il y a "<< n << "régions statistiquement indépendantes" <<std::endl;

    TCanvas* c2 = new TCanvas("c2", "Segmentation de l'axe plan", 800, 400);

    // Création d'un graphique vide (juste pour définir les axes)
    TH2F *frame = new TH2F("frame", "Statistically independent regions", 10, 5, 35, 10, 0, 1);
    frame->GetXaxis()->SetTitle("m_{X} (GeV)");
    frame->Draw();


    for (size_t i = 0; i < centres.size(); i++) {

        m = centres[i];

        double sigmaL = linearFitSigmaL->Eval(m);
        double sigmaR = linearFitSigmaR->Eval(m);

        double x_min = m - 1.5 * sigmaL;
        double x_max = m + 1.5 * sigmaR;

        TBox *box = new TBox(x_min, 0, x_max, 1);
        if (i%2==0){box->SetFillColorAlpha(kBlue, 0.3);  // Bleu avec 30% d'opacité
        box->Draw("same");}
        else{box->SetFillColorAlpha(kGreen, 0.3);  // Bleu avec 30% d'opacité
        box->Draw("same");}


    }

    frame->SetStats(0);


    c2->SaveAs("regions_statistiques_plan.png");




}