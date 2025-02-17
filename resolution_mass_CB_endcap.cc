#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "TCanvas.h"
#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooExponential.h"
#include "RooGaussian.h"
#include "RooCrystalBall.h"
#include "RooAddPdf.h"
#include "RooFitResult.h"
#include "RooDataSet.h"
#include "RooPlot.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TGraphErrors.h"


using namespace std;

using namespace RooFit ;

void resolution_mass_CB_endcap() {

    // Récupération des données
    std::unique_ptr<TFile> myFile( TFile::Open("diphoton_ntuple_collimated_yy.root") );
    auto Tree = myFile->Get<TTree>("diphotonTree");

    
    // Création de la dataset
    RooRealVar mass("mass", "mass", 0, 100, "GeV");
    RooRealVar diff("diff", "diff", -50, 50, "GeV");
    RooArgSet vars(mass, diff);
    RooArgSet vars_diff(diff);
    RooDataSet data("data", "data", vars);

    // Définir les valeurs initiales 

	const Float_t m_init = 0  ;
	const Float_t sigmaL_init = 0.4  ;
    const Float_t sigmaR_init = 0.4  ;
    const Float_t alphaL_init = 2 ;
    const Float_t alphaR_init = 2 ;
    const Float_t nL_init = 1.5 ;
    const Float_t nR_init = 1.5 ;

    // Définir la PDF
	
	RooRealVar m( "m" , "m" ,  0 , -2, 2 , "GeV" ) ;
	RooRealVar sigmaL( "sigmaL" , "sigmaL" ,  0.01 ,  10 , "GeV" ) ;
    RooRealVar sigmaR( "sigmaR" , "sigmaR" ,  0.01 ,  10 , "GeV" ) ;
    RooRealVar alphaL( "alphaL" , "alphaL" ,  0.01 ,  50 , "GeV" ) ;
    RooRealVar alphaR( "alphaR" , "alphaR" ,  0.01 ,  50, "GeV" ) ;
    RooRealVar nL( "nL" , "nL" ,  0.01 ,  100 ) ;
    RooRealVar nR( "nR" , "nR" ,  0.01 ,  100 ) ;
	RooCrystalBall PDF( "PDF" , "PDF" , diff , m , sigmaL, sigmaR, alphaL, nL, alphaR, nR ) ;


    // Création d'un TGraph pour stocker l'évolution des écarts-types
    TGraphErrors *sigmaLGraph = new TGraphErrors();
    TGraphErrors *sigmaRGraph = new TGraphErrors();

    // Création d'un TGraph pour stocker l'évolution des moyennes
    TGraphErrors *mGraph = new TGraphErrors();

    // Création d'un TGraph pour stocker l'évolution des alpha
    TGraphErrors *alphaLGraph = new TGraphErrors();
    TGraphErrors *alphaRGraph = new TGraphErrors();

    // Création d'un TGraph pour stocker l'évolution des n
    TGraphErrors *nLGraph = new TGraphErrors();
    TGraphErrors *nRGraph = new TGraphErrors();

    // Dessiner les histogrammes/pdf sur un canvas
    TCanvas* canvas = new TCanvas("canvas", "Canvas", 2000, 400);
    canvas->Divide(5, 1, 0.01, 0.1);


    // Ajouter un titre
    TPaveText *title = new TPaveText(0.1, 0.94, 0.9, 0.98, "NDC"); // Coordonnées en NDC (Normalized Device Coordinates)
    title->SetFillColor(0);
    title->SetTextAlign(22);  // Centré horizontalement et verticalement
    title->AddText("Mass resolution !(UU)-endcap Crystal Ball [5;10]");
    canvas->cd(0);
    title->Draw();  


    // Déclaration des variables pour stocker les valeurs des branches
    Float_t y_pt, yPrime_pt, y_eta, yPrime_eta, y_phi, yPrime_phi, y_e, yPrime_e,y_truth_eta, yPrime_truth_eta, y_truth_phi, yPrime_truth_phi, y_truth_pt, yPrime_truth_pt, y_truth_e, yPrime_truth_e ;
    Int_t HLT_g140_loose_L1EM22VHI, y_convType, yPrime_convType;
    Tree->SetBranchAddress("y_pt", &y_pt);
    Tree->SetBranchAddress("yPrime_pt", &yPrime_pt);
    Tree->SetBranchAddress("y_eta", &y_eta);
    Tree->SetBranchAddress("yPrime_eta", &yPrime_eta);
    Tree->SetBranchAddress("y_phi", &y_phi);
    Tree->SetBranchAddress("yPrime_phi", &yPrime_phi);
    Tree->SetBranchAddress("y_e", &y_e);
    Tree->SetBranchAddress("yPrime_e", &yPrime_e);
    Tree->SetBranchAddress("HLT_g140_loose_L1EM22VHI", &HLT_g140_loose_L1EM22VHI);
    Tree->SetBranchAddress("y_truth_eta", &y_truth_eta);
    Tree->SetBranchAddress("yPrime_truth_eta", &yPrime_truth_eta);
    Tree->SetBranchAddress("y_truth_phi", &y_truth_phi);
    Tree->SetBranchAddress("yPrime_truth_phi", &yPrime_truth_phi);

    Tree->SetBranchAddress("y_truth_pt", &y_truth_pt);
    Tree->SetBranchAddress("yPrime_truth_pt", &yPrime_truth_pt);
    Tree->SetBranchAddress("y_truth_eta", &y_truth_eta);
    Tree->SetBranchAddress("yPrime_truth_eta", &yPrime_truth_eta);
    Tree->SetBranchAddress("y_truth_phi", &y_truth_phi);
    Tree->SetBranchAddress("yPrime_truth_phi", &yPrime_truth_phi);
    Tree->SetBranchAddress("y_truth_e", &y_truth_e);
    Tree->SetBranchAddress("yPrime_truth_e", &yPrime_truth_e);

    Tree->SetBranchAddress("y_convType", &y_convType);
    Tree->SetBranchAddress("yPrime_convType", &yPrime_convType);



    // Boucle sur toutes les entrées de l'arbre
    int nEntries = Tree->GetEntries();
    cout << "le nombre d entries est "  << nEntries << endl;
    for (int iEntry = 0; iEntry < nEntries; ++iEntry) {
        // Charger l'entrée courante
        Tree->GetEntry(iEntry);

        if ((y_convType!=0 || yPrime_convType!=0) && y_pt>0 && yPrime_pt>0 && HLT_g140_loose_L1EM22VHI==1 && TMath::Abs(y_eta-y_truth_eta)<0.1 && TMath::Abs(yPrime_eta-yPrime_truth_eta)<0.1 && TMath::Abs(y_phi-y_truth_phi)<0.1 && TMath::Abs(yPrime_phi-yPrime_truth_phi)<0.1 && TMath::Abs(y_eta)<2.37 && TMath::Abs(yPrime_eta)<2.37 && TMath::Abs(y_eta)>1.52 && TMath::Abs(yPrime_eta)>1.52)  {

        //Vecteurs des deux photons
        TLorentzVector v1, v2, v3, v4;
        v1.SetPtEtaPhiE(y_pt,y_eta,y_phi,y_e);
        v2.SetPtEtaPhiE(yPrime_pt,yPrime_eta,yPrime_phi,yPrime_e);
        v3.SetPtEtaPhiE(y_truth_pt,y_truth_eta,y_truth_phi,y_truth_e);
        v4.SetPtEtaPhiE(yPrime_truth_pt,yPrime_truth_eta,yPrime_truth_phi,yPrime_truth_e);

        // Vecteur de l'action
        TLorentzVector axion = v1 + v2;
        TLorentzVector axion_truth = v3 + v4;

        // Calcul de la masse invariante
        mass.setVal(axion_truth.M());
        diff.setVal(axion.M()-axion_truth.M());

        // Remplir le dataset
        data.add(vars);

        }
    }

    data.Print("V");

    TH1D *h_mass = new TH1D("Résolution de la masse", "Résolution de la masse", 100, -6, 6);
    h_mass->SetXTitle("m_reco-m_truth (GeV)");  
    h_mass->SetTitle("Resolution de la masse !(uu)-endcap pour l'intervalle [5 , 35] ");
    TH1D *h = new TH1D("Mass distribution", "Mass Distribution", 100, 0, 45);
    h->SetXTitle("m_truth (GeV)");  
    h->SetTitle("Mass distribution !(uu)-endcap ");
    // Garder les points de bonne masse
    for (int j = 0; j < data.numEntries(); ++j) {
        mass.setVal(data.get(j)->getRealValue("mass"));
        diff.setVal(data.get(j)->getRealValue("diff"));
        double mass_val = mass.getVal();
        double diff_val = diff.getVal();
        h->Fill(mass_val);

        if (mass_val >= 5 && mass_val <= 35  ) {
            h_mass->Fill(diff_val);
        }
    }






    int loop_count = 0; 

    // Boucle pour créer 50 histogrammes
    for (int i = 5; i < 10; ++i) {

        // Réinitialiser les valeurs des paramètres
	    m.setVal( m_init ) ;
	    sigmaL.setVal(sigmaL_init) ;
        sigmaR.setVal(sigmaR_init);
        alphaL.setVal(alphaL_init);
        alphaR.setVal(alphaR_init);
        nL.setVal(nL_init) ;
        nR.setVal(nR_init) ;


        char data_name[30];  // Tableau pour construire un nom unique
        snprintf(data_name, sizeof(data_name), "data_diff_%d", i); 
        RooDataSet *data_diff = new RooDataSet(data_name, "data_diff", vars_diff);

        // Création d'un histogramme 
        char hist_name[30];  // Tableau pour construire un nom unique
        snprintf(hist_name, sizeof(hist_name), "h_mass_%d", i); 
        TH1D *h_mass = new TH1D(hist_name, "Résolution de la masse", 100, -6, 6);
        h_mass->SetXTitle("m_reco-m_truth (GeV)");  
        h_mass->SetTitle(Form("Resolution de la masse pour l'intervalle [%d , %d] ", i, i+1));

        //int loop_count = 0; 

        // Garder les points de bonne masse
        for (int j = 0; j < data.numEntries(); ++j) {
            mass.setVal(data.get(j)->getRealValue("mass"));
            diff.setVal(data.get(j)->getRealValue("diff"));
            double mass_val = mass.getVal();
            double diff_val = diff.getVal();

            if (mass_val >= i && mass_val <= i + 1  ) {
                //std::cout << "mass " << mass << " diff " << diff << std::endl;
                h_mass->Fill(diff_val);
                data_diff->add(vars_diff);
                //++loop_count;
            }

        }



        // Fit 
        m.setConstant(kFALSE);
        sigmaL.setConstant(kFALSE);
        sigmaR.setConstant(kFALSE);
        alphaL.setConstant(kFALSE);
        alphaR.setConstant(kFALSE);
        nL.setConstant(kFALSE);
        nR.setConstant(kFALSE);

        RooFitResult *myFit = PDF.fitTo(*data_diff, Save());

        myFit->Print( "v" ) ;

        // Extraire les valeurs ajustées
        float sigmaL_val = sigmaL.getVal();
        float sigmaR_val = sigmaR.getVal();
        float alphaL_val = alphaL.getVal();
        float alphaR_val = alphaR.getVal();
        float nL_val = nL.getVal();
        float nR_val = nR.getVal();
        float m_val = m.getVal();
        float sigmaL_val_error = sigmaL.getError();
        float sigmaR_val_error = sigmaR.getError();
        float alphaL_val_error = alphaL.getError();
        float alphaR_val_error = alphaR.getError();
        float nL_val_error = nL.getError();
        float nR_val_error = nR.getError();
        float m_val_error = m.getError();

        // Ajouter l'écart-type à notre graphique
        sigmaLGraph->SetPoint(i, i, sigmaL_val); 
        sigmaLGraph->SetPointError(i, 0, sigmaL_val_error); 

        sigmaRGraph->SetPoint(i, i, sigmaR_val); 
        sigmaRGraph->SetPointError(i, 0, sigmaR_val_error); 
        
        alphaLGraph->SetPoint(i, i, alphaL_val); 
        alphaLGraph->SetPointError(i, 0, alphaL_val_error); 

        alphaRGraph->SetPoint(i, i, alphaR_val); 
        alphaRGraph->SetPointError(i, 0, alphaR_val_error); 

        nLGraph->SetPoint(i, i, nL_val); 
        nLGraph->SetPointError(i, 0, nL_val_error); 

        nRGraph->SetPoint(i, i, nR_val); 
        nRGraph->SetPointError(i, 0, nR_val_error); 

        mGraph->SetPoint(i, i, m_val); 
        mGraph->SetPointError(i, 0, m_val_error); 

        //std::cout << "La boucle a été exécutée " << loop_count << " fois." << std::endl;

        RooPlot * myFrame = diff.frame( 50 ) ;
        myFrame->GetXaxis()->SetLimits(-15,15);
        //myFrame->GetXaxis()->SetLimits(m_val-5*TMath::Max(s1_val, s2_val), m_val+5*TMath::Max(s1_val, s2_val));
        data_diff->plotOn( myFrame ) ;
        PDF.plotOn( myFrame ) ;
        myFrame->SetXTitle("m_reco-m_truth (GeV)");  
        myFrame->SetTitle(Form("Resolution de la masse pour l'intervalle [%d , %d] ", i, i+1));

        // Dessiner l'histogramme sur le canvas
        canvas->cd(i - 4); 
        gPad->SetLogy();
        myFrame->Draw();
    }

    // Dessiner le Tgraph des écarts-types sur un autre Canvas 
    TCanvas* canvas_sGraph = new TCanvas("canvas_sGraph", "Canvas_sGraph", 800, 400);
    canvas_sGraph -> Divide(2,1);

    canvas_sGraph->cd(1);
    sigmaLGraph->GetXaxis()->SetLimits(5, 35);
    sigmaLGraph->SetTitle("Left standard deviation evolution");
    sigmaLGraph->GetXaxis()->SetTitle("Mass (GeV)");
    sigmaLGraph->GetYaxis()->SetTitle(" Standard deviation ");
    sigmaLGraph->Draw("AP");

    canvas_sGraph->cd(2);
    sigmaRGraph->GetXaxis()->SetLimits(5, 35);
    sigmaRGraph->SetTitle("Right standard deviation evolution");
    sigmaRGraph->GetXaxis()->SetTitle("Mass (GeV)");
    sigmaRGraph->GetYaxis()->SetTitle(" Standard deviation ");
    sigmaRGraph->Draw("AP");

    // Dessiner le Tgraph des alpha sur un autre Canvas 
    TCanvas* canvas_aGraph = new TCanvas("canvas_aGraph", "Canvas_aGraph", 800, 400);
    canvas_aGraph -> Divide(2,1);

    canvas_aGraph->cd(1);
    alphaLGraph->GetXaxis()->SetLimits(5, 35);
    alphaLGraph->SetTitle("Location of the left transition evolution");
    alphaLGraph->GetXaxis()->SetTitle("Mass (GeV)");
    alphaLGraph->GetYaxis()->SetTitle(" Location of transition to a power law on the left ");
    alphaLGraph->Draw("AP");

    canvas_aGraph->cd(2);
    alphaRGraph->GetXaxis()->SetLimits(5, 35);
    alphaRGraph->SetTitle("Location of the right transition evolution");
    alphaRGraph->GetXaxis()->SetTitle("Mass (GeV)");
    alphaRGraph->GetYaxis()->SetTitle(" Location of transition to a power law on the right ");
    alphaRGraph->Draw("AP");


    // Dessiner le Tgraph des n sur un autre Canvas 
    TCanvas* canvas_nGraph = new TCanvas("canvas_nGraph", "Canvas_nGraph", 800, 400);
    canvas_nGraph -> Divide(2,1);

    canvas_nGraph->cd(1);
    nLGraph->GetXaxis()->SetLimits(5, 35);
    nLGraph->SetTitle("Exponent of power-law tail on the left evolution");
    nLGraph->GetXaxis()->SetTitle("Mass (GeV)");
    nLGraph->GetYaxis()->SetTitle(" Exponent of power-law tail on the left ");
    nLGraph->Draw("AP");

    canvas_nGraph->cd(2);
    nRGraph->GetXaxis()->SetLimits(5, 35);
    nRGraph->SetTitle("Exponent of power-law tail on the right evolution");
    nRGraph->GetXaxis()->SetTitle("Mass (GeV)");
    nRGraph->GetYaxis()->SetTitle(" Exponent of power-law tail on the right ");
    nRGraph->Draw("AP");


    // Dessiner le Tgraph de la moyenne sur un autre Canvas 
    TCanvas* canvas_mGraph = new TCanvas("canvas_mGraph", "Canvas_mGraph", 400, 400);
    mGraph->GetXaxis()->SetLimits(5, 35);
    mGraph->SetTitle("Mean evolution");
    mGraph->GetXaxis()->SetTitle("Mass (GeV)");
    mGraph->GetYaxis()->SetTitle(" Mean ");
    mGraph->Draw("AP");

    // Dessiner le diagramme de résolution en masse
    TCanvas* canvas_h_mass = new TCanvas("canvas_h_mass", "Canvas_h_mass", 400, 400);
    h_mass->Draw();

    TCanvas* canvas_h = new TCanvas("canvas_h", "Canvas_h", 400, 400);
    h->Draw();


    // Fit
	
	//sig_m.setConstant( kFALSE ) ;
	//sig_s.setConstant( kFALSE ) ;

	//RooFitResult * myFit = sig_PDF.fitTo( data , Save() ) ;

	//myFit->Print( "v" ) ;


    // Créer un canvas pour afficher l'histogramme et la PDF

    //RooPlot * myFrame = mass.frame( 100 ) ;
	//data.plotOn( myFrame ) ;
	//sig_PDF.plotOn( myFrame ) ;


    canvas -> SaveAs("Resolution_mass_!(uu)-endcap_CB_log_5_10.pdf");

    //canvas_sGraph -> SaveAs("Evolution_standard-deviations_!(uu)-endcap_CB.pdf");

    //canvas_aGraph -> SaveAs("Evolution_alpha_!(uu)-endcap_CB.pdf");

    //canvas_nGraph -> SaveAs("Evolution_n_!(uu)-endcap_CB.pdf");

    //canvas_mGraph -> SaveAs("Evolution_mean_!(uu)-endcap_CB.pdf");

    //canvas_h_mass -> SaveAs("Resolution_mass_!(uu)-endcap.pdf");

    //canvas_h -> SaveAs("Mass_Distribution_!(uu)-endcap.pdf");

    myFile->Close();


}
