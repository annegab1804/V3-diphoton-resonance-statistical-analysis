
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
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TGraphErrors.h"


using namespace std;

using namespace RooFit ;

void resolution_mass() {

    // Récupération des données
    std::unique_ptr<TFile> myFile( TFile::Open("diphoton_ntuple_collimated_yy.root") );
    auto Tree = myFile->Get<TTree>("diphotonTree");

    // Création de la dataset
    RooRealVar mass("mass", "mass", 0, 100, "GeV");
    RooRealVar diff("diff", "diff", -10, 10, "GeV");
    RooArgSet vars(mass, diff);
    RooArgSet vars_diff(diff);
    RooDataSet data("data", "data", vars);
    RooDataSet data_barrel("data_barrel", "data_barrel", vars);
    RooDataSet data_endcap("data_endcap", "data_endcap", vars);


    // Dessiner les histogrammes/pdf sur un canvas
    TCanvas* canvas = new TCanvas("canvas", "Canvas", 2000, 400);
    canvas->Divide(5, 1, 0.01, 0.1);


    // Ajouter un titre
    TPaveText *title = new TPaveText(0.1, 0.94, 0.9, 0.98, "NDC"); // Coordonnées en NDC (Normalized Device Coordinates)
    title->SetFillColor(0);
    title->SetTextAlign(22);  // Centré horizontalement et verticalement
    title->AddText("Mass resolution UU-barrel [30;35]");
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

        if (y_pt>5 && yPrime_pt>20 && y_pt+ yPrime_pt > 140 && HLT_g140_loose_L1EM22VHI==1 && TMath::Abs(y_eta-y_truth_eta)<0.1 && TMath::Abs(yPrime_eta-yPrime_truth_eta)<0.1 && TMath::Abs(y_phi-y_truth_phi)<0.1 && TMath::Abs(yPrime_phi-yPrime_truth_phi)<0.1 )  {
        //if (y_pt>0 && yPrime_pt>0 && HLT_g140_loose_L1EM22VHI==1 && TMath::Abs(y_eta-y_truth_eta)<0.1 && TMath::Abs(yPrime_eta-yPrime_truth_eta)<0.1 && TMath::Abs(y_phi-y_truth_phi)<0.1 && TMath::Abs(yPrime_phi-yPrime_truth_phi)<0.1 )  {

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

        if (y_convType==0 && yPrime_convType==0 && TMath::Abs(y_eta)<1.37 && TMath::Abs(yPrime_eta)<1.37){
            data_barrel.add(vars);
        }

        if ((y_convType!=0 || yPrime_convType!=0) && TMath::Abs(y_eta)<2.37 && TMath::Abs(yPrime_eta)<2.37 && TMath::Abs(y_eta)>1.52 && TMath::Abs(yPrime_eta)>1.52){
            data_endcap.add(vars);
        }


        }
    }

    data.Print("V");

    TH1D *h_mass = new TH1D("Résolution de la masse", "Résolution de la masse", 100, -6, 6);
    TH1D *h_mass_barrel = new TH1D("Résolution de la masse uu-barrel", "Résolution de la masse uu-barrel", 100, -6, 6);
    TH1D *h_mass_endcap = new TH1D("Résolution de la masse !(uu)-endcap", "Résolution de la masse !(uu)-endcap", 100, -6, 6);
    h_mass->SetXTitle("m_reco-m_truth (GeV)");  
    h_mass->SetYTitle("events / GeV");  
    h_mass->SetTitle("Mass resolution for the range [5 , 35] ");
    h_mass->GetXaxis()->SetRangeUser(-4, 4);

    TH1D *h = new TH1D("Mass distribution", "Mass Distribution", 100, 0, 45);
    TH1D *h_barrel = new TH1D("Mass distribution uu-barrel", "Mass Distribution uu-barrel", 100, 0, 45);
    TH1D *h_endcap = new TH1D("Mass distribution !(uu)-endcap", "Mass Distribution !(uu)-endcap", 100, 0, 45);
    h->SetXTitle("m_truth (GeV)");  
    h->SetYTitle("events / GeV");  
    h->SetTitle("Mass distribution ");

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

    for (int j = 0; j < data_barrel.numEntries(); ++j) {
        mass.setVal(data_barrel.get(j)->getRealValue("mass"));
        diff.setVal(data_barrel.get(j)->getRealValue("diff"));
        double mass_val = mass.getVal();
        double diff_val = diff.getVal();
        h_barrel->Fill(mass_val);

        if (mass_val >= 5 && mass_val <= 35  ) {
            h_mass_barrel->Fill(diff_val);
        }
    }

    for (int j = 0; j < data_endcap.numEntries(); ++j) {
        mass.setVal(data_endcap.get(j)->getRealValue("mass"));
        diff.setVal(data_endcap.get(j)->getRealValue("diff"));
        double mass_val = mass.getVal();
        double diff_val = diff.getVal();
        h_endcap->Fill(mass_val);

        if (mass_val >= 5 && mass_val <= 35  ) {
            h_mass_endcap->Fill(diff_val);
        }
    }





    int loop_count = 0; 

    // Boucle pour créer 5 histogrammes
    for (int i = 30; i < 35; ++i) {

        
        char data_name[35];  // Tableau pour construire un nom unique
        snprintf(data_name, sizeof(data_name), "data_diff_%d", i); 
        RooDataSet *data_diff = new RooDataSet(data_name, "data_diff", vars_diff);

        // Création d'un histogramme 
        char hist_name[35];  // Tableau pour construire un nom unique
        snprintf(hist_name, sizeof(hist_name), "h_mass_%d", i); 
        TH1D *h_mass = new TH1D(hist_name, "Résolution de la masse", 100, -6, 6);

        h_mass->SetXTitle("m_reco-m_truth (GeV)");  
        h_mass->SetTitle(Form("Resolution de la masse pour l'intervalle [%d , %d] ", i, i+1));

        //int loop_count = 0; 

        // Garder les points de bonne masse
        for (int j = 0; j < data_barrel.numEntries(); ++j) {
            mass.setVal(data_barrel.get(j)->getRealValue("mass"));
            diff.setVal(data_barrel.get(j)->getRealValue("diff"));
            double mass_val = mass.getVal();
            double diff_val = diff.getVal();

            if (mass_val >= i && mass_val <= i + 1 ) {
                //std::cout << "mass " << mass << " diff " << diff << std::endl;
                h_mass->Fill(diff_val);
                data_diff->add(vars_diff);
                //++loop_count;
            }

        }


        RooPlot * myFrame = diff.frame( 50 ) ;
        myFrame->GetXaxis()->SetLimits(-10, 10);
        data_diff->plotOn( myFrame ) ;
        myFrame->SetXTitle("m_reco-m_truth (GeV)");  
        myFrame->SetTitle(Form("Resolution de la masse pour l'intervalle [%d , %d] ", i, i+1));

        // Dessiner l'histogramme sur le canvas
        canvas->cd(i - 29); 
        gPad->SetLogy();
        myFrame->Draw();
    }


    // Dessiner le diagramme de résolution en masse
    TCanvas* canvas_h_mass = new TCanvas("canvas_h_mass", "Canvas_h_mass", 400, 400);
    canvas_h_mass->SetLeftMargin(0.15);
    canvas_h_mass->SetRightMargin(0.05);
    
    h_mass->SetLineColor(kBlack);
    h_mass->SetStats(0);
    h_mass->GetXaxis()->SetTitle("m_{reco}-m_{truth}  [GeV]");
    h_mass->GetYaxis()->SetTitle("events / GeV");

    h_mass_barrel->SetLineColor(kRed);
    h_mass_barrel->SetStats(0);
    h_mass_endcap->SetLineColor(kBlue);
    h_mass_endcap->SetStats(0);

    h_mass_barrel->SetXTitle("m_{reco}-m_{truth} (GeV)");  
    h_mass_barrel->SetYTitle("#frac{1}{N} #frac{dN}{d m}");  
    h_mass_barrel->SetTitle("Normalised mass resolution for the range [5 , 35] ");
    h_mass_barrel->GetXaxis()->SetRangeUser(-4, 4);

    h_mass_barrel->DrawNormalized();
    h_mass->DrawNormalized("Same");
    h_mass_endcap->DrawNormalized("Same");


    TLegend* legend_mass = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend_mass->AddEntry(h_mass, "All", "l");
    legend_mass->AddEntry(h_mass_barrel, "uu-barrel", "l");
    legend_mass->AddEntry(h_mass_endcap, "!(uu)-endcap", "l");
    legend_mass->Draw();

    TCanvas* canvas_h = new TCanvas("canvas_h", "Canvas_h", 400, 400);
    canvas_h->SetLeftMargin(0.15);
    canvas_h->SetRightMargin(0.05);

    h->SetStats(0);
    h->GetXaxis()->SetTitle("m_{#gamma#gamma}  [GeV]");
    h->GetYaxis()->SetTitle("events / GeV");
    h->SetLineColor(kBlack);
    h_barrel->SetLineColor(kRed);
    h_endcap->SetLineColor(kBlue);

    
    h->Draw();
    h_barrel->Draw("Same");
    h_endcap->Draw("Same");

    TLegend* legend = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend->AddEntry(h, "All", "l");
    legend->AddEntry(h_barrel, "uu-barrel", "l");
    legend->AddEntry(h_endcap, "!(uu)-endcap", "l");
    legend->Draw();



    //canvas -> SaveAs("Resolution_mass_uu-barrel_log_30-35_Cut.pdf");

    canvas_h_mass -> SaveAs("Resolutions_mass_Cut_norm.pdf");

    canvas_h -> SaveAs("Mass_Distributions_Cut.pdf");

    myFile->Close();


}
