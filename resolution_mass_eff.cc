
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

void resolution_mass_eff() {

    // Récupération des données
    std::unique_ptr<TFile> myFile( TFile::Open("diphoton_ntuple_collimated_yy.root") );
    auto Tree = myFile->Get<TTree>("diphotonTree");

    //TFile *f_filteredTree = new TFile("f_filteredTree.root", "RECREATE");
    //TTree *filteredTree = Tree->CopyTree("y_pt>0 && yPrime_pt>0 && HLT_g140_loose_L1EM22VHI==1 && TMath::Abs(y_eta-y_truth_eta)<0.1 && TMath::Abs(yPrime_eta-yPrime_truth_eta)<0.1 && TMath::Abs(y_phi-y_truth_phi)<0.1 && TMath::Abs(yPrime_phi-yPrime_truth_phi)<0.1");

    // Création de la dataset
    RooRealVar mass("mass", "mass", 0, 100, "GeV");
    RooRealVar diff("diff", "diff", -10, 10, "GeV");
    RooArgSet vars(mass, diff);
    RooArgSet vars_diff(diff);
    RooDataSet data("data", "data", vars);
    RooDataSet data_barrel("data_barrel", "data_barrel", vars);
    RooDataSet data_barrel_Cut10("data_barrel_Cut10", "data_barrel_Cut10", vars);
    RooDataSet data_barrel_Cut15("data_barrel_Cut15", "data_barrel_Cut15", vars);
    RooDataSet data_barrel_Cut20("data_barrel_Cut20", "data_barrel_Cut20", vars);
    RooDataSet data_barrel_Cut25("data_barrel_Cut25", "data_barrel_Cut25", vars);
    RooDataSet data_barrel_Cut30("data_barrel_Cut30", "data_barrel_Cut30", vars);
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

        if (y_pt>5 && yPrime_pt>5 && y_pt+ yPrime_pt > 140 && HLT_g140_loose_L1EM22VHI==1 && TMath::Abs(y_eta-y_truth_eta)<0.1 && TMath::Abs(yPrime_eta-yPrime_truth_eta)<0.1 && TMath::Abs(y_phi-y_truth_phi)<0.1 && TMath::Abs(yPrime_phi-yPrime_truth_phi)<0.1 )  {
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

        if (y_pt>5 && yPrime_pt>10 && y_pt+ yPrime_pt > 140 && HLT_g140_loose_L1EM22VHI==1 && TMath::Abs(y_eta-y_truth_eta)<0.1 && TMath::Abs(yPrime_eta-yPrime_truth_eta)<0.1 && TMath::Abs(y_phi-y_truth_phi)<0.1 && TMath::Abs(yPrime_phi-yPrime_truth_phi)<0.1 )  {
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


        if (y_convType==0 && yPrime_convType==0 && TMath::Abs(y_eta)<1.37 && TMath::Abs(yPrime_eta)<1.37){
            data_barrel_Cut10.add(vars);
        }


        }

        if (y_pt>5 && yPrime_pt>15 && y_pt+ yPrime_pt > 140 && HLT_g140_loose_L1EM22VHI==1 && TMath::Abs(y_eta-y_truth_eta)<0.1 && TMath::Abs(yPrime_eta-yPrime_truth_eta)<0.1 && TMath::Abs(y_phi-y_truth_phi)<0.1 && TMath::Abs(yPrime_phi-yPrime_truth_phi)<0.1 )  {
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

        

        if (y_convType==0 && yPrime_convType==0 && TMath::Abs(y_eta)<1.37 && TMath::Abs(yPrime_eta)<1.37){
            data_barrel_Cut15.add(vars);
        }


        }

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


        if (y_convType==0 && yPrime_convType==0 && TMath::Abs(y_eta)<1.37 && TMath::Abs(yPrime_eta)<1.37){
            data_barrel_Cut20.add(vars);
        }



        }

        if (y_pt>5 && yPrime_pt>25 && y_pt+ yPrime_pt > 140 && HLT_g140_loose_L1EM22VHI==1 && TMath::Abs(y_eta-y_truth_eta)<0.1 && TMath::Abs(yPrime_eta-yPrime_truth_eta)<0.1 && TMath::Abs(y_phi-y_truth_phi)<0.1 && TMath::Abs(yPrime_phi-yPrime_truth_phi)<0.1 )  {
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



        if (y_convType==0 && yPrime_convType==0 && TMath::Abs(y_eta)<1.37 && TMath::Abs(yPrime_eta)<1.37){
            data_barrel_Cut25.add(vars);
        }



        }

        if (y_pt>5 && yPrime_pt>30 && y_pt+ yPrime_pt > 140 && HLT_g140_loose_L1EM22VHI==1 && TMath::Abs(y_eta-y_truth_eta)<0.1 && TMath::Abs(yPrime_eta-yPrime_truth_eta)<0.1 && TMath::Abs(y_phi-y_truth_phi)<0.1 && TMath::Abs(yPrime_phi-yPrime_truth_phi)<0.1 )  {
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


        if (y_convType==0 && yPrime_convType==0 && TMath::Abs(y_eta)<1.37 && TMath::Abs(yPrime_eta)<1.37){
            data_barrel_Cut30.add(vars);
        }



        }
    }

    data.Print("V");

    TH1D *h_mass = new TH1D("Résolution de la masse", "Résolution de la masse", 100, -6, 6);
    TH1D *h_mass_barrel = new TH1D("Résolution de la masse uu-barrel ", "Résolution de la masse uu-barrel", 100, -6, 6);
    TH1D *h_mass_barrel_Cut10 = new TH1D("Résolution de la masse uu-barrel Cut 10", "Résolution de la masse uu-barrel Cut 10", 100, -6, 6);
    TH1D *h_mass_barrel_Cut15 = new TH1D("Résolution de la masse uu-barrel Cut 15", "Résolution de la masse uu-barrel Cut 15", 100, -6, 6);
    TH1D *h_mass_barrel_Cut20 = new TH1D("Résolution de la masse uu-barrel Cut 20", "Résolution de la masse uu-barrel Cut 20", 100, -6, 6);
    TH1D *h_mass_barrel_Cut25 = new TH1D("Résolution de la masse uu-barrel Cut 25", "Résolution de la masse uu-barrel Cut 25", 100, -6, 6);
    TH1D *h_mass_barrel_Cut30 = new TH1D("Résolution de la masse uu-barrel Cut 30", "Résolution de la masse uu-barrel Cut 30", 100, -6, 6);
    TH1D *h_mass_endcap = new TH1D("Résolution de la masse !(uu)-endcap", "Résolution de la masse !(uu)-endcap", 100, -6, 6);
    h_mass->SetXTitle("m_reco-m_truth (GeV)");  
    h_mass->SetYTitle("events / GeV");  
    h_mass->SetTitle("Mass resolution for the range [5 , 35] ");
    h_mass->GetXaxis()->SetRangeUser(-4, 4);

    TH1D *h = new TH1D("Mass distribution", "Mass Distribution", 100, 0, 45);
    TH1D *h_barrel = new TH1D("Signal mass distribution", "Signal mass distribution", 100, 0, 45);
    TH1D *h_barrel_Cut10 = new TH1D("Mass distribution uu-barrel Cut 10", "Mass distribution uu-barrel Cut 10", 100, 0, 45);
    TH1D *h_barrel_Cut15 = new TH1D("Mass distribution uu-barrel Cut 15", "Mass distribution uu-barrel Cut 15", 100, 0, 45);
    TH1D *h_barrel_Cut20 = new TH1D("Mass distribution uu-barrel Cut 20", "Mass distribution uu-barrel Cut 20", 100, 0, 45);
    TH1D *h_barrel_Cut25 = new TH1D("Mass distribution uu-barrel Cut 25", "Mass distribution uu-barrel Cut 25", 100, 0, 45);
    TH1D *h_barrel_Cut30 = new TH1D("Mass distribution uu-barrel Cut 30", "Mass distribution uu-barrel Cut 30", 100, 0, 45);
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

    for (int j = 0; j < data_barrel_Cut10.numEntries(); ++j) {
        mass.setVal(data_barrel_Cut10.get(j)->getRealValue("mass"));
        diff.setVal(data_barrel_Cut10.get(j)->getRealValue("diff"));
        double mass_val = mass.getVal();
        double diff_val = diff.getVal();
        h_barrel_Cut10->Fill(mass_val);

        if (mass_val >= 5 && mass_val <= 35  ) {
            h_mass_barrel_Cut10->Fill(diff_val);
        }
    }

    for (int j = 0; j < data_barrel_Cut15.numEntries(); ++j) {
        mass.setVal(data_barrel_Cut15.get(j)->getRealValue("mass"));
        diff.setVal(data_barrel_Cut15.get(j)->getRealValue("diff"));
        double mass_val = mass.getVal();
        double diff_val = diff.getVal();
        h_barrel_Cut15->Fill(mass_val);

        if (mass_val >= 5 && mass_val <= 35  ) {
            h_mass_barrel_Cut15->Fill(diff_val);
        }
    }

    for (int j = 0; j < data_barrel_Cut20.numEntries(); ++j) {
        mass.setVal(data_barrel_Cut20.get(j)->getRealValue("mass"));
        diff.setVal(data_barrel_Cut20.get(j)->getRealValue("diff"));
        double mass_val = mass.getVal();
        double diff_val = diff.getVal();
        h_barrel_Cut20->Fill(mass_val);

        if (mass_val >= 5 && mass_val <= 35  ) {
            h_mass_barrel_Cut20->Fill(diff_val);
        }
    }

    for (int j = 0; j < data_barrel_Cut25.numEntries(); ++j) {
        mass.setVal(data_barrel_Cut25.get(j)->getRealValue("mass"));
        diff.setVal(data_barrel_Cut25.get(j)->getRealValue("diff"));
        double mass_val = mass.getVal();
        double diff_val = diff.getVal();
        h_barrel_Cut25->Fill(mass_val);

        if (mass_val >= 5 && mass_val <= 35  ) {
            h_mass_barrel_Cut25->Fill(diff_val);
        }
    }

    for (int j = 0; j < data_barrel_Cut30.numEntries(); ++j) {
        mass.setVal(data_barrel_Cut30.get(j)->getRealValue("mass"));
        diff.setVal(data_barrel_Cut30.get(j)->getRealValue("diff"));
        double mass_val = mass.getVal();
        double diff_val = diff.getVal();
        h_barrel_Cut30->Fill(mass_val);

        if (mass_val >= 5 && mass_val <= 35  ) {
            h_mass_barrel_Cut30->Fill(diff_val);
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

    h_mass_barrel_Cut10->SetLineColor(51 );
    h_mass_barrel_Cut15->SetLineColor(51 + 6*1);
    h_mass_barrel_Cut20->SetLineColor(51 + 6*2);
    h_mass_barrel_Cut25->SetLineColor(51 + 6*3);
    h_mass_barrel_Cut30->SetLineColor(51 + 6*4);

    TCanvas * canvas_h_mass_Cuts = new TCanvas("canvas_h_mass_Cuts", "canvas_h_mass_Cuts", 400, 400);
    canvas_h_mass_Cuts->SetRightMargin(0.12);

    auto legend_mass_Cuts = new TLegend(0.8,0.5,0.88,0.9);
    legend_mass_Cuts->SetHeader("Slice","C");

    h_mass_barrel->SetStats(0);
    h_mass_barrel->GetXaxis()->SetTitle("m_{reco}-m_{truth}  [GeV]");
    h_mass_barrel->GetYaxis()->SetTitle("events / GeV");
    h_mass_barrel->Draw("");
    legend_mass_Cuts->AddEntry(h_mass_barrel,"pT > 5 GeV");

    h_mass_barrel_Cut10->SetStats(0);
    h_mass_barrel_Cut10->GetXaxis()->SetTitle("m_{reco}-m_{truth}  [GeV]");
    h_mass_barrel_Cut10->GetYaxis()->SetTitle("events / GeV");
    h_mass_barrel_Cut10->Draw("Same");
    legend_mass_Cuts->AddEntry(h_mass_barrel_Cut10,"pT > 10 GeV");

    h_mass_barrel_Cut15->SetStats(0);
    h_mass_barrel_Cut15->GetXaxis()->SetTitle("m_{reco}-m_{truth}  [GeV]");
    h_mass_barrel_Cut15->GetYaxis()->SetTitle("events / GeV");
    h_mass_barrel_Cut15->Draw("Same");
    legend_mass_Cuts->AddEntry(h_mass_barrel_Cut15,"pT > 15 GeV");

    h_mass_barrel_Cut20->SetStats(0);
    h_mass_barrel_Cut20->GetXaxis()->SetTitle("m_{reco}-m_{truth}  [GeV]");
    h_mass_barrel_Cut20->GetYaxis()->SetTitle("events / GeV");
    h_mass_barrel_Cut20->Draw("Same");
    legend_mass_Cuts->AddEntry(h_mass_barrel_Cut20,"pT > 20 GeV");

    h_mass_barrel_Cut25->SetStats(0);
    h_mass_barrel_Cut25->GetXaxis()->SetTitle("m_{reco}-m_{truth}  [GeV]");
    h_mass_barrel_Cut25->GetYaxis()->SetTitle("events / GeV");
    h_mass_barrel_Cut25->Draw("Same");
    legend_mass_Cuts->AddEntry(h_mass_barrel_Cut25,"pT > 25 GeV");

    h_mass_barrel_Cut30->SetStats(0);
    h_mass_barrel->GetXaxis()->SetTitle("m_{reco}-m_{truth}  [GeV]");
    h_mass_barrel->GetYaxis()->SetTitle("events / GeV");
    h_mass_barrel->Draw("Same");
    legend_mass_Cuts->AddEntry(h_mass_barrel_Cut30,"pT > 30 GeV");

    legend_mass_Cuts->Draw();

    //canvas_h_mass_Cuts->SaveAs("canvas_h_mass_Cuts.pdf");

    h_barrel_Cut10->SetLineColor(51 );
    h_barrel_Cut15->SetLineColor(51 + 6*1);
    h_barrel_Cut20->SetLineColor(51 + 6*2);
    h_barrel_Cut25->SetLineColor(51 + 6*3);
    h_barrel_Cut30->SetLineColor(51 + 6*4);

    TCanvas * canvas_h_Cuts = new TCanvas("canvas_h_Cuts", "canvas_h_Cuts", 400, 400);
    canvas_h_Cuts->SetLeftMargin(0.15);
    canvas_h_Cuts->SetRightMargin(0.05);

    auto legend_Cuts = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend_Cuts->SetHeader("Slice","C");

    h_barrel->SetStats(0);
    h_barrel->GetXaxis()->SetTitle("m_{#gamma#gamma} [GeV]");
    h_barrel->GetYaxis()->SetTitle("events / GeV");
    h_barrel->Draw("");
    legend_Cuts->AddEntry(h_barrel,"pT > 5 GeV");

    h_barrel_Cut10->GetXaxis()->SetTitle("m_{#gamma#gamma} [GeV]");
    h_barrel_Cut10->SetStats(0);
    h_barrel_Cut10->GetYaxis()->SetTitle("events / GeV");
    h_barrel_Cut10->Draw("Same");
    legend_Cuts->AddEntry(h_barrel_Cut10,"pT > 10 GeV");

    h_barrel_Cut15->SetStats(0);
    h_barrel_Cut15->GetXaxis()->SetTitle("m_{#gamma#gamma}[GeV]");
    h_barrel_Cut15->GetYaxis()->SetTitle("events / GeV");
    h_barrel_Cut15->Draw("Same");
    legend_Cuts->AddEntry(h_barrel_Cut15,"pT > 15 GeV");

    h_barrel_Cut20->SetStats(0);
    h_barrel_Cut20->GetXaxis()->SetTitle("m_{#gamma#gamma}  [GeV]");
    h_barrel_Cut20->GetYaxis()->SetTitle("events / GeV");
    h_barrel_Cut20->Draw("Same");
    legend_Cuts->AddEntry(h_barrel_Cut20,"pT > 20 GeV");

    h_barrel_Cut25->SetStats(0);
    h_barrel_Cut25->GetXaxis()->SetTitle("m_{#gamma#gamma}GeV]");
    h_barrel_Cut25->GetYaxis()->SetTitle("events / GeV");
    h_barrel_Cut25->Draw("Same");
    legend_Cuts->AddEntry(h_barrel_Cut25,"pT > 25 GeV");

    h_barrel_Cut30->SetStats(0);
    h_barrel->GetXaxis()->SetTitle("m_{#gamma#gamma} [GeV]");
    h_barrel->GetYaxis()->SetTitle("events / GeV");
    h_barrel->Draw("Same");
    legend_Cuts->AddEntry(h_barrel_Cut30,"pT > 30 GeV");

    legend_Cuts->Draw();

    //canvas_h_Cuts->SaveAs("canvas_h_Cuts.pdf");

    // Créer le graphe
    TGraphErrors* graph_eff30 = new TGraphErrors();
    //graph_eff30->SetTitle("Efficiency evolution when pTsubl > 30");
    graph_eff30->GetXaxis()->SetTitle("m_{#gamma#gamma} [GeV]");  
    graph_eff30->GetYaxis()->SetTitle("Efficiency");
    graph_eff30->SetMarkerStyle(20);       
    graph_eff30->SetMarkerSize(0.2);

    TGraphErrors* graph_eff10 = new TGraphErrors();
    //graph_eff10->SetTitle("Efficiency evolution when pTsubl > 10");
    graph_eff10->SetTitle("Signal efficiency evolution");
    graph_eff10->GetXaxis()->SetTitle("m_{#gamma#gamma} [GeV]");  
    graph_eff10->GetYaxis()->SetTitle("Efficiency");
    graph_eff10->SetMarkerStyle(20);       
    graph_eff10->SetMarkerSize(0.2);

    TGraphErrors* graph_eff15 = new TGraphErrors();
    //graph_eff15->SetTitle("Efficiency evolution when pTsubl > 15");
    graph_eff15->GetXaxis()->SetTitle("m_{#gamma#gamma} [GeV]");  
    graph_eff15->GetYaxis()->SetTitle("Efficiency");
    graph_eff15->SetMarkerStyle(20);       
    graph_eff15->SetMarkerSize(0.2);

    TGraphErrors* graph_eff20 = new TGraphErrors();
    //graph_eff20->SetTitle("Efficiency evolution when pTsubl > 20");
    graph_eff20->GetXaxis()->SetTitle("m_{#gamma#gamma} [GeV]");  
    graph_eff20->GetYaxis()->SetTitle("Efficiency");
    graph_eff20->SetMarkerStyle(20);       
    graph_eff20->SetMarkerSize(0.2);

    TGraphErrors* graph_eff25 = new TGraphErrors();
    //graph_eff25->SetTitle("Efficiency evolution when pTsubl > 25");
    graph_eff25->GetXaxis()->SetTitle("m_{#gamma#gamma} [GeV]");  
    graph_eff25->GetYaxis()->SetTitle("Efficiency");
    graph_eff25->SetMarkerStyle(20);       
    graph_eff25->SetMarkerSize(0.2);

    // Afficher le contenu de chaque bin
    for (int bin = 1; bin <= h_barrel->GetNbinsX(); ++bin) { 

        double binCenter = h_barrel->GetBinCenter(bin);     
        double binContent = h_barrel->GetBinContent(bin);   
        double binError = h_barrel->GetBinError(bin);  

        double binContentCut10 = h_barrel_Cut10->GetBinContent(bin);   
        double binErrorCut10 = h_barrel_Cut10->GetBinError(bin);  
        double binContentCut15 = h_barrel_Cut15->GetBinContent(bin);   
        double binErrorCut15 = h_barrel_Cut15->GetBinError(bin); 
        double binContentCut20 = h_barrel_Cut20->GetBinContent(bin);   
        double binErrorCut20 = h_barrel_Cut20->GetBinError(bin); 
        double binContentCut25 = h_barrel_Cut25->GetBinContent(bin);   
        double binErrorCut25 = h_barrel_Cut25->GetBinError(bin); 
        double binContentCut30 = h_barrel_Cut30->GetBinContent(bin);   
        double binErrorCut30 = h_barrel_Cut30->GetBinError(bin); 

        double efficacite10 = binContentCut10/binContent;
        double incertitude10 = efficacite10* (1- efficacite10) / binContent;
        double efficacite15 = binContentCut15/binContent;
        double incertitude15 = efficacite15* (1- efficacite15) / binContent;
        double efficacite20 = binContentCut20/binContent;
        double incertitude20 = efficacite20* (1- efficacite20) / binContent;
        double efficacite25 = binContentCut25/binContent;
        double incertitude25 = efficacite25* (1- efficacite25) / binContent;
        double efficacite30 = binContentCut30/binContent;
        double incertitude30 = efficacite30* (1- efficacite30) / binContent;

        graph_eff10->SetPoint(bin-1, binCenter, efficacite10);
        graph_eff10->SetPointError(bin-1, 2.5*1e-3, incertitude10);
        graph_eff15->SetPoint(bin-1, binCenter, efficacite15);
        graph_eff15->SetPointError(bin-1, 2.5*1e-3, incertitude15);
        graph_eff20->SetPoint(bin-1, binCenter, efficacite20);
        graph_eff20->SetPointError(bin-1, 2.5*1e-3, incertitude20);
        graph_eff25->SetPoint(bin-1, binCenter, efficacite25);
        graph_eff25->SetPointError(bin-1, 2.5*1e-3, incertitude25);
        graph_eff30->SetPoint(bin-1, binCenter, efficacite30);
        graph_eff30->SetPointError(bin-1, 2.5*1e-3, incertitude30);


        
    }

    TCanvas *canvas_eff = new TCanvas("canvas_eff", "canvase_eff", 400, 400);   

    canvas_eff->SetTitle("Signal efficiency evolution");

    graph_eff10->GetXaxis()->SetRangeUser(8, 45);

    graph_eff10->GetYaxis()->SetRangeUser(0, 1.5);
    graph_eff10->SetMarkerColor(51);
    graph_eff10->SetLineColor(51);
    graph_eff10->Draw("AP");

    graph_eff15->GetYaxis()->SetRangeUser(0, 1.5);
    graph_eff15->SetMarkerColor(51+6);
    graph_eff15->SetLineColor(51+6);
    graph_eff15->Draw("SameP");

    graph_eff20->GetYaxis()->SetRangeUser(0, 1.5);
    graph_eff20->SetMarkerColor(51+12);
    graph_eff20->SetLineColor(51+12);
    graph_eff20->Draw("SameP");

    graph_eff25->GetYaxis()->SetRangeUser(0, 1.5);
    graph_eff25->SetMarkerColor(51+18);
    graph_eff25->SetLineColor(51+18);
    graph_eff25->Draw("SameP");

    graph_eff30->GetYaxis()->SetRangeUser(0, 1.5);
    graph_eff30->SetMarkerColor(51+24);
    graph_eff30->SetLineColor(51+24);
    graph_eff30->Draw("SameP");

    TLegend *legend3 = new TLegend(0.6, 0.7, 0.9, 0.9);  
    legend3->SetHeader("pT Cut","C");
    legend3->AddEntry(graph_eff10, "pT>10", "p"); 
    legend3->AddEntry(graph_eff15, "pT>15", "p");
    legend3->AddEntry(graph_eff20, "pT>20", "p");
    legend3->AddEntry(graph_eff25, "pT>25", "p");
    legend3->AddEntry(graph_eff30, "pT>30", "p");

    legend3->Draw();

    canvas_eff -> SaveAs("efficacité_barrel_signal.pdf");


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
    h_mass->Draw();
    h_mass_barrel->Draw("Same");
    h_mass_endcap->Draw("Same");


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

    //canvas_h_mass -> SaveAs("Resolutions_mass_Cut.pdf");

    //canvas_h -> SaveAs("Mass_Distributions_Cut.pdf");

    //graph_eff20->SetName("eff_20");

    //auto *outputFile = TFile::Open("eff_signal.root", "RECREATE");
    //graph_eff20->Write();
    //outputFile->Close();


    myFile->Close();


}
