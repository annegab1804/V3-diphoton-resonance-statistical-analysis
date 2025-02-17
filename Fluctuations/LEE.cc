#include "TFile.h"
#include "TGraph.h"
#include "iostream"
#include "vector"
#include "TH1F.h"
#include "TMath.h"

void LEE() {

    std::vector<std::string> runs = {"run1", "run2", "run3", "run4", "run5", "run6"};
    std::vector<std::string> fileNames = {
        "fluctuations.root", "fluctuations_a.root", "fluctuations_b.root",
        "fluctuations_c.root", "fluctuations_d.root"
    };

    std::vector<std::string> fileNames51 = { "fluctuations_b.root",
        "fluctuations_c.root", "fluctuations_d.root"
    };

    TH1F *h_qmax = new TH1F("h_qmax", "Qmax distribution", 250, 0, 25);

    for (const auto& run : runs) {
        for (const auto& fileName : fileNames) {
            std::string filePath = run + "/" + fileName;
            TFile *file = TFile::Open(filePath.c_str(), "READ");

            if (!file || file->IsZombie()) {
                std::cerr << "Erreur lors de l'ouverture du fichier ROOT " << filePath << std::endl;
                continue;
            }

            for (int i = 1; i <= 20; i++) {
                TString graphName = Form("graphq_%d", i);
                TGraph *graph = (TGraph*)file->Get(graphName);

                if (graph) {
                    double* yValues = graph->GetY();
                    double* xValues = graph->GetX();
                    int nPoints = graph->GetN();
                    double maxVal = -1e9;

                    for (int j = 0; j < nPoints; ++j) {
                        if (yValues[j] > maxVal) {
                            maxVal = yValues[j];
                        }
                    }

                    h_qmax->Fill(maxVal);

                    
                } else {
                    std::cout << "Erreur : TGraph " << graphName 
                              << " non trouvé dans " << filePath << std::endl;
                }
            }
            file->Close();
        }
    }


    for (const auto& fileName : fileNames) {
        std::string filePath = "run50/" + fileName;
        TFile *file = TFile::Open(filePath.c_str(), "READ");
        if (!file || file->IsZombie()) {
            std::cerr << "Erreur lors de l'ouverture du fichier ROOT " << filePath << std::endl;
            continue;
        }
        for (int i = 1; i <= 50; i++) {
            TString graphName = Form("graphq_%d", i);
            TGraph *graph = (TGraph*)file->Get(graphName);
            if (graph) {
                double* yValues = graph->GetY();
                double* xValues = graph->GetX();
                int nPoints = graph->GetN();
                double maxVal = -1e9;
                for (int j = 0; j < nPoints; ++j) {
                    if (yValues[j] > maxVal) {
                        maxVal = yValues[j];
                    }
                }
                h_qmax->Fill(maxVal);
            } else {
                std::cout << "Erreur : TGraph " << graphName 
                          << " non trouvé dans " << filePath << std::endl;
            }
        }

        file->Close();

    }

    for (const auto& fileName : fileNames51) {
        std::string filePath = "run51/" + fileName;
        TFile *file = TFile::Open(filePath.c_str(), "READ");
        if (!file || file->IsZombie()) {
            std::cerr << "Erreur lors de l'ouverture du fichier ROOT " << filePath << std::endl;
            continue;
        }
        for (int i = 1; i <= 50; i++) {
            TString graphName = Form("graphq_%d", i);
            TGraph *graph = (TGraph*)file->Get(graphName);
            if (graph) {
                double* yValues = graph->GetY();
                double* xValues = graph->GetX();
                int nPoints = graph->GetN();
                double maxVal = -1e9;
                for (int j = 0; j < nPoints; ++j) {
                    if (yValues[j] > maxVal) {
                        maxVal = yValues[j];
                    }
                }
                h_qmax->Fill(maxVal);
            } else {
                std::cout << "Erreur : TGraph " << graphName 
                          << " non trouvé dans " << filePath << std::endl;
            }
        }

        file->Close();

    }

    double median;
    double prob[1] = {0.5};  // 0.5 = 50% = médiane
    h_qmax->GetQuantiles(1, &median, prob);
    double median_chi2 = 0.455;  



    std::cout << "La médiane de qmax est : " << median << std::endl;

    double N_trials = exp(median - 0.455);
    std::cout << "Nombre de trials estimé : " << N_trials << std::endl;


    TCanvas *c1 = new TCanvas("c1", "Distribution de Qmax", 800, 600);

    // Définition de la fonction Chi2 avec correction pour éviter division par 0
    TF1 *chi2 = new TF1("chi2", "[0] * (1.0 / sqrt(2.0 * TMath::Pi() * TMath::Max(x, 1e-6))) * exp(-x / 2.0)", 0, 25);

    // Normalisation pour que l'aire de la courbe soit égale à celle de l'histogramme
    double binWidth = h_qmax->GetBinWidth(1);  // Largeur des bins
    double histIntegral = h_qmax->Integral() * binWidth; // Aire totale

    chi2->SetParameter(0, histIntegral); // Ajuste la normalisation

    // Dessiner l'histogramme et superposer la courbe Chi2
    h_qmax->SetStats(0);
    h_qmax->GetXaxis()->SetTitle("q_{max}");
    h_qmax->GetYaxis()->SetTitle("Number of dataset");
    h_qmax->Draw();
    //chi2->SetLineColor(kRed);
    //chi2->Draw("same");


    // Trouver les bornes de l'axe Y pour que les lignes soient bien visibles
    double ymin = 0;
    double ymax = h_qmax->GetMaximum()*1.05;  // 10% au-dessus du max pour l'esthétique

    // Tracer la médiane de l'histogramme
    TLine *line_hist = new TLine(median, ymin, median, ymax);
    line_hist->SetLineColor(kBlue);   // Couleur bleue pour l'histogramme
    line_hist->SetLineStyle(2);       // Style pointillé
    line_hist->SetLineWidth(2);
    line_hist->Draw("same");

    // Tracer la médiane de la loi Chi2_1
    //TLine *line_chi2 = new TLine(median_chi2, ymin, median_chi2, ymax);
    //line_chi2->SetLineColor(kRed); // Couleur verte pour Chi2_1
    //line_chi2->SetLineStyle(2);        // Style pointillé
    //line_chi2->SetLineWidth(2);
    //line_chi2->Draw("same");

    // Ajouter une légende
    TLegend *leg = new TLegend(0.6, 0.7, 0.85, 0.85);
    leg->AddEntry(h_qmax, "q_{max} distribution", "l");
    //leg->AddEntry(chi2, "#chi^{2} distribution", "l");
    leg->AddEntry(line_hist,  "q_{max} median", "l");
    //leg->AddEntry(line_chi2, "#chi^{2} median", "l");
    leg->Draw();

    // Sauvegarde du graphe
    c1->SaveAs("distribution_qmax_chi2.png");

    auto *outputFile = TFile::Open("qmax.root", "RECREATE");
    h_qmax->Write();
    outputFile->Close();

    }

