#include <TFile.h>
#include <TGraph.h>
#include <iostream>
#include <vector>

void compter_pic() {
    int peakCount = 0;

    //std::vector<std::string> fileNames = {
    //    "fluctuations.root", "fluctuations_a.root", "fluctuations_b.root",
    //    "fluctuations_c.root", "fluctuations_d.root"
    //};

    std::vector<std::string> fileNames = { "fluctuations_b.root",
        "fluctuations_c.root", "fluctuations_d.root"
    };

    for (const auto& fileName : fileNames) {
        TFile *file = TFile::Open(fileName.c_str(), "READ");

        if (!file || file->IsZombie()) {
            std::cerr << "Erreur lors de l'ouverture du fichier ROOT " << fileName << std::endl;
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


                if (maxVal >  9) {
                    peakCount++;

                    std::cout << "Dans " << fileName << ", le TGraph " << graphName 
                              << " a une fluctuation à plus de 3 sigma, valeur maximale : " << maxVal << std::endl;

                    // Créer des vecteurs pour les nouvelles valeurs p = TMath::Prob(q, 1)
                    std::vector<double> probValues;

                    // Calculer p = TMath::Prob(q, 1) pour chaque valeur de q
                    for (int l = 0; l < nPoints; ++l) {
                        double p = TMath::Prob(yValues[l], 1); // Calculer la probabilité pour q
                        probValues.push_back(p);
                    }

                    TGraph *graphProb = new TGraph(nPoints, xValues, &probValues[0]);
                    
                    //Créer un TCanvas pour afficher le graphique
                    TCanvas *c1 = new TCanvas(Form("c1_%s_%d", fileName.c_str(), i+1), Form("Graph %s %d", fileName.c_str(), i+1), 800, 600);
                    c1->SetLogy();
                    graphProb->GetYaxis()->SetRangeUser(1e-7, 1);
                    graphProb->Draw("AP");

                    double sigmaLevels[5] = {0.3173, 0.0455, 0.0027, 0.000063, 2.87e-7};
                    int colors[5] = {kBlack, kRed, kGreen+2, kMagenta, kCyan+2};

                    const char* sigmaLabels[5] = {"1#sigma", "2#sigma", "3#sigma", "4#sigma", "5#sigma"};

                    TLine *lines[5];
                    TLatex *labels[5];

                    for (int k = 0; k < 5; k++) {
                        // Tracer la ligne horizontale
                        lines[k] = new TLine(4, sigmaLevels[k], 36, sigmaLevels[k]);
                        lines[k]->SetLineColor(colors[k]);
                        lines[k]->SetLineStyle(2); // Pointillé
                        lines[k]->Draw("Same");

                        // Ajouter un label juste à côté de l'axe des ordonnées
                        labels[k] = new TLatex(3 - 0.2, sigmaLevels[k], sigmaLabels[k]); 
                        labels[k]->SetTextSize(0.035);  // Taille du texte
                        labels[k]->SetTextColor(colors[k]); // Couleur du texte
                        labels[k]->SetTextAlign(32);  // Alignement à droite
                        labels[k]->Draw();
                    }

                    TString outputFileName = Form("graph_%s_%d.png", fileName.c_str(), i);
        
                    // Optionnel : Enregistrer chaque graphique dans un fichier PDF
                    c1->SaveAs(outputFileName);

                    // Nettoyage de la mémoire
                    for (auto line : lines) delete line;
                    for (auto label : labels) delete label;
                    delete c1;


                }

            } else {
                std::cout << "Erreur : TGraph " << graphName 
                          << " non trouvé dans " << fileName << std::endl;
            }
        }
        file->Close();
    }

    std::cout << "Il y a " << peakCount << " graphiques avec une fluctuation a plus de 3 sigmas." << std::endl;
}
