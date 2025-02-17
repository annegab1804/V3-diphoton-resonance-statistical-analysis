#include <TGraph.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <iostream>

void resolution_relative(){

    const int n = 300;  // Nombre de points
    double x[n], y[n];

    double a = 0.0108157 ;  // Pente
    double b = 0.157893 ;  // Ordonnée à l'origine
    
    int index = 0;

    // Remplir les tableaux
    for (int i = 5000; i < 35000; i+=100)
    {
        x[index] = i*1e-3;  // Éviter x=0 pour éviter la division par zéro
        double f_x = a * x[index] + b;
        y[index] = f_x / x[index];  // Calcul du rapport f(x)/x
        index ++;
    }

    // Création du graphe
    TGraph *graph = new TGraph(index, x, y);
    TCanvas *c1 = new TCanvas("c1", "Ratio f(x)/x vs x", 800, 600);

    graph->SetTitle("#sigma_{R}/m ; m (GeV); #sigma_{R}/m");
    graph->SetMarkerStyle(20);
    graph->SetMarkerSize(0.8);
    graph->Draw("APL");  // A: axes, P: points, L: lignes

    c1->Draw();
    c1->SaveAs("resolution_relative.pdf");
}