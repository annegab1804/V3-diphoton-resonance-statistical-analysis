#include <TMath.h>
#include <iostream>

void find_sigma() {

    double prob = 0.885;  
    double p  = (prob+1)/2;

    // Calcul de la valeur quantile (inverse de la CDF)
    double quantile = TMath::NormQuantile(p);

    std::cout << "La valeur de k pour une probabilité de " << prob * 100 << "% est : " << quantile << std::endl;
}
