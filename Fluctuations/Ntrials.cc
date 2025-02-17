#include "TH1F.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TRandom3.h"
#include "TLegend.h"
#include <iostream>
#include "TMath.h"

void Ntrials() {


    TRandom3 randGen(0);

    int n = 1;
    int n_sigma = 0;
    double prob =0;
    double max_val = -1e9;


    //TH1F *h_max = new TH1F("h_max", "Distribution du max de Gauss(0,1)", 100, 0, 25);

    while(prob <0.115){

        for (int j=0; j<1000; ++j){

            max_val = -1e9;

            for (int i=0 ; i< n ; ++i){

                double val = randGen.Gaus(0, 1);

                if (TMath::Abs(val)>max_val){
                    max_val = TMath::Abs(val);
                }

            }

            if (max_val >= 3){
                n_sigma++;
            }

            //h_max->Fill(max_val);


        }

        prob = static_cast<double>(n_sigma) / 1000.0;

        std::cout << "Proba de dépasser 3 sigmas pour n = "<< n << " : " << prob << std::endl;

        n+=1;
        n_sigma=0;

    }


}