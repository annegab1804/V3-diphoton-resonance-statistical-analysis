
bkg_PDF.cc : crée la PDF du background sans bruit de fond irréductible pour le cut à 20 (il faut le document h_m_yy_JZ/h_m_yy_JZ_barrel_20.root )

resolution_mass.cc : crée la distribution et la résolution en masse du signal et comparant les filtres « uu-barrel » et « !(uu)-endcaps ». 

Resolution_mass_CB_uu_barrel.cc : sert à trouver les fonctions décrivants les paramètres de la DSCB Besoin du fichier "diphoton_ntuple_collimated_yy.root ».

resolution_mass_CB_endcaps.cc : pareil mais pour les « !(uu)-endcaps »

h_m_yy_JZ_fit.cc : trace l’évolution de q, p-value, Ns, sigma Ns pour différentes masses de résonances de 6 à 34GeV tous les 1MeV. C’est les pdf des uu-barrel

h_m_yy_JZ_fit1.cc : trace le résultat du fit bkg + signal a une masse donnée

Nevents.cc : calcule le nombre d’évènements du background réductible et irréductible du run 2 des documents dans Run2data.root

resolution_mass_eff.cc : trace l’efficacité de mes différents cuts sur la valeur du pt du subleading pour les JZ. Utilise "diphoton_ntuple_collimated_yy.root"

generateAsimov.cc: trace l’incetitude sur Ns pour l’asimov dataset. 

bkg_tot.cc : Utilise "h_m_yy_JZ/h_m_yy_JZ_barrel_20.root" et "run2data.root" pour créer le dossier "bkg_tot.root" qui contient l’histogramme du bkg avec en plus la partie irréductible

bkg_PDF_tot.cc : crée la PDF du bruit de fond mais cette fois ci avec le dossier "bkg_tot.root"  pour prendre en compte la partie irréductible 

DSCB.cc: trace la représentation générale d’une double sided crystal ball

resolution_relative.cc : trace sigmaR/m 

Run2data.root: dossier root avec trois histogrammes: bruit de fond réductible du run 2, bruit de fond irréductible du run 2, ratio des deux 

sigmaNs.root : contient l’évolution de l’incertitude de Ns pour l’asimov 




Dans le dossier fluctuations: 

h_m_yy_JZ_q_20.cc : crée deux fichiers root "fluctuations.root" et « p-value.root » qui contiennent 20 graphs de q et de la p-value en fonction de la masse de résonance pour 20 datasets différents 
Pareil pour les fichier a,b,c,d
Quand il y a écrit 50 c’est qu’il en crée 50 et 100 c’est qu’il en crée 100

LEE.cc : crée la distribution de qmax pour mes 1000 datasets du bruit de fond. Utilise les dossiers run1 à run 6 et run50 et run51. Sauvegarde dans le fichier root "qmax.root"

Ntrials.cc: trouve le nombre de régions indépendantes si on suppose que c’est la somme de N gaussiennes. 

N.cc: trouve le nombre de régions indépendantes si on suppose que c’est le nombre de signal qui ne se recouvre pas 

find_sigma.cc: trouve à quoi correspond une fluctuation qui apparait 11.5 pourcent des cas. 

Qmax.root: distribution de qmax 

Les dossiers run 1 à 6 et run 50 et 51 contiennent les résultats de mes 1000 datasets 
