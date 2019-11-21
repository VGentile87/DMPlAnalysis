#include "DMPlsMtAnalyzer.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

using namespace std;

double DMPlsMtAnalyzer::mtrk_phi(int mt_ngr, double *x_mt, double *y_mt)
{

  float phi_mt;

  if(mt_ngr==2) phi_mt=TMath::ATan((y_mt[1]-y_mt[0])/(x_mt[1]-x_mt[0]));
  
  if(mt_ngr>2){   /// fit lineare se il numero di grani è maggiore di 2
    TGraph *grmt = new TGraph(mt_ngr,x_mt,y_mt);
    grmt->Fit("pol1","Q");
    TF1 *pol1 = grmt->GetFunction("pol1");
    phi_mt = TMath::ATan(pol1->GetParameter(1));
    pol1->ReleaseParameter(1);
    delete grmt;
  }
  return phi_mt;
}






