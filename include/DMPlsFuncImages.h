// point2d.h
#ifndef DMPlsFuncImages_H
#define DMPlsFuncImages_H

#include <TVectorD.h>
#include <TGraph.h>
#include <iostream>
#include <vector>
#include <TMath.h>
#include <TRandom.h>
#include <TRandom3.h>
#include <TROOT.h>
#include <TEllipse.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TH2F.h>
#include <TF1.h>
#include <TF2.h>
#include <TStyle.h>
#include <TSystem.h>

#include "/home/vale/Dropbox/DMDS/include/DMRViewHeader.h"
#include "/home/vale/Dropbox/DMDS/include/DMRAffine2D.h"
#include "/home/vale/Dropbox/DMDS/include/DMRCluster.h"
#include "/home/vale/Dropbox/DMDS/include/DMRGrain.h"
#include "/home/vale/Dropbox/DMDS/include/DMRMicrotrack.h"
#include "/home/vale/Dropbox/DMDS/include/DMRImage.h"
#include "/home/vale/Dropbox/DMDS/include/DMRRun.h"
#include "/home/vale/Dropbox/DMDS/include/DMRRunHeader.h"


using namespace std;

/*namespace DMPlsAnalyzer {
	class Plasmon;
}*/
class DMPlsFuncImages
{
public:
TVectorD fit_ellipse(TGraph *g);
TVectorD ConicToParametric(const TVectorD &conic);
TGraph *TestGraphDLSF(Bool_t randomize);
void fitEllipseTGraphDLSF(TGraph *g, double *fit_info);
TH2F * merged_histo(DMRRun *aRun, DMRView *view, int ihd, int igr, int *clset);
void gr_spectrum(TH2F *h, double *gr_scale_result);
void scaling(DMRRun *aRun, int ihd, int icl, int ipol, double *scale_result);
float maxpeak(DMRRun *aRun, int ihd, int icl, int ipol);
void images_info(DMRRun *aRun,int ihd, int icl, int ipol, double *im_info);
int otsu_method(float *histogram, long int total_pixels);
void images_ellipse(DMRRun *aRun,DMRView* view, int ihd, int icl, int ipol, double *fit_info, TCanvas *c1);

private:
double x;
double y;
};
#endif // DMPlsAnalyzer_H
