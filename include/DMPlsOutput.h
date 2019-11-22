// point2d.h
#ifndef DMPlsOutput_H
#define DMPlsOutput_H


#include <iostream>
#include <fstream>
#include <vector>
#include <TMath.h>

#include <TGraph.h>
#include <TCanvas.h>
#include <TH2F.h>
#include <TH1F.h>
#include <TF1.h>
#include <TF2.h>
#include <TTree.h>
#include <TMath.h>
#include <TFile.h>
#include <TTree.h>

using namespace std;

TH2F * hcutarea;
TH2F * hradius;
TH2F * hxmap;
TH2F * hymap;
TH2F * hrphimap;
TH2F * hphimap;
TH1F * hphicell;
TH1F * hrdist;
TH1F * hrdist2;
TH1F * hrdist_ldust;
//TH1F * hmtrk;
TH1F * hclx;
TH1F * hcly;

TH1F * hrclx[7];
TH1F * hrcly[7];
TH1F * hrclx_ev[7];
TH1F * hrcly_ev[7];
TH1F * hzrclx[7];
TH1F * hzrcly[7];

ofstream mybfcl;   /// lista bfcl per ogni grano
ofstream bfcl8;    /// lista bfcl (solo 8pol) per funzioni immagini animate
ofstream yandex;   /// lista bfcl (solo 8pol) per funzioni immagini animate
ofstream cut8;     /// lista bfcl (solo 8pol) per funzioni immagini animate con tagli
ofstream log_col;  /// lista single grain per color matching
ofstream sig;     /// lista id grani segnale
ofstream bkg;     /// lista id grani fondo

TGraph * grNcl;
TCanvas *c1;
TFile * f_out;
TTree * Tree_out;

/// TREE VARIABLES
Bool_t eLargeDust;
Bool_t eGoodZone;
Bool_t eShadow;
Double_t eClustTypeDust;
Double_t eCleanPar;
Int_t eChannel;
Int_t eGrainID;
Int_t ePolID;
Int_t eHeaderID;
Int_t eViewID;
Int_t eBfcPolID;
Int_t eBfcGap;
Int_t eNcl;
Int_t eNgr;
Int_t eNclFr;
Int_t ePuls;
Int_t eBfcID;
Int_t eNumIm;
Int_t eImCols;
Int_t eImRows;
Int_t eFlag;
Double_t eImArea;
Double_t eGrainPol;
Double_t eGrainx;
Double_t eGrainy;
Double_t ePreGrainx;
Double_t ePreGrainy;
Double_t eGrainz;
Double_t eClx;
Double_t eCly;
Double_t eClz;
Double_t eGrainrRms;
Double_t eGrainxRms;
Double_t eGrainyRms;
Double_t eGrainzRms;
Double_t eGrainMinRms;
Double_t eGrainMajRms;
Double_t eGrainPhiRms;
Double_t eGrainMin;
Double_t eGrainMaj;
Double_t eGrainEll;
Double_t eGrainPhi;
Double_t eGrainTheta;
Double_t eEllPrjX;
Double_t eEllPrjY;
Double_t eClustx;
Double_t eClusty;
Double_t eMClustx;
Double_t eMClusty;
Double_t eClustz;
Double_t eClustMin;
Double_t eClustMaj;
Double_t eClustEll;
Double_t eClustPhi;
Double_t eClustMaxPeak;
Double_t eClustMeanBkg;
Double_t eVolume;
Double_t eBfcVolume;
Double_t eFrBfVolume;
Double_t eArea;
Double_t eAreaRms;
Double_t eGrainArea;
Double_t eBfcArea;
Double_t eFrBfArea;
Double_t eBfcBorderFrame;
Double_t eXView;
Double_t eYView;
Double_t eZ;
Double_t eZlen;
Int_t eChain;
Int_t eSameFrame;
Int_t eNFrame;
Int_t eFrame;
Int_t eSetFrame;
Int_t eDeltaFrame;
Double_t eClDist;
Int_t eSetID;
Int_t eSetNCopy;
Int_t eSetNStatic;
Double_t eSetGapZ;
Double_t eSetGapZ2;
Double_t eSetGap;
Double_t eSetPath;
Double_t eSetMeanPath;
Double_t eSetRmsPath;
Double_t eSetMaxPath;
Double_t eSetMaxDist;
Double_t eSetNpeaks;
Double_t eSetNpeaksMax;
Double_t eSetNpeaksNcopy;
Double_t eSetNpeaksDist;
Double_t eSetNpeaksDistRms;
Double_t eSetNpeaksPhi;
Double_t eSetNpeaksDVol;
Double_t eSetNpeaksDNpx;
Double_t eSetNpeaksDBri;
Double_t eSetNpeaksMaxPhiAmp;
Double_t eSetNpeaksMeanPhi;
Double_t eSetMaxPol1;
Double_t eSetMaxPol2;
Double_t eSetMaxPixDist;
Double_t eSetMaxBar;
Double_t eSetBrAmp;
Double_t eSetBkgAmp;
Double_t eSetPeakAmp;
Double_t eSetMaxAmp;
Double_t eSetRmsAmp;
Double_t eSetMeanAmp;
Double_t eSetBrMaxPol;
Double_t eSetBrMinPol;
Double_t eSetXRms;
Double_t eSetYRms;
Double_t eSetXBar;
Double_t eSetYBar;
Double_t eSetXMaxBar;
Double_t eSetYMaxBar;
Double_t eSetXMinBar;
Double_t eSetYMinBar;
Double_t eSetPhiRms;
Double_t eSetPhiMean;
Double_t eSetChi2;
Double_t eMinDistGrain;
Double_t eSetPhi;
Double_t eSetPhiBar;
Double_t eSetPhiMaxAmp;
Double_t eSetTheBar;
Double_t eSetPhiFit;
Double_t eSetMeanBright;
Double_t eIsolated;
Double_t eSetNpxRms;
Double_t eSetVolRms;
Double_t eSetVolRatio;
Double_t eBfcMeanBkg;
Double_t eBfcSigPeak;
Double_t eBfcWeigth;
Double_t eSclMeanBkg;
Double_t eSclSigPeak;
Double_t eGrMeanBkg;
Double_t eGrSigPeak;
Double_t eGrSclMeanBkg;
Double_t eGrSclSigPeak;
Double_t eGrWeigth;
Double_t eDeltaPhi;
Double_t eDeltaPol;
/////// Microtracks /////////
Int_t eMTrk;
Int_t eMTID;
Int_t eMTGr;
Int_t eMTNfr;
Double_t eMTPhi;
Double_t eMTThe;
Double_t eMTLen;
Double_t eMTChi2;
Double_t eMTBrDif;
Double_t eMTNpxDif;
Double_t eMTZDif;
////////////////////////////////
Double_t eClFitMin;
Double_t eClFitMaj;
Double_t eClFitEll;
Double_t eClFitPhi;
Double_t eClFitx;
Double_t eClFity;
Double_t eSetFitBar;
Double_t eSetFitPhi;
Double_t eSetFitx;
Double_t eSetFity;     
Double_t eSetFitPhiMean;
/////////////////////////

/*namespace DMPlsAnalyzer {
	class Plasmon;
}*/
class DMPlsOutput
{
 public:

  std::tuple <TH2F*,TH1F*,TH1F*,TH1F*,TH1F*,TH1F*> createHistos (double len_view_x, double len_view_y, double totXbin, double totYbin);
  std::tuple <ofstream,ofstream,ofstream,ofstream,ofstream,ofstream,ofstream> createLogs ();
  std::tuple <TGraph*> createGraphs ();
  std::tuple <TCanvas*> createCanvas ();
  std::tuple <TFile*> createFile ();
  std::tuple <TTree*> createTree ();

    
 private:
  double x;
  double y;
};
#endif // DMPlsOutput_H
