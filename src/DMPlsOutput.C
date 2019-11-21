#include "../include/DMPlsOutput.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
using namespace std;


std::tuple <TH2F*,TH1F*,TH1F*,TH1F*,TH1F*,TH1F*> DMPlsOutput::createHistos (double len_view_x, double len_view_y, double totXbin, double totYbin)
{

  const double pi=TMath::Pi();
  
 /// HISTOGRAMS
  hcutarea = new TH2F("cut_area","",len_view_x*100,-len_view_x/2.,len_view_x/2.,len_view_y*100,-len_view_y/2.,len_view_y/2.);  // mappa area_tagliata
  hrdist = new TH1F("rdist","",256,0,20);  // distanza mutua tra grani dopo i cuts primari (ldust, nearby ldust, no_ell_fit, minor_cut, long_chains)
  hrdist2 = new TH1F("rdist2","",128,0,1); // distanza tra i picchi dei grani npeaks
  hrdist_ldust = new TH1F("rdist_ldust","",10000,0,500); // distanza di ogni grano no large dust dai grani di large dust (espressa in unità di raggio del large dust)
  //hmtrk = new TH1F("dist_mtrk","",10000,0,500); // distanza minima tra grani appartenenti a due microtracce differenti
  
  hclx = new TH1F("clx","",800,-10,10);  // distanza mutua tra cluster di una collezione da otto [nm]
  hcly = new TH1F("cly","",800,-10,10);  // distanza mutua tra cluster di una collezione da otto [nm]


  return std::make_tuple(hcutarea,hrdist,hrdist2,hrdist_ldust,hclx,hcly);
}

std::tuple <ofstream,ofstream,ofstream,ofstream,ofstream,ofstream,ofstream> DMPlsOutput::createLogs ()
{
  mybfcl.open("pred_bfcl.txt");    /// lista bfcl per ogni grano
  bfcl8.open("bfcl8copy.txt");     /// lista bfcl (solo 8pol) per funzioni immagini animate
  yandex.open("yandex_bfcl.txt");     /// lista bfcl (solo 8pol) per funzioni immagini animate
  cut8.open("bfcl8_with_cuts.txt");     /// lista bfcl (solo 8pol) per funzioni immagini animate con tagli
  log_col.open("grains_with_cuts.txt");     /// lista single grains per color matching
  sig.open("sig_grains.txt");          /// lista id grani di segnale
  bkg.open("bkg_grains.txt");          /// lista id grani di fondo

  return std::make_tuple(move(mybfcl),move(bfcl8),move(yandex),move(cut8),move(log_col),move(sig),move(bkg));
}

std::tuple <TGraph*> DMPlsOutput::createGraphs ()
{
 grNcl = new TGraph(); // ncl per view
 grNcl->SetName("ncl_per_view");

 return std::make_tuple(grNcl);
}

std::tuple <TCanvas*> DMPlsOutput::createCanvas ()
{
  c1 = new TCanvas();
  return std::make_tuple(c1);
}

std::tuple <TFile*> DMPlsOutput::createFile ()
{
  f_out = new TFile("debug8_test_grain.root","RECREATE");
  return std::make_tuple(f_out);
}

std::tuple <TTree*> DMPlsOutput::createTree ()
{
 
  Tree_out = new TTree("tree1","data");
  
  Tree_out->Branch("eHeaderID",&eHeaderID,"eHeaderID/I");
  Tree_out->Branch("eViewID",&eViewID,"eViewID/I");
  Tree_out->Branch("eFlag",&eFlag,"eFlag/I");
  Tree_out->Branch("eGrainID",&eGrainID,"eGrainID/I");
  Tree_out->Branch("ePolID",&ePolID,"ePolID/I");
  Tree_out->Branch("eBfcID",&eBfcID,"eBfcID/I");
  Tree_out->Branch("eBfcPolID",&eBfcPolID,"eBfcPolID/I");
  Tree_out->Branch("eBfcGap",&eBfcGap,"eBfcGap/I");
  Tree_out->Branch("eNcl",&eNcl,"eNcl/I");
  Tree_out->Branch("eNgr",&eNgr,"eNgr/I");
  Tree_out->Branch("eNclFr",&eNclFr,"eNclFr/I");
  Tree_out->Branch("eGrainx",&eGrainx,"eGrainx/D");
  Tree_out->Branch("eGrainy",&eGrainy,"eGrainy/D");
  Tree_out->Branch("eGrainz",&eGrainz,"eGrainz/D");
  Tree_out->Branch("eClustx",&eClustx,"eClustx/D");
  Tree_out->Branch("eClusty",&eClusty,"eClusty/D");
  Tree_out->Branch("eClustz",&eClustz,"eClustz/D");
  Tree_out->Branch("eClustMin",&eClustMin,"eClustMin/D");
  Tree_out->Branch("eClustMaj",&eClustMaj,"eClustMaj/D");
  Tree_out->Branch("eClustEll",&eClustEll,"eClustEll/D");
  Tree_out->Branch("eClustPhi",&eClustPhi,"eClustPhi/D");
  Tree_out->Branch("ePuls",&ePuls,"ePuls/I");
  Tree_out->Branch("eGrainMin",&eGrainMin,"eGrainMin/D");
  Tree_out->Branch("eGrainMaj",&eGrainMaj,"eGrainMaj/D");
  Tree_out->Branch("eGrainEll",&eGrainEll,"eGrainEll/D");
  Tree_out->Branch("eGrainPhi",&eGrainPhi,"eGrainPhi/D");
  Tree_out->Branch("eGrainTheta",&eGrainTheta,"eGrainTheta/D");
  Tree_out->Branch("eEllPrjX",&eEllPrjX,"eEllPrjX/D");
  Tree_out->Branch("eEllPrjY",&eEllPrjY,"eEllPrjY/D");
  Tree_out->Branch("eGoodZone",&eGoodZone,"eGoodZone/B");
  Tree_out->Branch("eCleanPar",&eCleanPar,"eCleanPar/D");
  Tree_out->Branch("eLargeDust",&eLargeDust,"eLargeDust/B");
  Tree_out->Branch("eClustTypeDust",&eClustTypeDust,"eClustTypeDust/D");
  Tree_out->Branch("eVolume",&eVolume,"eVolume/D");
  Tree_out->Branch("eBfcVolume",&eBfcVolume,"eBfcVolume/D");
  Tree_out->Branch("eFrBfVolume",&eFrBfVolume,"eFrBfVolume/D");
  Tree_out->Branch("eArea",&eArea,"eArea/D");
  Tree_out->Branch("eBfcArea",&eBfcArea,"eBfcArea/D");
  Tree_out->Branch("eFrBfArea",&eFrBfArea,"eFrBfArea/D");
  Tree_out->Branch("eGrainArea",&eGrainArea,"eGrainArea/D");
  Tree_out->Branch("eXView",&eXView,"eXView/D");
  Tree_out->Branch("eYView",&eYView,"eYView/D");
  Tree_out->Branch("eZ",&eZ,"eZ/D");
  Tree_out->Branch("eZlen",&eZlen,"eZlen/D");
  Tree_out->Branch("eNFrame",&eNFrame,"eNFrame/I");
  Tree_out->Branch("eSetFrame",&eSetFrame,"eSetFrame/I");
  Tree_out->Branch("eSetNCopy",&eSetNCopy,"eSetNCopy/I");
  Tree_out->Branch("eSetNStatic",&eSetNStatic,"eSetNStatic/I");
  Tree_out->Branch("eSetXRms",&eSetXRms,"eSetXRms/D");
  Tree_out->Branch("eSetYRms",&eSetYRms,"eSetYRms/D");
  Tree_out->Branch("eSetXBar",&eSetXBar,"eSetXBar/D");
  Tree_out->Branch("eSetYBar",&eSetYBar,"eSetYBar/D");
  Tree_out->Branch("eSetXMaxBar",&eSetXMaxBar,"eSetXMaxBar/D");
  Tree_out->Branch("eSetYMaxBar",&eSetYMaxBar,"eSetYMaxBar/D");
  Tree_out->Branch("eSetXMinBar",&eSetXMinBar,"eSetXMinBar/D");
  Tree_out->Branch("eSetYMinBar",&eSetYMinBar,"eSetYMinBar/D");
  Tree_out->Branch("eSetNpeaks",&eSetNpeaks,"eSetNpeaks/D");
  Tree_out->Branch("eSetNpeaksMax",&eSetNpeaksMax,"eSetNpeaksMax/D");
  Tree_out->Branch("eSetNpeaksNcopy",&eSetNpeaksNcopy,"eSetNpeaksNcopy/D");
  Tree_out->Branch("eSetNpeaksDist",&eSetNpeaksDist,"eSetNpeaksDist/D");
  Tree_out->Branch("eSetNpeaksDistRms",&eSetNpeaksDistRms,"eSetNpeaksDistRms/D");
  Tree_out->Branch("eSetNpeaksPhi",&eSetNpeaksPhi,"eSetNpeaksPhi/D");
  Tree_out->Branch("eSetNpeaksDVol",&eSetNpeaksDVol,"eSetNpeaksDVol/D");
  Tree_out->Branch("eSetNpeaksDNpx",&eSetNpeaksDNpx,"eSetNpeaksDNpx/D");
  Tree_out->Branch("eSetNpeaksDBri",&eSetNpeaksDBri,"eSetNpeaksDBri/D");
  Tree_out->Branch("eSetNpeaksMaxPhiAmp",&eSetNpeaksMaxPhiAmp,"eSetNpeaksMaxPhiAmp/D");
  Tree_out->Branch("eSetNpeaksMeanPhi",&eSetNpeaksMeanPhi,"eSetNpeaksMeanPhi/D");
  Tree_out->Branch("eSetGapZ",&eSetGapZ,"eSetGapZ/D");
  Tree_out->Branch("eSetGapZ2",&eSetGapZ2,"eSetGapZ2/D");
  Tree_out->Branch("eSetGap",&eSetGap,"eSetGap/D");
  Tree_out->Branch("eSetPath",&eSetPath,"eSetPath/D");
  Tree_out->Branch("eSetMeanPath",&eSetMeanPath,"eSetMeanPath/D");
  Tree_out->Branch("eSetRmsPath",&eSetRmsPath,"eSetRmsPath/D");
  Tree_out->Branch("eSetMaxPath",&eSetMaxPath,"eSetMaxPath/D");
  Tree_out->Branch("eSetMaxDist",&eSetMaxDist,"eSetMaxDist/D");
  Tree_out->Branch("eSetMaxPol1",&eSetMaxPol1,"eSetMaxPol1/D");
  Tree_out->Branch("eSetMaxPol2",&eSetMaxPol2,"eSetMaxPol2/D");
  Tree_out->Branch("eSetMaxBar",&eSetMaxBar,"eSetMaxBar/D");
  Tree_out->Branch("eSetMaxAmp",&eSetMaxAmp,"eSetMaxAmp/D");
  Tree_out->Branch("eSetRmsAmp",&eSetRmsAmp,"eSetRmsAmp/D");
  Tree_out->Branch("eSetMeanAmp",&eSetMeanAmp,"eSetMeanAmp/D");
  Tree_out->Branch("eSetPhi",&eSetPhi,"eSetPhi/D");
  Tree_out->Branch("eSetPhiRms",&eSetPhiRms,"eSetPhiRms/D");
  Tree_out->Branch("eSetPhiMean",&eSetPhiMean,"eSetPhiMean/D");
  Tree_out->Branch("eSetPhiBar",&eSetPhiBar,"eSetPhiBar/D");
  Tree_out->Branch("eSetPhiMaxAmp",&eSetPhiMaxAmp,"eSetPhiMaxAmp/D");
  Tree_out->Branch("eSetTheBar",&eSetTheBar,"eSetTheBar/D");
  Tree_out->Branch("eSetMeanBright",&eSetMeanBright,"eSetMeanBright/D");
  Tree_out->Branch("eSetBrAmp",&eSetBrAmp,"eSetBrAmp/D");
  Tree_out->Branch("eSetBrMaxPol",&eSetBrMaxPol,"eSetBrMaxPol/D");
  Tree_out->Branch("eSetBrMinPol",&eSetBrMinPol,"eSetBrMinPol/D");
  Tree_out->Branch("eSetBkgAmp",&eSetBkgAmp,"eSetBkgAmp/D");
  Tree_out->Branch("eSetPeakAmp",&eSetPeakAmp,"eSetPeakAmp/D");
  Tree_out->Branch("eIsolated",&eIsolated,"eIsolated/D");
  Tree_out->Branch("eSetNpxRms",&eSetNpxRms,"eSetNpxRms/D");
  Tree_out->Branch("eSetVolRms",&eSetVolRms,"eSetVolRms/D");
  Tree_out->Branch("eSetVolRatio",&eSetVolRatio,"eSetVolRatio/D");
  Tree_out->Branch("eMTrk",&eMTrk,"eMTrk/I");
  Tree_out->Branch("eMTID",&eMTID,"eMTID/I");
  Tree_out->Branch("eMTGr",&eMTGr,"eMTGr/I");
  Tree_out->Branch("eMTNfr",&eMTNfr,"eMTNfr/I");
  Tree_out->Branch("eMTPhi",&eMTPhi,"eMTPhi/D");
  Tree_out->Branch("eMTThe",&eMTThe,"eMThe/D"); 
  Tree_out->Branch("eMTLen",&eMTLen,"eMTLen/D");
  Tree_out->Branch("eMTBrDif",&eMTBrDif,"eMTBrDif/D");
  Tree_out->Branch("eMTNpxDif",&eMTNpxDif,"eMTNpxDif/D");
  Tree_out->Branch("eMTZDif",&eMTZDif,"eMTZDif/D");
  Tree_out->Branch("eChannel",&eChannel ,"eChannel/I");
  Tree_out->Branch("eDeltaPhi",&eDeltaPhi ,"eDeltaPhi/D");
  Tree_out->Branch("eDeltaPol",&eDeltaPol ,"eDeltaPol/D");

  return std::make_tuple(Tree_out);
}

