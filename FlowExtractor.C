#include <iostream>
#include "Riostream.h"
#include <cstdlib>
#include <string>
#include <vector>
#include <math.h>
#include <set>
#include <map>

//---
//need this stuff to compile in linux:
#include <TROOT.h>
#include <TApplication.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TStreamerElement.h>
#include <TStyle.h>
#include "TSystemDirectory.h"
//---

#include "TFile.h"
#include "TNtuple.h"
#include "TTree.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TMath.h"
#include "TChain.h"
#include "TLegend.h"
#include "TFractionFitter.h"
#include "TVirtualFitter.h"
#include "TCut.h"
#include "TObject.h"
#include "TLine.h"
#include "TLeaf.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TGaxis.h"
#include "TSpectrum.h"
#include "TMath.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TMultiGraph.h"
#include "TChain.h"
#include "TLegend.h"
#include "TFractionFitter.h"
#include "TVirtualFitter.h"
#include "TCut.h"
#include "TObject.h"
#include "TLine.h"
#include "TSpline.h"
#include "TPaveText.h"
#include "Fit/FitResult.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TRandom3.h"
#include "TDatime.h"
#include "Math/MinimizerOptions.h"

using namespace std;

void FlowExtractor()
{
  int FlowOrder = 2;
  // Get # of events from Kaon TTree file
  // TFile * tf_evt_in  = new TFile("/star/data01/pwg/dchen/Ana/7p2GeV_FXT_2018/KKinvM/kaonTree/merged_merged_D9AAA11F292FF0BA1047384EC836F91E.picoDst.result.root","READ");
  // TH1D * h_evt   = (TH1D*) tf_evt_in -> Get("h_evt");
  // double d_N_evt    = h_evt -> Integral();

  // normal/mixed event invM
  TFile * tf_evt_nm_in = new TFile("/star/data01/pwg/dchen/Ana/7p2GeV_FXT_2018/KKinvM/phiInvM/merged_KKinvMOutput.root","READ");
  TFile * tf_evt_mx_in = new TFile("/star/data01/pwg/dchen/Ana/7p2GeV_FXT_2018/KKinvM/phiInvM/merged_KKinvMOutputMixedEvents.root","READ");

  TH1D *mHistKKInvMpT[12];
  for(int pt=0; pt<12; pt++)
  {
    mHistKKInvMpT[pt] = (TH1D*) tf_evt_nm_in->Get(Form("histKKInvMpT%d",pt));
    mHistKKInvMpT[pt]   -> SetXTitle("M_{K+K-}");
    mHistKKInvMpT[pt]   -> SetYTitle("Counts");
    //mHistKKInvMpT[pt]->Rebin();
  }

  TH1D *mHistKKInvMpTMixed[12];
  for(int pt=0; pt<12; pt++)
  {
    mHistKKInvMpTMixed[pt] = (TH1D*) tf_evt_mx_in->Get(Form("histKKInvMpT%d_Mixed",pt));
    mHistKKInvMpTMixed[pt]   -> SetXTitle("M_{K+K-}");
    mHistKKInvMpTMixed[pt]   -> SetYTitle("Counts");
    //mHistKKInvMpTMixed[pt]->Rebin();
  }

  //--- Normalization range
  double a_d_int_range[4] =
    {
      0.99,
      1.008,
      1.04,
      1.09
    };

  // Get the bin of the Normalization range
  int a_iBin_range[4];
  for(int ijk = 0; ijk < 4; ijk++) a_iBin_range[ijk] =  mHistKKInvMpT[0] -> FindFixBin(a_d_int_range[ijk]);

  gStyle->SetOptDate(0);

  // TCanvas To draw nomal/mixed invariant mass and flow vs invariant mass
  TCanvas * TC_invM = new TCanvas("TC_invM","TC_invM",1920,1080);
  TC_invM->Divide(4,3);
  TC_invM->SaveAs(Form("PhiInvMassFlowV%d.pdf[",FlowOrder));

  // Draw same and normalized-mx events together
  for(int pt=0; pt<12; pt++)
  {
   TC_invM->cd(++pt);
   mHistKKInvMpT[pt]->SetTitle(Form("%d <= p_{T} <= %d (GeV/c)",pt*0.3,(pt+1.)*0.3));
   // mHistKKInvMpTMixed[pt]->SetTitle(Form("%d <= p_{T} <= %d (GeV/c)",pt*0.3,(pt+1.)*0.3));
   mHistKKInvMpT[pt]->GetXaxis()->SetTitleSize(0.08);
   mHistKKInvMpT[pt]->GetXaxis()->SetTitleOffset(1);
   mHistKKInvMpT[pt]->GetXaxis()->SetRangeUser(0.98,1.1);
   // mHistKKInvMpTMixed[pt]->GetXaxis()->SetTitleSize(0.08);
   // mHistKKInvMpTMixed[pt]->GetXaxis()->SetTitleOffset(1);
   // mHistKKInvMpTMixed[pt]->GetXaxis()->SetRangeUser(0.98,1.1);
   // mHistKKInvMpTMixed[pt]-> SetFillColor(kYellow);
   // mHistKKInvMpTMixed[pt]-> SetFillStyle(3001);
   mHistKKInvMpT[pt]->Draw();
   // mHistKKInvMpT[pt]->Draw("same");



  }

  TC_invM->SaveAs(Form("PhiInvMassFlowV%d.pdf]",FlowOrder));

  return;
}
