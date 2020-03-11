#include <iostream>
#include <fstream>
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
  gStyle->SetOptStat(0);
  int FlowOrder = 2;
  // Get # of events from Kaon TTree file
  // TFile * tf_evt_in  = new TFile("/star/data01/pwg/dchen/Ana/7p2GeV_FXT_2018/KKinvM/kaonTree/merged_merged_D9AAA11F292FF0BA1047384EC836F91E.picoDst.result.root","READ");
  // TH1D * h_evt   = (TH1D*) tf_evt_in -> Get("h_evt");
  // double d_N_evt    = h_evt -> Integral();

  // normal/mixed event invM
  TFile * tf_evt_nm_in = new TFile("/mnt/c/Users/pjska/github/FlowExtractor/input/merged_KKinvMOutput.root","READ");
  TFile * tf_evt_mx_in = new TFile("/mnt/c/Users/pjska/github/FlowExtractor/input/merged_KKinvMOutputMixedEvents.root","READ");

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


  // gStyle->SetOptDate(0);

  // TCanvas To draw nomal/mixed invariant mass and flow vs invariant mass
  TCanvas * TC_invM = new TCanvas("TC_invM","TC_invM",800,600);
  TC_invM->Divide();
  TC_invM->SaveAs(Form("PhiInvMassFlowV%d.pdf[",FlowOrder));

  // Draw same and normalized-mx events together
  for(int pt=0; pt<12; pt++)
  {
    TPad* thePad = (TPad*)TC_invM->cd();

    thePad->cd();
    TH1D* sameEventInvM = (TH1D*) tf_evt_nm_in->Get(Form("histKKInvMpT%d",pt));
    sameEventInvM->SetTitle(Form("%1.1f <= p_{T} <= %1.1f (GeV/c)",pt*0.3,(pt+1.)*0.3));
    sameEventInvM   -> SetXTitle("M_{K+K-} (GeV/c^{2})");
    sameEventInvM   -> SetYTitle("Counts");
    sameEventInvM->GetXaxis()->SetRangeUser(0.98,1.1);
    TH1D* mixedEventInvM = (TH1D*) tf_evt_mx_in->Get(Form("histKKInvMpT%d_Mixed",pt));
    mixedEventInvM->SetTitle(Form("%1.1f <= p_{T} <= %1.1f (GeV/c)",pt*0.3,(pt+1.)*0.3));
    mixedEventInvM   -> SetXTitle("M_{K+K-} (GeV/c^{2})");
    mixedEventInvM   -> SetYTitle("Counts");
    mixedEventInvM->GetXaxis()->SetRangeUser(0.98,1.1);
    mixedEventInvM -> SetFillColor(kYellow);
    mixedEventInvM -> SetFillStyle(3001);

    // Get the bin of the Normalization range
    int a_iBin_range[4];
    for(int ijk = 0; ijk < 4; ijk++) a_iBin_range[ijk] =  sameEventInvM -> FindFixBin(a_d_int_range[ijk]);
    //right bg Normalization
    double d_r_area     = a_d_int_range[3]-a_d_int_range[2];
    double d_r_same_int = sameEventInvM -> Integral(a_iBin_range[2],a_iBin_range[3]);
    double d_r_mx_int   = mixedEventInvM -> Integral(a_iBin_range[2],a_iBin_range[3]);
    double d_r_norm     = (d_r_mx_int != 0.0) ? d_r_same_int/d_r_mx_int : 1.0;
    cout<<" R: "<<d_r_area<<" : "<<d_r_same_int<<" : "<<d_r_mx_int<<endl;
    //left bg Normalization
    double d_l_area =  a_d_int_range[1]-a_d_int_range[0];
    double d_l_same_int = sameEventInvM -> Integral(a_iBin_range[0],a_iBin_range[1]);
    double d_l_mx_int   = mixedEventInvM -> Integral(a_iBin_range[0],a_iBin_range[1]);
    double d_l_norm     = (d_l_mx_int != 0.0) ? d_l_same_int/d_l_mx_int : 1.0;
    cout<<" L: "<<d_l_area<<" : "<<d_l_same_int<<" : "<<d_l_mx_int<<endl;
    cout<<" l: "<<d_l_norm<<" r: "<<d_r_norm<<endl;
    double d_norm = ((d_r_norm/d_r_area) + (d_l_norm/d_l_area))/( (1.0/d_r_area) + (1.0/d_l_area));
    cout<<" d_norm = "<<d_norm<<endl;

    // Normalize mixed event invariant mass
    mixedEventInvM -> Sumw2();
    mixedEventInvM -> Scale(d_norm);

    sameEventInvM->Draw();
    mixedEventInvM->Draw("HISTsames");

    // normalized areas
    TLine * a_TL_int[4];
    for(int k = 0; k < 4; k++)
    {
      a_TL_int[k] = new TLine(a_d_int_range[k],0.0001,a_d_int_range[k],1.*sameEventInvM -> GetMaximum());
      a_TL_int[k] -> SetLineColor(kRed);
      a_TL_int[k] -> Draw();
    }
    TC_invM->Update();
    TC_invM->SaveAs(Form("PhiInvMassFlowV%d.pdf",FlowOrder));

    // Substract normalized mixedEventInvM from sameEventInvM to get Signal
    TH1D * HistSignal = (TH1D*) sameEventInvM -> Clone("HistSignal");
    HistSignal -> Reset();
    HistSignal -> Sumw2();

    for(int ijk = 1; ijk < HistSignal->GetNbinsX()+1; ijk++)
    {
      double d_center   = sameEventInvM -> GetBinCenter(ijk);
      double d_same     = sameEventInvM -> GetBinContent(ijk);
      double d_same_err = sameEventInvM -> GetBinError(ijk);
      double d_mx       = mixedEventInvM -> GetBinContent(ijk);
      double d_mx_err   = mixedEventInvM -> GetBinError(ijk);
      double d_sig      = d_same - d_mx;
      double d_sig_err  = sqrt(d_same_err*d_same_err+d_mx_err*d_mx_err);

      HistSignal -> SetBinContent(ijk,d_sig);
      HistSignal -> SetBinError(ijk,d_sig_err);
    }

    HistSignal -> Draw("E");
    gStyle->SetOptFit(1111);

    //Fit function
    //fit to Gauss plus constant
    TFormula * GausPlus = new TFormula("GausPlus","gaus(0)+[3]");
    TF1 * tf1_phi = new TF1("polygauss_single",GausPlus->GetExpFormula(),1.,1.06);

    //fit to a simple gauss first to get seed
    TF1 * tf1_gauss = new TF1("tf1_gauss","gaus",0.9,1.1);
    HistSignal -> Fit(tf1_gauss,"0","R",1.01,1.03);
    // seeds
    double d_seed_p0    = tf1_gauss -> GetParameter(0);
    double d_seed_mean  = tf1_gauss -> GetParameter(1);
    double d_seed_sigma = tf1_gauss -> GetParameter(2);

    tf1_phi -> SetParameter(1,d_seed_mean);
    tf1_phi -> SetParameter(2,d_seed_sigma);
    tf1_phi -> SetParLimits(1,d_seed_mean-d_seed_sigma,d_seed_mean+d_seed_sigma);
    tf1_phi -> SetParLimits(2,0.66*d_seed_sigma,1.5*d_seed_sigma);
    tf1_phi -> SetLineColor(kBlue);
    int FitStatus = HistSignal   -> Fit(tf1_phi,"E+","R",0.98,1.1);
    cout << "FitStatus= " << FitStatus << endl;
    HistSignal   -> SetTitle(Form("%s After BG Subtraction",HistSignal->GetTitle() ));

    TC_invM->Update();
    TC_invM->SaveAs(Form("PhiInvMassFlowV%d.pdf",FlowOrder));

  }
  tf_evt_nm_in->Close();
  tf_evt_mx_in->Close();
  TC_invM->SaveAs(Form("PhiInvMassFlowV%d.pdf]",FlowOrder));

  // return;
}
