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
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TMultiGraph.h"
#include "TChain.h"
#include "TLegend.h"
#include "TFractionFitter.h"
#include "TVirtualFitter.h"
#include "TCut.h"
#include "TObject.h"
#include "TSpline.h"
#include "TPaveText.h"
#include "Fit/FitResult.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TRandom3.h"
#include "TDatime.h"
#include "Math/MinimizerOptions.h"

using namespace std;

double sigmaRange = 5.;
double dParBg[3];
double dParSig[4];
double BackgroundFitting(double *x, double *p);
double TotalFitting(double *x, double *p);

void FlowExtractor()
{
  gStyle->SetOptStat(0);
  int FlowOrder = 2;
  double d_Flow[12];
  double d_Flow_err[12];
  TH1D *HistFlowVsPt = new TH1D("HistFlowVsPt",Form("#phi meson v_{%d} VS p_{T}",FlowOrder),14,-0.3,3.9);
  HistFlowVsPt   -> SetXTitle("p_{T} (GeV/c)");
  HistFlowVsPt   -> SetYTitle("v_{2}");
  HistFlowVsPt   -> GetYaxis()->SetRangeUser(-0.3,0.5);
  HistFlowVsPt   -> GetXaxis()->SetRangeUser(-0.1,2.1);

  // normal/mixed event invM input
  TFile * tf_evt_nm_in = new TFile("/mnt/c/Users/pjska/github/FlowExtractor/input/merged_KKinvMOutput.root","READ");
  TFile * tf_evt_mx_in = new TFile("/mnt/c/Users/pjska/github/FlowExtractor/input/merged_KKinvMOutputMixedEvents.root","READ");
  // flow VS Invariant Mass input
  TFile * tf_flow_invM_in = new TFile("/mnt/c/Users/pjska/github/FlowExtractor/input/phiV2VsInvMOutput.picoDst.result.root","READ");
  if( !tf_flow_invM_in->IsOpen() ) std::cout<<"No flow input!"<<std::endl;
  if(  tf_flow_invM_in->IsOpen() ) {
      std::cout<<"flow file loaded successfully!"<<std::endl;
  }

  //--- Normalization range
  double a_d_int_range[4] =
    {
      0.99,
      1.008,
      1.04,
      1.09
    };

  // TCanvas To draw nomal/mixed invariant mass and flow vs invariant mass
  TCanvas * TC_invM = new TCanvas("TC_invM","TC_invM",800,600);
  TC_invM->SaveAs(Form("PhiInvMassFlowV%d.pdf[",FlowOrder));
  TPaveText* label = new TPaveText(0.2,0.3,0.8,0.9);
  label->AddText(Form("#phi Meson Elliptic Flow v%d Analysis",FlowOrder));
  label->AddText("2018 Fixed-Target 26.5 (7.2) GeV ");
  label->Draw();
  TC_invM->Update();
  TC_invM->SaveAs(Form("PhiInvMassFlowV%d.pdf",FlowOrder));

  // Draw same and normalized-mx events together
  for(int pt=0; pt<7; pt++)
  {
    TPad* thePad = (TPad*)TC_invM->cd();

    thePad->cd();
    TH1D* sameEventInvM = (TH1D*) tf_evt_nm_in->Get(Form("histKKInvMpT%d",pt));
    sameEventInvM->SetTitle(Form("%1.1f <= p_{T} < %1.1f (GeV/c)",pt*0.3,(pt+1.)*0.3));
    sameEventInvM   -> SetXTitle("M_{K+K-} (GeV/c^{2})");
    sameEventInvM   -> SetYTitle("Counts");
    sameEventInvM->GetXaxis()->SetRangeUser(0.98,1.1);
    TH1D* mixedEventInvM = (TH1D*) tf_evt_mx_in->Get(Form("histKKInvMpT%d_Mixed",pt));
    mixedEventInvM->SetTitle(Form("%1.1f <= p_{T} < %1.1f (GeV/c)",pt*0.3,(pt+1.)*0.3));
    mixedEventInvM   -> SetXTitle("M_{K+K-} (GeV/c^{2})");
    mixedEventInvM   -> SetYTitle("Counts");
    mixedEventInvM->GetXaxis()->SetRangeUser(0.98,1.1);
    mixedEventInvM -> SetFillColor(kYellow);
    mixedEventInvM -> SetFillStyle(3001);
    TProfile* flowVsInvMass = (TProfile*) tf_flow_invM_in->Get(Form("histPhiV2VsInvMpT%d_pfx",pt));
    flowVsInvMass->SetTitle(Form("%1.1f <= p_{T} < %1.1f (GeV/c)",pt*0.3,(pt+1.)*0.3));
    flowVsInvMass   -> SetXTitle("M_{K+K-} (GeV/c^{2})");
    flowVsInvMass   -> SetYTitle("v_{2}");
    flowVsInvMass->GetXaxis()->SetRangeUser(0.98,1.1);
    flowVsInvMass->GetYaxis()->SetRangeUser(-1.,1.);

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
    //fit Signal with Gauss plus constant
    TFormula * GausPlus = new TFormula("GausPlus","gaus(0)+[3]");
    TF1 * tf1_Signal = new TF1("polygauss_single",GausPlus->GetExpFormula(),0.99,1.1);

    //fit to a simple gauss first to get seed
    TF1 * tf1_gauss = new TF1("tf1_gauss","gaus",0.9,1.1);
    HistSignal -> Fit(tf1_gauss,"0","R",1.01,1.03);
    // seeds
    double d_seeds_p0    = tf1_gauss -> GetParameter(0);
    double d_seeds_mean  = tf1_gauss -> GetParameter(1);
    double d_seeds_sigma = tf1_gauss -> GetParameter(2);

    tf1_Signal -> SetParameter(1,d_seeds_mean);
    tf1_Signal -> SetParameter(2,d_seeds_sigma);
    tf1_Signal -> SetParLimits(1,d_seeds_mean-d_seeds_sigma,d_seeds_mean+d_seeds_sigma);
    tf1_Signal -> SetParLimits(2,0.66*d_seeds_sigma,1.5*d_seeds_sigma);
    tf1_Signal -> SetLineColor(kBlue);

    int FitStatus = HistSignal   -> Fit(tf1_Signal,"E+","R",0.99,1.1);
    cout << "FitStatus= " << FitStatus << endl;
    HistSignal   -> SetTitle(Form("%s Signal fit with Gaussian",HistSignal->GetTitle() ));

    // To count how many #phi mesons HistSignal has
    dParSig[0]    = tf1_Signal -> GetParameter(0);
    dParSig[1]    = tf1_Signal -> GetParameter(1);
    dParSig[2]    = tf1_Signal -> GetParameter(2);
    dParSig[3]    = tf1_Signal -> GetParameter(3);
    int iBin_3sigint_low = HistSignal -> FindFixBin(dParSig[1] - (3*dParSig[2]));
    int iBin_3sigint_hi  = HistSignal -> FindFixBin(dParSig[1] + (3*dParSig[2]));
    double d_3sig_integral_error;
    double d_3sig_integral = HistSignal -> IntegralAndError(iBin_3sigint_low,iBin_3sigint_hi,d_3sig_integral_error,"");

    TPaveText * ptext = new TPaveText(0.1,0.65,0.30,0.9,"NDCARC");
    ptext -> AddText(Form("Mean: %.4f",dParSig[1]));
    ptext -> AddText(Form("Sigma: %.4f",dParSig[2]));
    ptext -> AddText(Form("3#sigma Int: %.4f",d_3sig_integral));
    ptext -> AddText(Form("3#sigma Int Err: %.4f",d_3sig_integral_error));
    ptext -> Draw("same");

    TC_invM->Update();
    TC_invM->SaveAs(Form("PhiInvMassFlowV%d.pdf",FlowOrder));

    // Fitting the background
    TH1D * HistBackground = (TH1D*) mixedEventInvM -> Clone("HistBackground");
    HistBackground   -> SetTitle(Form("%s Background fit with Pol2",HistBackground->GetTitle() ));

    TF1 * tf1_Background = new TF1("tf1_pol","[0] + [1]*x + [2]*x**2",(dParSig[1] - (sigmaRange*dParSig[2])),(dParSig[1] + (sigmaRange*dParSig[2])));
    HistBackground -> SetFillColor(kYellow);
    HistBackground -> SetFillStyle(3001);
    HistBackground->SetMaximum(HistBackground->GetBinContent(HistBackground->GetMaximumBin())*1.4);
    HistBackground->Draw();
    HistBackground->Fit(tf1_Background,"E","R",(dParSig[1] - (sigmaRange*dParSig[2])),(dParSig[1] + (sigmaRange*dParSig[2])));
    sameEventInvM->Draw("same");

    TC_invM->Update();
    TC_invM->SaveAs(Form("PhiInvMassFlowV%d.pdf",FlowOrder));

    dParBg[0]=tf1_Background->GetParameter(0);
    dParBg[1]=tf1_Background->GetParameter(1);
    dParBg[2]=tf1_Background->GetParameter(2);

    flowVsInvMass->Draw();
    TF1 * tf1_backgroundFlow = new TF1("tf1_backgroundFlow",BackgroundFitting,0.99,1.1,/*1*//*2*/3/*4*/);
    TF1 * tf1_totalFlow = new TF1("tf1_totalFlow",TotalFitting,0.99,1.1,/*2*//*3*/4/*5*/);

    flowVsInvMass->Fit(tf1_backgroundFlow,"E","R",0.99,1.1);
    double d_V2_bg_p0 = tf1_backgroundFlow->GetParameter(0);
    double d_V2_bg_p1 = tf1_backgroundFlow->GetParameter(1);
    double d_V2_bg_p2 = tf1_backgroundFlow->GetParameter(2);

    tf1_totalFlow->SetParameter(0,d_V2_bg_p0);
    tf1_totalFlow->SetParameter(1,d_V2_bg_p1);
    tf1_totalFlow->SetParameter(2,d_V2_bg_p2);

    tf1_totalFlow -> SetLineColor(kBlue);
    flowVsInvMass->Fit(tf1_totalFlow,"E+","R",0.99,1.1);

    d_Flow[pt]     = tf1_totalFlow->GetParameter(/*1*//*2*/3/*4*/);
    d_Flow_err[pt] = tf1_totalFlow->GetParError(/*1*//*2*/3/*4*/);
    HistFlowVsPt->SetBinContent(pt+2,d_Flow[pt]);
    HistFlowVsPt->SetBinError(pt+2,d_Flow_err[pt]);

    TPaveText * ptextFlow = new TPaveText(0.2,0.65,0.40,0.9,"NDCARC");
    ptextFlow -> AddText(Form("v2^{sig}: %.4f",d_Flow[pt]));
    ptextFlow -> AddText(Form("v2^{sig} Error: %.4f",d_Flow_err[pt]));
    ptextFlow->Draw("same");

    TC_invM->Update();
    TC_invM->SaveAs(Form("PhiInvMassFlowV%d.pdf",FlowOrder));

  }

  HistFlowVsPt->Draw();
  TLine * lineZero = new TLine(-0.3,0.,2.7,0.);
  lineZero->SetLineColor(1);
  lineZero->SetLineStyle(2);
  lineZero->Draw("same");
  TC_invM->Update();
  TC_invM->SaveAs(Form("PhiInvMassFlowV%d.pdf",FlowOrder));


  tf_evt_nm_in->Close();
  tf_evt_mx_in->Close();
  TC_invM->SaveAs(Form("PhiInvMassFlowV%d.pdf]",FlowOrder));

  // return;
}

// Fitting functions for Flow VS Invariant Mass
double BackgroundFitting(double *x, double *p)
{
  if( x[0] >= (dParSig[1] - (sigmaRange*dParSig[2])) && x[0] <= (dParSig[1] + (sigmaRange*dParSig[2])))
   {
     return (((dParBg[0] + dParBg[1] * x[0] + dParBg[2] * pow(x[0],2))
     /((dParSig[0]*exp(-0.5*pow(((x[0]-dParSig[1])/dParSig[2]),2))+dParSig[3])+(dParBg[0] + dParBg[1] * x[0] + dParBg[2] * pow(x[0],2))))
     */*p[0]*//*(p[0] + p[1] * x[0])*/(p[0] + p[1] * x[0] + p[2] * pow(x[0],2))/*(p[0] + p[1] * x[0] + p[2] * pow(x[0],2) + p[3]*pow(x[0],3)*/);
   }
   else
   {
     return /*p[0]*//*(p[0] + p[1] * x[0])*/(p[0] + p[1] * x[0] + p[2] * pow(x[0],2))/*(p[0] + p[1] * x[0] + p[2] * pow(x[0],2) + p[3]*pow(x[0],3)*/;
   }
}
double TotalFitting(double *x, double *p)
{
  if( x[0] >= (dParSig[1] - (sigmaRange*dParSig[2])) && x[0] <= (dParSig[1] + (sigmaRange*dParSig[2])))
   {
     return (
       (p[/*1*//*2*/3/*4*/]*((dParSig[0]*exp(-0.5*pow(((x[0]-dParSig[1])/dParSig[2]),2))+dParSig[3])
     /((dParSig[0]*exp(-0.5*pow(((x[0]-dParSig[1])/dParSig[2]),2))+dParSig[3])+(dParBg[0] + dParBg[1] * x[0] + dParBg[2] * pow(x[0],2)))))
     +
       (((dParBg[0] + dParBg[1] * x[0] + dParBg[2] * pow(x[0],2))
     /((dParSig[0]*exp(-0.5*pow(((x[0]-dParSig[1])/dParSig[2]),2))+dParSig[3])+(dParBg[0] + dParBg[1] * x[0] + dParBg[2] * pow(x[0],2))))
     */*p[0]*//*(p[0] + p[1] * x[0])*/(p[0] + p[1] * x[0] + p[2] * pow(x[0],2))/*(p[0] + p[1] * x[0] + p[2] * pow(x[0],2) + p[3]*pow(x[0],3))*/));
   }
   else
   {
     return /*p[0]*//*(p[0] + p[1] * x[0])*/(p[0] + p[1] * x[0] + p[2] * pow(x[0],2))/*(p[0] + p[1] * x[0] + p[2] * pow(x[0],2) + p[3]*pow(x[0],3))*/;
   }
}
