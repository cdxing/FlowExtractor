// C++ headers
#include <iostream>
#include <fstream>
#include "Riostream.h"
#include <cstdlib>
#include <string>
#include <vector>
#include <math.h>
#include <set>
#include <map>
#include <array>
//need this stuff to compile in linux:
#include <TROOT.h>
#include <TApplication.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TStreamerElement.h>
#include <TStyle.h>
#include "TSystemDirectory.h"

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
#include "TGaxis.h"
#include "TSpectrum.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TMultiGraph.h"
#include "TSpline.h"
#include "TPaveText.h"
#include "Fit/FitResult.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TRandom3.h"
#include "TDatime.h"
#include "Math/MinimizerOptions.h"
/* *****************************************************************************
 * This macro draws the QA plots for 7.2 GeV phi-meson v1 analysis
 *
 * Author: Ding Chen
 * Date: April 7, 2021
 *
 *******************************************************************************
*/
using namespace std;
const Double_t _y_CM = -2.03;

const Double_t _massPion     = 0.13957039;
const Double_t _massKaon     = 0.493677;
const Double_t _massProton   = 0.938272081;

const Double_t _massPhi      = 1.019461;
const Double_t _massKaonShort  = 0.497611;

const Double_t _massLambda  = 1.115683;
const Double_t _massOmegaMinus  = 1.67245;
const Double_t _massXiMinus  = 1.32171;

const Double_t _sigmaRange = 5.; // Sigma of the Fitting range
Double_t dParBg[3]; // Bkg fitting parameters
Double_t dParSig[4]; // Sig + Bkg fitting parameters
Double_t proportion(Double_t *x, Double_t *p);
Double_t d_v2nq(Double_t d_v2, Int_t ncq );
Double_t d_mTm0nq(Double_t d_pT, Int_t ncq, Double_t mass );
Double_t BackgroundFitting(Double_t *x, Double_t *p);
Double_t TotalFitting(Double_t *x, Double_t *p);

void plot7p2QA(){
  gStyle->SetOptStat(0);
  gStyle->SetOptDate(0);
  gStyle->SetEndErrorSize(6);
  gStyle->SetOptTitle(0);

  TFile * file_PID_QA = new TFile("/mnt/c/Users/pjska/github/FlowExtractor/res_sys/result_sys_flow/merged_sys_primary_var0_iter3_.root","READ");
  if( !file_PID_QA->IsOpen() ) std::cout<<"No PID QA input!"<<std::endl;
  if(  file_PID_QA->IsOpen() ) {
      std::cout<<"PID QA file loaded successfully!"<<std::endl;
  }
  TFile * file_phi_QA = new TFile("/mnt/c/Users/pjska/github/FlowExtractor/res_sys/result_sys_flow/hadd_PhiMesonAna_OUTPUT_sys_primary_var0_iter3_.root","READ");
  if( !file_phi_QA->IsOpen() ) std::cout<<"No phi QA input!"<<std::endl;
  if(  file_phi_QA->IsOpen() ) {
      std::cout<<"#phi QA file loaded successfully!"<<std::endl;
  }
  TH2D* h2_dedx_KP = (TH2D*) file_PID_QA->Get("hist_dEdx_kaonPlus");
  TH2D* h2_dedx_KM = (TH2D*) file_PID_QA->Get("hist_dEdx_kaonMinus");
  TH2D* h2_mass_KP = (TH2D*) file_PID_QA->Get("hist_mass_kaonPlus");
  TH2D* h2_mass_KM = (TH2D*) file_PID_QA->Get("hist_mass_kaonMinus");
  TH2D* h2_K_plus_pT_vs_y = (TH2D*) file_PID_QA->Get("hist_pt_y_kaonPlus");
  TH2D* h2_K_minus_pT_vs_y = (TH2D*) file_PID_QA->Get("hist_pt_y_kaonMinus");
  TH2D* h2_phi_pT_vs_y = (TH2D*) file_phi_QA->Get("hist_SE_pt_y_PhiMeson");
  TCanvas* c1_1 = new TCanvas("c1_1","dedx QA",200,0,1024,768);
  c1_1->SetLogz();
  c1_1->cd();
  h2_dedx_KP->GetZaxis()->SetRangeUser(1.,200000);
  h2_dedx_KP->GetYaxis()->SetRangeUser(-2.,10.5);
  h2_dedx_KP->GetXaxis()->SetRangeUser(-3.,3.);
  h2_dedx_KP->Draw("colz");
  h2_dedx_KM->GetZaxis()->SetRangeUser(1.,50000);
  h2_dedx_KM->Draw("Samecolz");
  // dedx distribution
  TPaveText * ptxt_prelim_dedx = new TPaveText(0.12,0.6,0.45,0.85,"NDCARC");
  ptxt_prelim_dedx->SetFillStyle(0);
  ptxt_prelim_dedx->SetBorderSize(0);
  ptxt_prelim_dedx->SetFillColor(0);
  ptxt_prelim_dedx -> AddText("STAR Au+Au #sqrt{s_{NN}}=7.2 GeV FXT ");
  // ptxt_prelim_dedx -> AddText("STAR preliminary");
  ptxt_prelim_dedx->Draw("same");

  TCanvas* c1_2 = new TCanvas("c1_2","mass2 QA",200,0,1024,768);
  c1_2->SetLogz();
  c1_2->cd();
  h2_mass_KP->GetZaxis()->SetRangeUser(1.,200000);
  h2_mass_KP->GetYaxis()->SetRangeUser(0.,0.55);
  h2_mass_KP->GetXaxis()->SetRangeUser(-3.,3.);
  h2_mass_KP->Draw("colz");
  h2_mass_KM->GetZaxis()->SetRangeUser(1.,200000);
  h2_mass_KM->Draw("Samecolz");
  TPaveText * ptxt_prelim_mass2 = new TPaveText(0.15,0.6,0.48,0.85,"NDCARC");
  ptxt_prelim_mass2->SetFillStyle(0);
  ptxt_prelim_mass2->SetBorderSize(0);
  ptxt_prelim_mass2->SetFillColor(0);
  ptxt_prelim_mass2 -> AddText("STAR Au+Au #sqrt{s_{NN}}=7.2 GeV FXT ");
  // ptxt_prelim_mass2 -> AddText("STAR preliminary");
  ptxt_prelim_mass2->Draw("same");
  // mass2 distribution

  TCanvas* c1 = new TCanvas("c1","K^{+} QA",200,0,1024,768);
  c1->SetLogz();
  c1->cd();
  h2_K_plus_pT_vs_y->GetYaxis()->SetRangeUser(0,3.5);
  h2_K_plus_pT_vs_y->GetXaxis()->SetLimits(-0.98,2.52);
  h2_K_plus_pT_vs_y->GetXaxis()->SetRangeUser(-0.08,2.02);
  h2_K_plus_pT_vs_y->Draw("colz");
  // K+ distribution

  TLine* tl = new TLine(0,0,0,3.5);
  tl->SetLineStyle(7);
  tl->Draw("same");

  TPaveText * ptxt_y_CM = new TPaveText(0.12,0.1,0.15,0.15,"NDCARC");
  ptxt_y_CM->SetFillColor(0);
  ptxt_y_CM -> AddText("y_{CM}");
  // ptxt_y_CM->Draw("same");
  TPaveText * ptxt_Kplus = new TPaveText(0.15,0.15,0.38,0.25,"NDCARC");
  ptxt_Kplus->SetFillColor(0);
  ptxt_Kplus -> AddText("K^{+}");
  ptxt_Kplus->Draw("same");
  TPaveText * ptxt_preli = new TPaveText(0.15,0.7,0.35,0.85,"NDCARC");
  ptxt_preli->SetFillColor(0);
  ptxt_preli -> AddText("Au+Au 7.2 GeV FXT");
  ptxt_preli -> AddText("y_{target} = 2.02");
  ptxt_preli -> AddText("STAR preliminary");
  ptxt_preli->Draw("same");

  TCanvas* c2 = new TCanvas("c2","#phi QA",200,0,1024,768);
  c2->SetLogz();
  c2->cd();

  h2_K_minus_pT_vs_y->GetYaxis()->SetRangeUser(0,3.5);
  h2_K_minus_pT_vs_y->GetXaxis()->SetLimits(-0.98,2.52);
  h2_K_minus_pT_vs_y->GetXaxis()->SetRangeUser(-0.08,2.02);
  h2_K_minus_pT_vs_y->Draw("colz");
  // #phi distribution
  tl->Draw("same");
  // ptxt_y_CM->Draw("same");
  TPaveText * ptxt_Kminus = new TPaveText(0.15,0.15,0.38,0.25,"NDCARC");
  ptxt_Kminus->SetFillColor(0);
  ptxt_Kminus -> AddText("K^{-}");
  ptxt_Kminus->Draw("same");
  ptxt_preli->Draw("same");

  TCanvas* c3 = new TCanvas("c3","#phi QA",200,0,1024,768);
  c3->SetLogz();
  c3->cd();
  h2_phi_pT_vs_y->GetXaxis()->SetLimits(-0.98,2.52);
  h2_phi_pT_vs_y->GetXaxis()->SetRangeUser(-0.08,2.02);
  h2_phi_pT_vs_y->GetYaxis()->SetRangeUser(0,3.5);
  h2_phi_pT_vs_y->GetYaxis()->SetTitle("p_{T} [GeV/c]");
  h2_phi_pT_vs_y->GetXaxis()->SetTitle("y");
  h2_phi_pT_vs_y->Draw("colz");
  // #phi distribution
  tl->Draw("same");
  // ptxt_y_CM->Draw("same");
  TPaveText * ptxt_phi = new TPaveText(0.15,0.15,0.38,0.25,"NDCARC");
  ptxt_phi->SetFillColor(0);
  ptxt_phi -> AddText("K^{+}K^{-} pair");
  ptxt_phi->Draw("same");
  ptxt_preli->Draw("same");

  // invM method
  TFile * file_KK_InvM_Input = new TFile("/mnt/c/Users/pjska/github/FlowExtractor/res_sys/result_sys_invM/merged_merged_sys_primary_var0_iter1_.root","READ");
  if( !file_KK_InvM_Input->IsOpen() ) std::cout<<"No SE/ME input!"<<std::endl;
  if(  file_KK_InvM_Input->IsOpen() ) {
      std::cout<<"#phi InvM loaded successfully!"<<std::endl;
  }
  // TFile * file_flow_invM_Input = new TFile("/mnt/c/Users/pjska/github/FlowExtractor/res_sys/result_sys_flow/hadd_PhiMesonAna_OUTPUT_sys_primary_var0_iter3_.root","READ");
  TFile * file_flow_invM_Input = new TFile("/mnt/c/Users/pjska/github/FlowExtractor/res_v2_7p2/0_EvtPlnR1/merged_merged_PhiMesonAna_OUTPUT_sys_primary_var0_iter2_652813E64F931A3A7865DC8AA3CF9F7E_.root","READ");
  if( !file_flow_invM_Input->IsOpen() ) std::cout<<"No flow input!"<<std::endl;
  if(  file_flow_invM_Input->IsOpen() ) {
      std::cout<<"flow file loaded successfully!"<<std::endl;
  }
  TH1D *Hist_Input_SE_InvM_rapSetA_centSetA;
  TH1D *Hist_Input_ME_InvM_rapSetA_centSetA;
  TProfile *Profile_Input_v1_reso_rapSetA_centSetA;
  Hist_Input_SE_InvM_rapSetA_centSetA = (TH1D*) file_KK_InvM_Input->Get("Hist_SE_InvM_rapSetA2_centSetA1");
  // Hist_Input_SE_InvM_rapSetA_centSetA = (TH1D*) file_flow_invM_Input->Get("Hist_SE_InvM_rapSetA2_centSetA1");
  Hist_Input_ME_InvM_rapSetA_centSetA = (TH1D*) file_KK_InvM_Input->Get("Hist_ME_InvM_rapSetA2_centSetA1");
  // Hist_Input_ME_InvM_rapSetA_centSetA = (TH1D*) file_flow_invM_Input->Get("Hist_rotation_InvM_rapSetA2_centSetA1");
  Profile_Input_v1_reso_rapSetA_centSetA = (TProfile*) file_flow_invM_Input->Get("Hist_v1_reso_rapSetA2_centSetA1_pfx");

  TCanvas* c4 = new TCanvas("c4","#phi invM",200,0,1024,768);
  c4->cd();
  Hist_Input_SE_InvM_rapSetA_centSetA->GetYaxis()->SetRangeUser(-0.1*(Double_t)Hist_Input_SE_InvM_rapSetA_centSetA->GetMaximum(),1.1*(Double_t)Hist_Input_SE_InvM_rapSetA_centSetA->GetMaximum());
  Hist_Input_SE_InvM_rapSetA_centSetA->GetXaxis()->SetRangeUser(0.99,1.09);
  Hist_Input_SE_InvM_rapSetA_centSetA->SetTitle("");
  Hist_Input_SE_InvM_rapSetA_centSetA->GetXaxis()->SetTitle("m_{inv} [GeV/c^{2}]");
  Hist_Input_SE_InvM_rapSetA_centSetA->GetXaxis()->SetTitleOffset(1);
  Hist_Input_ME_InvM_rapSetA_centSetA->GetXaxis()->SetRangeUser(0.99,1.09);
  // Get the bin of the Normalization range
  // ----------------------- Normalization range -------------------------------
  Double_t a_d_int_range[4] ={
    0.99,
    1.008,
    1.04,
    1.09
  };
  int a_iBin_range[4];
  for(int ijk = 0; ijk < 4; ijk++) a_iBin_range[ijk] =  Hist_Input_SE_InvM_rapSetA_centSetA -> FindFixBin(a_d_int_range[ijk]);
  //right bg Normalization
  Double_t d_r_area     = a_d_int_range[3]-a_d_int_range[2];
  Double_t d_r_same_int = Hist_Input_SE_InvM_rapSetA_centSetA -> Integral(a_iBin_range[2],a_iBin_range[3]);
  Double_t d_r_mx_int   = Hist_Input_ME_InvM_rapSetA_centSetA -> Integral(a_iBin_range[2],a_iBin_range[3]);
  Double_t d_r_norm     = (d_r_mx_int != 0.0) ? d_r_same_int/d_r_mx_int : 1.0;
  cout<<" R: "<<d_r_area<<" : "<<d_r_same_int<<" : "<<d_r_mx_int<<endl;
  //left bg Normalization
  Double_t d_l_area =  a_d_int_range[1]-a_d_int_range[0];
  Double_t d_l_same_int = Hist_Input_SE_InvM_rapSetA_centSetA -> Integral(a_iBin_range[0],a_iBin_range[1]);
  Double_t d_l_mx_int   = Hist_Input_ME_InvM_rapSetA_centSetA -> Integral(a_iBin_range[0],a_iBin_range[1]);
  Double_t d_l_norm     = (d_l_mx_int != 0.0) ? d_l_same_int/d_l_mx_int : 1.0;
  cout<<" L: "<<d_l_area<<" : "<<d_l_same_int<<" : "<<d_l_mx_int<<endl;
  cout<<" l: "<<d_l_norm<<" r: "<<d_r_norm<<endl;
  Double_t d_norm = ((d_r_norm/d_r_area) + (d_l_norm/d_l_area))/( (1.0/d_r_area) + (1.0/d_l_area));
  cout<<" d_norm = "<<d_norm<<endl;

  // Normalize mixed event invariant mass
  Hist_Input_SE_InvM_rapSetA_centSetA -> SetMarkerStyle(1);
  Hist_Input_ME_InvM_rapSetA_centSetA -> Sumw2();
  Hist_Input_ME_InvM_rapSetA_centSetA -> Scale(d_norm);
  Hist_Input_SE_InvM_rapSetA_centSetA->Draw();
  Hist_Input_ME_InvM_rapSetA_centSetA -> SetLineColor(kRed);
  Hist_Input_ME_InvM_rapSetA_centSetA -> SetFillColor(kRed);
  Hist_Input_ME_InvM_rapSetA_centSetA -> SetFillStyle(3002);
  Hist_Input_ME_InvM_rapSetA_centSetA->Draw("HISTsames");

  // Substract normalized ME from SE to get Signal
  TH1D * HistSignal = (TH1D*)Hist_Input_SE_InvM_rapSetA_centSetA -> Clone("HistSignal");
  // HistSignal->SetLineColor(kRed);
  HistSignal -> Reset();
  HistSignal -> Sumw2();
  for(int ijk = 1; ijk < HistSignal->GetNbinsX()+1; ijk++)
  {
    Double_t d_center   = Hist_Input_SE_InvM_rapSetA_centSetA -> GetBinCenter(ijk);
    Double_t d_same     = Hist_Input_SE_InvM_rapSetA_centSetA -> GetBinContent(ijk);
    Double_t d_same_err = Hist_Input_SE_InvM_rapSetA_centSetA -> GetBinError(ijk);
    Double_t d_mx       = Hist_Input_ME_InvM_rapSetA_centSetA -> GetBinContent(ijk);
    Double_t d_mx_err   = Hist_Input_ME_InvM_rapSetA_centSetA -> GetBinError(ijk);
    Double_t d_sig      = d_same - d_mx;
    Double_t d_sig_err  = sqrt(d_same_err*d_same_err+d_mx_err*d_mx_err);

    HistSignal -> SetBinContent(ijk,d_sig);
    HistSignal -> SetBinError(ijk,d_sig_err);
  }
  HistSignal->SetMarkerStyle(2);
  HistSignal->SetMarkerColor(kBlue);
  HistSignal->SetLineColor(kBlue);
  HistSignal -> Draw("HISTsamesE");
  // gStyle->SetOptFit(1111);

  //Fit function
  //fit Signal with Gauss plus constant
  TFormula * GausPlus = new TFormula("GausPlus","gaus(0)+[3]");
  TF1 * tf1_Signal = new TF1("polygauss_single",GausPlus->GetExpFormula(),0.99,1.09);
  //fit to a simple gauss first to get seed
  TF1 * tf1_gauss = new TF1("tf1_gauss","gaus",0.9,1.1);
  HistSignal -> Fit(tf1_gauss,"0","R",1.01,1.03);
  // seeds
  Double_t d_seeds_p0    = tf1_gauss -> GetParameter(0);
  Double_t d_seeds_mean  = tf1_gauss -> GetParameter(1);
  Double_t d_seeds_sigma = tf1_gauss -> GetParameter(2);

  tf1_Signal -> SetParameter(1,d_seeds_mean);
  tf1_Signal -> SetParameter(2,d_seeds_sigma);
  tf1_Signal -> SetParLimits(1,d_seeds_mean-d_seeds_sigma,d_seeds_mean+d_seeds_sigma);
  tf1_Signal -> SetParLimits(2,0.66*d_seeds_sigma,1.5*d_seeds_sigma);
  tf1_Signal -> SetLineColor(kBlue);

  int FitStatus = HistSignal   -> Fit(tf1_Signal,"E+","R",0.99,1.09);
  tf1_Signal->Draw("same");
  c4->Update();
  cout << "FitStatus= " << FitStatus << endl;
  // To count how many #phi mesons HistSignal has
  dParSig[0]    = tf1_Signal -> GetParameter(0);
  dParSig[1]    = tf1_Signal -> GetParameter(1);
  dParSig[2]    = tf1_Signal -> GetParameter(2);
  dParSig[3]    = tf1_Signal -> GetParameter(3);
  int iBin_3sigint_low = HistSignal -> FindFixBin(dParSig[1] - (3*dParSig[2]));
  int iBin_3sigint_hi  = HistSignal -> FindFixBin(dParSig[1] + (3*dParSig[2]));
  Double_t d_3sig_integral_error;
  Double_t d_3sig_integral = HistSignal -> IntegralAndError(iBin_3sigint_low,iBin_3sigint_hi,d_3sig_integral_error,"");
  TPaveText * ptext = new TPaveText(0.15,0.7,0.30,0.85,"NDCARC");
  ptext->SetFillColor(0);
  ptext -> AddText(Form("Mean: %.4f",dParSig[1]));
  ptext -> AddText(Form("Sigma: %.4f",dParSig[2]));
  ptext -> AddText(Form("#chi^{2}/NDF : %.2f / %d",(Double_t)tf1_Signal->GetChisquare(),(Int_t)tf1_Signal->GetNDF()));
  // ptext -> AddText(Form("3#sigma Int: %.4f",d_3sig_integral));
  // ptext -> AddText(Form("3#sigma Int Err: %.4f",d_3sig_integral_error));
  ptext -> Draw("same");
  TPaveText * ptxt_prelim_1 = new TPaveText(0.6,0.7,0.88,0.85,"NDCARC");
  ptxt_prelim_1->SetFillColor(0);
  ptxt_prelim_1 -> AddText("Au+Au 7.2 GeV FXT 10-40%");
  ptxt_prelim_1 -> AddText("1.02 < y < 1.52, y_{target} = 2.02");
  // ptxt_prelim_1 -> AddText("STAR preliminary");
  ptxt_prelim_1->Draw("same");



  // Fitting the background
  TF1 * tf1_Background = new TF1("tf1_pol","[0] + [1]*x + [2]*x**2",(dParSig[1] - (_sigmaRange*dParSig[2])),(dParSig[1] + (_sigmaRange*dParSig[2])));
  Hist_Input_ME_InvM_rapSetA_centSetA -> SetFillColor(kRed);
  Hist_Input_ME_InvM_rapSetA_centSetA -> SetFillStyle(3002);
  Hist_Input_ME_InvM_rapSetA_centSetA->Fit(tf1_Background,"E+","R",(dParSig[1] - (_sigmaRange*dParSig[2])),(dParSig[1] + (_sigmaRange*dParSig[2])));
  tf1_Background->SetLineColor(kRed);
  tf1_Background->Draw("same");
  TLegend *leg_invM = new TLegend(0.4,0.7,0.6,0.88);
  leg_invM->SetBorderSize(0);
  leg_invM->AddEntry(Hist_Input_SE_InvM_rapSetA_centSetA, "Data (Sig+Bkg)", "lp");
  leg_invM->AddEntry(HistSignal, "Data (Sig)", "lep");
  leg_invM->AddEntry(tf1_Signal, "Fit (Sig)", "l");
  leg_invM->AddEntry(Hist_Input_ME_InvM_rapSetA_centSetA, "Data (Bkg)", "f");
  leg_invM->AddEntry(tf1_Background,"Fit (Bkg) ","l");
  leg_invM->Draw("same");

  dParBg[0]=tf1_Background->GetParameter(0);
  dParBg[1]=tf1_Background->GetParameter(1);
  dParBg[2]=tf1_Background->GetParameter(2);
  TF1 * tf1_backgroundFlow = new TF1("tf1_backgroundFlow",BackgroundFitting,0.99,1.09,/*1*//*2*/3/*4*/);
  TF1 * tf1_totalFlow = new TF1("tf1_totalFlow",TotalFitting,0.99,1.09,/*2*//*3*/4/*5*/);

  TCanvas* c5 = new TCanvas("c5","#phi v1",200,0,1024,768);
  c5->cd();
  Profile_Input_v1_reso_rapSetA_centSetA->GetYaxis()->SetTitle("v_{1}(m_{inv})");
  Profile_Input_v1_reso_rapSetA_centSetA->GetYaxis()->SetTitleOffset(1.2);
  Profile_Input_v1_reso_rapSetA_centSetA->GetXaxis()->SetRangeUser(0.99,1.09);
  Profile_Input_v1_reso_rapSetA_centSetA->GetXaxis()->SetTitleOffset(1.);
  Profile_Input_v1_reso_rapSetA_centSetA->Draw();
  Profile_Input_v1_reso_rapSetA_centSetA->Fit(tf1_backgroundFlow,"E+","R",0.99,1.09);
  Double_t d_V2_bg_p0 = tf1_backgroundFlow->GetParameter(0);
  Double_t d_V2_bg_p1 = tf1_backgroundFlow->GetParameter(1);
  Double_t d_V2_bg_p2 = tf1_backgroundFlow->GetParameter(2);
  tf1_totalFlow->SetParameter(0,d_V2_bg_p0);
  tf1_totalFlow->SetParameter(1,d_V2_bg_p1);
  tf1_totalFlow->SetParameter(2,d_V2_bg_p2);
  tf1_totalFlow -> SetLineColor(kBlue);
  Profile_Input_v1_reso_rapSetA_centSetA->Fit(tf1_totalFlow,"E+","R",0.99,1.09);
  Double_t d_FLow_rapSetA_centSetA = tf1_totalFlow->GetParameter(/*1*//*2*/3/*4*/);
  Double_t d_Flow_err_rapSetA_centSetA = tf1_totalFlow->GetParError(/*1*//*2*/3/*4*/);
  TPaveText * ptxt_phi_v1 = new TPaveText(0.2,0.15,0.4,0.25,"NDCARC");
  ptxt_phi_v1->SetFillColor(0);
  ptxt_phi_v1 -> AddText(Form("v_{1}^{sig}: %.4f %c %.4f",d_FLow_rapSetA_centSetA,177,d_Flow_err_rapSetA_centSetA));
  ptxt_phi_v1 -> AddText(Form("#chi^{2}/NDF : %.2f / %d",(Double_t)tf1_totalFlow->GetChisquare(),(Int_t)tf1_totalFlow->GetNDF()));
  ptxt_phi_v1->Draw("same");
  TLegend *leg_flowExtract = new TLegend(0.4,0.75,0.68,0.88);
  leg_flowExtract->SetBorderSize(0);
  leg_flowExtract->AddEntry(Profile_Input_v1_reso_rapSetA_centSetA, "Data (Sig+Bkg)", "lep");
  leg_flowExtract->AddEntry(tf1_totalFlow, "Fit (Sig + Bkg)", "l");
  leg_flowExtract->AddEntry(tf1_backgroundFlow,"Fit (Bkg) ","l");
  leg_flowExtract->Draw("same");
  TPaveText * ptxt_prelim_2 = new TPaveText(0.55,0.15,0.85,0.3,"NDCARC");
  ptxt_prelim_2->SetFillColor(0);
  ptxt_prelim_2 -> AddText("Au+Au 7.2 GeV FXT 10-40%");
  ptxt_prelim_2 -> AddText("1.02 < y < 1.52, y_{target} = 2.02");
  ptxt_prelim_2 -> AddText("STAR preliminary");
  ptxt_prelim_2->Draw("same");
}

Double_t proportion(Double_t *x, Double_t *p)
{
  return ( p[0] * x[0]);
}

Double_t d_v2nq(Double_t d_v2, Int_t ncq ){
  Double_t v2nq = d_v2/(Double_t)ncq;
  return v2nq;
}
Double_t d_mTm0nq(Double_t d_pT, Int_t ncq, Double_t mass ){ // pT to (mT - m0)/nq
  Double_t mTm0nq =  (sqrt(mass * mass + d_pT * d_pT) - mass) / (Double_t)ncq;
  return mTm0nq;
}

// Fitting functions for Flow VS Invariant Mass
Double_t BackgroundFitting(Double_t *x, Double_t *p)
{
  if( x[0] >= (dParSig[1] - (_sigmaRange*dParSig[2])) && x[0] <= (dParSig[1] + (_sigmaRange*dParSig[2])))
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
Double_t TotalFitting(Double_t *x, Double_t *p)
{
  if( x[0] >= (dParSig[1] - (_sigmaRange*dParSig[2])) && x[0] <= (dParSig[1] + (_sigmaRange*dParSig[2])))
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
