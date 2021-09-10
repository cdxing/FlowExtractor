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
#include "TProfile.h"
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
#include "FlowExtractor.h"

// 27 GeV
using namespace std;
// -------------------------- set Some fitting prerequsites --------------------

Double_t dParBg[3]; // Bkg fitting parameters
Double_t dParSig[4]; // Sig + Bkg fitting parameters
Double_t proportion(Double_t *x, Double_t *p);
Double_t BackgroundFitting(Double_t *x, Double_t *p);
Double_t TotalFitting(Double_t *x, Double_t *p);

// ======================== (1) Analysis Start =================================
void FlowExtractor( /*TString invMFileName = "./res_sys/result_sys_invM/merged_merged_sys_primary_var0_iter1_.root",*/
                   // TString FlowFileName =  "./res_sys/result_sys_flow/hadd_PhiMesonAna_OUTPUT_sys_primary_var0_iter3_.root" ,
                   TString FlowFileName =  "../out_FlowHists.root" ,
                    // double inputParameter1 = 0.
                    Int_t   inputp2 = 0, // sysErr cut Indexes 0-15
                    Int_t   inputp3 = 0, // sysErr cut variations, each systematic check has 2 or 3 vertions
                    Int_t   inputp4 = 3 // Iteration of the analysis is. In this analysis, 2 iterations is enough
){
  Int_t sys_cutN = inputp2; // sysErr cut Indexes 0-15
  Int_t sys_varN = inputp3; // sysErr cut variations, each systematic check has 2 or 3 vertions
  Int_t sys_iterN = inputp4; // Iteration of the analysis is. In this analysis, 2 iterations is enough
  string sys_object[17]  = {"primary", "etaGap", "etaRange",
                            "vz", "vr", "dedx", "dca",
                            "nHitsFit", "ratio", "nSigK", "mass2",
                            "pT", "dipAngle", "vtxDiff", "mthdDiff",
                            "binning",
                            "TPCpid"};
  std::cout << "sys_cutN == "<< sys_cutN <<": "<< sys_object[sys_cutN] << std::endl;
  TString outTxt = "./out/out_flow_invM_";
  TString outHead = outTxt;
  outTxt.Append(sys_object[sys_cutN]);
  outTxt.Append(Form("_var%d_iter%d_", sys_varN, sys_iterN));
  outTxt.Append(".txt");
  std::ofstream flowFile(outTxt,ofstream::out);
  // ---------------------- Analysis Setup -------------------------------------
  // gStyle->SetOptStat(0);
  gROOT->ForceStyle();
  // ----- InvMass plots in different centrality and pT or y bins --------------

  // directed and elliptic flow. Indexes: 0: v1, 1: v2; raw, reso; pT/y SetA; cent SetA; jkk
  Double_t d_FLow_typeID_cent[2][2][8][4]; // pt SetA, cent SetA
  Double_t d_Flow_err_typeID_cent[2][2][8][4]; // pt SetA, cent SetA
  // ---------------------- Input files and plots ------------------------------
  // SE/ME invM input
  TFile * file_flow_invM_Input = new TFile(FlowFileName,"READ");
  if( !file_flow_invM_Input->IsOpen() ) std::cout<<"No flow input!"<<std::endl;
  if(  file_flow_invM_Input->IsOpen() ) {
      std::cout<<"flow file loaded successfully!"<<std::endl;
  }
  TH1F *mHist_Input_SE_InvM_typeID_cent[8][4];
  TH1F *mHist_Input_ME_InvM_typeID_cent[8][4];
  TProfile *mProfile_Input_flow_reso_typeID_cent[8][4];
  TString Centrality_01[4] = {"0080","0010","1040","4080"};

  for(Int_t cent = 0; cent < Bin_Centrality_01; cent++)
  {
      for(Int_t rap_bin = 0; rap_bin < Bin_rap; rap_bin++)
      {
        mHist_Input_SE_InvM_typeID_cent[rap_bin][cent] = (TH1F*) file_flow_invM_Input->Get(Form("InvMass_SE_rapbin%d_cent%s",rap_bin+1,Centrality_01[cent].Data()));
        mHist_Input_ME_InvM_typeID_cent[rap_bin][cent] = (TH1F*) file_flow_invM_Input->Get(Form("InvMass_rot_rapbin%d_cent%s",rap_bin+1,Centrality_01[cent].Data()));
        mProfile_Input_flow_reso_typeID_cent[rap_bin][cent] = (TProfile*) file_flow_invM_Input->Get(Form("flow_InvMass_rapbin%d_cent%s",rap_bin+1,Centrality_01[cent].Data()));
  }
  }
  // ---------------------- Output files and plots -----------------------------
  TString outFile = ".phiflow.result.root";
  outFile.Prepend(Form("_var%d_iter%d_", sys_varN, sys_iterN));
  outFile.Prepend(sys_object[sys_cutN]);
  outFile.Prepend(outHead);
  // outFile.Prepend("./out_sys/out_sys_Crosscheck/");
  TFile *outputFile = new TFile(outFile,"recreate");
  TCanvas *canvas_InvM_typeID_cent = new TCanvas("canvas_InvM_typeID_cent","canvas_InvM_typeID_cent",1920,1080);
  TCanvas *canvas_flow_reso_typeID_cent = new TCanvas("canvas_flow_reso_typeID_cent","canvas_flow_reso_typeID_cent",1920,1080);
  canvas_InvM_typeID_cent->Divide(8,4);
  canvas_flow_reso_typeID_cent->Divide(8,4);
  TCanvas *canvas_flow_reso_dist_typeID_cent = new TCanvas("canvas_flow_reso_dist_typeID_cent","canvas_flow_reso_dist_typeID_cent",1920,1080);
  canvas_flow_reso_dist_typeID_cent->Divide(2,2);
  TGraphErrors *mTGE_flow_reso_dist_typeID_cent[4];

  // ======================== (2) Fit SE and ME InvM plots =====================
  // ----------------------- Normalization range -------------------------------

  for(int type_id=0; type_id<Bin_rap; type_id++)
  {
    for(int cent=0; cent<Bin_Centrality_01;cent++){
      canvas_InvM_typeID_cent->cd((type_id+1)+8*cent);
      mHist_Input_SE_InvM_typeID_cent[type_id][cent]->GetYaxis()->SetRangeUser(-0.1*(Double_t)mHist_Input_SE_InvM_typeID_cent[type_id][cent]->GetMaximum(),1.1*(Double_t)mHist_Input_SE_InvM_typeID_cent[type_id][cent]->GetMaximum());
      mHist_Input_SE_InvM_typeID_cent[type_id][cent]->GetXaxis()->SetRangeUser(0.99,1.08);
      mHist_Input_ME_InvM_typeID_cent[type_id][cent]->GetXaxis()->SetRangeUser(0.99,1.08);
      // Get the bin of the Normalization range
      int a_iBin_range[4];
      for(int ijk = 0; ijk < 4; ijk++) a_iBin_range[ijk] =  mHist_Input_SE_InvM_typeID_cent[type_id][cent] -> FindFixBin(a_d_int_range[ijk]);
      //right bg Normalization
      Double_t d_r_area     = a_d_int_range[3]-a_d_int_range[2];
      Double_t d_r_same_int = mHist_Input_SE_InvM_typeID_cent[type_id][cent] -> Integral(a_iBin_range[2],a_iBin_range[3]);
      Double_t d_r_mx_int   = mHist_Input_ME_InvM_typeID_cent[type_id][cent] -> Integral(a_iBin_range[2],a_iBin_range[3]);
      Double_t d_r_norm     = (d_r_mx_int != 0.0) ? d_r_same_int/d_r_mx_int : 1.0;
      cout<<" R: "<<d_r_area<<" : "<<d_r_same_int<<" : "<<d_r_mx_int<<endl;
      //left bg Normalization
      Double_t d_l_area =  a_d_int_range[1]-a_d_int_range[0];
      Double_t d_l_same_int = mHist_Input_SE_InvM_typeID_cent[type_id][cent] -> Integral(a_iBin_range[0],a_iBin_range[1]);
      Double_t d_l_mx_int   = mHist_Input_ME_InvM_typeID_cent[type_id][cent] -> Integral(a_iBin_range[0],a_iBin_range[1]);
      Double_t d_l_norm     = (d_l_mx_int != 0.0) ? d_l_same_int/d_l_mx_int : 1.0;
      cout<<" L: "<<d_l_area<<" : "<<d_l_same_int<<" : "<<d_l_mx_int<<endl;
      cout<<" l: "<<d_l_norm<<" r: "<<d_r_norm<<endl;
      Double_t d_norm = ((d_r_norm/d_r_area) + (d_l_norm/d_l_area))/( (1.0/d_r_area) + (1.0/d_l_area));
      cout<<" d_norm = "<<d_norm<<endl;

      // Normalize mixed event invariant mass
      mHist_Input_SE_InvM_typeID_cent[type_id][cent] -> SetMarkerStyle(1);
      mHist_Input_ME_InvM_typeID_cent[type_id][cent] -> Sumw2();
      mHist_Input_ME_InvM_typeID_cent[type_id][cent] -> Scale(d_norm);
      mHist_Input_SE_InvM_typeID_cent[type_id][cent]->Draw();
      mHist_Input_ME_InvM_typeID_cent[type_id][cent] -> SetLineColor(kRed);
      mHist_Input_ME_InvM_typeID_cent[type_id][cent] -> SetFillColor(kRed);
      mHist_Input_ME_InvM_typeID_cent[type_id][cent] -> SetFillStyle(3002);
      mHist_Input_ME_InvM_typeID_cent[type_id][cent]->Draw("HISTsames");
      // Substract normalized ME from SE to get Signal
      TH1D * HistSignal = (TH1D*) mHist_Input_SE_InvM_typeID_cent[type_id][cent] -> Clone("HistSignal");
      // HistSignal->SetLineColor(kRed);
      HistSignal -> Reset();
      HistSignal -> Sumw2();

      for(int ijk = 1; ijk < HistSignal->GetNbinsX()+1; ijk++)
      {
        Double_t d_center   = mHist_Input_SE_InvM_typeID_cent[type_id][cent] -> GetBinCenter(ijk);
        Double_t d_same     = mHist_Input_SE_InvM_typeID_cent[type_id][cent] -> GetBinContent(ijk);
        Double_t d_same_err = mHist_Input_SE_InvM_typeID_cent[type_id][cent] -> GetBinError(ijk);
        Double_t d_mx       = mHist_Input_ME_InvM_typeID_cent[type_id][cent] -> GetBinContent(ijk);
        Double_t d_mx_err   = mHist_Input_ME_InvM_typeID_cent[type_id][cent] -> GetBinError(ijk);
        Double_t d_sig      = d_same - d_mx;
        Double_t d_sig_err  = sqrt(d_same_err*d_same_err+d_mx_err*d_mx_err);

        HistSignal -> SetBinContent(ijk,d_sig);
        HistSignal -> SetBinError(ijk,d_sig_err);
      }
      HistSignal->SetMarkerStyle(2);
      HistSignal->SetMarkerColor(kBlue);
      HistSignal->SetLineColor(kBlue);
      HistSignal -> Draw("HISTsamesE");
      gStyle->SetOptFit(1111);

      //Fit function
      //fit Signal with Gauss plus constant
      TFormula * GausPlus = new TFormula("GausPlus","gaus(0)+[3]");
      TF1 * tf1_Signal = new TF1("polygauss_single",GausPlus->GetExpFormula(),0.99,1.08);
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

      int FitStatus = HistSignal   -> Fit(tf1_Signal,"E+","R",0.99,1.08);
      tf1_Signal->Draw("same");
      canvas_InvM_typeID_cent->cd((type_id+1)+8*cent)->Update();
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
      TPaveText * ptext = new TPaveText(0.1,0.65,0.30,0.9,"NDCARC");
      ptext -> AddText(Form("Mean: %.4f",dParSig[1]));
      ptext -> AddText(Form("Sigma: %.4f",dParSig[2]));
      ptext -> AddText(Form("3#sigma Int: %.4f",d_3sig_integral));
      ptext -> AddText(Form("3#sigma Int Err: %.4f",d_3sig_integral_error));
      ptext -> Draw("same");
      // Fitting the background
      TF1 * tf1_Background = new TF1("tf1_pol","[0] + [1]*x + [2]*x**2",(dParSig[1] - (_sigmaRange*dParSig[2])),(dParSig[1] + (_sigmaRange*dParSig[2])));
      mHist_Input_ME_InvM_typeID_cent[type_id][cent] -> SetFillColor(kRed);
      mHist_Input_ME_InvM_typeID_cent[type_id][cent] -> SetFillStyle(3002);
      mHist_Input_ME_InvM_typeID_cent[type_id][cent]->Fit(tf1_Background,"E+","R",(dParSig[1] - (_sigmaRange*dParSig[2])),(dParSig[1] + (_sigmaRange*dParSig[2])));
      tf1_Background->SetLineColor(kRed);
      tf1_Background->Draw("same");
      dParBg[0]=tf1_Background->GetParameter(0);
      dParBg[1]=tf1_Background->GetParameter(1);
      dParBg[2]=tf1_Background->GetParameter(2);
      TF1 * tf1_backgroundFlow = new TF1("tf1_backgroundFlow",BackgroundFitting,0.99,1.08,/*1*//*2*/3/*4*/);
      TF1 * tf1_totalFlow = new TF1("tf1_totalFlow",TotalFitting,0.99,1.08,/*2*//*3*/4/*5*/);

      Double_t d_V2_bg_p0;// = tf1_backgroundFlow->GetParameter(0);
      Double_t d_V2_bg_p1;// = tf1_backgroundFlow->GetParameter(1);
      Double_t d_V2_bg_p2;// = tf1_backgroundFlow->GetParameter(2);

      canvas_flow_reso_typeID_cent->cd((type_id+1)+8*cent);
      mProfile_Input_flow_reso_typeID_cent[type_id][cent]->GetXaxis()->SetRangeUser(0.99,1.08);
      mProfile_Input_flow_reso_typeID_cent[type_id][cent]->Draw();
      mProfile_Input_flow_reso_typeID_cent[type_id][cent]->Fit(tf1_backgroundFlow,"E+","R",0.99,1.08);
      d_V2_bg_p0 = tf1_backgroundFlow->GetParameter(0);
      d_V2_bg_p1 = tf1_backgroundFlow->GetParameter(1);
      d_V2_bg_p2 = tf1_backgroundFlow->GetParameter(2);
      tf1_totalFlow->SetParameter(0,d_V2_bg_p0);
      tf1_totalFlow->SetParameter(1,d_V2_bg_p1);
      tf1_totalFlow->SetParameter(2,d_V2_bg_p2);
      tf1_totalFlow -> SetLineColor(kBlue);
      mProfile_Input_flow_reso_typeID_cent[type_id][cent]->Fit(tf1_totalFlow,"E+","R",0.99,1.08);
      d_FLow_typeID_cent[0][1][type_id][cent] = tf1_totalFlow->GetParameter(/*1*//*2*/3/*4*/);
      d_Flow_err_typeID_cent[0][1][type_id][cent] = tf1_totalFlow->GetParError(/*1*//*2*/3/*4*/);
      TPaveText * ptextFlow_v2_reso_typeID_cent = new TPaveText(0.2,0.8,0.6,0.9,"NDCARC");
      ptextFlow_v2_reso_typeID_cent -> AddText(Form("v_{1}^{sig}: %.4f %c %.4f",d_FLow_typeID_cent[0][1][type_id][cent],177,d_Flow_err_typeID_cent[0][1][type_id][cent]));
      ptextFlow_v2_reso_typeID_cent -> AddText(Form("#chi^{2}/NDF : %.2f / %d",(Double_t)tf1_totalFlow->GetChisquare(),(Int_t)tf1_totalFlow->GetNDF()));
      ptextFlow_v2_reso_typeID_cent->Draw("same");
    }
  }
  const int n_typeID_cent = 8;
  for(int cent=0; cent<Bin_Centrality_01;cent++){
    TLine *l1_typeID_cent = new TLine(-2,0,2,0);
    l1_typeID_cent->SetLineStyle(2);
    Double_t x[n_typeID_cent] = {-0.8, -0.45, -0.2, -0.05, 0.05, 0.2, 0.45, 0.8};
    Double_t ex[n_typeID_cent] = {0.2, 0.15, 0.1, 0.05, 0.05, 0.1, 0.15, 0.2};

    Double_t y_v2_reso[n_typeID_cent] = {d_FLow_typeID_cent[0][1][0][cent], d_FLow_typeID_cent[0][1][1][cent],
      d_FLow_typeID_cent[0][1][2][cent],d_FLow_typeID_cent[0][1][3][cent],
    d_FLow_typeID_cent[0][1][4][cent],d_FLow_typeID_cent[0][1][5][cent],
  d_FLow_typeID_cent[0][1][6][cent],d_FLow_typeID_cent[0][1][7][cent]};
    Double_t ey_v2_reso[n_typeID_cent] = {d_Flow_err_typeID_cent[0][1][0][cent], d_Flow_err_typeID_cent[0][1][1][cent],
      d_Flow_err_typeID_cent[0][1][2][cent], d_Flow_err_typeID_cent[0][1][3][cent]
    , d_Flow_err_typeID_cent[0][1][4][cent], d_Flow_err_typeID_cent[0][1][5][cent]
  , d_Flow_err_typeID_cent[0][1][6][cent], d_Flow_err_typeID_cent[0][1][7][cent]};
    canvas_flow_reso_dist_typeID_cent->cd(cent+1);
    mTGE_flow_reso_dist_typeID_cent[cent] = new TGraphErrors(n_typeID_cent,x,y_v2_reso,ex,ey_v2_reso);
    mTGE_flow_reso_dist_typeID_cent[cent]->GetXaxis()->SetTitle("y");
    mTGE_flow_reso_dist_typeID_cent[cent]->GetYaxis()->SetTitle("v_{1}");
    if(cent<3){
      mTGE_flow_reso_dist_typeID_cent[cent]->GetXaxis()->SetRangeUser(-1.1,1.1);
      mTGE_flow_reso_dist_typeID_cent[cent]->GetYaxis()->SetRangeUser(-0.15,0.15);
    } else{
      mTGE_flow_reso_dist_typeID_cent[cent]->GetXaxis()->SetRangeUser(-1.1,1.1);
      mTGE_flow_reso_dist_typeID_cent[cent]->GetYaxis()->SetRangeUser(-0.15,0.15);
    }
    mTGE_flow_reso_dist_typeID_cent[cent]->SetMarkerColor(4);
    mTGE_flow_reso_dist_typeID_cent[cent]->SetMarkerStyle(24);
    mTGE_flow_reso_dist_typeID_cent[cent]->Draw("AP");
    l1_typeID_cent->Draw("same");
  }
  mTGE_flow_reso_dist_typeID_cent[0]->SetTitle(Form("v_{1}, %s",Centrality_01[0].Data()));
  mTGE_flow_reso_dist_typeID_cent[1]->SetTitle(Form("v_{1}, %s",Centrality_01[1].Data()));
  mTGE_flow_reso_dist_typeID_cent[2]->SetTitle(Form("v_{1}, %s",Centrality_01[2].Data()));
  mTGE_flow_reso_dist_typeID_cent[3]->SetTitle(Form("v_{1}, %s",Centrality_01[3].Data()));
  gStyle->SetOptFit(1);
  TF1 * tf1_dv1dy[4];
  TPaveText * ptext_flow_type[4];
  for(int i = 0; i<4; i++){
    tf1_dv1dy[i] = new TF1(Form("myfunction_%d",i),proportion,-1.,1.,1);
    mTGE_flow_reso_dist_typeID_cent[i]->Fit(tf1_dv1dy[i],"E+","R",-0.6,0.6);
    ptext_flow_type[i] = new TPaveText(0.2,0.8,0.6,0.9,"NDCARC");
    ptext_flow_type[i] -> AddText(Form("dv_{1}/dy: %.4f %c %.4f",(Double_t)tf1_dv1dy[i]->GetParameter(0),177,(Double_t)tf1_dv1dy[i]->GetParError(0)));
    ptext_flow_type[i] -> AddText(Form("#chi^{2}/NDF : %.2f / %d",(Double_t)tf1_dv1dy[i]->GetChisquare(),(Int_t)tf1_dv1dy[i]->GetNDF()));
    canvas_flow_reso_dist_typeID_cent->cd(i+1);
    ptext_flow_type[i]->Draw("same");
  }

  cout << endl;
  flowFile << Centrality_01[0].Data()<<endl;
  for(Int_t rap_bin = 0; rap_bin < Bin_rap; rap_bin++){
    flowFile << d_FLow_typeID_cent[0][1][rap_bin][0] <<", ";
  }
  flowFile<<endl;
  for(Int_t rap_bin = 0; rap_bin < Bin_rap; rap_bin++){
    flowFile << d_Flow_err_typeID_cent[0][1][rap_bin][0] <<", ";
  }
  flowFile<<endl;
  flowFile << Centrality_01[1].Data()<<endl;
  for(Int_t rap_bin = 0; rap_bin < Bin_rap; rap_bin++){
    flowFile << d_FLow_typeID_cent[0][1][rap_bin][1] <<", ";
  }
  flowFile<<endl;
  for(Int_t rap_bin = 0; rap_bin < Bin_rap; rap_bin++){
    flowFile << d_Flow_err_typeID_cent[0][1][rap_bin][1] <<", ";
  }
  flowFile<<endl;
  flowFile << Centrality_01[2].Data()<<endl;
  for(Int_t rap_bin = 0; rap_bin < Bin_rap; rap_bin++){
    flowFile << d_FLow_typeID_cent[0][1][rap_bin][2] <<", ";
  }
  flowFile<<endl;
  for(Int_t rap_bin = 0; rap_bin < Bin_rap; rap_bin++){
    flowFile << d_Flow_err_typeID_cent[0][1][rap_bin][2] <<", ";
  }
  flowFile<<endl;
  flowFile << Centrality_01[3].Data()<<endl;
  for(Int_t rap_bin = 0; rap_bin < Bin_rap; rap_bin++){
    flowFile << d_FLow_typeID_cent[0][1][rap_bin][3] <<", ";
  }
  flowFile<<endl;
  for(Int_t rap_bin = 0; rap_bin < Bin_rap; rap_bin++){
    flowFile << d_Flow_err_typeID_cent[0][1][rap_bin][3] <<", ";
  }
  flowFile<<endl;
  cout << endl;
  outputFile->cd();
  canvas_InvM_typeID_cent->Write();
  canvas_flow_reso_typeID_cent->Write();
  canvas_flow_reso_dist_typeID_cent->Write();
  flowFile.close();
}

Double_t proportion(Double_t *x, Double_t *p)
{
  return ( p[0] * x[0]);
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
