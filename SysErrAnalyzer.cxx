#include <iostream>
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

#include "Riostream.h"
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
#include "TGaxis.h"
#include "TSpectrum.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TMultiGraph.h"
#include "TSpline.h"
#include "TPaveText.h"
#include "TProfile.h"
#include "Fit/FitResult.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TRandom3.h"
#include "TDatime.h"
#include "Math/MinimizerOptions.h"

using namespace std;
const Double_t _y_CM = -2.03;
const Int_t _N_sysCut = 13;
const Double_t _massPion     = 0.13957039;
const Double_t _massKaon     = 0.493677;
const Double_t _massProton   = 0.938272081;

const Double_t _massPhi      = 1.019461;
const Double_t _massKaonShort  = 0.497611;

const Double_t _massLambda  = 1.115683;
const Double_t _massOmegaMinus  = 1.67245;
const Double_t _massXiMinus  = 1.32171;
void PaintBin (TH1D *h, Int_t bin, Int_t color);
void SysErrAnalyzer()
{
  // (1) ================ Input:  v1 and dv1/dy of different variations  =====================================
  string sys_object[17]  = {"primary", "etaGap", "etaRange",
                            "vz", "vr", "dedx", "dca",
                            "nHitsFit", "ratio", "nSigK", "mass2",
                            "pT", "dipAngle", "vtxDiff", "mthdDiff",
                            "binning",
                            "TPCpid"};
  std::ifstream inputReso0("/mnt/c/Users/pjska/github/FlowExtractor/out_sys/phi_v1_y_sys_primary_var0_iter3_.txt");
  Double_t d_v1_dv1dy_statErr0[8];
  // primary
  std::ifstream inputReso1_1("/mnt/c/Users/pjska/github/FlowExtractor/out_sys/phi_v1_y_sys_etaGap_var1_iter3_.txt");
  std::ifstream inputReso1_2("/mnt/c/Users/pjska/github/FlowExtractor/out_sys/phi_v1_y_sys_etaGap_var2_iter3_.txt");
  Double_t d_v1_dv1dy_statErr1_1[8];
  Double_t d_v1_dv1dy_statErr1_2[8];
  // etaGap
  std::ifstream inputReso2_1("/mnt/c/Users/pjska/github/FlowExtractor/out_sys/phi_v1_y_sys_etaRange_var2_iter3_.txt");
  Double_t d_v1_dv1dy_statErr2_1[8];
  // etaRange
  std::ifstream inputReso3_1("/mnt/c/Users/pjska/github/FlowExtractor/out_sys/phi_v1_y_sys_vz_var1_iter3_.txt");
  std::ifstream inputReso3_2("/mnt/c/Users/pjska/github/FlowExtractor/out_sys/phi_v1_y_sys_vz_var2_iter3_.txt");
  Double_t d_v1_dv1dy_statErr3_1[8];
  Double_t d_v1_dv1dy_statErr3_2[8];
  // vz
  std::ifstream inputReso4_1("/mnt/c/Users/pjska/github/FlowExtractor/out_sys/phi_v1_y_sys_vr_var1_iter3_.txt");
  std::ifstream inputReso4_2("/mnt/c/Users/pjska/github/FlowExtractor/out_sys/phi_v1_y_sys_vr_var2_iter3_.txt");
  Double_t d_v1_dv1dy_statErr4_1[8];
  Double_t d_v1_dv1dy_statErr4_2[8];
  // vr
  std::ifstream inputReso5_1("/mnt/c/Users/pjska/github/FlowExtractor/out_sys/phi_v1_y_sys_dedx_var1_iter3_.txt");
  std::ifstream inputReso5_2("/mnt/c/Users/pjska/github/FlowExtractor/out_sys/phi_v1_y_sys_dedx_var2_iter3_.txt");
  Double_t d_v1_dv1dy_statErr5_1[8];
  Double_t d_v1_dv1dy_statErr5_2[8];
  // dedx
  std::ifstream inputReso6_1("/mnt/c/Users/pjska/github/FlowExtractor/out_sys/phi_v1_y_sys_dca_var1_iter3_.txt");
  Double_t d_v1_dv1dy_statErr6_1[8];
  // dca
  std::ifstream inputReso7_1("/mnt/c/Users/pjska/github/FlowExtractor/out_sys/phi_v1_y_sys_nHitsFit_var1_iter3_.txt");
  std::ifstream inputReso7_2("/mnt/c/Users/pjska/github/FlowExtractor/out_sys/phi_v1_y_sys_nHitsFit_var2_iter3_.txt");
  Double_t d_v1_dv1dy_statErr7_1[8];
  Double_t d_v1_dv1dy_statErr7_2[8];
  // nHitsFit
  std::ifstream inputReso8_1("/mnt/c/Users/pjska/github/FlowExtractor/out_sys/phi_v1_y_sys_ratio_var1_iter3_.txt");
  std::ifstream inputReso8_2("/mnt/c/Users/pjska/github/FlowExtractor/out_sys/phi_v1_y_sys_ratio_var2_iter3_.txt");
  Double_t d_v1_dv1dy_statErr8_1[8];
  Double_t d_v1_dv1dy_statErr8_2[8];
  // ratio (nHitsFit/nHitsPoss)
  std::ifstream inputReso9_1("/mnt/c/Users/pjska/github/FlowExtractor/out_sys/phi_v1_y_sys_nSigK_var1_iter3_.txt");
  std::ifstream inputReso9_2("/mnt/c/Users/pjska/github/FlowExtractor/out_sys/phi_v1_y_sys_nSigK_var2_iter3_.txt");
  Double_t d_v1_dv1dy_statErr9_1[8];
  Double_t d_v1_dv1dy_statErr9_2[8];
  // nSigK
  std::ifstream inputReso10_1("/mnt/c/Users/pjska/github/FlowExtractor/out_sys/phi_v1_y_sys_mass2_var1_iter3_.txt");
  std::ifstream inputReso10_2("/mnt/c/Users/pjska/github/FlowExtractor/out_sys/phi_v1_y_sys_mass2_var2_iter3_.txt");
  Double_t d_v1_dv1dy_statErr10_1[8];
  Double_t d_v1_dv1dy_statErr10_2[8];
  // mass2
  std::ifstream inputReso11_1("/mnt/c/Users/pjska/github/FlowExtractor/out_sys/phi_v1_y_sys_pT_var1_iter3_.txt");
  std::ifstream inputReso11_2("/mnt/c/Users/pjska/github/FlowExtractor/out_sys/phi_v1_y_sys_pT_var2_iter3_.txt");
  Double_t d_v1_dv1dy_statErr11_1[8];
  Double_t d_v1_dv1dy_statErr11_2[8];
  // pT
  std::ifstream inputReso12_1("/mnt/c/Users/pjska/github/FlowExtractor/out_sys/phi_v1_y_sys_dipAngle_var1_iter3_.txt");
  Double_t d_v1_dv1dy_statErr12_1[8];
  // dipAngle
  std::ifstream inputReso13_1("/mnt/c/Users/pjska/github/FlowExtractor/out_sys/phi_v1_y_sys_vtxDiff_var1_iter3_.txt");
  std::ifstream inputReso13_2("/mnt/c/Users/pjska/github/FlowExtractor/out_sys/phi_v1_y_sys_vtxDiff_var2_iter3_.txt");
  Double_t d_v1_dv1dy_statErr13_1[8];
  Double_t d_v1_dv1dy_statErr13_2[8];
  // vtxdiff


  for(int i=0;i<8;i++){
    inputReso0   >> d_v1_dv1dy_statErr0[i];

    inputReso1_1 >> d_v1_dv1dy_statErr1_1[i];
    inputReso1_2 >> d_v1_dv1dy_statErr1_2[i];

    inputReso2_1 >> d_v1_dv1dy_statErr2_1[i];

    inputReso3_1 >> d_v1_dv1dy_statErr3_1[i];
    inputReso3_2 >> d_v1_dv1dy_statErr3_2[i];

    inputReso4_1 >> d_v1_dv1dy_statErr4_1[i];
    inputReso4_2 >> d_v1_dv1dy_statErr4_2[i];

    inputReso5_1 >> d_v1_dv1dy_statErr5_1[i];
    inputReso5_2 >> d_v1_dv1dy_statErr5_2[i];

    inputReso6_1 >> d_v1_dv1dy_statErr6_1[i];

    inputReso7_1 >> d_v1_dv1dy_statErr7_1[i];
    inputReso7_2 >> d_v1_dv1dy_statErr7_2[i];

    inputReso8_1 >> d_v1_dv1dy_statErr8_1[i];
    inputReso8_2 >> d_v1_dv1dy_statErr8_2[i];

    inputReso9_1 >> d_v1_dv1dy_statErr9_1[i];
    inputReso9_2 >> d_v1_dv1dy_statErr9_2[i];

    inputReso10_1 >> d_v1_dv1dy_statErr10_1[i];
    inputReso10_2 >> d_v1_dv1dy_statErr10_2[i];

    inputReso11_1 >> d_v1_dv1dy_statErr11_1[i];
    inputReso11_2 >> d_v1_dv1dy_statErr11_2[i];

    inputReso12_1 >> d_v1_dv1dy_statErr12_1[i];

    inputReso13_1 >> d_v1_dv1dy_statErr13_1[i];
    inputReso13_2 >> d_v1_dv1dy_statErr13_2[i];
  }


  // (2) =========== Processor:  TProfile/ TH1D to calculate Std Dev of each selection =======================
  TFile * outFile = new TFile("systematic_error_output.result.root","RECREATE" );
  TProfile  *TP_profile[4][_N_sysCut]; // 13 selections
  TProfile  *TP_profile2[_N_sysCut]; // 13 selections
  TProfile  *TP_profile3[_N_sysCut]; // 13 selections
  TProfile  *TP_profile4[_N_sysCut]; // 13 selections
  TH1D      *h_bin[4][_N_sysCut];

  for(int i = 0 ; i<_N_sysCut;i++){
    string  s_sysObj = sys_object[i+1];
    char * cstr = new char [s_sysObj.length()+1];
    std::strcpy (cstr, s_sysObj.c_str());

    for(int j = 0; j<3;j++){
      TP_profile[j][i]= new TProfile(Form("TProfile_v1_bin%d_%s",j+1,cstr),Form("Variations of v_{1} in y bin%d, %s cut", j+1 ,cstr),1,0,1);
      TP_profile[j][i]->GetXaxis()->SetTitle("v_{1}");
      h_bin[j][i]= new TH1D(Form("h_v1_bin%d_%s",j,cstr),Form("Variations of v_{1} in y bin%d, %s cut",j,cstr),10000000,-0.1,0.1);
      h_bin[j][i]->GetXaxis()->SetTitle("v_{1}");
    }
    TP_profile[3][i]= new TProfile(Form("TProfile_dv1dy_%s",cstr),Form("Variations of dv_{1}/dy, %s cut",cstr),1,0,1);
    TP_profile[3][i]->GetXaxis()->SetTitle("dv_{1}/dy");
    h_bin[3][i]= new TH1D(Form("h_dv1dy_%s",cstr),Form("Variations of dv_{1}/dy, %s cut",cstr),10000000,0.0,0.02);
    h_bin[3][i]->GetXaxis()->SetTitle("dv_{1}/dy");
  }

  // Fill default value to all the TProfile and TH1D
  for(int j =0 ; j<4; j++){
    for(int i = 0 ; i<_N_sysCut;i++){
      TP_profile[j][i]->Fill(0.5, d_v1_dv1dy_statErr0[j]);
      h_bin[j][i]->Fill(d_v1_dv1dy_statErr0[j]);
    }
  }

  // Fill variations to each TProfile and TH1D
  for(int j =0 ; j<4; j++){
    TP_profile[j][0]->Fill(0.5, d_v1_dv1dy_statErr1_1[j]);
    h_bin[j][0]->Fill(d_v1_dv1dy_statErr1_1[j]);
    TP_profile[j][0]->Fill(0.5, d_v1_dv1dy_statErr1_2[j]);
    h_bin[j][0]->Fill(d_v1_dv1dy_statErr1_2[j]);
    // etaGap
    TP_profile[j][1]->Fill(0.5, d_v1_dv1dy_statErr2_1[j]);
    h_bin[j][1]->Fill(d_v1_dv1dy_statErr2_1[j]);
    // etaRange
    TP_profile[j][2]->Fill(0.5, d_v1_dv1dy_statErr3_1[j]);
    h_bin[j][2]->Fill(d_v1_dv1dy_statErr3_1[j]);
    TP_profile[j][2]->Fill(0.5, d_v1_dv1dy_statErr3_2[j]);
    h_bin[j][2]->Fill(d_v1_dv1dy_statErr3_2[j]);
    // vz
    TP_profile[j][3]->Fill(0.5, d_v1_dv1dy_statErr4_1[j]);
    h_bin[j][3]->Fill(d_v1_dv1dy_statErr4_1[j]);
    TP_profile[j][3]->Fill(0.5, d_v1_dv1dy_statErr4_2[j]);
    h_bin[j][3]->Fill(d_v1_dv1dy_statErr4_2[j]);
    // vr
    TP_profile[j][4]->Fill(0.5, d_v1_dv1dy_statErr5_1[j]);
    h_bin[j][4]->Fill(d_v1_dv1dy_statErr5_1[j]);
    TP_profile[j][4]->Fill(0.5, d_v1_dv1dy_statErr5_2[j]);
    h_bin[j][4]->Fill(d_v1_dv1dy_statErr5_2[j]);
    // dedx
    TP_profile[j][5]->Fill(0.5, d_v1_dv1dy_statErr6_1[j]);
    h_bin[j][5]->Fill(d_v1_dv1dy_statErr6_1[j]);
    // dca
    TP_profile[j][6]->Fill(0.5, d_v1_dv1dy_statErr7_1[j]);
    h_bin[j][6]->Fill(d_v1_dv1dy_statErr7_1[j]);
    TP_profile[j][6]->Fill(0.5, d_v1_dv1dy_statErr7_2[j]);
    h_bin[j][6]->Fill(d_v1_dv1dy_statErr7_2[j]);
    // nHitsFit
    TP_profile[j][7]->Fill(0.5, d_v1_dv1dy_statErr8_1[j]);
    h_bin[j][7]->Fill(d_v1_dv1dy_statErr8_1[j]);
    TP_profile[j][7]->Fill(0.5, d_v1_dv1dy_statErr8_2[j]);
    h_bin[j][7]->Fill(d_v1_dv1dy_statErr8_2[j]);
    // ratio (nHitsFit/nHitsPoss)
    TP_profile[j][8]->Fill(0.5, d_v1_dv1dy_statErr9_1[j]);
    h_bin[j][8]->Fill(d_v1_dv1dy_statErr9_1[j]);
    TP_profile[j][8]->Fill(0.5, d_v1_dv1dy_statErr9_2[j]);
    h_bin[j][8]->Fill(d_v1_dv1dy_statErr9_2[j]);
    // nSigK
    TP_profile[j][9]->Fill(0.5, d_v1_dv1dy_statErr10_1[j]);
    h_bin[j][9]->Fill(d_v1_dv1dy_statErr10_1[j]);
    TP_profile[j][9]->Fill(0.5, d_v1_dv1dy_statErr10_2[j]);
    h_bin[j][9]->Fill(d_v1_dv1dy_statErr10_2[j]);
    // mass2
    TP_profile[j][10]->Fill(0.5, d_v1_dv1dy_statErr11_1[j]);
    h_bin[j][10]->Fill(d_v1_dv1dy_statErr11_1[j]);
    TP_profile[j][10]->Fill(0.5, d_v1_dv1dy_statErr11_2[j]);
    h_bin[j][10]->Fill(d_v1_dv1dy_statErr11_2[j]);
    // pT
    TP_profile[j][11]->Fill(0.5, d_v1_dv1dy_statErr12_1[j]);
    h_bin[j][11]->Fill(d_v1_dv1dy_statErr12_1[j]);
    // dipAngle
    TP_profile[j][12]->Fill(0.5, d_v1_dv1dy_statErr13_1[j]);
    h_bin[j][12]->Fill(d_v1_dv1dy_statErr13_1[j]);
    TP_profile[j][12]->Fill(0.5, d_v1_dv1dy_statErr13_2[j]);
    h_bin[j][12]->Fill(d_v1_dv1dy_statErr13_2[j]);
    // vtxdiff
  }
  // (4) =========== Output:  Sys Err of v1 and dv1/dy by Sqrt(Sum(Std Dev)) =======================
  TString outTxt = "systematic_error_output";
  outTxt.Append(".txt");
  std::ofstream sysFile(outTxt,ofstream::out);
  Double_t d_sys_err_bin[4] = {0.0,0.0,0.0,0.0};
  for(int j =0 ; j<4; j++){
    for(int i = 0 ; i<_N_sysCut;i++){
      Double_t d_stdDev = TP_profile[j][i]->GetStdDev(2);
      // Double_t d_stdDev = TP_profile[j][i]->GetBinError(1);
      cout << "bin " << j << ", sys Object: "<< sys_object[i+1] << ", Std Dev:" << d_stdDev << endl;
      sysFile << "bin " << j << ", sys Object: "<< sys_object[i+1] << ", Std Dev:" << d_stdDev << endl;
      d_sys_err_bin[j] += d_stdDev * d_stdDev;
    }
    cout << "bin " << j <<  ", Std Err:" << sqrt(d_sys_err_bin[j]) << endl;
    cout <<  " "<<endl;
    sysFile << "bin " << j <<  ", Std Err:" << sqrt(d_sys_err_bin[j]) << endl;
    sysFile <<  " "<<endl;
  }
// TPaveText on the final result
  // TPaveText * ptext1 = new TPaveText(0.2,0.65,0.40,0.9,"NDCARC");
  // Double_t d_mean_bin4 = TP_profile1        -> GetBinContent(1);
  // Double_t d_err_on_mean_bin4 = TP_profile1 -> GetBinError(1);
  // ptext1 -> AddText(Form("Mean: %.4f",d_mean_bin4));
  // ptext1 -> AddText(Form("Error on mean: %.4f",d_err_on_mean_bin4));
  //
  // TPaveText * ptext4 = new TPaveText(0.2,0.65,0.40,0.9,"NDCARC");
  // Double_t d_mean_slope = TP_profile4        -> GetBinContent(1);
  // Double_t d_err_on_mean_slope = TP_profile4 -> GetBinError(1);
  // ptext4 -> AddText(Form("Mean: %.4f",d_mean_slope));
  // ptext4 -> AddText(Form("Error on mean: %.4f",d_err_on_mean_slope));
  //
  // TCanvas * TC_inv1 = new TCanvas("TC_inv1","TC_inv1",1280,720);
  // TC_inv1->cd();
  // h_bin1->Draw();
  // ptext1->Draw("same");
  TCanvas * TC_inv1[_N_sysCut];
  Int_t i_bin = h_bin[3][0]->FindBin(0.0142375);
  for(int i = 0 ; i<_N_sysCut;i++){
    TC_inv1[i] = new TCanvas(Form("TC_inv%d",i),Form("TC_inv%d",i),800,600);
    TC_inv1[i]->cd();
    h_bin[3][i]->Draw();
    PaintBin (h_bin[3][i], i_bin, kRed);

    TC_inv1[i]->SaveAs(Form("/mnt/c/Users/pjska/Desktop/h_dv1dy_%d.png",i));
  }


  outFile->Write();
  sysFile.close();
}

void PaintBin (TH1D *h, Int_t bin, Int_t color)
{
   printf("%d %d %d\n", bin, color, 1 /*h->GetBinContent(bin)*/);
   TBox *b = new TBox(h->GetBinLowEdge(bin),
                      h->GetMinimum(),
                      h->GetBinWidth(bin)+h->GetBinLowEdge(bin),
                    1/*h->GetBinContent(bin)*/);
   b->SetFillColor(color);
   b->Draw();

}
