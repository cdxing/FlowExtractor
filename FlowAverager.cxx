#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TGraph.h"
#include "TMath.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TProfile3D.h"
/* *****************************************************************************
 * This code is used to Average flow in certain centrality bins from TProfile3D
 * TProfile3D here stores flow in centrality, pt and rapidity bin
 *
 * Author: Ding Chen
 * Date: Oct 9, 2020
 *
 * Update: Ding Chen
 * To Average v2 of KP, KM and Kaon pair backgrounds, in order to answer question
 * from Fuqiang
 * Date: April, 2021
 *******************************************************************************
*/
// Define Global contstant
const Int_t _Ncentralities = 9;
const Int_t rapidityBins   = 15;
const Double_t _y_CM = -2.03;

void FlowAverager(const char* inFile = "/mnt/c/Users/pjska/github/FlowExtractor/res_v2_7p2/pwg/merged_merged_sys_primary_var0_iter4_22F13B3805886BE920AC171E7672141C_NoKaonTpcEp.root")
{
  TFile* F = TFile::Open(inFile);
  if (!F || F->IsZombie()){
    std::cout<<"Cannot open the file"<<inFile<<std::endl;
    return;
  }
  Int_t ptBins = 15; Double_t ptLow = 0.0; Double_t ptHigh = 3.0;

  // TProfile3D *tp3d_proton_v1 = (TProfile3D*)F->Get("profile3D_Proton_v1");
  TProfile3D *tp3d_proton_v2 = (TProfile3D*)F->Get("profile3D_Phi_bkg_v2");
  // TProfile3D *tp3d_KaonPlus_v1 = (TProfile3D*)F->Get("profile3D_KP_v1");
  TProfile3D *tp3d_KaonPlus_v2 = (TProfile3D*)F->Get("profile3D_KP_v2");
  // TProfile3D *tp3d_KaonMinus_v1 = (TProfile3D*)F->Get("profile3D_KM_v1");
  TProfile3D *tp3d_KaonMinus_v2 = (TProfile3D*)F->Get("profile3D_KM_v2");
  // TProfile3D *tp3d_pionPlus = (TProfile3D*)F->Get("profile3D_pionPlus_v1");
  // TProfile3D *tp3d_pionMinus = (TProfile3D*)F->Get("profile3D_pionMinus_v1");
  // TProfile3D *tp3d_kaonPlus = (TProfile3D*)F->Get("profile3D_KP_v1");
  // TProfile3D *tp3d_kaonMinus = (TProfile3D*)F->Get("profile3D_KM_v1");

  // TProfile *profile_correlation_tpc_east_tpc_west_input = (TProfile*)F->Get("profile_correlation_tpc_east_tpc_west");
  // TProfile *profile_correlation_tpc_east_bbc_east_input = (TProfile*)F->Get("profile_correlation_tpc_east_bbc_east");
  // TProfile *profile_correlation_tpc_west_bbc_east_input = (TProfile*)F->Get("profile_correlation_tpc_west_bbc_east");

  // TH1D *h_correlation_tpc_east_tpc_west_input = new TH1D("h_correlation_tpc_east_tpc_west_input","<cos(#psi^{TPC east}_{1} #minus #psi^{TPC west}_{1})>",7,0.5,7+0.5);
  // TH1D *h_correlation_tpc_east_bbc_east_input = new TH1D("h_correlation_tpc_east_bbc_east_input","<cos(#psi^{TPC east}_{1} #minus #psi^{BBC east}_{1})>",7,0.5,7+0.5);
  // TH1D *h_correlation_tpc_west_bbc_east_input = new TH1D("h_correlation_tpc_west_bbc_east_input","<cos(#psi^{TPC west}_{1} #minus #psi^{BBC east}_{1})>",7,0.5,7+0.5);

  // for(int i = 1;i<=7;i++)
  // {
  //   h_correlation_tpc_east_tpc_west_input->SetBinContent(i,profile_correlation_tpc_east_tpc_west_input->GetBinContent(i));
  //   h_correlation_tpc_east_tpc_west_input->SetBinError(i,profile_correlation_tpc_east_tpc_west_input->GetBinError(i));
  //
  //   h_correlation_tpc_east_bbc_east_input->SetBinContent(i,profile_correlation_tpc_east_bbc_east_input->GetBinContent(i));
  //   h_correlation_tpc_east_bbc_east_input->SetBinError(i,profile_correlation_tpc_east_bbc_east_input->GetBinError(i));
  //
  //   h_correlation_tpc_west_bbc_east_input->SetBinContent(i,profile_correlation_tpc_west_bbc_east_input->GetBinContent(i));
  //   h_correlation_tpc_west_bbc_east_input->SetBinError(i,profile_correlation_tpc_west_bbc_east_input->GetBinError(i));
  //
  //
  // }

  // TProfile2D *tp2d_proton_v1_vs_y = tp3d_proton_v1->Project3DProfile("xz");
  // TProfile2D *tp2d_proton_v1_vs_pt = tp3d_proton_v1->Project3DProfile("xy");

  TProfile2D *tp2d_proton_v2_vs_y = tp3d_proton_v2->Project3DProfile("xz");
  TProfile2D *tp2d_proton_v2_vs_pt = tp3d_proton_v2->Project3DProfile("xy");

  TProfile2D *tp2d_KaonPlus_v2_vs_y = tp3d_KaonPlus_v2->Project3DProfile("xz");
  TProfile2D *tp2d_KaonPlus_v2_vs_pt = tp3d_KaonPlus_v2->Project3DProfile("xy");

  TProfile2D *tp2d_KaonMinus_v2_vs_y = tp3d_KaonMinus_v2->Project3DProfile("xz");
  TProfile2D *tp2d_KaonMinus_v2_vs_pt = tp3d_KaonMinus_v2->Project3DProfile("xy");
  // TProfile2D *tp2d_pionPlus = tp3d_pionPlus->Project3DProfile("xz");
  // TProfile2D *tp2d_pionMinus = tp3d_pionMinus->Project3DProfile("xz");
  // TProfile2D *tp2d_kaonPlus = tp3d_kaonPlus->Project3DProfile("xz");
  // TProfile2D *tp2d_kaonMinus = tp3d_kaonMinus->Project3DProfile("xz");

  // TProfile *profile_resolution2 = new TProfile("profile_resolution2","<cos(#psi^{TPC west}_{1} #minus #psi^{BBC east}_{1})>*<cos(#psi^{TPC east}_{1} #minus #psi^{BBC east}_{1})>/<cos(#psi^{TPC east}_{1} #minus #psi^{TPC west}_{1})>",7,0.5,7+0.5,-1.0,1.0,"");
  // profile_resolution2->GetXaxis()->SetTitle("Centrality (%)");
  // profile_resolution2->GetYaxis()->SetTitle("Resolution Squared");
  //
  // profile_resolution2->Multiply(profile_correlation_tpc_east_bbc_east_input,profile_correlation_tpc_west_bbc_east_input,1,1,"");
  // profile_resolution2->Divide(profile_correlation_tpc_east_tpc_west_input);

  // TH1D *h_resolution2 = new TH1D("h_resolution2","#frac{<cos(#psi^{TPC west}_{1} #minus #psi^{BBC east}_{1})> #times <cos(#psi^{TPC east}_{1} #minus #psi^{BBC east}_{1})>}{<cos(#psi^{TPC east}_{1} #minus #psi^{TPC west}_{1})>}",7,0.5,7+0.5);
  // h_resolution2->GetXaxis()->SetTitle("Centrality (%)");
  // h_resolution2->GetYaxis()->SetTitle("Resolution Squared");

  // TH1D *h_resolution = new TH1D("h_resolution","#sqrt{#frac{<cos(#psi^{TPC west}_{1} #minus #psi^{BBC east}_{1})> #times <cos(#psi^{TPC east}_{1} #minus #psi^{BBC east}_{1})>}{<cos(#psi^{TPC east}_{1} #minus #psi^{TPC west}_{1})>}}",7,0.5,7+0.5);
  // h_resolution->GetXaxis()->SetTitle("Centrality (%)");
  // h_resolution->GetYaxis()->SetTitle("Resolution");

  // h_resolution2->Multiply(h_correlation_tpc_east_bbc_east_input,h_correlation_tpc_west_bbc_east_input,1,1,"");
  // h_resolution2->Divide(h_correlation_tpc_east_tpc_west_input);



  // h_resolution2->Draw();
  TH1D * h_PRO_v1_vs_y_10_40_7p7 = new TH1D("h_PRO_v1_vs_y_10_40_7p7","h_PRO_v1_vs_y_10_40_7p7",15,-2.9-_y_CM,0.1-_y_CM);
  h_PRO_v1_vs_y_10_40_7p7->GetYaxis()->SetTitle("v_{1}");
  h_PRO_v1_vs_y_10_40_7p7->GetXaxis()->SetTitle("y-y_{CM}");
  //10-40%

  double Rapidity[8]= {-0.7,-0.5,-0.3,-0.1,0.1,0.3,0.5,0.7};
  double Rapidity_Error[8]={0.0,};
  double y2[8]={-0.0218351,-0.00944496,-0.00447499,-0.00132103,0.001304121,0.004387375,0.0093865102,0.02143059};
  double ey2[8]={0.0006705,0.0006228,0.0005701,0.0005629,0.0005672,0.000569,0.0005728, 0.0006806};
  TCanvas *TCanv_v1_vs_y = new TCanvas("TCanv_v1_vs_y","dv1/dy vs. ",200,10,1024,768);
  TCanv_v1_vs_y->DrawFrame(-2.15, -0.04, 2.15, 0.15);
  TCanv_v1_vs_y->cd();
  TGraphErrors *TGE_proton_v1_7p7 = new TGraphErrors(8, Rapidity, y2, Rapidity_Error, ey2);
  TGE_proton_v1_7p7->SetMarkerStyle(8);
  TGE_proton_v1_7p7->SetMarkerColor(kRed);
  TGE_proton_v1_7p7->SetLineColor(kRed);
  TGE_proton_v1_7p7->SetMarkerSize(1);
  TGE_proton_v1_7p7->Draw("P");
  TH1D * h_PRO_v1_vs_y_10_40 = new TH1D("h_PRO_v1_vs_y_10_40","h_PRO_v1_vs_y_10_40",15,-2.9-_y_CM,0.1-_y_CM);
  h_PRO_v1_vs_y_10_40->GetYaxis()->SetTitle("v_{1}");
  h_PRO_v1_vs_y_10_40->GetXaxis()->SetTitle("y-y_{CM}");

  Int_t i_cent_section[10] = {0, 5, 10, 20, 30, 40, 50 ,60, 70, 80};
  Int_t i_cent_section_BES_II[4] = {0, 10, 40, 80};

  TH1D * h_proton_v2_vs_pt_7p2[9];
  TH1D * h_proton_v2_vs_pt_7p2_BES_II[3], * h_KaonPlus_v2_vs_pt_7p2_BES_II[3], * h_KaonMinus_v2_vs_pt_7p2_BES_II[3];
  for(int cent = 0; cent<9;cent++){
    h_proton_v2_vs_pt_7p2[cent] = new TH1D(Form("h_proton_v2_vs_pt_7p2_cent%d",cent),
    Form("Au+Au 7.2 GeV FXT, %d - %d %%, #phi^{bkg} M_{inv}: [1.04, 1.09]",i_cent_section[cent],i_cent_section[cent+1]),
    ptBins,ptLow,ptHigh);
    h_proton_v2_vs_pt_7p2[cent]->GetYaxis()->SetTitle("v_{2}");
    h_proton_v2_vs_pt_7p2[cent]->GetYaxis()->SetTitleOffset(1);
    h_proton_v2_vs_pt_7p2[cent]->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  }
  for(int cent=0;cent <3;cent++){
    h_proton_v2_vs_pt_7p2_BES_II[cent] = new TH1D(Form("h_proton_v2_vs_pt_7p2_BES_II_cent%d",cent),
    Form("Au+Au 7.2 GeV FXT, %d - %d %%, #phi^{bkg} M_{inv}: [1.04, 1.09]",i_cent_section_BES_II[cent],i_cent_section_BES_II[cent+1]),
    ptBins,ptLow,ptHigh);
    h_proton_v2_vs_pt_7p2_BES_II[cent]->GetYaxis()->SetTitle("v_{2}");
    h_proton_v2_vs_pt_7p2_BES_II[cent]->GetYaxis()->SetTitleOffset(1);
    h_proton_v2_vs_pt_7p2_BES_II[cent]->GetXaxis()->SetTitle("p_{T} [GeV/c]");
    h_KaonPlus_v2_vs_pt_7p2_BES_II[cent] = new TH1D(Form("h_KaonPlus_v2_vs_pt_7p2_BES_II_cent%d",cent),
    Form("Au+Au 7.2 GeV FXT, %d - %d %%, K^{+}",i_cent_section_BES_II[cent],i_cent_section_BES_II[cent+1]),
    ptBins,ptLow,ptHigh);
    h_KaonPlus_v2_vs_pt_7p2_BES_II[cent]->GetYaxis()->SetTitle("v_{2}");
    h_KaonPlus_v2_vs_pt_7p2_BES_II[cent]->GetYaxis()->SetTitleOffset(1);
    h_KaonPlus_v2_vs_pt_7p2_BES_II[cent]->GetXaxis()->SetTitle("p_{T} [GeV/c]");
    h_KaonMinus_v2_vs_pt_7p2_BES_II[cent] = new TH1D(Form("h_KaonMinus_v2_vs_pt_7p2_BES_II_cent%d",cent),
    Form("Au+Au 7.2 GeV FXT, %d - %d %%, K^{-}",i_cent_section_BES_II[cent],i_cent_section_BES_II[cent+1]),
    ptBins,ptLow,ptHigh);
    h_KaonMinus_v2_vs_pt_7p2_BES_II[cent]->GetYaxis()->SetTitle("v_{2}");
    h_KaonMinus_v2_vs_pt_7p2_BES_II[cent]->GetYaxis()->SetTitleOffset(1);
    h_KaonMinus_v2_vs_pt_7p2_BES_II[cent]->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  }
  // TH1D * h_pionPlus_v1_vs_y_10_30 = new TH1D("h_pionPlus_v1_vs_y_10_30","h_pionPlus_v1_vs_y_10_30",15,-2.9,0.1);
  // TH1D * h_pionMinus_v1_vs_y_10_30 = new TH1D("h_pionMinus_v1_vs_y_10_30","h_pionMinus_v1_vs_y_10_30",15,-2.9,0.1);
  // TH1D * h_kaonPlus_v1_vs_y_10_30 = new TH1D("h_kaonPlus_v1_vs_y_10_30","h_kaonPlus_v1_vs_y_10_30",15,-2.9,0.1);
  // TH1D * h_kaonMinus_v1_vs_y_10_30 = new TH1D("h_kaonMinus_v1_vs_y_10_30","h_kaonMinus_v1_vs_y_10_30",15,-2.9,0.1);


  double  v_1_vs_y_proton[_Ncentralities][rapidityBins] = {0};
  double  v_1_vs_y_proton_err[_Ncentralities][rapidityBins] = {0};
  double  v_1_vs_y_proton_avrg[rapidityBins]  = {0};
  double  v_1_vs_y_proton_wt[rapidityBins]     = {0};
  double  v_1_vs_y_proton_fnl[rapidityBins] = {0};
  double  v_1_vs_y_proton_errbar[rapidityBins] = {0};

  // Proton
  double  v_2_vs_pt_proton[9][15] = {0};
  double  v_2_vs_pt_proton_err[9][15] = {0};

  double  v_2_vs_pt_proton_avrg[9][15]  = {0};
  double  v_2_vs_pt_proton_wt[9][15]     = {0};
  double  v_2_vs_pt_proton_fnl[9][15] = {0};
  double  v_2_vs_pt_proton_errbar[9][15] = {0};
  // 9 centrality bins, 15 pT bins

  double  v_2_vs_pt_proton_avrg_BES_II[3][15]  = {0};
  double  v_2_vs_pt_proton_wt_BES_II[3][15]     = {0};
  double  v_2_vs_pt_proton_fnl_BES_II[3][15] = {0};
  double  v_2_vs_pt_proton_errbar_BES_II[3][15] = {0};
  // 3 centrality bins, 0-10%, 10 - 40%, 40 -80%,  15 pT bins

  // KaonPlus
  double  v_2_vs_pt_KaonPlus[9][15] = {0};
  double  v_2_vs_pt_KaonPlus_err[9][15] = {0};

  double  v_2_vs_pt_KaonPlus_avrg[9][15]  = {0};
  double  v_2_vs_pt_KaonPlus_wt[9][15]     = {0};
  double  v_2_vs_pt_KaonPlus_fnl[9][15] = {0};
  double  v_2_vs_pt_KaonPlus_errbar[9][15] = {0};
  // 9 centrality bins, 15 pT bins

  double  v_2_vs_pt_KaonPlus_avrg_BES_II[3][15]  = {0};
  double  v_2_vs_pt_KaonPlus_wt_BES_II[3][15]     = {0};
  double  v_2_vs_pt_KaonPlus_fnl_BES_II[3][15] = {0};
  double  v_2_vs_pt_KaonPlus_errbar_BES_II[3][15] = {0};
  // 3 centrality bins, 0-10%, 10 - 40%, 40 -80%,  15 pT bins

  // KaonMinus
  double  v_2_vs_pt_KaonMinus[9][15] = {0};
  double  v_2_vs_pt_KaonMinus_err[9][15] = {0};

  double  v_2_vs_pt_KaonMinus_avrg[9][15]  = {0};
  double  v_2_vs_pt_KaonMinus_wt[9][15]     = {0};
  double  v_2_vs_pt_KaonMinus_fnl[9][15] = {0};
  double  v_2_vs_pt_KaonMinus_errbar[9][15] = {0};
  // 9 centrality bins, 15 pT bins

  double  v_2_vs_pt_KaonMinus_avrg_BES_II[3][15]  = {0};
  double  v_2_vs_pt_KaonMinus_wt_BES_II[3][15]     = {0};
  double  v_2_vs_pt_KaonMinus_fnl_BES_II[3][15] = {0};
  double  v_2_vs_pt_KaonMinus_errbar_BES_II[3][15] = {0};
  // 3 centrality bins, 0-10%, 10 - 40%, 40 -80%,  15 pT bins

  double  value4[_Ncentralities][rapidityBins] = {0};
  double  error4[_Ncentralities][rapidityBins] = {0};
  double  value4_avg[rapidityBins]  = {0};
  double  weight4[rapidityBins]     = {0};
  double  value4_fnal[rapidityBins] = {0};
  double  errorbar4[rapidityBins] = {0};

  double  value5[_Ncentralities][rapidityBins] = {0};
  double  error5[_Ncentralities][rapidityBins] = {0};
  double  value5_avg[rapidityBins]  = {0};
  double  weight5[rapidityBins]     = {0};
  double  value5_fnal[rapidityBins] = {0};
  double  errorbar5[rapidityBins] = {0};

  double  d_resolution2[_Ncentralities] ={0.137326,0.220693,0.24983,0.208653,0.132604,1,1,1,1};
  double  d_resolution2_err[_Ncentralities] ={0};

  double  d_resolution[_Ncentralities] ={0};
  double  d_resolution_err[_Ncentralities] ={0};

  // for(int i=0;i<7;i++)
  // {
  //   d_resolution2[i] = h_resolution2->GetBinContent(i+1);
  //   d_resolution2_err[i] = h_resolution2->GetBinError(i+1);
  //   if(d_resolution2[i]>0)
  //   {
  //     d_resolution[i] = TMath::Sqrt(d_resolution2[i]);
  //     d_resolution_err[i] = 0.5 * d_resolution2_err[i]/TMath::Sqrt(d_resolution2[i]);
  //     h_resolution->SetBinContent(i+1,d_resolution[i]);
  //     h_resolution->SetBinError(i+1,d_resolution_err[i]);
  //   }
  //   std::cout<< Form("resolution %d = ",i+1) << d_resolution[i]<<std::endl;
  //
  // }
  // h_resolution->Draw();

  d_resolution[0] = 0.314819;
  d_resolution[1] = 0.380233;
  d_resolution[2] = 0.43834;
  d_resolution[3] = 0.473238;
  d_resolution[4] = 0.492292;
  d_resolution[5] = 0.415132;
  d_resolution[6] = 0.329119;
  d_resolution[7] = 0.34041;
  d_resolution[8] = 0.438046;



  for(int i=0;i<_Ncentralities;i++)
  {
    /*
    for(int j=0;j<rapidityBins;j++) // vs. y
    {
      double d_v1_vs_y_raw_proton           = tp2d_proton_v1_vs_y->GetBinContent(j+1,i+1);
      double d_v1_vs_y_raw_proton_err       = tp2d_proton_v1_vs_y->GetBinError(j+1,i+1);
      // double d_v1_raw_pionPlus      = tp2d_pionPlus->GetBinContent(j+1,i+3);
      // double d_v1_raw_pionPlus_err  = tp2d_pionPlus->GetBinError(j+1,i+3);
      // double d_v1_raw_pionMinus     = tp2d_pionMinus->GetBinContent(j+1,i+3);
      // double d_v1_raw_pionMinus_err = tp2d_pionMinus->GetBinError(j+1,i+3);
      // double d_v1_raw_kaonPlus      = tp2d_kaonPlus->GetBinContent(j+1,i+3);
      // double d_v1_raw_kaonPlus_err  = tp2d_kaonPlus->GetBinError(j+1,i+3);
      // double d_v1_raw_kaonMinus     = tp2d_kaonMinus->GetBinContent(j+1,i+3);
      // double d_v1_raw_kaonMinus_err = tp2d_kaonMinus->GetBinError(j+1,i+3);

      // v_1_vs_y_proton[i][j] = d_v1_vs_y_raw_proton;
      // v_1_vs_y_proton_err[i][j] = d_v1_vs_y_raw_proton_err;
      v_1_vs_y_proton[i][j] = d_v1_vs_y_raw_proton/d_resolution[i];
      v_1_vs_y_proton_err[i][j] = d_v1_vs_y_raw_proton_err/d_resolution[i];
      // v_1_vs_y_proton_err[i][j] = TMath::Sqrt((d_v1_vs_y_raw_proton_err/d_resolution[i])*(d_v1_vs_y_raw_proton_err/d_resolution[i])
      //              + (0.5*d_v1_vs_y_raw_proton*d_resolution2_err[i]/(pow(d_resolution[i],3)))*(0.5*d_v1_vs_y_raw_proton*d_resolution2_err[i]/(pow(d_resolution[i],3))));

      // value2[i][j] = d_v1_raw_pionPlus/d_resolution[i+2];
      // error2[i][j] = TMath::Sqrt((d_v1_raw_pionPlus_err/d_resolution[i+2])*(d_v1_raw_pionPlus_err/d_resolution[i+2])
      //              + (0.5*d_v1_raw_pionPlus*d_resolution2_err[i+2]/(pow(d_resolution[i+2],3)))*(0.5*d_v1_raw_pionPlus*d_resolution2_err[i+2]/(pow(d_resolution[i+2],3))));
      //
      // value3[i][j] = d_v1_raw_pionMinus/d_resolution[i+2];
      // error3[i][j] = TMath::Sqrt((d_v1_raw_pionMinus_err/d_resolution[i+2])*(d_v1_raw_pionMinus_err/d_resolution[i+2])
      //              + (0.5*d_v1_raw_pionMinus*d_resolution2_err[i+2]/(pow(d_resolution[i+2],3)))*(0.5*d_v1_raw_pionMinus*d_resolution2_err[i+2]/(pow(d_resolution[i+2],3))));
      //
      // value4[i][j] = d_v1_raw_kaonPlus/d_resolution[i+2];
      // error4[i][j] = TMath::Sqrt((d_v1_raw_kaonPlus_err/d_resolution[i+2])*(d_v1_raw_kaonPlus_err/d_resolution[i+2])
      //              + (0.5*d_v1_raw_kaonPlus*d_resolution2_err[i+2]/(pow(d_resolution[i+2],3)))*(0.5*d_v1_raw_kaonPlus*d_resolution2_err[i+2]/(pow(d_resolution[i+2],3))));
      //
      // value5[i][j] = d_v1_raw_kaonMinus/d_resolution[i+2];
      // error5[i][j] = TMath::Sqrt((d_v1_raw_kaonMinus_err/d_resolution[i+2])*(d_v1_raw_kaonMinus_err/d_resolution[i+2])
      //              + (0.5*d_v1_raw_kaonMinus*d_resolution2_err[i+2]/(pow(d_resolution[i+2],3)))*(0.5*d_v1_raw_kaonMinus*d_resolution2_err[i+2]/(pow(d_resolution[i+2],3))));

    }
    */
    for(int k=0;k<ptBins;k++){ // vs. pt
      double d_v2_vs_pt_raw_proton           = tp2d_proton_v2_vs_pt->GetBinContent(k+1,i+1);
      double d_v2_vs_pt_raw_proton_err       = tp2d_proton_v2_vs_pt->GetBinError(k+1,i+1);
      v_2_vs_pt_proton[i][k] = d_v2_vs_pt_raw_proton/d_resolution2[i];
      v_2_vs_pt_proton_err[i][k] = d_v2_vs_pt_raw_proton_err/d_resolution2[i];
      cout << "v2 phi bkg: "<< v_2_vs_pt_proton[i][k]<<endl;
      double d_v2_vs_pt_raw_KaonPlus           = tp2d_KaonPlus_v2_vs_pt->GetBinContent(k+1,i+1);
      double d_v2_vs_pt_raw_KaonPlus_err       = tp2d_KaonPlus_v2_vs_pt->GetBinError(k+1,i+1);
      double d_v2_vs_pt_raw_KaonMinus           = tp2d_KaonMinus_v2_vs_pt->GetBinContent(k+1,i+1);
      double d_v2_vs_pt_raw_KaonMinus_err       = tp2d_KaonMinus_v2_vs_pt->GetBinError(k+1,i+1);

      v_2_vs_pt_KaonPlus[i][k] = d_v2_vs_pt_raw_KaonPlus/d_resolution2[i];
      v_2_vs_pt_KaonPlus_err[i][k] = d_v2_vs_pt_raw_KaonPlus_err/d_resolution2[i];
      cout << "v2 KP : "<< v_2_vs_pt_KaonPlus[i][k]<<endl;
      v_2_vs_pt_KaonMinus[i][k] = d_v2_vs_pt_raw_KaonMinus/d_resolution2[i];
      v_2_vs_pt_KaonMinus_err[i][k] = d_v2_vs_pt_raw_KaonMinus_err/d_resolution2[i];
      cout << "v2 KM : "<< v_2_vs_pt_KaonMinus[i][k]<<endl;
    }
  }
  /*
  for(int j=0;j<rapidityBins;j++) // vs. y
  {
    for(int i=2;i<5;i++) // centrality = 2, 10-20%; centrality = 5, 40-50%
    {
      if(v_1_vs_y_proton_err[i][j] != 0)
      {
        v_1_vs_y_proton_avrg[j] += (v_1_vs_y_proton[i][j]) / ((v_1_vs_y_proton_err[i][j]) * (v_1_vs_y_proton_err[i][j]));
        v_1_vs_y_proton_wt[j]    += 1/((v_1_vs_y_proton_err[i][j]) * (v_1_vs_y_proton_err[i][j]));
      }

      // if(error2[i][j] != 0)
      // {
      //   value2_avg[j] += (value2[i][j]) / ((error2[i][j]) * (error2[i][j]));
      //   weight2[j]    += 1/((error2[i][j]) * (error2[i][j]));
      // }
      //
      // if(error3[i][j] != 0)
      // {
      //   value3_avg[j] += (value3[i][j]) / ((error3[i][j]) * (error3[i][j]));
      //   weight3[j]    += 1/((error3[i][j]) * (error3[i][j]));
      // }
      //
      // if(error4[i][j] != 0)
      // {
      //   value4_avg[j] += (value4[i][j]) / ((error4[i][j]) * (error4[i][j]));
      //   weight4[j]    += 1/((error4[i][j]) * (error4[i][j]));
      // }
      //
      // if(error5[i][j] != 0)
      // {
      //   value5_avg[j] += (value5[i][j]) / ((error5[i][j]) * (error5[i][j]));
      //   weight5[j]    += 1/((error5[i][j]) * (error5[i][j]));
      // }

    }
    if(v_1_vs_y_proton_wt[j]!=0)
    {
      v_1_vs_y_proton_fnl[j] = v_1_vs_y_proton_avrg[j]/v_1_vs_y_proton_wt[j];
      v_1_vs_y_proton_errbar[j] = TMath::Sqrt(1/v_1_vs_y_proton_wt[j]);
    }

    // if(weight4[j]!=0)
    // {
    //   value4_fnal[j] = value4_avg[j]/weight4[j];
    //   errorbar4[j] = TMath::Sqrt(1/weight4[j]);
    // }
    //
    // if(weight5[j]!=0)
    // {
    //   value5_fnal[j] = value5_avg[j]/weight5[j];
    //   errorbar5[j] = TMath::Sqrt(1/weight5[j]);
    // }

  }
  */

  for(int k=0;k<ptBins;k++){ // vs. pt

    for(int i = 0; i <2;i++){ // 0 - 10 %
      if(v_2_vs_pt_proton_err[i][k] != 0)
      {
        v_2_vs_pt_proton_avrg_BES_II[0][k] += (v_2_vs_pt_proton[i][k]) / ((v_2_vs_pt_proton_err[i][k]) * (v_2_vs_pt_proton_err[i][k]));
        v_2_vs_pt_proton_wt_BES_II[0][k]    += 1/((v_2_vs_pt_proton_err[i][k]) * (v_2_vs_pt_proton_err[i][k]));
      }
      if(v_2_vs_pt_KaonPlus_err[i][k] != 0)
      {
        v_2_vs_pt_KaonPlus_avrg_BES_II[0][k] += (v_2_vs_pt_KaonPlus[i][k]) / ((v_2_vs_pt_KaonPlus_err[i][k]) * (v_2_vs_pt_KaonPlus_err[i][k]));
        v_2_vs_pt_KaonPlus_wt_BES_II[0][k]    += 1/((v_2_vs_pt_KaonPlus_err[i][k]) * (v_2_vs_pt_KaonPlus_err[i][k]));
      }
      if(v_2_vs_pt_KaonMinus_err[i][k] != 0)
      {
        v_2_vs_pt_KaonMinus_avrg_BES_II[0][k] += (v_2_vs_pt_KaonMinus[i][k]) / ((v_2_vs_pt_KaonMinus_err[i][k]) * (v_2_vs_pt_KaonMinus_err[i][k]));
        v_2_vs_pt_KaonMinus_wt_BES_II[0][k]    += 1/((v_2_vs_pt_KaonMinus_err[i][k]) * (v_2_vs_pt_KaonMinus_err[i][k]));
      }
    }

    for(int i = 2; i <5;i++){ // 10 - 40 %
      if(v_2_vs_pt_proton_err[i][k] != 0)
      {
        v_2_vs_pt_proton_avrg_BES_II[1][k] += (v_2_vs_pt_proton[i][k]) / ((v_2_vs_pt_proton_err[i][k]) * (v_2_vs_pt_proton_err[i][k]));
        v_2_vs_pt_proton_wt_BES_II[1][k]    += 1/((v_2_vs_pt_proton_err[i][k]) * (v_2_vs_pt_proton_err[i][k]));
      }
      if(v_2_vs_pt_KaonPlus_err[i][k] != 0)
      {
        v_2_vs_pt_KaonPlus_avrg_BES_II[1][k] += (v_2_vs_pt_KaonPlus[i][k]) / ((v_2_vs_pt_KaonPlus_err[i][k]) * (v_2_vs_pt_KaonPlus_err[i][k]));
        v_2_vs_pt_KaonPlus_wt_BES_II[1][k]    += 1/((v_2_vs_pt_KaonPlus_err[i][k]) * (v_2_vs_pt_KaonPlus_err[i][k]));
      }
      if(v_2_vs_pt_KaonMinus_err[i][k] != 0)
      {
        v_2_vs_pt_KaonMinus_avrg_BES_II[1][k] += (v_2_vs_pt_KaonMinus[i][k]) / ((v_2_vs_pt_KaonMinus_err[i][k]) * (v_2_vs_pt_KaonMinus_err[i][k]));
        v_2_vs_pt_KaonMinus_wt_BES_II[1][k]    += 1/((v_2_vs_pt_KaonMinus_err[i][k]) * (v_2_vs_pt_KaonMinus_err[i][k]));
      }
    }

    for(int i = 5; i <9;i++){ // 40 - 80 %
      if(v_2_vs_pt_proton_err[i][k] != 0)
      {
        v_2_vs_pt_proton_avrg_BES_II[2][k] += (v_2_vs_pt_proton[i][k]) / ((v_2_vs_pt_proton_err[i][k]) * (v_2_vs_pt_proton_err[i][k]));
        v_2_vs_pt_proton_wt_BES_II[2][k]    += 1/((v_2_vs_pt_proton_err[i][k]) * (v_2_vs_pt_proton_err[i][k]));
      }
      if(v_2_vs_pt_KaonPlus_err[i][k] != 0)
      {
        v_2_vs_pt_KaonPlus_avrg_BES_II[2][k] += (v_2_vs_pt_KaonPlus[i][k]) / ((v_2_vs_pt_KaonPlus_err[i][k]) * (v_2_vs_pt_KaonPlus_err[i][k]));
        v_2_vs_pt_KaonPlus_wt_BES_II[2][k]    += 1/((v_2_vs_pt_KaonPlus_err[i][k]) * (v_2_vs_pt_KaonPlus_err[i][k]));
      }
      if(v_2_vs_pt_KaonMinus_err[i][k] != 0)
      {
        v_2_vs_pt_KaonMinus_avrg_BES_II[2][k] += (v_2_vs_pt_KaonMinus[i][k]) / ((v_2_vs_pt_KaonMinus_err[i][k]) * (v_2_vs_pt_KaonMinus_err[i][k]));
        v_2_vs_pt_KaonMinus_wt_BES_II[2][k]    += 1/((v_2_vs_pt_KaonMinus_err[i][k]) * (v_2_vs_pt_KaonMinus_err[i][k]));
      }
    }

    for(int i = 0; i< 3; i++ ){ // 0 - 10%, 10 - 40 %, 40 - 80%
      if(v_2_vs_pt_proton_wt_BES_II[i][k]!=0){
        v_2_vs_pt_proton_fnl_BES_II[i][k] = v_2_vs_pt_proton_avrg_BES_II[i][k]/v_2_vs_pt_proton_wt_BES_II[i][k];
        v_2_vs_pt_proton_errbar_BES_II[i][k] = TMath::Sqrt(1/v_2_vs_pt_proton_wt_BES_II[i][k]);
      }
      if(v_2_vs_pt_KaonPlus_wt_BES_II[i][k]!=0){
        v_2_vs_pt_KaonPlus_fnl_BES_II[i][k] = v_2_vs_pt_KaonPlus_avrg_BES_II[i][k]/v_2_vs_pt_KaonPlus_wt_BES_II[i][k];
        v_2_vs_pt_KaonPlus_errbar_BES_II[i][k] = TMath::Sqrt(1/v_2_vs_pt_KaonPlus_wt_BES_II[i][k]);
      }
      if(v_2_vs_pt_KaonMinus_wt_BES_II[i][k]!=0){
        v_2_vs_pt_KaonMinus_fnl_BES_II[i][k] = v_2_vs_pt_KaonMinus_avrg_BES_II[i][k]/v_2_vs_pt_KaonMinus_wt_BES_II[i][k];
        v_2_vs_pt_KaonMinus_errbar_BES_II[i][k] = TMath::Sqrt(1/v_2_vs_pt_KaonMinus_wt_BES_II[i][k]);
      }
    }

    for(int i = 0; i <9;i++){ // each and every centralties
      if(v_2_vs_pt_proton_err[i][k] != 0)
      {
        v_2_vs_pt_proton_avrg[i][k] += (v_2_vs_pt_proton[i][k]) / ((v_2_vs_pt_proton_err[i][k]) * (v_2_vs_pt_proton_err[i][k]));
        v_2_vs_pt_proton_wt[i][k]    += 1/((v_2_vs_pt_proton_err[i][k]) * (v_2_vs_pt_proton_err[i][k]));
      }
    }
    for(int i = 0; i <9;i++){ // each and every centralties
      if(v_2_vs_pt_proton_wt[i][k]!=0)
      {
        v_2_vs_pt_proton_fnl[i][k] = v_2_vs_pt_proton_avrg[i][k]/v_2_vs_pt_proton_wt[i][k];
        v_2_vs_pt_proton_errbar[i][k] = TMath::Sqrt(1/v_2_vs_pt_proton_wt[i][k]);
      }
    }
  }

  /*
  for(int j=0;j<rapidityBins;j++) // vs. y
  {
    h_PRO_v1_vs_y_10_40->SetBinContent((j+1),v_1_vs_y_proton_fnl[j]);
    h_PRO_v1_vs_y_10_40->SetBinError((j+1),v_1_vs_y_proton_errbar[j]);

    // h_pionPlus_v1_vs_y_10_30->SetBinContent((j+1),value2_fnal[j]);
    // h_pionPlus_v1_vs_y_10_30->SetBinError((j+1),errorbar2[j]);
    //
    // h_pionMinus_v1_vs_y_10_30->SetBinContent((j+1),value3_fnal[j]);
    // h_pionMinus_v1_vs_y_10_30->SetBinError((j+1),errorbar3[j]);
    //
    // h_kaonPlus_v1_vs_y_10_30->SetBinContent((j+1),value4_fnal[j]);
    // h_kaonPlus_v1_vs_y_10_30->SetBinError((j+1),errorbar4[j]);
    //
    // h_kaonMinus_v1_vs_y_10_30->SetBinContent((j+1),value5_fnal[j]);
    // h_kaonMinus_v1_vs_y_10_30->SetBinError((j+1),errorbar5[j]);
    //
  }
  */

  for(int k=0;k<ptBins;k++) // vs. pt Fill histograms
  {
    for(int cent=0;cent<9;cent++){
      h_proton_v2_vs_pt_7p2[cent]->SetBinContent((k+1),v_2_vs_pt_proton_fnl[cent][k]);
      h_proton_v2_vs_pt_7p2[cent]->SetBinError((k+1),v_2_vs_pt_proton_errbar[cent][k]);
    }
    for(int centBES = 0; centBES<3;centBES++){
      h_proton_v2_vs_pt_7p2_BES_II[centBES]->SetBinContent((k+1),v_2_vs_pt_proton_fnl_BES_II[centBES][k]);
      h_proton_v2_vs_pt_7p2_BES_II[centBES]->SetBinError((k+1),v_2_vs_pt_proton_errbar_BES_II[centBES][k]);

      h_KaonPlus_v2_vs_pt_7p2_BES_II[centBES]->SetBinContent((k+1),v_2_vs_pt_KaonPlus_fnl_BES_II[centBES][k]);
      h_KaonPlus_v2_vs_pt_7p2_BES_II[centBES]->SetBinError((k+1),v_2_vs_pt_KaonPlus_errbar_BES_II[centBES][k]);

      h_KaonMinus_v2_vs_pt_7p2_BES_II[centBES]->SetBinContent((k+1),v_2_vs_pt_KaonMinus_fnl_BES_II[centBES][k]);
      h_KaonMinus_v2_vs_pt_7p2_BES_II[centBES]->SetBinError((k+1),v_2_vs_pt_KaonMinus_errbar_BES_II[centBES][k]);
    }
  }
  h_PRO_v1_vs_y_10_40->Draw("same");
  TLegend *legend1 = new TLegend(0.1,0.65,0.6,0.9);
  legend1->AddEntry(h_PRO_v1_vs_y_10_40,"7.2 GeV FXT 10-40% Proton v_{1}","lep");
  legend1->AddEntry(TGE_proton_v1_7p7,"7.7 GeV COL 10-40% Proton v_{1}","lep");
  legend1->Draw("same");
  TCanvas *TCanv_v2_vs_pt = new TCanvas("TCanv_v2_vs_pt","v_{2} vs. p_{T}",200,10,1024,768);
  // TCanv_v2_vs_pt->DrawFrame(0, -0.15, 3, 0.15);
  TCanv_v2_vs_pt->Divide(3,3);
  TLine *tl_0 = new TLine(0,0,2,0);
  tl_0->SetLineStyle(2);
  for(int cent=0;cent<9;cent++){
    TCanv_v2_vs_pt->cd(cent+1);
    h_proton_v2_vs_pt_7p2[cent]->GetXaxis()->SetRangeUser(0,2.0);
    h_proton_v2_vs_pt_7p2[cent]->Draw();
    tl_0->Draw("same");
  }

  TCanvas *TCanv_v2_vs_pt_BES_II = new TCanvas("TCanv_v2_vs_pt_BES_II","v_{2} vs. p_{T}",200,10,1024,768);
  TCanv_v2_vs_pt_BES_II->Divide(3,3);
  for(int cent=0;cent<3;cent++){
    TCanv_v2_vs_pt_BES_II->cd(cent+1);
    h_proton_v2_vs_pt_7p2_BES_II[cent]->GetXaxis()->SetRangeUser(0,2.0);
    // h_proton_v2_vs_pt_7p2_BES_II[cent]->GetYaxis()->SetRangeUser(-0.001,0.002);
    h_proton_v2_vs_pt_7p2_BES_II[cent]->GetYaxis()->SetRangeUser(-0.04,0.15);
    h_proton_v2_vs_pt_7p2_BES_II[cent]->Draw();
    tl_0->Draw("same");
    TCanv_v2_vs_pt_BES_II->cd(cent+4);
    h_KaonPlus_v2_vs_pt_7p2_BES_II[cent]->GetXaxis()->SetRangeUser(0,2.0);
    // h_KaonPlus_v2_vs_pt_7p2_BES_II[cent]->GetYaxis()->SetRangeUser(-0.002,0.008);
    h_KaonPlus_v2_vs_pt_7p2_BES_II[cent]->GetYaxis()->SetRangeUser(-0.04,0.15);
    h_KaonPlus_v2_vs_pt_7p2_BES_II[cent]->Draw();
    tl_0->Draw("same");
    TCanv_v2_vs_pt_BES_II->cd(cent+7);
    h_KaonMinus_v2_vs_pt_7p2_BES_II[cent]->GetXaxis()->SetRangeUser(0,2.0);
    // h_KaonMinus_v2_vs_pt_7p2_BES_II[cent]->GetYaxis()->SetRangeUser(-0.002,0.01);
    h_KaonMinus_v2_vs_pt_7p2_BES_II[cent]->GetYaxis()->SetRangeUser(-0.04,0.15);
    h_KaonMinus_v2_vs_pt_7p2_BES_II[cent]->Draw();
    tl_0->Draw("same");
  }

  TFile* o1 = new TFile("averageflow_pwg_bkg.root","RECREATE");
  h_PRO_v1_vs_y_10_40->Write();
  for(int cent=0;cent<9;cent++){
    h_proton_v2_vs_pt_7p2[cent]->Write();
  }
  for(int cent=0;cent<3;cent++){
    h_proton_v2_vs_pt_7p2_BES_II[cent]->Write();
    h_KaonPlus_v2_vs_pt_7p2_BES_II[cent]->Write();
    h_KaonMinus_v2_vs_pt_7p2_BES_II[cent]->Write();
  }
  // h_pionPlus_v1_vs_y_10_30->Write();
  // h_pionMinus_v1_vs_y_10_30->Write();
  // h_kaonPlus_v1_vs_y_10_30->Write();
  // h_kaonMinus_v1_vs_y_10_30->Write();
  // profile_resolution2->Write();
  // profile_correlation_tpc_east_tpc_west_input->Write();
  // profile_correlation_tpc_east_bbc_east_input->Write();
  // profile_correlation_tpc_west_bbc_east_input->Write();
  // h_correlation_tpc_east_tpc_west_input->Write();
  // h_correlation_tpc_east_bbc_east_input->Write();
  // h_correlation_tpc_west_bbc_east_input->Write();
  // h_resolution2->Write();
  // h_resolution->Write();


  return;

}
