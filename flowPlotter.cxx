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

Double_t proportion(Double_t *x, Double_t *p);
Double_t Pol3plus1(Double_t *x, Double_t *p);
Double_t d_v2nq(Double_t d_v2, Int_t ncq );
Double_t d_mTm0nq(Double_t d_pT, Int_t ncq, Double_t mass );

void flowPlotter(){
  // ========================= (1) dv_{1}/dy of phi-meson =========================
  TCanvas *TCanv_dv1dy = new TCanvas("TCanv_dv1dy","dv_{1}/dy vs. ",200,10,1024,768);
  TCanv_dv1dy->DrawFrame(2., -0.06, 250., 0.2);
  // 4.5 GeV phi
  double x_phi_4p5[1]    = {4.5};
  double xErr_phi_4p5[1] = {0};
  double y_phi_4p5[1]      = {-0.01301};
  double yErr_stat_phi_4p5[1] = {0.01566};
         yErr_stat_phi_4p5[0] *= sqrt(2); // 2 points are flipped, only two real data points for fitting, stat err * sqrt(2)
  double yErr_sys_phi_4p5[1]  = {0.013012};
         yErr_sys_phi_4p5[0]  *= sqrt(2);

  // BES-I phi
  double x_phi_BES[6]    = { 11.5, 14.5, 19.6, 27, 39, 200};
  double xErr_phi_BES[6] = {0};
  double y_phi_BES[6]      = {0.0199, -0.0311, -0.0304, -0.02,   -0.013,   -0.00366944};
  double yErr_stat_phi_BES[6] = {0.0183, 0.00938, 0.00763, 0.00644, 0.006142, 0.00203012};
  double yErr_sys_phi_BES[6]  = {0.03,   0.021,   0.014,   0.015,   0.012,    0.000533};

  // 3 GeV from Guannan Xie
  double x_phi_3[1]    = {3};
  double xErr_phi_3[1] = {0};
  double y_phi_3[1]      = {0.1205};
  double yErr_stat_phi_3[1] = {0.0755};
  double yErr_sys_phi_3[1]  = {0.01345};
  // 7.7 GeV BES-I data, credit Guannan
  double x_phi_7p7[1]    = {7.7};
  double xErr_phi_7p7[1] = {0};
  double y_phi_7p7[1]      = {0.0153485};//{0.0159611};
  double yErr_stat_phi_7p7[1] = {0.0391923};//{0.0403535};
  // 7.2 GeV
  double x_phi_7p2[1]    = {7.2};
  double xErr_phi_7p2[1] = {0};
  double y_phi_7p2[1]      = {0.0138904};
  double yErr_stat_phi_7p2[1] = {0.00916586};
  double yErr_sys_phi_7p2[1] = {0.00647123};
  // data set (1) with stat and sys errors

  // Now draw data set (1)
  // ------------- We first have to draw it only with the stat errors ------------
  gStyle->SetOptStat(0);
  gStyle->SetOptDate(0);
  gStyle->SetEndErrorSize(6);
  gStyle->SetOptTitle(0);
  TCanv_dv1dy->SetGridx(0);
  TCanv_dv1dy->SetGridy(0);
  TCanv_dv1dy->SetLogx();
  TH2D * histTemp = new TH2D("histTemp","histTemp",1000,2.,250,1000,-0.06,0.2);
  histTemp->GetYaxis()->SetTitle("Directed Flow Slope dv_{1} /dy|_{y=0}");
  histTemp->GetYaxis()->SetTitleOffset(1.2);
  histTemp->GetXaxis()->SetTitleOffset(1.3);
  histTemp->GetXaxis()->SetTitle("Collision Energy #sqrt{s_{NN}} (GeV)");
  // histTemp->SetStats(FALSE);
  histTemp->Draw();
  //4.5 GeV phi
  TGraphErrors *graph_4p5 = new TGraphErrors(1, x_phi_4p5, y_phi_4p5, xErr_phi_4p5, yErr_stat_phi_4p5);
  graph_4p5->SetMarkerStyle(34);
  graph_4p5->SetMarkerColor(kBlack);
  graph_4p5->SetLineColor(kBlack);
  graph_4p5->SetMarkerSize(2);
  graph_4p5->Draw("sameP");
  // BES-I phi
  TGraphErrors *graph_BES = new TGraphErrors(6, x_phi_BES, y_phi_BES, xErr_phi_BES, yErr_stat_phi_BES);
  graph_BES->SetMarkerStyle(28);
  graph_BES->SetMarkerColor(kBlack);
  graph_BES->SetLineColor(kBlack);
  graph_BES->SetMarkerSize(2);
  graph_BES->Draw("P");
  // 3 GeV from Guannan Xie
  TGraphErrors *graph_3 = new TGraphErrors(1, x_phi_3, y_phi_3, xErr_phi_3, yErr_stat_phi_3);
  graph_3->SetMarkerStyle(20);
  graph_3->SetMarkerColor(kRed);
  graph_3->SetLineColor(kRed);
  graph_3->SetMarkerSize(2);
  graph_3->Draw("sameP");
  // 7.7 GeV BES-I data, credit Guannan
  TGraphErrors *graph_7p7 = new TGraphErrors(1, x_phi_7p7, y_phi_7p7, xErr_phi_7p7, yErr_stat_phi_7p7);
  graph_7p7->SetMarkerStyle(28);
  graph_7p7->SetMarkerColor(kBlack);
  graph_7p7->SetLineColor(kBlack);
  graph_7p7->SetMarkerSize(2);
  graph_7p7->Draw("P");
  // 7.2 GeV
  TGraphErrors *graph_7p2 = new TGraphErrors(1, x_phi_7p2, y_phi_7p2, xErr_phi_7p2, yErr_stat_phi_7p2);
  graph_7p2->SetMarkerStyle(29);
  graph_7p2->SetMarkerColor(kBlue);
  graph_7p2->SetLineColor(kBlue);
  graph_7p2->SetMarkerSize(2);
  graph_7p2->Draw("P");
  TGraphErrors *graph_7p2_sys = new TGraphErrors(1, x_phi_7p2, y_phi_7p2, xErr_phi_7p2, yErr_sys_phi_7p2);
  graph_7p2_sys->SetMarkerColor(kBlue);
  graph_7p2_sys->SetLineColor(kBlue);
  graph_7p2_sys->Draw("[]");
  // TLgend
  TLegend *legend = new TLegend(0.3,0.65,0.75,0.85);
  legend->SetBorderSize(0);
  legend->AddEntry((TObject*)0, "#phi-meson", "");
  legend->AddEntry(graph_3,"3 GeV   FXT 10-40% ","p");
  legend->AddEntry(graph_4p5,"4.5 GeV FXT 0-30% ","p");
  legend->AddEntry(graph_7p2,"7.2 GeV FXT 10-40% ","p");
  // legend->AddEntry(graph_7p7,"7.7 GeV 10-40%","p");
  legend->AddEntry(graph_BES,"BES-I 10-40%","p");
  legend->Draw("same");
  TPaveText * ptxt_3_200 = new TPaveText(0.65,0.7,0.88,0.85,"NDCARC");
  ptxt_3_200->SetFillColor(0);
  ptxt_3_200 -> AddText("Au+Au 3 - 7.2 GeV FXT");
  ptxt_3_200 -> AddText("STAR preliminary");
  ptxt_3_200->Draw("same");
  // ---------------- Now we have to somehow depict the sys errors ---------------
  //4.5 GeV phi
  TGraphErrors *graph_4p5_sys = new TGraphErrors(1, x_phi_4p5, y_phi_4p5, xErr_phi_4p5, yErr_sys_phi_4p5);
  graph_4p5_sys->SetLineColor(kBlack);
  graph_4p5_sys->Draw("[]");
  // BES-I phi
  TGraphErrors *graph_BES_sys = new TGraphErrors(6, x_phi_BES, y_phi_BES, xErr_phi_BES, yErr_sys_phi_BES);
  graph_BES_sys->SetLineColor(kBlack);
  graph_BES_sys->Draw("[]");
  // 3 GeV from Guannan Xie
  TGraphErrors *graph_3_sys = new TGraphErrors(1, x_phi_3, y_phi_3, xErr_phi_3, yErr_sys_phi_3);
  graph_3_sys->SetLineColor(kRed);
  graph_3_sys->Draw("[]");
  TLine *line = new TLine(2, 0, 250, 0);
  line->SetLineStyle(7);
  line->Draw("same");
  // ========================== (2) v1 vs. y at 7.2 GeV ========================
  TCanvas *TCanv_v1_vs_y = new TCanvas("TCanv_v1_vs_y","v_{1} vs. y at 7.2 GeV",200,10,1024,768);
  TCanv_v1_vs_y->DrawFrame(0, -0.04, 2.1, 0.1);
  // 7.2 GeV v1 vs. y
  double x_rap[3]    = { -1.25-_y_CM, -0.75-_y_CM, -0.25-_y_CM};
  double x_rap_flip[3]    = { -(-1.25-_y_CM), -(-0.75-_y_CM), -(-0.25-_y_CM)};
  double x_rap_jkk[3]    = { -1.25-_y_CM+0.02, -0.75-_y_CM+0.02, -0.25-_y_CM+0.02};
  double xErr_rap_newProd[3] = {0.25, 0.25, 0.25};
  double xErr_rap[3] = {0};
  double y_v1[3]      = {0.00678684, 0.0191681, 0.0203576};
  double y_v1_newProd[3] = {5.3034e-05, 0.0294281, 0.085239};
  double y_v1_newProd_flip[3] = {-5.3034e-05, -0.0294281, -0.085239};
  double y_v1_newProd_jkk[3] = {5.794e-05, 0.02952, 0.08509};
  double y_0[3]      = {0,0,0};
  double yErr_stat_v1[3] = {0.0220404, 0.0129847, 0.0461443};
  double yErr_stat_v1_newProd[3] = { 0.0255991, 0.0138418, 0.0488292};
  double yErr_stat_v1_newProd_jkk[3] = { 0.0077*3, 0.005338*3, 0.02125*3};
  double yErr_sys_v1[3] = {0.00881589, 0.00917299, 0.0286432};
  TCanv_v1_vs_y->cd();
  TH2D * histTemp1 = new TH2D("histTemp1","histTemp1",1000,0,2.1,1000,-0.04,0.15);
  histTemp1->GetYaxis()->SetTitle("v_{1}");
  histTemp1->GetYaxis()->SetTitleOffset(1);
  histTemp1->GetXaxis()->SetTitle("y (COM frame)");
  histTemp1->Draw();
  // 7.2 GeV
  TGraphErrors *graph_v1_vs_y_7p2 = new TGraphErrors(3, x_rap, y_v1, xErr_rap, yErr_stat_v1);
  graph_v1_vs_y_7p2->SetMarkerStyle(30);
  graph_v1_vs_y_7p2->SetMarkerColor(kBlue);
  graph_v1_vs_y_7p2->SetLineColor(kBlue);
  graph_v1_vs_y_7p2->SetMarkerSize(2);
  graph_v1_vs_y_7p2->Draw("P");
  TGraphErrors *graph_v1_vs_y_7p2_sys = new TGraphErrors(3, x_rap, y_v1, 0, yErr_sys_v1);
  graph_v1_vs_y_7p2_sys->SetMarkerColor(kBlue);
  // graph_v1_vs_y_7p2_sys->SetFillColor(4);
  // graph_v1_vs_y_7p2_sys->SetFillStyle(3001);
  graph_v1_vs_y_7p2_sys->SetLineColor(kBlue);
  graph_v1_vs_y_7p2_sys->Draw("[]");
  TLine *line1 = new TLine(0, 0, 2.1, 0);
  line1->SetLineStyle(7);
  line1->Draw("same");
  TLine *line2 = new TLine(-_y_CM, -0.04, -_y_CM, 0.1);
  line2->SetLineStyle(7);
  // line2->SetTitle("y_{target}")
  // line2->Draw("same");
  TF1 * tf1_dv1dy = new TF1("tf1_dv1dy",proportion,0.,2.,1);
  graph_v1_vs_y_7p2->Fit(tf1_dv1dy,"E+","R",-1.42-_y_CM,/*2.*/0.-_y_CM);
  TPaveText * ptext = new TPaveText(0.12,0.72,0.4,0.85,"NDCARC");
  ptext->SetFillColor(0);
  ptext -> AddText("Au+Au 7.2 GeV FXT 10-40%");
  ptext -> AddText("y_{target} = 2.02");
  ptext -> AddText("STAR preliminary");
  TPaveText * ptxt_y_target = new TPaveText(0.82,0.84,0.87,0.88,"NDCARC");
  ptxt_y_target->SetFillColor(0);
  ptxt_y_target -> AddText("y_{target} = 2.02");
  // ptext -> AddText(Form("dv_{1}/dy: %.4f %c %.4f %c %.4f",0.0138904/*(Double_t)tf1_dv1dy->GetParameter(0)*/,177,0.00916586/*(Double_t)tf1_dv1dy->GetParError(0)*/,177,0.00647123));
  // ptext -> AddText(Form("#chi^{2}/ndf = %.2f / %d",(Double_t)tf1_dv1dy->GetChisquare(),(Int_t)tf1_dv1dy->GetNDF()));
  TLegend *leg_7p2_v1_phi = new TLegend(0.25,0.5,0.75,0.8);
  leg_7p2_v1_phi->SetFillStyle(0);
  leg_7p2_v1_phi->SetBorderSize(0);
  leg_7p2_v1_phi->AddEntry((TObject*)0,Form("dv_{1} /dy|_{y=0} = %.4f %c %.4f %c %.4f",(Double_t)y_phi_7p2[0]/*(Double_t)tf1_dv1dy->GetParameter(0)*/,177,(Double_t)yErr_stat_phi_7p2[0]/*(Double_t)tf1_dv1dy->GetParError(0)*/,177,0.00647123),"");
  leg_7p2_v1_phi->Draw("same");
  ptext->Draw("same");
  // ptxt_y_target->Draw("same");
  // ========================== (2.1) v1 vs. y at 7.2 GeV  7.7 GeV ========================
  double yy[4] = { -0.45, -0.15, 0.15, 0.45};
  double yyerr[4] = {0.15, 0.15, 0.15, 0.15};
  double v11[4] = { -0.0173913, 0.00217391, -0.00797101, 0.0057971};
  double v11err[4] = {0.0275362, 0.0210145, 0.0195652 , 0.0271739};

  TGraphErrors* gr7p7Cen_v11_slop;
  // gr7p7Cen_v11_slop = new TGraphErrors(4, yy, v11, yyerr, v11err);
  gr7p7Cen_v11_slop = new TGraphErrors(4, yy, v11, 0, v11err);
  TCanvas *TCanv_7p2_7p7 = new TCanvas("TCanv_7p2_7p7","7.2 7.7 combined ",200,10,1024,768);
  // TCanv_7p2_7p7->DrawFrame(-1., -0.05, 2.1, 0.15);
  TCanv_7p2_7p7->DrawFrame(-2.1, -0.15, 2.1, 0.15);
  TCanv_7p2_7p7->cd();
  // TH2D * histTemp_7p2_7p7 = new TH2D("histTemp_7p2_7p7","histTemp_7p2_7p7",1000,-1,2.1,1000,-0.05,0.15);
  TH2D * histTemp_7p2_7p7 = new TH2D("histTemp_7p2_7p7","histTemp_7p2_7p7",1000,-2.1,2.1,1000,-0.15,0.15);
  histTemp_7p2_7p7->GetYaxis()->SetTitle("v_{1}");
  histTemp_7p2_7p7->GetYaxis()->SetTitleOffset(1);
  histTemp_7p2_7p7->GetXaxis()->SetTitle("y");
  histTemp_7p2_7p7->Draw();
  TGraphErrors *graph_v1_vs_y_7p2_newProd_flip = new TGraphErrors(3, x_rap_flip, y_v1_newProd_flip, xErr_rap_newProd, yErr_stat_v1_newProd);
  graph_v1_vs_y_7p2_newProd_flip->SetMarkerStyle(26);
  graph_v1_vs_y_7p2_newProd_flip->SetMarkerColor(kRed);
  graph_v1_vs_y_7p2_newProd_flip->SetLineColor(kRed);
  graph_v1_vs_y_7p2_newProd_flip->SetMarkerSize(2);
  // graph_v1_vs_y_7p2_newProd_flip->SetMarkerStyle(73);
  graph_v1_vs_y_7p2_newProd_flip->Draw("P");
  TGraphErrors *graph_v1_vs_y_7p2_newProd = new TGraphErrors(3, x_rap, y_v1_newProd, xErr_rap_newProd, yErr_stat_v1_newProd);
  graph_v1_vs_y_7p2_newProd->SetMarkerStyle(23);
  graph_v1_vs_y_7p2_newProd->SetMarkerColor(kRed);
  graph_v1_vs_y_7p2_newProd->SetLineColor(kRed);
  graph_v1_vs_y_7p2_newProd->SetMarkerSize(2);
  graph_v1_vs_y_7p2_newProd->Draw("P");

  TGraphErrors *graph_v1_vs_y_7p2_newProd_jkk = new TGraphErrors(3, x_rap_jkk, y_v1_newProd_jkk, /*xErr_rap_newProd*/0, yErr_stat_v1_newProd_jkk);
  graph_v1_vs_y_7p2_newProd_jkk->SetMarkerStyle(32);
  graph_v1_vs_y_7p2_newProd_jkk->SetMarkerColor(kRed+2);
  graph_v1_vs_y_7p2_newProd_jkk->SetLineColor(kRed+2);
  graph_v1_vs_y_7p2_newProd_jkk->SetMarkerSize(2);
  graph_v1_vs_y_7p2_newProd_jkk->Draw("P");
  gr7p7Cen_v11_slop->SetMarkerStyle(20);
  gr7p7Cen_v11_slop->SetMarkerSize(2);
  gr7p7Cen_v11_slop->SetMarkerColor(kGreen+3);
  gr7p7Cen_v11_slop->SetLineWidth(2);
  gr7p7Cen_v11_slop->SetLineColor(kGreen+3);
  gr7p7Cen_v11_slop->Draw("Psame");
  graph_v1_vs_y_7p2->Draw("Psame");
  double yy_7p2_7p7[7] = { -0.45, -0.15, 0.15, 0.45, -1.25-_y_CM, -0.75-_y_CM, -0.25-_y_CM};
  double yyerr_7p2_7p7[7] = {0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15};
  double v11_7p2_7p7[7] = { -0.0173913, 0.00217391, -0.00797101, 0.0057971, 0.00678684, 0.0191681, 0.0203576};
  double v11err_7p2_7p7[7] = {0.0275362, 0.0210145, 0.0195652 , 0.0271739, 0.0220404, 0.0129847, 0.0461443};
  double v11_7p2_7p7_newProd[7] = { -0.0173913, 0.00217391, -0.00797101, 0.0057971, 5.3034e-05, 0.0294281, 0.085239};
  double v11err_7p2_7p7_newProd[7] = {0.0275362, 0.0210145, 0.0195652 , 0.0271739, 0.0255991, 0.0138418, 0.0488292};
  TGraphErrors *graph_v1_vs_y_7p2_7p7 = new TGraphErrors(7, yy_7p2_7p7, v11_7p2_7p7, 0, v11err_7p2_7p7);
  graph_v1_vs_y_7p2_7p7->SetMarkerStyle(31);
  graph_v1_vs_y_7p2_7p7->SetMarkerColor(kBlack);
  graph_v1_vs_y_7p2_7p7->SetLineColor(kBlack);
  graph_v1_vs_y_7p2_7p7->SetMarkerSize(2);
  // graph_v1_vs_y_7p2_7p7->Draw("Psame");
  TGraphErrors *graph_v1_vs_y_7p2_7p7_newProd = new TGraphErrors(7, yy_7p2_7p7, v11_7p2_7p7_newProd, 0, v11err_7p2_7p7_newProd);
  graph_v1_vs_y_7p2_7p7_newProd->SetMarkerStyle(31);
  graph_v1_vs_y_7p2_7p7_newProd->SetMarkerColor(kBlack);
  graph_v1_vs_y_7p2_7p7_newProd->SetLineColor(kBlack);
  graph_v1_vs_y_7p2_7p7_newProd->SetMarkerSize(2);
  // graph_v1_vs_y_7p2_7p7_newProd->Draw("Psame");
  TF1 * tf1_dv1dy_7p2_7p7 = new TF1("tf1_dv1dy_7p2_7p7",proportion,0.,2.,1);
  graph_v1_vs_y_7p2_7p7->Fit(tf1_dv1dy_7p2_7p7,"E+","R",-0.6,/*2.*/0.-_y_CM);
  graph_v1_vs_y_7p2_7p7->GetFunction("tf1_dv1dy_7p2_7p7")->Draw("same");
  TF1 * tf1_dv1dy_7p2_7p7_pol3 = new TF1("tf1_dv1dy_7p2_7p7_pol3",Pol3plus1,0.,2.,2);
  graph_v1_vs_y_7p2_7p7_newProd->Fit(tf1_dv1dy_7p2_7p7_pol3,"E+","R",-0.6,/*2.*/0.-_y_CM);
  graph_v1_vs_y_7p2_7p7_newProd->GetFunction("tf1_dv1dy_7p2_7p7_pol3")->SetLineColor(kBlue);
  graph_v1_vs_y_7p2_7p7_newProd->GetFunction("tf1_dv1dy_7p2_7p7_pol3")->Draw("same");
  TPaveText * ptext_7p2_7p7 = new TPaveText(0.12,0.75,0.35,0.85,"NDCARC");
  ptext_7p2_7p7->SetFillColor(0);
  ptext_7p2_7p7 -> AddText("Au+Au 7.2 GeV FXT");
  ptext_7p2_7p7 -> AddText("STAR preliminary");
  // ptext_7p2_7p7 -> AddText(Form("dv_{1}/dy: %.4f %c %.4f",(Double_t)tf1_dv1dy_7p2_7p7->GetParameter(0),177,(Double_t)tf1_dv1dy_7p2_7p7->GetParError(0)));
  // ptext_7p2_7p7 -> AddText(Form("#chi^{2}/ndf = %.2f / %d",(Double_t)tf1_dv1dy_7p2_7p7->GetChisquare(),(Int_t)tf1_dv1dy_7p2_7p7->GetNDF()));
  ptext_7p2_7p7->Draw("same");
  TPaveText * ptext_7p2_7p7_target = new TPaveText(0.78,0.12,0.88,0.22,"NDCARC");
  ptext_7p2_7p7_target->SetFillColor(0);
  ptext_7p2_7p7_target -> AddText("y_{target} = 2.02");
  ptext_7p2_7p7_target->Draw("same");
  graph_v1_vs_y_7p2_sys->Draw("[]");
  TLine *TL_7p2_7p7 = new TLine(-2.1, 0, 2.1, 0);
  TL_7p2_7p7->SetLineStyle(7);
  TL_7p2_7p7->Draw("same");
  TLine *TL_7p2_7p7_1 = new TLine(0, -0.15, 0, 0.15);
  TL_7p2_7p7_1->SetLineStyle(7);
  TL_7p2_7p7_1->Draw("same");
  TLegend *legend_7p2_7p7 = new TLegend(0.4,0.6,0.8,0.8);
  legend_7p2_7p7->SetFillStyle(0);
  legend_7p2_7p7->SetBorderSize(0);
  legend_7p2_7p7->AddEntry(gr7p7Cen_v11_slop,"#phi 7.7 GeV COL 10-40% ","p");
  legend_7p2_7p7->AddEntry(graph_v1_vs_y_7p2,"#phi 7.2 GeV FXT 10-40% ","p");
  legend_7p2_7p7->AddEntry(graph_v1_vs_y_7p2_newProd,"#phi 7.2 GeV FXT 10-40% new Prod ","p");
  legend_7p2_7p7->AddEntry(graph_v1_vs_y_7p2_newProd_jkk,"#phi 7.2 GeV FXT 10-40% new Prod - Jackknife ","p");
  legend_7p2_7p7->AddEntry(graph_v1_vs_y_7p2_newProd_flip,"#phi 7.2 GeV FXT 10-40% new Prod - flipped ","p");
  legend_7p2_7p7->AddEntry((TObject*)0,Form("dv_{1} /dy|_{y=0} = %.4f %c %.4f, linear",(Double_t)tf1_dv1dy_7p2_7p7->GetParameter(0),177,(Double_t)tf1_dv1dy_7p2_7p7->GetParError(0)),"");
  legend_7p2_7p7->AddEntry((TObject*)0,Form("p1 = %.4f %c %.4f, newProd, Pol3",(Double_t)tf1_dv1dy_7p2_7p7_pol3->GetParameter(0),177,(Double_t)tf1_dv1dy_7p2_7p7_pol3->GetParError(0)),"");
  legend_7p2_7p7->AddEntry((TObject*)0,Form("p3 = %.4f %c %.4f, newProd, Pol3",(Double_t)tf1_dv1dy_7p2_7p7_pol3->GetParameter(1),177,(Double_t)tf1_dv1dy_7p2_7p7_pol3->GetParError(1)),"");
  legend_7p2_7p7->Draw("same");
  // ========================== (3) v2 vs. pT at 7.2 GeV =======================
  TCanvas *TCanv_v2_vs_pt_10_40 = new TCanvas("TCanv_v2_vs_pt_10_40","v2 vs. pT 10-40%",200,10,1024,768);
  TCanv_v2_vs_pt_10_40->DrawFrame(0, -0.1, 3, 0.2);
  // 7.2 GeV v2 vs. pT 10-40%
  double x_ptSetA[2]    = { 0.9,1.8};
  double x_ptSetA_jkk[2]    = { 0.9+0.02,1.8+0.02};
  double xErr_ptSetA[2] = {0.3,0.6};
  double x_ptSetB_3bin[3]    = { 0.85,1.55,2.5};
  double x_ptSetB_3bin_jkk[3]    = { 0.85+0.02,1.55+0.02,2.5+0.02};
  double xErr_ptSetB_3bin[3] = {0.25,0.45,0.5};
  // double y_v2_10_40[2]      = {0.0256657, 0.0120596}; // EPD
  // double yErr_stat_v2_10_40[2] = {0.023679, 0.0240889};

  // double y_v2_10_40[2]      = {-0.0149938, -0.0150853}; // w bin
  // double yErr_stat_v2_10_40[2] = {0.0277206, 0.028806};
  double y_v2_10_40[2]      = {-0.0097, -0.0205}; // w bin
  double yErr_stat_v2_10_40[2] = {0.0281, 0.0277};

  double y_v2_10_40_jkk[2]      = {-0.00966106, -0.0205263}; // w bin
  double yErr_stat_v2_10_40_jkk[2] = {0.014085270558323343*3, 0.009257925680356804*3};

  // double y_v2_10_40_3bin[3]      = { -0.0285367, -0.0089344, 0.112995}; // w bin
  // double yErr_stat_v2_10_40_3bin[3] = { 0.0322494, 0.0270015, 0.16299};
  // double y_v2_10_40_3bin_jkk[3]      = {-0.02838, -0.008862, 0.1135}; // w bin
  // double yErr_stat_v2_10_40_3bin_jkk[3] = {0.008845*3, 0.008841*3, 0.07765*3};

  double y_v2_10_40_3bin[3]      = { -0.0276968, -0.00838587, 0.0991779}; // w bin
  double yErr_stat_v2_10_40_3bin[3] = { 0.031929, 0.0260998, 0.150601};
  double y_v2_10_40_3bin_jkk[3]      = {-0.027653, -0.0083664, 0.0985826}; // w bin
  double yErr_stat_v2_10_40_3bin_jkk[3] = {0.0333665, 0.0272353, 0.13186};
  // double yErr_stat_v2_10_40_3bin_jkk[3] = {0.009892, 0.00487,0.01565};
  // v2 10-40% ptSetB_centSetA: -0.0276872, -0.0100532, 0.0504745
  // v2 Err 10-40% ptSetB_centSetA: 0.0310942, 0.0257785, 0.178318
  double x_ptSetB[4]    = { 0.55,0.85,1.2,1.7};
  double xErr_ptSetB[4] = {0.15,0.15,0.2,0.3};

  // double y_v2_10_40[4]      = { 0.0120157, -0.0221078, -0.0647598, 0.0635793};
  // double yErr_stat_v2_10_40[4] = {0.0487792, 0.0416835, 0.0301251, 0.0388423};

  // 7.2 GeV v2 vs. pT 10-40% Shaowei
  double y_v2_10_40_Shaowei[2]      = {-0.0191553, -0.030874};
  double yErr_stat_v2_10_40_Shaowei[2] = {0.0193706, 0.0224251};
  // 7.2 GeV v2 vs. pT 10-40% Guannan etasub
  float pt_7p2Cen_10_40[2] = {0.85, 1.85};
  float pterr_7p2Cen_10_40[2]={0.35,0.65};
  float v2_7p2Cen_10_40[2] = {0.0683646 , -0.000300613 };
  float v2stat_7p2Cen_10_40[2] = {0.0584518, 0.0710356};
  // 7.2 GeV v2 vs. pT 10-40% Guannan invM
  float pt_7p2Cen_10_40_invM[2] = {0.85 + 0.05, 1.85 + 0.05};
  float pterr_7p2Cen_10_40_invM[2]={0.35,0.65};
  float v2_7p2Cen_10_40_invM[2] = {0.0704069 , -0.00410226 };
  float v2stat_7p2Cen_10_40_invM[2] = {0.0365758, 0.0502626};
  // centrality: 10-40%
  // beam energy: 7.7 GeV
  // Event plane method: eta-sub
  // ---------------- Particle species: Phi ----------------
  Double_t pt_bin_center_10_40[1] = {0.9945};
  Double_t pt_bin_center_errCen_10_40[1] = {0};
  Double_t v2_values_10_40[1] = {0.0399507};
  Double_t v2_stat_error_10_40[1] = {0.0193815};
  // ****************** beam energy: 11.5 GeV ******************
  // Event plane method: eta-sub
  // ---------------- Particle species: Phi ----------------
  Double_t pt_bin_center_11p5[3] = {0.7705,1.2325,1.8925};
  Double_t v2_values_11p5[3] = {0.0357841,0.0761367,0.0778874};
  Double_t v2_stat_error_11p5[3] = {0.0120707,0.0127943,0.0261273};
  Double_t v2_syst_low_error_11p5[3] = {0.00683794,0.00634618,0.00700603};
  Double_t v2_syst_high_error_11p5[3] = {0.00540595,0.00642469,0.00644049};
  // ****************** beam energy: 19.6 GeV ******************
  // Event plane method: eta-sub
  // ---------------- Particle species: Phi ----------------
  Double_t pt_bin_center_19p6[5] = {0.6515,0.9985,1.3735,1.7775,2.2645};
  Double_t v2_values_19p6[5] = {0.014229,0.0481282,0.0734469,0.0872473,0.113327};
  Double_t v2_stat_error_19p6[5] = {0.00862292,0.00630574,0.00811914,0.0124024,0.017897};
  Double_t v2_syst_low_error_19p6[5] = {0.00529471,0.00438684,0.00180629,0.00541846,0.00201318};
  Double_t v2_syst_high_error_19p6[5] = {0.00512403,0.00249747,0.00339875,0.00440116,0.00222428};
  // ****************** beam energy: 27 GeV ******************
  // Event plane method: eta-sub
  // ---------------- Particle species: Phi ----------------
  Double_t pt_bin_center_27[5] = {0.6505,0.9995,1.3765,1.8445,2.4725};
  Double_t v2_values_27[5] = {0.0201674,0.0486692,0.0869156,0.0978724,0.116611};
  Double_t v2_stat_error_27[5] = {0.00582714,0.00414427,0.00513633,0.00669162,0.0122933};
  Double_t v2_syst_low_error_27[5] = {0.00238713,0.00415245,0.00439713,0.00501242,0.00545236};
  Double_t v2_syst_high_error_27[5] = {0.00204073,0.00189032,0.00548138,0.00402443,0.0044512};
  // ****************** beam energy: 39 GeV ******************
  // Event plane method: eta-sub
  // ---------------- Particle species: Phi ----------------
  Double_t pt_bin_center_39[9] = {0.4905,0.7055,0.9055,1.1055,1.2955,1.5825,1.9805,2.3805,2.8715};
  Double_t v2_values_39[9] = {0.0251251,0.0190003,0.0436195,0.06456,0.0827717,0.0947783,0.116937,0.135436,0.111607};
  Double_t v2_stat_error_39[9] = {0.00805145,0.00455958,0.00398961,0.00399311,0.00446543,0.00402194,0.00582456,0.00900157,0.0125194};
  Double_t v2_syst_low_error_39[9] = {0.00479174,0.003753,0.00214025,0.00134574,0.00119334,0.00199355,0.00184564,0.00328231,0.000280756};
  Double_t v2_syst_high_error_39[9] = {0.0036222,0.002499,0.0012339,0.00152178,0.00401162,0.00157706,0.00226953,0.00396367,0.000235831};
  // ****************** beam energy: 62.4 GeV ******************
  // Event plane method: eta-sub
  // ---------------- Particle species: Phi ----------------
  Double_t pt_bin_center_62p4[10] = {0.4925,0.7055,0.9055,1.1055,1.2955,1.5055,1.7045,1.9055,2.2355,2.8625};
  Double_t v2_values_62p4[10] = {0.0365446,0.0203302,0.0487051,0.0690831,0.0882974,0.0955015,0.108206,0.132425,0.11579,0.153856};
  Double_t v2_stat_error_62p4[10] = {0.0109052,0.00607439,0.00524557,0.00517082,0.005725,0.00657996,0.00782458,0.00937029,0.00803008,0.0156828};
  Double_t v2_syst_low_error_62p4[10] = {0.00403066,0.00221206,0.00182418,0.00198662,0.00124976,0.00159041,0.00144111,0.00255449,0.00563747,0.00114042};
  Double_t v2_syst_high_error_62p4[10] = {0.00511719,0.00241609,0.00236243,0.00191797,0.00180674,0.00222915,0.00184768,0.00316373,0.00464342,0.00135326};

  TCanv_v2_vs_pt_10_40->cd();
  TH2D * histTemp2 = new TH2D("histTemp2","histTemp2",1000,0,3,1000,-0.08,0.26);
  histTemp2->GetYaxis()->SetTitle("v_{2}");
  histTemp2->GetYaxis()->SetTitleOffset(1.2);
  histTemp2->GetXaxis()->SetTitle("pT [GeV/c]");
  histTemp2->Draw();
  TGraphErrors* gr7p7Cen_10_40 = new TGraphErrors(1, pt_bin_center_10_40, v2_values_10_40, pt_bin_center_errCen_10_40, v2_stat_error_10_40);
  gr7p7Cen_10_40->SetMarkerStyle(23);
  gr7p7Cen_10_40->SetMarkerColor(kBlack);
  gr7p7Cen_10_40->SetLineColor(kBlack);
  gr7p7Cen_10_40->SetMarkerSize(2);
  gr7p7Cen_10_40->Draw("P");
  TGraphErrors* gr11p5Cen_10_40 = new TGraphErrors(3, pt_bin_center_11p5, v2_values_11p5, 0, v2_stat_error_11p5);
  gr11p5Cen_10_40->SetMarkerStyle(23);
  gr11p5Cen_10_40->SetMarkerColor(3);
  gr11p5Cen_10_40->SetLineColor(3);
  gr11p5Cen_10_40->SetMarkerSize(2);
  gr11p5Cen_10_40->Draw("P");
  TGraphErrors* gr19p6Cen_10_40 = new TGraphErrors(5, pt_bin_center_19p6, v2_values_19p6, 0, v2_stat_error_19p6);
  gr19p6Cen_10_40->SetMarkerStyle(23);
  gr19p6Cen_10_40->SetMarkerColor(4);
  gr19p6Cen_10_40->SetLineColor(4);
  gr19p6Cen_10_40->SetMarkerSize(2);
  gr19p6Cen_10_40->Draw("P");
  TGraphErrors* gr27Cen_10_40 = new TGraphErrors(5, pt_bin_center_27, v2_values_27, 0, v2_stat_error_27);
  gr27Cen_10_40->SetMarkerStyle(23);
  gr27Cen_10_40->SetMarkerColor(6);
  gr27Cen_10_40->SetLineColor(6);
  gr27Cen_10_40->SetMarkerSize(2);
  gr27Cen_10_40->Draw("P");
  TGraphErrors* gr39Cen_10_40 = new TGraphErrors(9, pt_bin_center_39, v2_values_39, 0, v2_stat_error_39);
  gr39Cen_10_40->SetMarkerStyle(23);
  gr39Cen_10_40->SetMarkerColor(7);
  gr39Cen_10_40->SetLineColor(7);
  gr39Cen_10_40->SetMarkerSize(2);
  gr39Cen_10_40->Draw("P");
  TGraphErrors* gr62p4Cen_10_40 = new TGraphErrors(10, pt_bin_center_62p4, v2_values_62p4, 0, v2_stat_error_62p4);
  gr62p4Cen_10_40->SetMarkerStyle(23);
  gr62p4Cen_10_40->SetMarkerColor(11);
  gr62p4Cen_10_40->SetLineColor(11);
  gr62p4Cen_10_40->SetMarkerSize(2);
  gr62p4Cen_10_40->Draw("P");
  TGraphErrors* gr7p2Cen_10_40_invM = new TGraphErrors(2, pt_7p2Cen_10_40_invM, v2_7p2Cen_10_40_invM, pterr_7p2Cen_10_40_invM, v2stat_7p2Cen_10_40_invM);
  gr7p2Cen_10_40_invM->SetMarkerStyle(4);
  gr7p2Cen_10_40_invM->SetMarkerColor(kRed);
  gr7p2Cen_10_40_invM->SetLineColor(kRed);
  gr7p2Cen_10_40_invM->SetMarkerSize(2);
  // gr7p2Cen_10_40_invM->Draw("P");

  TGraphErrors* gr7p2Cen_10_40 = new TGraphErrors(2, pt_7p2Cen_10_40, v2_7p2Cen_10_40, pterr_7p2Cen_10_40, v2stat_7p2Cen_10_40);
  gr7p2Cen_10_40->SetMarkerStyle(8);
  gr7p2Cen_10_40->SetMarkerColor(kRed);
  gr7p2Cen_10_40->SetLineColor(kRed);
  gr7p2Cen_10_40->SetMarkerSize(2);
  // gr7p2Cen_10_40->Draw("P");
  // 7.2 GeV v2 vs. pT 10-40%
  TGraphErrors *graph_v2_vs_pT_10_40_7p2 = new TGraphErrors(2, x_ptSetA, y_v2_10_40, /*xErr_ptSetA*/0, yErr_stat_v2_10_40);
  graph_v2_vs_pT_10_40_7p2->SetMarkerStyle(29);
  graph_v2_vs_pT_10_40_7p2->SetMarkerColor(kBlue+2);
  graph_v2_vs_pT_10_40_7p2->SetLineColor(kBlue+2);
  graph_v2_vs_pT_10_40_7p2->SetMarkerSize(2);
  graph_v2_vs_pT_10_40_7p2->Draw("P");
  // 7.2 GeV v2 vs. pT 10-40%
  TGraphErrors *graph_v2_vs_pT_10_40_7p2_jkk = new TGraphErrors(2, x_ptSetA_jkk, y_v2_10_40_jkk, /*xErr_ptSetA*/0, yErr_stat_v2_10_40_jkk);
  graph_v2_vs_pT_10_40_7p2_jkk->SetMarkerStyle(30);
  graph_v2_vs_pT_10_40_7p2_jkk->SetMarkerColor(kBlue);
  graph_v2_vs_pT_10_40_7p2_jkk->SetLineColor(kBlue);
  graph_v2_vs_pT_10_40_7p2_jkk->SetMarkerSize(2);
  graph_v2_vs_pT_10_40_7p2_jkk->Draw("P");
  // 3bin
  TGraphErrors *graph_v2_vs_pT_10_40_7p2_3bin = new TGraphErrors(3, x_ptSetB_3bin, y_v2_10_40_3bin, xErr_ptSetB_3bin, yErr_stat_v2_10_40_3bin);
  graph_v2_vs_pT_10_40_7p2_3bin->SetMarkerStyle(23);
  graph_v2_vs_pT_10_40_7p2_3bin->SetMarkerColor(kRed);
  graph_v2_vs_pT_10_40_7p2_3bin->SetLineColor(kRed);
  graph_v2_vs_pT_10_40_7p2_3bin->SetMarkerSize(2);
  graph_v2_vs_pT_10_40_7p2_3bin->Draw("P");
  // 3bin jkk
  TGraphErrors *graph_v2_vs_pT_10_40_7p2_3bin_jkk = new TGraphErrors(3, x_ptSetB_3bin_jkk, y_v2_10_40_3bin_jkk, /*xErr_ptSetB_3bin*/0, yErr_stat_v2_10_40_3bin_jkk);
  graph_v2_vs_pT_10_40_7p2_3bin_jkk->SetMarkerStyle(32);
  graph_v2_vs_pT_10_40_7p2_3bin_jkk->SetMarkerColor(kRed+2);
  graph_v2_vs_pT_10_40_7p2_3bin_jkk->SetLineColor(kRed+2);
  graph_v2_vs_pT_10_40_7p2_3bin_jkk->SetMarkerSize(2);
  graph_v2_vs_pT_10_40_7p2_3bin_jkk->Draw("P");
  // 7.2 GeV v2 vs. pT 10-40% Shaowei
  TGraphErrors *graph_v2_vs_pT_10_40_7p2_Shaowei = new TGraphErrors(2, x_ptSetA, y_v2_10_40_Shaowei, xErr_ptSetA, yErr_stat_v2_10_40_Shaowei);
  graph_v2_vs_pT_10_40_7p2_Shaowei->SetMarkerStyle(21);
  graph_v2_vs_pT_10_40_7p2_Shaowei->SetMarkerColor(kGreen+3);
  graph_v2_vs_pT_10_40_7p2_Shaowei->SetLineColor(kGreen+3);
  graph_v2_vs_pT_10_40_7p2_Shaowei->SetMarkerSize(2);
  graph_v2_vs_pT_10_40_7p2_Shaowei->Draw("P");
  TLine *line1_1 = new TLine(0, 0, 3, 0);
  line1_1->SetLineStyle(7);
  line1_1->Draw("same");
  TLegend *legend1 = new TLegend(0.4,0.65,0.9,0.9);
  legend1->AddEntry(graph_v2_vs_pT_10_40_7p2,"7.2 GeV 10-40% invM","p");
  legend1->AddEntry(graph_v2_vs_pT_10_40_7p2_jkk,"7.2 GeV 10-40% invM - Jackknife","p");
  legend1->AddEntry(graph_v2_vs_pT_10_40_7p2_3bin,"7.2 GeV 10-40% invM - 3 pT bin","p");
  legend1->AddEntry(graph_v2_vs_pT_10_40_7p2_3bin_jkk,"7.2 GeV 10-40% invM - 3 pT bin - Jackknife","p");
  // legend1->AddEntry(gr7p2Cen_10_40_invM,"7.2 GeV 10-40% invM - Guannan","p");
  // legend1->AddEntry(gr7p2Cen_10_40,"7.2 GeV 10-40% etasub - Guannan","p");
  legend1->AddEntry(graph_v2_vs_pT_10_40_7p2_Shaowei,"7.2 GeV 10-40% etasub - Shaowei","p");
  legend1->AddEntry(gr7p7Cen_10_40,"7.7 GeV 10-40% etasub - BES-I","p");
  legend1->AddEntry(gr11p5Cen_10_40,"11.5 GeV 10-40% etasub - BES-I","p");
  legend1->AddEntry(gr19p6Cen_10_40,"19.6 GeV 10-40% etasub - BES-I","p");
  legend1->AddEntry(gr27Cen_10_40,"27 GeV 10-40% etasub - BES-I","p");
  legend1->AddEntry(gr39Cen_10_40,"39 GeV 10-40% etasub - BES-I","p");
  legend1->AddEntry(gr62p4Cen_10_40,"62.4 GeV 10-40% etasub - BES-I","p");
  legend1->Draw("same");

  // -------------------------- 0 - 60% / 80% ----------------------------------
  TCanvas *TCanv_v2_vs_pt_0_60 = new TCanvas("TCanv_v2_vs_pt_0_60","v2 vs. pT 0-60%",200,10,1024,768);
  TCanv_v2_vs_pt_0_60->DrawFrame(0, -0.1, 3, 0.2);

  // 7.2 GeV v2 vs. pT 0-60%
  // double x_ptSetA[2]    = { 0.9,1.8};
  // double xErr_ptSetA[2] = {0.3,0.6};
  double y_v2_0_60[2]      = {0.0108515, -0.0102652};
  double yErr_stat_v2_0_60[2] = {0.0250194, 0.0256442};

  // 7.2 GeV v2 vs. pT 0-70% Shaowei
  double x_0_70[3]    = {0.812924, 1.19676, 1.6985};
  double xErr_0_70[3] = {0};
  double y_v2_0_70_Shaowei[3]  = {-0.0246507,-0.017738,-0.0521873};
  double yErr_stat_v2_0_70_Shaowei[3] = {0.01228,0.0112587,0.0169663};
  TCanv_v2_vs_pt_0_60->cd();

  // 7.2 GeV v2 vs. pT 0-60% Guannan
  //strict eta sub method
  float pt_7p2Cen_0_60[2] = {0.85, 1.85};
  float pterr_7p2Cen_0_60[2]={0.35,0.65};
  float v2_7p2Cen_0_60[2] = {0.086357, 0.0136125 };
  float v2stat_7p2Cen_0_60[2] = {0.0647058, 0.0837308};
  //invMass sub method
  float pt_7p2Cen_0_60_invM[2] = {0.85 + 0.05, 1.85 + 0.05};
  float pterr_7p2Cen_0_60_invM[2]={0.35,0.65};
  float v2_7p2Cen_0_60_invM[2] = {0.0412306, 0.0123695 };
  float v2stat_7p2Cen_0_60_invM[2] = {0.0420239, 0.0593082};
  // centrality: 0-80%
  // beam energy: 7.7 GeV
  // Event plane method: eta-sub
  // ---------------- Particle species: Phi ----------------
  Double_t pt_bin_center[2] = {0.9185,1.6685};
  Double_t pt_bin_center_err[2] = {0};
  Double_t v2_values[2] = {0.0478335,-0.049273};
  Double_t v2_stat_error[2] = {0.0247765,0.0605828};
  histTemp2->Draw();
  TGraphErrors* gr7p7Cen_0_80 = new TGraphErrors(2, pt_bin_center, v2_values, pt_bin_center_err, v2_stat_error);
  gr7p7Cen_0_80->SetMarkerStyle(23);
  gr7p7Cen_0_80->SetMarkerColor(kBlack);
  gr7p7Cen_0_80->SetLineColor(kBlack);
  gr7p7Cen_0_80->SetMarkerSize(2);
  gr7p7Cen_0_80->Draw("P");
  TGraphErrors* gr7p2Cen_0_60_invM = new TGraphErrors(2, pt_7p2Cen_0_60_invM, v2_7p2Cen_0_60_invM, pterr_7p2Cen_0_60_invM, v2stat_7p2Cen_0_60_invM);
  gr7p2Cen_0_60_invM->SetMarkerStyle(4);
  gr7p2Cen_0_60_invM->SetMarkerColor(kRed);
  gr7p2Cen_0_60_invM->SetLineColor(kRed);
  gr7p2Cen_0_60_invM->SetMarkerSize(2);
  gr7p2Cen_0_60_invM->Draw("P");
  TGraphErrors* gr7p2Cen_0_60 = new TGraphErrors(2, pt_7p2Cen_0_60, v2_7p2Cen_0_60, pterr_7p2Cen_0_60, v2stat_7p2Cen_0_60);
  gr7p2Cen_0_60->SetMarkerStyle(8);
  gr7p2Cen_0_60->SetMarkerColor(kRed);
  gr7p2Cen_0_60->SetLineColor(kRed);
  gr7p2Cen_0_60->SetMarkerSize(2);
  gr7p2Cen_0_60->Draw("P");
  // 7.2 GeV v2 vs. pT 0-60%
  TGraphErrors *graph_v2_vs_pT_0_60_7p2 = new TGraphErrors(2, x_ptSetA, y_v2_0_60, xErr_ptSetA, yErr_stat_v2_0_60);
  graph_v2_vs_pT_0_60_7p2->SetMarkerStyle(30);
  graph_v2_vs_pT_0_60_7p2->SetMarkerColor(kBlue);
  graph_v2_vs_pT_0_60_7p2->SetLineColor(kBlue);
  graph_v2_vs_pT_0_60_7p2->SetMarkerSize(2);
  graph_v2_vs_pT_0_60_7p2->Draw("P");
  // 7.2 GeV v2 vs. pT 0-70% Shaowei
  TGraphErrors *graph_v2_vs_pT_0_70_7p2_SL = new TGraphErrors(3, x_0_70, y_v2_0_70_Shaowei, xErr_0_70, yErr_stat_v2_0_70_Shaowei);
  graph_v2_vs_pT_0_70_7p2_SL->SetMarkerStyle(21);
  graph_v2_vs_pT_0_70_7p2_SL->SetMarkerColor(kGreen+3);
  graph_v2_vs_pT_0_70_7p2_SL->SetLineColor(kGreen+3);
  graph_v2_vs_pT_0_70_7p2_SL->SetMarkerSize(2);
  graph_v2_vs_pT_0_70_7p2_SL->Draw("P");
  TLine *line1_2 = new TLine(0, 0, 3, 0);
  line1_2->SetLineStyle(7);
  line1_2->Draw("same");
  TLegend *legend1_1 = new TLegend(0.4,0.65,0.9,0.9);
  legend1_1->AddEntry(graph_v2_vs_pT_0_60_7p2,"7.2 GeV 0-60% invM - This analysis","p");
  legend1_1->AddEntry(gr7p2Cen_0_60_invM,"7.2 GeV 0-60% invM - Guannan","p");
  legend1_1->AddEntry(gr7p2Cen_0_60,"7.2 GeV 0-60% etasub - Guannan","p");
  legend1_1->AddEntry(graph_v2_vs_pT_0_70_7p2_SL,"7.2 GeV 0-70% etasub - Shaowei","p");
  legend1_1->AddEntry(gr7p7Cen_0_80,"7.7 GeV 0-80% etasub - BES-I","p");
  legend1_1->Draw("same");
  //   v1 10-40% rapSetA_centSetA: -0.00596924, 0.016996, 0.0421228
  // v1 err 10-40% rapSetA_centSetA: 0.0223695, 0.0134569, 0.0488063
  // v2 10-40% ptSetA_centSetA: 0.0256657, 0.0120596
  // v2 Err 10-40% ptSetA_centSetA: 0.023679, 0.0240889
  // v2 40-60% ptSetA_centSetA: 0.0108515, -0.0102652
  // v2 Err 40-60% ptSetA_centSetA: 0.0250194, 0.0256442
  // cout<< "size of array test" << sizeof(v2_stat_error)/sizeof(v2_stat_error[0]) <<endl;
  // ================== NCQ scaling =================================
  // 10 - 40%
  TCanvas *TCanv_v2_vs_mt_10_40 = new TCanvas("TCanv_v2_vs_mt_10_40","v2 vs. (mT - m0)/nq 0-60%",200,10,1024,768);
  TCanv_v2_vs_mt_10_40->DrawFrame(0, -0.04, 1.6, 0.13);
  TH2D * histTemp4 = new TH2D("histTemp4","histTemp4",1000,0,1.,1000,-0.04,0.13);
  histTemp4->GetYaxis()->SetTitle("v_{2}/n_{q}");
  histTemp4->GetYaxis()->SetTitleOffset(1.2);
  histTemp4->GetXaxis()->SetTitleOffset(1.2);
  histTemp4->GetXaxis()->SetTitle("(m_{T}-m_{0})/n_{q} (GeV/c^{2})");
  histTemp4->Draw();
  // centrality: 10%-40%

  // ************************************************************************
  // ****************** beam energy: 7.7 GeV ******************
  // Event plane method: eta-sub
  // ---------------- Particle species: Phi ----------------
  Double_t pt_bin_center_phi_7p2[2] = {0.9, 1.8};
  Double_t v2_values_phi_7p2[2] = {-0.0097, -0.0205};
  Double_t v2_stat_error_phi_7p2[2] = {0.0281, 0.0277};
  // ---------------- Particle species: Phi ----------------
  Double_t pt_bin_center_phi_7p2_3bin[3] = {0.85, 1.55, 2.5};
  Double_t v2_values_phi_7p2_3bin[3] = { -0.0276968, -0.00838587, 0.0991779};
  Double_t v2_stat_error_phi_7p2_3bin[3] = { 0.031929, 0.0260998, 0.150601};
  // ****************** beam energy: 7.7 GeV ******************
  // Event plane method: eta-sub
  // ---------------- Particle species: Phi ----------------
  Double_t pt_bin_center_phi_7p7[1] = {0.9945};
  Double_t v2_values_phi_7p7[1] = {0.0399507};
  Double_t v2_stat_error_phi_7p7[1] = {0.0193815};
  Double_t v2_syst_low_error_phi_7p7[1] = {0.0098416};
  Double_t v2_syst_high_error_phi_7p7[1] = {0.00575247};
  // ---------------- Particle species: Lambda ----------------
  Double_t pt_bin_center_Lambda_7p7[10] = {0.4935,0.7055,0.9055,1.1055,1.2955,1.5055,1.7055,1.9055,2.1695,2.6175};
  Double_t v2_values_Lambda_7p7[10] = {0.0056175,0.036717,0.0481013,0.0792141,0.0815696,0.0867856,0.108612,0.122262,0.133421,0.153319};
  Double_t v2_stat_error_Lambda_7p7[10] = {0.00771747,0.0048333,0.00414777,0.00425718,0.00485913,0.00600927,0.00776441,0.0104924,0.0120231,0.0227567};
  Double_t v2_syst_low_error_Lambda_7p7[10] = {0.00303084,0.00143998,0.00101268,0.000986658,0.00162709,0.00191824,0.00258821,0.00401433,0.00404447,0.0104834};
  Double_t v2_syst_high_error_Lambda_7p7[10] = {0.00399962,0.00110434,0.00167376,0.0018407,0.00107831,0.000745689,0.00170209,0.00187684,0.00493309,0.0104299};
  // ---------------- Particle species: K0S ----------------
  Double_t pt_bin_center_K0S_7p7[8] = {0.3055,0.5055,0.7055,0.9055,1.1055,1.2955,1.5675,2.0025};
  Double_t v2_values_K0S_7p7[8] = {0.0053047,0.0263517,0.037153,0.0669851,0.0704264,0.0828911,0.083941,0.127632};
  Double_t v2_stat_error_K0S_7p7[8] = {0.00689945,0.00437349,0.00392164,0.00440512,0.00562355,0.00775186,0.00927243,0.0198935};
  Double_t v2_syst_low_error_K0S_7p7[8] = {2.50309e-05,0.000185801,0.000318705,0.000232844,0.000408465,0.000216072,0.0027795,0.00365314};
  Double_t v2_syst_high_error_K0S_7p7[8] = {2.53773e-05,0.00018327,0.000312561,0.000239669,0.000408919,0.000206365,0.0026328,0.00343635};
  // ---------------- Particle species: Proton ----------------
  Double_t pt_bin_center_proton_7p7[12] = {0.315,0.495,0.705,0.885,1.095,1.275,1.485,1.695,1.875,2.085,2.265,2.625};
  Double_t v2_values_proton_7p7[12] = {0.00643856,0.0196543,0.0362226,0.0541206,0.0781074,0.090889,0.10682,0.121606,0.131171,0.13338,0.14221,0.152974};
  Double_t v2_stat_error_proton_7p7[12] = {0.0028076,0.00108471,0.00093993,0.000987168,0.00115896,0.00144923,0.00189495,0.00257513,0.00358486,0.00494914,0.0070575,0.00820795};
  Double_t v2_syst_low_error_proton_7p7[12] = {0.00276658,0.000576754,0.00109959,0.000118159,0.000515559,0.000388675,0.000212473,0.000415654,0.000904369,0.000963559,0.00125468,0.00391717};
  Double_t v2_syst_high_error_proton_7p7[12] = {0.000944912,0.000311883,0.000503305,0.000107656,0.000468809,0.000749411,0.00011288,0.000274677,0.000484163,0.000938662,0.000643771,0.00192167};
  // ---------------- Particle species: PiP ----------------
  Double_t pt_bin_center_PiP_7p7[7] = {0.285,0.495,0.675,0.885,1.155,1.545,1.935};
  Double_t v2_values_PiP_7p7[7] = {0.0225062,0.0415791,0.0560643,0.0678896,0.0819491,0.100776,0.120508};
  Double_t v2_stat_error_PiP_7p7[7] = {0.000479722,0.000616763,0.000883004,0.00131197,0.00165343,0.00380281,0.00915605};
  Double_t v2_syst_low_error_PiP_7p7[7] = {3.61284e-10,1.42286e-05,0.000107432,2.12591e-06,1.33104e-05,0.000353813,0.00177874};
  Double_t v2_syst_high_error_PiP_7p7[7] = {4.7195e-11,2.84574e-05,0.000217747,4.25174e-06,6.65502e-06,0.000186354,0.00513822};
  // ---------------- Particle species: KP ----------------
  Double_t pt_bin_center_KP_7p7[7] = {0.315,0.495,0.705,0.885,1.155,1.545,1.965};
  Double_t v2_values_KP_7p7[7] = {0.0123349,0.0284567,0.0493235,0.0674997,0.080886,0.100798,0.0878592};
  Double_t v2_stat_error_KP_7p7[7] = {0.00259688,0.00169075,0.00177886,0.00220443,0.00239161,0.00478715,0.0106033};
  Double_t v2_syst_low_error_KP_7p7[7] = {0.000213675,1.68026e-05,8.83694e-05,9.47739e-05,5.02883e-05,0.000454596,0.00393434};
  Double_t v2_syst_high_error_KP_7p7[7] = {0.000181583,3.12802e-05,5.80217e-05,0.000169471,8.27577e-05,0.000264862,0.00378213};  // centrality: 0-80%
  // ======== Loop that convert v2 vs. pT to v2/nq vs. (mT - m0)/nq ============
  for(int i=0; i<sizeof(pt_bin_center_phi_7p2)/sizeof(pt_bin_center_phi_7p2[0]); i++){
    v2_values_phi_7p2[i] = d_v2nq(v2_values_phi_7p2[i], 2);
    v2_stat_error_phi_7p2[i] = d_v2nq(v2_stat_error_phi_7p2[i], 2);
    pt_bin_center_phi_7p2[i] = d_mTm0nq(pt_bin_center_phi_7p2[i],  2,  _massPhi );
  }
  for(int i=0; i<sizeof(pt_bin_center_phi_7p2_3bin)/sizeof(pt_bin_center_phi_7p2_3bin[0]); i++){
    v2_values_phi_7p2_3bin[i] = d_v2nq(v2_values_phi_7p2_3bin[i], 2);
    v2_stat_error_phi_7p2_3bin[i] = d_v2nq(v2_stat_error_phi_7p2_3bin[i], 2);
    pt_bin_center_phi_7p2_3bin[i] = d_mTm0nq(pt_bin_center_phi_7p2_3bin[i],  2,  _massPhi );
  }
  for(int i=0; i<sizeof(pt_bin_center_phi_7p7)/sizeof(pt_bin_center_phi_7p7[0]); i++){
    v2_values_phi_7p7[i] = d_v2nq(v2_values_phi_7p7[i], 2);
    v2_stat_error_phi_7p7[i] = d_v2nq(v2_stat_error_phi_7p7[i], 2);
    pt_bin_center_phi_7p7[i] = d_mTm0nq(pt_bin_center_phi_7p7[i],  2,  _massPhi );
  }
  for(int i=0; i<sizeof(pt_bin_center_Lambda_7p7)/sizeof(pt_bin_center_Lambda_7p7[0]); i++){
    v2_values_Lambda_7p7[i] = d_v2nq(v2_values_Lambda_7p7[i], 2);
    v2_stat_error_Lambda_7p7[i] = d_v2nq(v2_stat_error_Lambda_7p7[i], 2);
    pt_bin_center_Lambda_7p7[i] = d_mTm0nq(pt_bin_center_Lambda_7p7[i],  2,  _massLambda );
  }
  for(int i=0; i<sizeof(pt_bin_center_K0S_7p7)/sizeof(pt_bin_center_K0S_7p7[0]); i++){
    v2_values_K0S_7p7[i] = d_v2nq(v2_values_K0S_7p7[i], 2);
    v2_stat_error_K0S_7p7[i] = d_v2nq(v2_stat_error_K0S_7p7[i], 2);
    pt_bin_center_K0S_7p7[i] = d_mTm0nq(pt_bin_center_K0S_7p7[i],  2,  _massKaonShort );
  }
  for(int i=0; i<sizeof(pt_bin_center_proton_7p7)/sizeof(pt_bin_center_proton_7p7[0]); i++){
    v2_values_proton_7p7[i] = d_v2nq(v2_values_proton_7p7[i], 3);
    v2_stat_error_proton_7p7[i] = d_v2nq(v2_stat_error_proton_7p7[i], 3);
    pt_bin_center_proton_7p7[i] = d_mTm0nq(pt_bin_center_proton_7p7[i],  3,  _massProton );
  }
  for(int i=0; i<sizeof(pt_bin_center_PiP_7p7)/sizeof(pt_bin_center_PiP_7p7[0]); i++){
    v2_values_PiP_7p7[i] = d_v2nq(v2_values_PiP_7p7[i], 2);
    v2_stat_error_PiP_7p7[i] = d_v2nq(v2_stat_error_PiP_7p7[i], 2);
    pt_bin_center_PiP_7p7[i] = d_mTm0nq(pt_bin_center_PiP_7p7[i],  2,  _massPion );
  }
  for(int i=0; i<sizeof(pt_bin_center_KP_7p7)/sizeof(pt_bin_center_KP_7p7[0]); i++){
    v2_values_KP_7p7[i] = d_v2nq(v2_values_KP_7p7[i], 2);
    v2_stat_error_KP_7p7[i] = d_v2nq(v2_stat_error_KP_7p7[i], 2);
    pt_bin_center_KP_7p7[i] = d_mTm0nq(pt_bin_center_KP_7p7[i],  2,  _massPion );
  }
  TGraphErrors* gr7p7CenMt_10_40 = new TGraphErrors(1, pt_bin_center_phi_7p7, v2_values_phi_7p7, 0, v2_stat_error_phi_7p7);
  gr7p7CenMt_10_40->SetMarkerStyle(23);
  gr7p7CenMt_10_40->SetMarkerColor(kBlack);
  gr7p7CenMt_10_40->SetLineWidth(2);
  gr7p7CenMt_10_40->SetLineColor(kBlack);
  gr7p7CenMt_10_40->SetMarkerSize(2);
  gr7p7CenMt_10_40->Draw("P");
  TGraphErrors* gr7p7CenMt_10_40_Lambda = new TGraphErrors(10, pt_bin_center_Lambda_7p7, v2_values_Lambda_7p7, 0, v2_stat_error_Lambda_7p7);
  gr7p7CenMt_10_40_Lambda->SetMarkerStyle(22);
  gr7p7CenMt_10_40_Lambda->SetMarkerColor(kGreen+2);
  gr7p7CenMt_10_40_Lambda->SetLineWidth(2);
  gr7p7CenMt_10_40_Lambda->SetLineColor(kGreen+2);
  gr7p7CenMt_10_40_Lambda->SetMarkerSize(2);
  gr7p7CenMt_10_40_Lambda->Draw("P");
  TGraphErrors* gr7p7CenMt_10_40_K0S = new TGraphErrors(8, pt_bin_center_K0S_7p7,v2_values_K0S_7p7 , 0, v2_stat_error_K0S_7p7);
  gr7p7CenMt_10_40_K0S->SetMarkerStyle(28);
  gr7p7CenMt_10_40_K0S->SetMarkerColor(kRed+2);
  gr7p7CenMt_10_40_K0S->SetLineWidth(2);
  gr7p7CenMt_10_40_K0S->SetLineColor(kRed+2);
  gr7p7CenMt_10_40_K0S->SetMarkerSize(2);
  gr7p7CenMt_10_40_K0S->Draw("P");
  TGraphErrors* gr7p7CenMt_10_40_Proton = new TGraphErrors(12, pt_bin_center_proton_7p7,v2_values_proton_7p7 , 0, v2_stat_error_proton_7p7);
  gr7p7CenMt_10_40_Proton->SetMarkerStyle(20);
  gr7p7CenMt_10_40_Proton->SetMarkerColor(kAzure+2);
  gr7p7CenMt_10_40_Proton->SetLineWidth(2);
  gr7p7CenMt_10_40_Proton->SetLineColor(kAzure+2);
  gr7p7CenMt_10_40_Proton->SetMarkerSize(2);
  gr7p7CenMt_10_40_Proton->Draw("P");
  TGraphErrors* gr7p7CenMt_10_40_PiP = new TGraphErrors(7, pt_bin_center_PiP_7p7,v2_values_PiP_7p7 , 0, v2_stat_error_PiP_7p7);
  gr7p7CenMt_10_40_PiP->SetMarkerStyle(30);
  gr7p7CenMt_10_40_PiP->SetMarkerColor(kMagenta+2);
  gr7p7CenMt_10_40_PiP->SetLineWidth(2);
  gr7p7CenMt_10_40_PiP->SetLineColor(kMagenta+2);
  gr7p7CenMt_10_40_PiP->SetMarkerSize(2);
  gr7p7CenMt_10_40_PiP->Draw("P");
  TGraphErrors* gr7p7CenMt_10_40_KP = new TGraphErrors(7, pt_bin_center_KP_7p7,v2_values_KP_7p7 , 0, v2_stat_error_KP_7p7);
  gr7p7CenMt_10_40_KP->SetMarkerStyle(20);
  gr7p7CenMt_10_40_KP->SetMarkerColor(kGray+2);
  gr7p7CenMt_10_40_KP->SetLineWidth(2);
  gr7p7CenMt_10_40_KP->SetLineColor(kGray+2);
  gr7p7CenMt_10_40_KP->SetMarkerSize(2);
  gr7p7CenMt_10_40_KP->Draw("P");
  TGraphErrors *graph_v2_vs_mT_10_40_7p2 = new TGraphErrors(2, pt_bin_center_phi_7p2, v2_values_phi_7p2, 0, v2_stat_error_phi_7p2);
  graph_v2_vs_mT_10_40_7p2->SetMarkerStyle(29);
  graph_v2_vs_mT_10_40_7p2->SetMarkerColor(kBlue);
  graph_v2_vs_mT_10_40_7p2->SetLineColor(kBlue);
  graph_v2_vs_mT_10_40_7p2->SetLineWidth(2);
  graph_v2_vs_mT_10_40_7p2->SetMarkerSize(2);
  graph_v2_vs_mT_10_40_7p2->Draw("P");
  TGraphErrors *graph_v2_vs_mT_10_40_7p2_3bin = new TGraphErrors(3, pt_bin_center_phi_7p2_3bin, v2_values_phi_7p2_3bin, 0, v2_stat_error_phi_7p2_3bin);
  graph_v2_vs_mT_10_40_7p2_3bin->SetMarkerStyle(23);
  graph_v2_vs_mT_10_40_7p2_3bin->SetMarkerColor(kRed);
  graph_v2_vs_mT_10_40_7p2_3bin->SetLineColor(kRed);
  graph_v2_vs_mT_10_40_7p2_3bin->SetLineWidth(2);
  graph_v2_vs_mT_10_40_7p2_3bin->SetMarkerSize(2);
  graph_v2_vs_mT_10_40_7p2_3bin->Draw("P");
  TLine *line1_4 = new TLine(0, 0, 1., 0);
  line1_4->SetLineStyle(7);
  line1_4->Draw("same");
  TLegend *legend1_4 = new TLegend(0.4,0.65,0.9,0.9);
  legend1_4->AddEntry(graph_v2_vs_mT_10_40_7p2,"#phi FXT 7.2 GeV 10-40%","p");
  legend1_4->AddEntry(graph_v2_vs_mT_10_40_7p2_3bin,"#phi FXT 7.2 GeV 10-40%, 3 bin","p");
  legend1_4->AddEntry(gr7p7CenMt_10_40,"#phi 7.7 GeV 10-40% - BES-I","p");
  // legend1_4->AddEntry(gr7p7CenMt_0_80_XiM,"#Xi^{-} 7.7 GeV 0-80% - BES-I","p");
  // legend1_4->AddEntry(gr7p7CenMt_0_80_OmegaM,"#Omega^{-} 7.7 GeV 0-80% - BES-I","p");
  legend1_4->AddEntry(gr7p7CenMt_10_40_Lambda,"#Lambda 7.7 GeV 10-40% - BES-I","p");
  legend1_4->AddEntry(gr7p7CenMt_10_40_K0S,"K_{s}^{0} 7.7 GeV 10-40% - BES-I","p");
  legend1_4->AddEntry(gr7p7CenMt_10_40_Proton,"p 7.7 GeV 10-40% - BES-I","p");
  legend1_4->AddEntry(gr7p7CenMt_10_40_PiP,"#pi^{+} 7.7 GeV 10-40% - BES-I","p");
  legend1_4->AddEntry(gr7p7CenMt_10_40_KP,"K^{+} 7.7 GeV 10-40% - BES-I","p");
  legend1_4->Draw("same");
  // 0 - 60/80%
  TCanvas *TCanv_v2_vs_mt_0_60 = new TCanvas("TCanv_v2_vs_mt_0_60","v2 vs. (mT - m0)/nq 0-60%",200,10,1024,768);
  TCanv_v2_vs_mt_0_60->DrawFrame(0, -0.06, 1.6, 0.1);
  TH2D * histTemp3 = new TH2D("histTemp3","histTemp3",1000,0,1.6,1000,-0.06,0.1);
  histTemp3->GetYaxis()->SetTitle("v_{2}/n_{q}");
  histTemp3->GetYaxis()->SetTitleOffset(1);
  histTemp3->GetXaxis()->SetTitle("(m_{T}-m_{0})/n_{q} (GeV/c^{2})");
  histTemp3->Draw();
  // beam energy: 7.7 GeV
  // Event plane method: eta-sub
  // ---------------- Particle species: Phi ----------------
  Double_t pt_bin_center_7p7[2] = {0.9185,1.6685};
  Double_t v2_values_7p7[2] = {0.0478335,-0.049273};
  Double_t v2_stat_error_7p7[2] = {0.0247765,0.0605828};
  // Double_t v2_syst_low_error[2] = {0.00964971,0.0232742};
  // Double_t v2_syst_high_error[2] = {0.00910412,0.0277136};
  // Double_t v2_syst_global_error = 0.00369459;
  // ---------------- Particle species: XiM ----------------
  Double_t pt_bin_center_XiM[3] = {0.9545,1.4595,2.0975};
  Double_t v2_values_XiM[3] = {-0.0263616,0.0561645,0.130441};
  Double_t v2_stat_error_XiM[3] = {0.0277499,0.0247706,0.0442127};
  // Double_t v2_syst_low_error[3] = {0.0148688,0.00998551,0.0207038};
  // Double_t v2_syst_high_error[3] = {0.0127196,0.0127153,0.0237823};
  // Double_t v2_syst_global_error = 0.00674074;
  // ---------------- Particle species: OmegaM ----------------
  Double_t pt_bin_center_OmegaM[1] = {1.5805};
  Double_t v2_values_OmegaM[1] = {0.0870084};
  Double_t v2_stat_error_OmegaM[1] = {0.213381};
  // ---------------- Particle species: _Lambda ----------------
  Double_t pt_bin_center_Lambda[7] = {0.4935,0.7055,0.9045,1.1055,1.2955,1.6335,2.2135};
  Double_t v2_values_Lambda[7] = {0.0143684,0.0350052,0.0465096,0.0638577,0.0662386,0.0722904,0.0820912};
  Double_t v2_stat_error_Lambda[7] = {0.00866406,0.00536647,0.00456596,0.00467365,0.00531188,0.00469397,0.0122371};
  // ---------------- Particle species: _K0S ----------------
  Double_t pt_bin_center_K0S[7] = {0.3055,0.5045,0.7055,0.9055,1.1045,1.3635,1.8275};
  Double_t v2_values_K0S[7] = {0.00080378,0.0220384,0.0332416,0.0540146,0.0557628,0.0674778,0.0669422};
  Double_t v2_stat_error_K0S[7] = {0.00760497,0.00477263,0.00425852,0.00475247,0.00607755,0.0069004,0.0138193};
  // ---------------- Particle species: _Proton ----------------
  Double_t pt_bin_center_Proton[13] = {0.315,0.495,0.705,0.885,1.095,1.275,1.485,1.695,1.875,2.085,2.265,2.475,2.835};
  Double_t v2_values_Proton[13] = {0.00691689,0.0161549,0.0312339,0.0462411,0.0616304,0.0705497,0.0824948,0.0931209,0.0985242,0.0993178,0.107344,0.105471,0.0974711};
  Double_t v2_stat_error_Proton[13] = {0.00306278,0.00118166,0.00102077,0.00106674,0.00124804,0.00155582,0.00203117,0.00275382,0.00383099,0.00528279,0.00751169,0.0110623,0.0138618};
  // ---------------- Particle species: _PiP ----------------
  Double_t pt_bin_center_PiP[7] = {0.285,0.495,0.675,0.885,1.155,1.545,1.935};
  Double_t v2_values_PiP[7] = {0.0190878,0.0361773,0.0478764,0.0578427,0.0667559,0.0828597,0.0895746};
  Double_t v2_stat_error_PiP[7] = {0.00052733,0.000676753,0.000966526,0.00143331,0.00180211,0.00413486,0.00979729};
  // ---------------- Particle species: _KP -----------------
  Double_t pt_bin_center_KP[8] = {0.315,0.495,0.705,0.885,1.155,1.545,1.935,2.325};
  Double_t v2_values_KP[8] = {0.010691,0.0263934,0.0397986,0.0537884,0.0622878,0.0765433,0.077197,0.0529762};
  Double_t v2_stat_error_KP[8] = {0.00280781,0.00182188,0.0019098,0.00236217,0.0025566,0.00510552,0.01125,0.0305002};
  // ======== Loop that convert v2 vs. pT to v2/nq vs. (mT - m0)/nq ============
  for(int i=0; i<sizeof(pt_bin_center_7p7)/sizeof(pt_bin_center_7p7[0]); i++){
    v2_values_7p7[i] = d_v2nq(v2_values_7p7[i], 2);
    v2_stat_error_7p7[i] = d_v2nq(v2_stat_error_7p7[i], 2);
    pt_bin_center_7p7[i] = d_mTm0nq(pt_bin_center_7p7[i],  2,  _massPhi );
  }
  for(int i=0; i<sizeof(pt_bin_center_XiM)/sizeof(pt_bin_center_XiM[0]); i++){
    v2_values_XiM[i] = d_v2nq(v2_values_XiM[i], 3);
    v2_stat_error_XiM[i] = d_v2nq(v2_stat_error_XiM[i], 3);
    pt_bin_center_XiM[i] = d_mTm0nq(pt_bin_center_XiM[i],  3,  _massXiMinus );
  }
  v2_values_OmegaM[0] = d_v2nq(v2_values_OmegaM[0], 3);
  v2_stat_error_OmegaM[0] = d_v2nq(v2_stat_error_OmegaM[0], 3);
  pt_bin_center_OmegaM[0] = d_mTm0nq(pt_bin_center_OmegaM[0],  3,  _massOmegaMinus );
  for(int i=0; i<sizeof(pt_bin_center_Lambda)/sizeof(pt_bin_center_Lambda[0]); i++){
    v2_values_Lambda[i] = d_v2nq(v2_values_Lambda[i], 3);
    v2_stat_error_Lambda[i] = d_v2nq(v2_stat_error_Lambda[i], 3);
    pt_bin_center_Lambda[i] = d_mTm0nq(pt_bin_center_Lambda[i],  3,  _massLambda );
  }
  for(int i=0; i<sizeof(pt_bin_center_K0S)/sizeof(pt_bin_center_K0S[0]); i++){
    v2_values_K0S[i] = d_v2nq(v2_values_K0S[i], 2);
    v2_stat_error_K0S[i] = d_v2nq(v2_stat_error_K0S[i], 2);
    pt_bin_center_K0S[i] = d_mTm0nq(pt_bin_center_K0S[i],  2,  _massKaonShort );
  }
  for(int i=0; i<sizeof(pt_bin_center_Proton)/sizeof(pt_bin_center_Proton[0]); i++){
    v2_values_Proton[i] = d_v2nq(v2_values_Proton[i], 3);
    v2_stat_error_Proton[i] = d_v2nq(v2_stat_error_Proton[i], 3);
    pt_bin_center_Proton[i] = d_mTm0nq(pt_bin_center_Proton[i],  3,  _massProton );
  }
  for(int i=0; i<sizeof(pt_bin_center_PiP)/sizeof(pt_bin_center_PiP[0]); i++){
    v2_values_PiP[i] = d_v2nq(v2_values_PiP[i], 2);
    v2_stat_error_PiP[i] = d_v2nq(v2_stat_error_PiP[i], 2);
    pt_bin_center_PiP[i] = d_mTm0nq(pt_bin_center_PiP[i],  2,  _massPion );
  }
  for(int i=0; i<sizeof(pt_bin_center_KP)/sizeof(pt_bin_center_KP[0]); i++){
    v2_values_KP[i] = d_v2nq(v2_values_KP[i], 2);
    v2_stat_error_KP[i] = d_v2nq(v2_stat_error_KP[i], 2);
    pt_bin_center_KP[i] = d_mTm0nq(pt_bin_center_KP[i],  2,  _massKaon );
  }
  for(int i=0; i<sizeof(x_ptSetA)/sizeof(x_ptSetA[0]); i++){
    y_v2_0_60[i] = d_v2nq(y_v2_0_60[i], 2);
    yErr_stat_v2_0_60[i] = d_v2nq(yErr_stat_v2_0_60[i], 2);
    x_ptSetA[i] = d_mTm0nq(x_ptSetA[i],  2,  _massPhi );
  }

  TGraphErrors* gr7p7CenMt_0_80 = new TGraphErrors(2, pt_bin_center_7p7, v2_values_7p7, 0, v2_stat_error_7p7);
  gr7p7CenMt_0_80->SetMarkerStyle(23);
  gr7p7CenMt_0_80->SetMarkerColor(kRed);
  gr7p7CenMt_0_80->SetLineWidth(2);
  gr7p7CenMt_0_80->SetLineColor(kRed);
  gr7p7CenMt_0_80->SetMarkerSize(2);
  gr7p7CenMt_0_80->Draw("P");
  TGraphErrors* gr7p7CenMt_0_80_XiM = new TGraphErrors(3, pt_bin_center_XiM, v2_values_XiM, 0, v2_stat_error_XiM);
  gr7p7CenMt_0_80_XiM->SetMarkerStyle(22);
  gr7p7CenMt_0_80_XiM->SetMarkerColor(16);
  gr7p7CenMt_0_80_XiM->SetLineWidth(2);
  gr7p7CenMt_0_80_XiM->SetLineColor(16);
  gr7p7CenMt_0_80_XiM->SetMarkerSize(2);
  // gr7p7CenMt_0_80_XiM->Draw("P");
  TGraphErrors* gr7p7CenMt_0_80_OmegaM = new TGraphErrors(1, pt_bin_center_OmegaM, v2_values_OmegaM, 0, v2_stat_error_OmegaM);
  gr7p7CenMt_0_80_OmegaM->SetMarkerStyle(21);
  gr7p7CenMt_0_80_OmegaM->SetMarkerColor(7);
  gr7p7CenMt_0_80_OmegaM->SetLineWidth(2);
  gr7p7CenMt_0_80_OmegaM->SetLineColor(7);
  gr7p7CenMt_0_80_OmegaM->SetMarkerSize(2);
  // gr7p7CenMt_0_80_OmegaM->Draw("P");
  TGraphErrors* gr7p7CenMt_0_80_Lambda = new TGraphErrors(7, pt_bin_center_Lambda, v2_values_Lambda, 0, v2_stat_error_Lambda);
  gr7p7CenMt_0_80_Lambda->SetMarkerStyle(22);
  gr7p7CenMt_0_80_Lambda->SetMarkerColor(kGreen+2);
  gr7p7CenMt_0_80_Lambda->SetLineWidth(2);
  gr7p7CenMt_0_80_Lambda->SetLineColor(kGreen+2);
  gr7p7CenMt_0_80_Lambda->SetMarkerSize(2);
  gr7p7CenMt_0_80_Lambda->Draw("P");
  TGraphErrors* gr7p7CenMt_0_80_K0S = new TGraphErrors(7, pt_bin_center_K0S,v2_values_K0S , 0, v2_stat_error_K0S);
  gr7p7CenMt_0_80_K0S->SetMarkerStyle(28);
  gr7p7CenMt_0_80_K0S->SetMarkerColor(kRed+2);
  gr7p7CenMt_0_80_K0S->SetLineWidth(2);
  gr7p7CenMt_0_80_K0S->SetLineColor(kRed+2);
  gr7p7CenMt_0_80_K0S->SetMarkerSize(2);
  gr7p7CenMt_0_80_K0S->Draw("P");
  TGraphErrors* gr7p7CenMt_0_80_Proton = new TGraphErrors(13, pt_bin_center_Proton,v2_values_Proton , 0, v2_stat_error_Proton);
  gr7p7CenMt_0_80_Proton->SetMarkerStyle(20);
  gr7p7CenMt_0_80_Proton->SetMarkerColor(kAzure+2);
  gr7p7CenMt_0_80_Proton->SetLineWidth(2);
  gr7p7CenMt_0_80_Proton->SetLineColor(kAzure+2);
  gr7p7CenMt_0_80_Proton->SetMarkerSize(2);
  gr7p7CenMt_0_80_Proton->Draw("P");
  TGraphErrors* gr7p7CenMt_0_80_PiP = new TGraphErrors(7, pt_bin_center_PiP,v2_values_PiP , 0, v2_stat_error_PiP);
  gr7p7CenMt_0_80_PiP->SetMarkerStyle(30);
  gr7p7CenMt_0_80_PiP->SetMarkerColor(kMagenta+2);
  gr7p7CenMt_0_80_PiP->SetLineWidth(2);
  gr7p7CenMt_0_80_PiP->SetLineColor(kMagenta+2);
  gr7p7CenMt_0_80_PiP->SetMarkerSize(2);
  gr7p7CenMt_0_80_PiP->Draw("P");
  TGraphErrors* gr7p7CenMt_0_80_KP = new TGraphErrors(8, pt_bin_center_KP,v2_values_KP , 0, v2_stat_error_KP);
  gr7p7CenMt_0_80_KP->SetMarkerStyle(20);
  gr7p7CenMt_0_80_KP->SetMarkerColor(kGray+2);
  gr7p7CenMt_0_80_KP->SetLineWidth(2);
  gr7p7CenMt_0_80_KP->SetLineColor(kGray+2);
  gr7p7CenMt_0_80_KP->SetMarkerSize(2);
  gr7p7CenMt_0_80_KP->Draw("P");
  TGraphErrors *graph_v2_vs_mT_0_60_7p2 = new TGraphErrors(2, x_ptSetA, y_v2_0_60, 0, yErr_stat_v2_0_60);
  graph_v2_vs_mT_0_60_7p2->SetMarkerStyle(32);
  graph_v2_vs_mT_0_60_7p2->SetMarkerColor(kBlue);
  graph_v2_vs_mT_0_60_7p2->SetLineColor(kBlue);
  graph_v2_vs_mT_0_60_7p2->SetLineWidth(2);
  graph_v2_vs_mT_0_60_7p2->SetMarkerSize(2);
  graph_v2_vs_mT_0_60_7p2->Draw("P");
  TLine *line1_3 = new TLine(0, 0, 1.6, 0);
  line1_3->SetLineStyle(7);
  line1_3->Draw("same");
  TLegend *legend1_2 = new TLegend(0.4,0.65,0.9,0.9);
  legend1_2->AddEntry(graph_v2_vs_mT_0_60_7p2,"#phi FXT 7.2 GeV 0-60%","p");
  legend1_2->AddEntry(gr7p7CenMt_0_80,"#phi 7.7 GeV 0-80% - BES-I","p");
  // legend1_2->AddEntry(gr7p7CenMt_0_80_XiM,"#Xi^{-} 7.7 GeV 0-80% - BES-I","p");
  // legend1_2->AddEntry(gr7p7CenMt_0_80_OmegaM,"#Omega^{-} 7.7 GeV 0-80% - BES-I","p");
  legend1_2->AddEntry(gr7p7CenMt_0_80_Lambda,"#Lambda 7.7 GeV 0-80% - BES-I","p");
  legend1_2->AddEntry(gr7p7CenMt_0_80_K0S,"K_{s}^{0} 7.7 GeV 0-80% - BES-I","p");
  legend1_2->AddEntry(gr7p7CenMt_0_80_Proton,"p 7.7 GeV 0-80% - BES-I","p");
  legend1_2->AddEntry(gr7p7CenMt_0_80_PiP,"#pi^{+} 7.7 GeV 0-80% - BES-I","p");
  legend1_2->AddEntry(gr7p7CenMt_0_80_KP,"K^{+} 7.7 GeV 0-80% - BES-I","p");
  legend1_2->Draw("same");
  // ========================== (4) 4.5 GeV and above =======================
  // TCanvas *c1 = new TCanvas("c1","c1",200,10,1024,768);
  // c1->DrawFrame(3., -0.05, 250., 0.05);
  // TH2D * histTemp_4p5 = new TH2D("histTemp_4p5","histTemp_4p5",1000,3.,250,1000,-0.05,0.05);
  // histTemp_4p5->GetYaxis()->SetTitle("dv_{1}/dy");
  // histTemp_4p5->GetYaxis()->SetTitleOffset(1.4);
  // histTemp_4p5->GetXaxis()->SetTitle("#sqrt{s_{NN}} (GeV)");
  // histTemp_4p5->Draw();
  // //4.5 GeV phi
  // double x0[1]    = {4.5};
  // double zero0[1] = {0};
  // double py0[1]      = {-0.01301};
  // double ey_stat0[1] = {0.01566};
  // ey_stat0[0] *= sqrt(2);
  // double ey_sys0[1]  = {0.013012};
  // ey_sys0[0] *= sqrt(2);

  // ========================== (5) v1 systematic Ana 7.2 GeV =======================
  string sys_object[17]  = {"primary", "etaGap", "etaRange",
                            "vz", "vr", "dedx", "dca",
                            "nHitsFit", "ratio", "nSigK", "mass2",
                            "pT", "dipAngle", "vtxDiff", "mthdDiff",
                            "binning",
                            "TPCpid"
                          };
  // default [3]
  double sys_x0[3] = {-1.25-_y_CM, -0.75-_y_CM, -0.25-_y_CM};
  double sys_xErr0[3] = {0, 0, 0};
  double sys_y0[3] = {0.00678684, 0.0191681, 0.0203576};
  double sys_yErr0[3] = {0.0220404, 0.0129847, 0.0461443 };
  // etaGap 1
  double sys_x1_1[3] = {-1.25-_y_CM+0.01, -0.75-_y_CM+0.01, -0.25-_y_CM+0.01};
  double sys_xErr1_1[3] = {0, 0, 0};
  double sys_y1_1[3] = {0.000480216, 0.0175015, 0.0308636};
  double sys_yErr1_1[3] = {0.0232041, 0.0136474, 0.0482818 };
  // etaGap 2
  double sys_x1_2[3] = {-1.25-_y_CM+0.02, -0.75-_y_CM+0.02, -0.25-_y_CM+0.02};
  double sys_xErr1_2[3] = {0, 0, 0};
  double sys_y1_2[3] = {-0.00350079, 0.0205319, 0.0328723};
  double sys_yErr1_2[3] = {0.0220764, 0.0129927, 0.046171 };
  // etaRange 1
  double sys_x2_1[3] = {-1.25-_y_CM+0.03, -0.75-_y_CM+0.03, -0.25-_y_CM+0.03};
  double sys_xErr2_1[3] = {0, 0, 0};
  double sys_y2_1[3] = {-0.00556124, 0.00787956, 0.00998329};
  double sys_yErr2_1[3] = {0.0231609, 0.0135959, 0.0491899 };
  // vz 1
  double sys_x3_1[3] = {-1.25-_y_CM+0.04, -0.75-_y_CM+0.04, -0.25-_y_CM+0.04};
  double sys_xErr3_1[3] = {0, 0, 0};
  double sys_y3_1[3] = {0.00854142, 0.0197135, 0.0191546};
  double sys_yErr3_1[3] = {0.0221096, 0.0130339, 0.0462281 };
  // vz 2
  double sys_x3_2[3] = {-1.25-_y_CM+0.05, -0.75-_y_CM+0.05, -0.25-_y_CM+0.05};
  double sys_xErr3_2[3] = {0, 0, 0};
  double sys_y3_2[3] = {0.00421395, 0.0197106, 0.014646};
  double sys_yErr3_2[3] = {0.0221157, 0.013038, 0.046256 };
  // vr 1
  double sys_x4_1[3] = {-1.25-_y_CM+0.06, -0.75-_y_CM+0.06, -0.25-_y_CM+0.06};
  double sys_xErr4_1[3] = {0, 0, 0};
  double sys_y4_1[3] = {0.0068826, 0.0190863, 0.0217084};
  double sys_yErr4_1[3] = {0.0220416, 0.0129902, 0.0460079 };
  // vr 2
  double sys_x4_2[3] = {-1.25-_y_CM+0.07, -0.75-_y_CM+0.07, -0.25-_y_CM+0.07};
  double sys_xErr4_2[3] = {0, 0, 0};
  double sys_y4_2[3] = {0.00646351, 0.0190505, 0.013418};
  double sys_yErr4_2[3] = {0.0221013, 0.0130185, 0.0462275 };
  // dedx 1
  double sys_x5_1[3] = {-1.25-_y_CM+0.08, -0.75-_y_CM+0.08, -0.25-_y_CM+0.08};
  double sys_xErr5_1[3] = {0, 0, 0};
  double sys_y5_1[3] = {0.00556072, 0.0198019, 0.0254163};
  double sys_yErr5_1[3] = {0.0220475, 0.012998, 0.0462534 };
  // dedx 2
  double sys_x5_2[3] = {-1.25-_y_CM+0.09, -0.75-_y_CM+0.09, -0.25-_y_CM+0.09};
  double sys_xErr5_2[3] = {0, 0, 0};
  double sys_y5_2[3] = {0.00719935, 0.0189774, 0.0299739};
  double sys_yErr5_2[3] = {0.022779, 0.0132255, 0.0466932 };
  // dca 1
  double sys_x6_1[3] = {-1.25-_y_CM+0.1, -0.75-_y_CM+0.1, -0.25-_y_CM+0.1};
  double sys_xErr6_1[3] = {0, 0, 0};
  double sys_y6_1[3] = {0.00111875, 0.023637, 0.0396517};
  double sys_yErr6_1[3] = {0.0229959, 0.013812, 0.0505166 };
  // nHitsFit 1
  double sys_x7_1[3] = {-1.25-_y_CM+0.11, -0.75-_y_CM+0.11, -0.25-_y_CM+0.11};
  double sys_xErr7_1[3] = {0, 0, 0};
  double sys_y7_1[3] = {0.00675443, 0.0201602, 0.0176864};
  double sys_yErr7_1[3] = {0.0220832, 0.0130061, 0.0462062 };
  // nHitsFit 2
  double sys_x7_2[3] = {-1.25-_y_CM+0.12, -0.75-_y_CM+0.12, -0.25-_y_CM+0.12};
  double sys_xErr7_2[3] = {0, 0, 0};
  double sys_y7_2[3] = {0.00866891, 0.0180645, 0.01943};
  double sys_yErr7_2[3] = {0.0220924, 0.013033, 0.0462338 };
  // ratio 1
  double sys_x8_1[3] = {-1.25-_y_CM+0.13, -0.75-_y_CM+0.13, -0.25-_y_CM+0.13};
  double sys_xErr8_1[3] = {0, 0, 0};
  double sys_y8_1[3] = {0.00794904, 0.0193234, 0.0141304};
  double sys_yErr8_1[3] = {0.0220748, 0.0130031, 0.0462072 };
  // ratio 2
  double sys_x8_2[3] = {-1.25-_y_CM+0.14, -0.75-_y_CM+0.14, -0.25-_y_CM+0.14};
  double sys_xErr8_2[3] = {0, 0, 0};
  double sys_y8_2[3] = {0.00935401, 0.020488, 0.0207583};
  double sys_yErr8_2[3] = {0.0221918, 0.0130297, 0.0463469 };
  // nSigK 1
  double sys_x9_1[3] = {-1.25-_y_CM+0.15, -0.75-_y_CM+0.15, -0.25-_y_CM+0.15};
  double sys_xErr9_1[3] = {0, 0, 0};
  double sys_y9_1[3] = {0.00427665, 0.0203281, 0.0197152};
  double sys_yErr9_1[3] = {0.0231993, 0.0135066, 0.0479187 };
  // nSigK 2
  double sys_x9_2[3] = {-1.25-_y_CM+0.16, -0.75-_y_CM+0.16, -0.25-_y_CM+0.16};
  double sys_xErr9_2[3] = {0, 0, 0};
  double sys_y9_2[3] = {0.00635179, 0.0171804, 0.0158751};
  double sys_yErr9_2[3] = {0.0218175, 0.0126818, 0.0457916 };
  // mass2 1
  double sys_x10_1[3] = {-1.25-_y_CM+0.17, -0.75-_y_CM+0.17, -0.25-_y_CM+0.17};
  double sys_xErr10_1[3] = {0, 0, 0};
  double sys_y10_1[3] = {0.00172858, 0.0228337, 0.0183178};
  double sys_yErr10_1[3] = {0.0226613, 0.0130783, 0.0464455 };
  // mass2 2
  double sys_x10_2[3] = {-1.25-_y_CM+0.18, -0.75-_y_CM+0.18, -0.25-_y_CM+0.18};
  double sys_xErr10_2[3] = {0, 0, 0};
  double sys_y10_2[3] = {0.00412049, 0.0200491, 0.0227649};
  double sys_yErr10_2[3] = {0.02206, 0.0130793, 0.0464965 };
  // pT 1
  double sys_x11_1[3] = {-1.25-_y_CM+0.19, -0.75-_y_CM+0.19, -0.25-_y_CM+0.19};
  double sys_xErr11_1[3] = {0, 0, 0};
  double sys_y11_1[3] = {0.00741527, 0.0202824, 0.0236148};
  double sys_yErr11_1[3] = {0.0220722, 0.0129888, 0.0456992 };
  // pT 2
  double sys_x11_2[3] = {-1.25-_y_CM+0.2, -0.75-_y_CM+0.2, -0.25-_y_CM+0.2};
  double sys_xErr11_2[3] = {0, 0, 0};
  double sys_y11_2[3] = {0.00819094, 0.0169838, 0.0715593};
  double sys_yErr11_2[3] = {0.0223794, 0.0167388, 0.0702266 };
  // dipAngle 1
  double sys_x12_1[3] = {-1.25-_y_CM+0.21, -0.75-_y_CM+0.21, -0.25-_y_CM+0.21};
  double sys_xErr12_1[3] = {0, 0, 0};
  double sys_y12_1[3] = {0.00978963, 0.00676206, 0.00360617};
  double sys_yErr12_1[3] = {0.0173277, 0.0117513, 0.0450136 };
  // vtxDiff 1
  double sys_x13_1[3] = {-1.25-_y_CM+0.22, -0.75-_y_CM+0.22, -0.25-_y_CM+0.22};
  double sys_xErr13_1[3] = {0, 0, 0};
  double sys_y13_1[3] = {0.00676052, 0.0191058, 0.0203598};
  double sys_yErr13_1[3] = {0.0219867, 0.0129556, 0.0460639 };
  // vtxDiff 2
  double sys_x13_2[3] = {-1.25-_y_CM+0.23, -0.75-_y_CM+0.23, -0.25-_y_CM+0.23};
  double sys_xErr13_2[3] = {0, 0, 0};
  double sys_y13_2[3] = {0.0067717, 0.0191691, 0.0203616};
  double sys_yErr13_2[3] = {0.0220418, 0.0129849, 0.0461483 };
  // ----------------- TGEs ----------------------------------------
  TCanvas *c_sys = new TCanvas("c_sys","c_sys",200,10,1024,768);
  c_sys->DrawFrame(0, -0.04, 2.1, 0.15);
  c_sys->cd();
  // default TGraphErrors(3
  TGraphErrors *graph0 = new TGraphErrors(3, sys_x0, sys_y0, sys_xErr0, sys_yErr0);
  graph0->SetMarkerStyle(8);
  graph0->SetMarkerColor(kBlack);
  graph0->SetLineColor(kBlack);
  graph0->SetMarkerSize(2);
  graph0->Draw("P");
  // etaGap 1
  TGraphErrors *graph1_1 = new TGraphErrors(3, sys_x1_1, sys_y1_1, sys_xErr1_1, sys_yErr1_1);
  graph1_1->SetMarkerStyle(4);
  graph1_1->SetMarkerColor(kBlack);
  graph1_1->SetLineColor(kBlack);
  graph1_1->SetMarkerSize(2);
  graph1_1->Draw("P");
  // etaGap 2
  TGraphErrors *graph1_2 = new TGraphErrors(3, sys_x1_2, sys_y1_2, sys_xErr1_2, sys_yErr1_2);
  graph1_2->SetMarkerStyle(32);
  graph1_2->SetMarkerColor(kBlack);
  graph1_2->SetLineColor(kBlack);
  graph1_2->SetMarkerSize(2);
  graph1_2->Draw("P");
  // etaRange 1
  TGraphErrors *graph2_1 = new TGraphErrors(3, sys_x2_1, sys_y2_1, sys_xErr2_1, sys_yErr2_1);
  graph2_1->SetMarkerStyle(30);
  graph2_1->SetMarkerColor(kBlack);
  graph2_1->SetLineColor(kBlack);
  graph2_1->SetMarkerSize(2);
  graph2_1->Draw("P");
  // vz 1
  TGraphErrors *graph3_1 = new TGraphErrors(3, sys_x3_1, sys_y3_1, sys_xErr3_1, sys_yErr3_1);
  graph3_1->SetMarkerStyle(28);
  graph3_1->SetMarkerColor(kBlack);
  graph3_1->SetLineColor(kBlack);
  graph3_1->SetMarkerSize(2);
  graph3_1->Draw("P");
  // vz 2
  TGraphErrors *graph3_2 = new TGraphErrors(3, sys_x3_2, sys_y3_2, sys_xErr3_2, sys_yErr3_2);
  graph3_2->SetMarkerStyle(27);
  graph3_2->SetMarkerColor(kBlack);
  graph3_2->SetLineColor(kBlack);
  graph3_2->SetMarkerSize(2);
  graph3_2->Draw("P");
  // vr 1
  TGraphErrors *graph4_1 = new TGraphErrors(3, sys_x4_1, sys_y4_1, sys_xErr4_1, sys_yErr4_1);
  graph4_1->SetMarkerStyle(26);
  graph4_1->SetMarkerColor(kBlack);
  graph4_1->SetLineColor(kBlack);
  graph4_1->SetMarkerSize(2);
  graph4_1->Draw("P");
  // vr 2
  TGraphErrors *graph4_2 = new TGraphErrors(3, sys_x4_2, sys_y4_2, sys_xErr4_2, sys_yErr4_2);
  graph4_2->SetMarkerStyle(25);
  graph4_2->SetMarkerColor(kBlack);
  graph4_2->SetLineColor(kBlack);
  graph4_2->SetMarkerSize(2);
  graph4_2->Draw("P");
  // dedx 1
  TGraphErrors *graph5_1 = new TGraphErrors(3, sys_x5_1, sys_y5_1, sys_xErr5_1, sys_yErr5_1);
  graph5_1->SetMarkerStyle(39);
  graph5_1->SetMarkerColor(kBlack);
  graph5_1->SetLineColor(kBlack);
  graph5_1->SetMarkerSize(2);
  graph5_1->Draw("P");
  // dedx 2
  TGraphErrors *graph5_2 = new TGraphErrors(3, sys_x5_2, sys_y5_2, sys_xErr5_2, sys_yErr5_2);
  graph5_2->SetMarkerStyle(46);
  graph5_2->SetMarkerColor(kBlack);
  graph5_2->SetLineColor(kBlack);
  graph5_2->SetMarkerSize(2);
  graph5_2->Draw("P");
  // dca 1
  TGraphErrors *graph6_1 = new TGraphErrors(3, sys_x6_1, sys_y6_1, sys_xErr6_1, sys_yErr6_1);
  graph6_1->SetMarkerStyle(44);
  graph6_1->SetMarkerColor(kBlack);
  graph6_1->SetLineColor(kBlack);
  graph6_1->SetMarkerSize(2);
  graph6_1->Draw("P");
  // nHitsFit 1
  TGraphErrors *graph7_1 = new TGraphErrors(3, sys_x7_1, sys_y7_1, sys_xErr7_1, sys_yErr7_1);
  graph7_1->SetMarkerStyle(42);
  graph7_1->SetMarkerColor(kBlack);
  graph7_1->SetLineColor(kBlack);
  graph7_1->SetMarkerSize(2);
  graph7_1->Draw("P");
  // nHitsFit 2
  TGraphErrors *graph7_2 = new TGraphErrors(3, sys_x7_2, sys_y7_2, sys_xErr7_2, sys_yErr7_2);
  graph7_2->SetMarkerStyle(40);
  graph7_2->SetMarkerColor(kBlack);
  graph7_2->SetLineColor(kBlack);
  graph7_2->SetMarkerSize(2);
  graph7_2->Draw("P");
  // ratio 1
  TGraphErrors *graph8_1 = new TGraphErrors(3, sys_x8_1, sys_y8_1, sys_xErr8_1, sys_yErr8_1);
  graph8_1->SetMarkerStyle(38);
  graph8_1->SetMarkerColor(kBlack);
  graph8_1->SetLineColor(kBlack);
  graph8_1->SetMarkerSize(2);
  graph8_1->Draw("P");
  // ratio 2
  TGraphErrors *graph8_2 = new TGraphErrors(3, sys_x8_2, sys_y8_2, sys_xErr8_2, sys_yErr8_2);
  graph8_2->SetMarkerStyle(37);
  graph8_2->SetMarkerColor(kBlack);
  graph8_2->SetLineColor(kBlack);
  graph8_2->SetMarkerSize(2);
  graph8_2->Draw("P");
  // nSigK 1
  TGraphErrors *graph9_1 = new TGraphErrors(3, sys_x9_1, sys_y9_1, sys_xErr9_1, sys_yErr9_1);
  graph9_1->SetMarkerStyle(36);
  graph9_1->SetMarkerColor(kBlack);
  graph9_1->SetLineColor(kBlack);
  graph9_1->SetMarkerSize(2);
  graph9_1->Draw("P");
  // nSigK 2
  TGraphErrors *graph9_2 = new TGraphErrors(3, sys_x9_2, sys_y9_2, sys_xErr9_2, sys_yErr9_2);
  graph9_2->SetMarkerStyle(35);
  graph9_2->SetMarkerColor(kBlack);
  graph9_2->SetLineColor(kBlack);
  graph9_2->SetMarkerSize(2);
  graph9_2->Draw("P");
  // mass2 1
  TGraphErrors *graph10_1 = new TGraphErrors(3, sys_x10_1, sys_y10_1, sys_xErr10_1, sys_yErr10_1);
  graph10_1->SetMarkerStyle(34);
  graph10_1->SetMarkerColor(kBlack);
  graph10_1->SetLineColor(kBlack);
  graph10_1->SetMarkerSize(2);
  graph10_1->Draw("P");
  // mass2 2
  TGraphErrors *graph10_2 = new TGraphErrors(3, sys_x10_2, sys_y10_2, sys_xErr10_2, sys_yErr10_2);
  graph10_2->SetMarkerStyle(33);
  graph10_2->SetMarkerColor(kBlack);
  graph10_2->SetLineColor(kBlack);
  graph10_2->SetMarkerSize(2);
  graph10_2->Draw("P");
  // pT 1
  TGraphErrors *graph11_1 = new TGraphErrors(3, sys_x11_1, sys_y11_1, sys_xErr11_1, sys_yErr11_1);
  graph11_1->SetMarkerStyle(23);
  graph11_1->SetMarkerColor(kBlack);
  graph11_1->SetLineColor(kBlack);
  graph11_1->SetMarkerSize(2);
  graph11_1->Draw("P");
  // pT 2
  TGraphErrors *graph11_2 = new TGraphErrors(3, sys_x11_2, sys_y11_2, sys_xErr11_2, sys_yErr11_2);
  graph11_2->SetMarkerStyle(22);
  graph11_2->SetMarkerColor(kBlack);
  graph11_2->SetLineColor(kBlack);
  graph11_2->SetMarkerSize(2);
  graph11_2->Draw("P");
  // dipAngle 1
  TGraphErrors *graph12_1 = new TGraphErrors(3, sys_x12_1, sys_y12_1, sys_xErr12_1, sys_yErr12_1);
  graph12_1->SetMarkerStyle(21);
  graph12_1->SetMarkerColor(kBlack);
  graph12_1->SetLineColor(kBlack);
  graph12_1->SetMarkerSize(2);
  graph12_1->Draw("P");
  // vtxDiff 1
  TGraphErrors *graph13_1 = new TGraphErrors(3, sys_x13_1, sys_y13_1, sys_xErr13_1, sys_yErr13_1);
  graph13_1->SetMarkerStyle(48);
  graph13_1->SetMarkerColor(kBlack);
  graph13_1->SetLineColor(kBlack);
  graph13_1->SetMarkerSize(2);
  graph13_1->Draw("P");
  // vtxDiff 2
  TGraphErrors *graph13_2 = new TGraphErrors(3, sys_x13_2, sys_y13_2, sys_xErr13_2, sys_yErr13_2);
  graph13_2->SetMarkerStyle(49);
  graph13_2->SetMarkerColor(kBlack);
  graph13_2->SetLineColor(kBlack);
  graph13_2->SetMarkerSize(2);
  graph13_2->Draw("P");

  TLegend *legend_sys = new TLegend(0.1,0.5,0.78,0.9);
  legend_sys->SetFillStyle(0);
  legend_sys->SetBorderSize(0);
  legend_sys->SetNColumns(3);
  legend_sys->AddEntry(graph0,Form("Default                   dv_{1}/dy = %.5f %c %.5f", 0.0142375,177,0.00983601),"p");
  legend_sys->AddEntry(graph2_1,Form("#eta range [-4.95,-4.35]  dv_{1}/dy = %.5f %c %.5f", 0.00465417,177,0.0100203),"p");
  legend_sys->AddEntry(graph12_1,Form("Dip-angle > 0            dv_{1}/dy = %.5f %c %.5f", 0.00636387,177,0.00858689),"p");
  legend_sys->AddEntry(graph1_1,Form("#eta gap 0.15             dv_{1}/dy = %.5f %c %.5f", 0.0121755,177,0.0102094),"p");
  legend_sys->AddEntry(graph1_2,Form("#eta gap 0.1              dv_{1}/dy = %.5f %c %.5f", 0.0137032,177,0.00976577),"p");
  legend_sys->AddEntry(graph3_1,Form("198.4 < vz < 201.6 cm     dv_{1}/dy = %.5f %c %.5f", 0.014874,177,0.0099082),"p");
  legend_sys->AddEntry(graph3_2,Form("197.6 < vz < 202.4 cm     dv_{1}/dy = %.5f %c %.5f", 0.0142309,177,0.00985455),"p");
  legend_sys->AddEntry(graph4_1,Form("#sqrt{vx^{2} + (vy+2)^{2}} < 1.6 cm dv_{1}/dy = %.5f %c %.5f", 0.014195,177,0.00983872),"p");
  legend_sys->AddEntry(graph4_2,Form("#sqrt{vx^{2} + (vy+2)^{2}} < 2.4 cm dv_{1}/dy = %.5f %c %.5f", 0.0141089,177,0.0098533),"p");
  legend_sys->AddEntry(graph5_1,Form("nHitsDedx > 10            dv_{1}/dy = %.5f %c %.5f", 0.0144912,177,0.00984469),"p");
  legend_sys->AddEntry(graph5_2,Form("nHitsDedx > 20            dv_{1}/dy = %.5f %c %.5f", 0.0141844,177,0.0100262),"p");
  legend_sys->AddEntry(graph6_1,Form("DCA < 1                   dv_{1}/dy = %.5f %c %.5f", 0.016428,177,0.010445),"p");
  legend_sys->AddEntry(graph7_1,Form("nHitsFit > 10             dv_{1}/dy = %.5f %c %.5f", 0.0149157,177,0.00987638),"p");
  legend_sys->AddEntry(graph7_2,Form("nHitsFit > 20             dv_{1}/dy = %.5f %c %.5f", 0.0137583,177,0.0098668),"p");
  legend_sys->AddEntry(graph8_1,Form("#frac{nHitsFit}{nHitsPoss} > 0.45 dv_{1}/dy = %.5f %c %.5f", 0.0145176,177,0.00986885),"p");
  legend_sys->AddEntry(graph8_2,Form("#frac{nHitsFit}{nHitsPoss} > 0.55 dv_{1}/dy = %.5f %c %.5f", 0.015532,177,0.00994258),"p");
  legend_sys->AddEntry(graph9_1,Form("|n_{#sigma_K}| < 1.8      dv_{1}/dy = %.5f %c %.5f", 0.0146922,177,0.0102212),"p");
  legend_sys->AddEntry(graph9_2,Form("|n_{#sigma_K}| < 2.2      dv_{1}/dy = %.5f %c %.5f", 0.0128178,177,0.0095815),"p");
  legend_sys->AddEntry(graph10_1,Form("0.17<m^{2}<0.31          dv_{1}/dy = %.5f %c %.5f", 0.0160833,177,0.00995992),"p");
  legend_sys->AddEntry(graph10_2,Form("0.15<m^{2}<0.33          dv_{1}/dy = %.5f %c %.5f", 0.0144387,177,0.00988421),"p");
  legend_sys->AddEntry(graph11_1,Form("p_{T} > 0                dv_{1}/dy = %.5f %c %.5f", 0.0150988,177,0.00987733),"p");
  legend_sys->AddEntry(graph11_2,Form("p_{T} > 0.4              dv_{1}/dy = %.5f %c %.5f", 0.0127883,177,0.0120877),"p");
  legend_sys->AddEntry(graph13_1,Form("#Delta vz_{ME} < 0.2     dv_{1}/dy = %.5f %c %.5f", 0.0141903,177,0.00981312),"p");
  legend_sys->AddEntry(graph13_2,Form("#Delta vz_{ME} < 1.0     dv_{1}/dy = %.5f %c %.5f", 0.0142359,177,0.00983609),"p");
  legend_sys->Draw("same");
  TLine *line3 = new TLine(0, 0, 2.1, 0);
  line3->SetLineStyle(7);
  line3->Draw("same");
  TLine *line4 = new TLine(2.03, -0.04, 2.03, 0.15);
  line4->SetLineStyle(7);
  line4->Draw("same");

  // ========================== (6) v1 vs. pT at 7.2 GeV ========================
  TCanvas *TCanv_v1_vs_pt = new TCanvas("TCanv_v1_vs_pt","v_{1} vs. pt at 7.2 GeV",200,10,1024,768);
  TCanv_v1_vs_pt->DrawFrame(0, -0.04, 2.1, 0.1);
  // 7.2 GeV v1 vs. pt
  double x_ptSetA_default[2]    = { 0.9, 1.8};
  double xErr_ptSetA_default[2] = {0.3, 0.6};
  double y_v1_ptSetA[2]      = {0.00434645, 0.026056};
  double yErr_stat_v1_ptSetA[2] = {0.0182416, 0.0184488};
  TCanv_v1_vs_pt->cd();
  TH2D * histTemp6 = new TH2D("histTemp6","histTemp6",1000,0.5,2.5,1000,-0.02,0.05);
  histTemp6->GetYaxis()->SetTitle("v_{1}");
  histTemp6->GetYaxis()->SetTitleOffset(1.4);
  histTemp6->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  histTemp6->Draw();
  // 7.2 GeV
  TGraphErrors *graph_v1_vs_ptSetA_7p2 = new TGraphErrors(2, x_ptSetA_default, y_v1_ptSetA, xErr_ptSetA_default, yErr_stat_v1_ptSetA);
  graph_v1_vs_ptSetA_7p2->SetMarkerStyle(30);
  graph_v1_vs_ptSetA_7p2->SetMarkerColor(kBlue);
  graph_v1_vs_ptSetA_7p2->SetLineColor(kBlue);
  graph_v1_vs_ptSetA_7p2->SetMarkerSize(2);
  graph_v1_vs_ptSetA_7p2->Draw("P");
  TLine *tl_6 = new TLine(0.5, 0, 2.5, 0);
  tl_6->SetLineStyle(7);
  tl_6->Draw("same");
  TPaveText * ptext_ptSetA = new TPaveText(0.12,0.65,0.5,0.85,"NDCARC");
  ptext_ptSetA->SetFillColor(0);
  ptext_ptSetA -> AddText("STAR Au+Au #sqrt{s_{NN}} 7.2 GeV FXT 10-40%");
  ptext_ptSetA -> AddText("#phi meson");
  ptext_ptSetA -> AddText("Integral all y");
  ptext_ptSetA -> AddText("y_{target} = 2.02  ");
  // ptext_ptSetA -> AddText("STAR preliminary");
  TLegend *leg_7p2_v1_phi_ptSetA = new TLegend(0.25,0.5,0.75,0.8);
  leg_7p2_v1_phi_ptSetA->SetFillStyle(0);
  leg_7p2_v1_phi_ptSetA->SetBorderSize(0);
  // leg_7p2_v1_phi_ptSetA->AddEntry((TObject*)0,Form("dv_{1}/dy = %.4f %c %.4f %c %.4f",0.0138904/*(Double_t)tf1_dv1dy->GetParameter(0)*/,177,0.00916586/*(Double_t)tf1_dv1dy->GetParError(0)*/,177,0.00647123),"");
  // leg_7p2_v1_phi_ptSetA->Draw("same");
  ptext_ptSetA->Draw("same");

  // ========================== (7) v2 2 bin jackknife 7.2 GeV ========================
  TCanvas *TCanv_jkk = new TCanvas("TCanv_jkk","TCanv_jkk",200,10,1024,768);
  TCanv_jkk->Divide(2,1);
  Double_t jkk0[10] = {-0.0224922, -0.0222856, -0.0099591, 0.000997738,
    -0.0153799, -0.0036787, -0.0278074, -0.00608307, -0.0134696, 0.0238224};
  Double_t jkk1[10] = {-0.0140789, -0.0389647, -0.0166912, -0.0300542,
    -0.0187818, -0.0310772, -0.0184702, -0.00846529, -0.0104781, -0.0182829};
  TH1D *h_v2_bin1 = new TH1D("h_v2_bin1","h_v2_bin1",8, -0.05,0.05);
  TH1D *h_v2_bin2 = new TH1D("h_v2_bin2","h_v2_bin2",8, -0.05,0.05);
  for(int i=0;i<10;i++){
    h_v2_bin1->Fill(jkk0[i]);
    h_v2_bin2->Fill(jkk1[i]);
  }
  TCanv_jkk->cd(1);
  h_v2_bin1->Draw();
  TCanv_jkk->cd(2);
  h_v2_bin2->Draw();

  // ========================== (8) v2 2 bin jackknife 7.2 GeV ========================
  TCanvas *TCanv_jkk_3bin = new TCanvas("TCanv_jkk_3bin","TCanv_jkk_3bin",200,10,1024,768);
  TCanv_jkk_3bin->Divide(3,1);
  Double_t jkk0_3bin[10] = {-0.0341356, -0.0164558, -0.0147401, -0.044426, -0.031164,
    -0.023199, -0.0320284, -0.0362466, -0.0301097, -0.0212607};
  Double_t jkk1_3bin[10] =  {-0.0111124, 0.003465, -0.00825907, -0.015967, -0.0155545,
    -0.00968208, 0.00260028, -0.0251693, -0.0074882, -0.00145464};
  Double_t jkk2_3bin[10] =   {0.0437703, 0.165346, 0.0167534, 0.0965503, 0.095658,
    0.297014, 0.0773467, 0.180491, 0.105279, 0.0566139};

  TH1D *h_v2_jkk_bin1 = new TH1D("h_v2_jkk_bin1","h_v2_jkk_bin1",8, -0.05,0.05);
  TH1D *h_v2_jkk_bin2 = new TH1D("h_v2_jkk_bin2","h_v2_jkk_bin2",8, -0.05,0.05);
  TH1D *h_v2_jkk_bin3 = new TH1D("h_v2_jkk_bin3","h_v2_jkk_bin3",16, -0.3,0.3);
  for(int i=0;i<10;i++){
    h_v2_jkk_bin1->Fill(jkk0_3bin[i]);
    h_v2_jkk_bin2->Fill(jkk1_3bin[i]);
    h_v2_jkk_bin3->Fill(jkk2_3bin[i]);
  }
  TCanv_jkk_3bin->cd(1);
  h_v2_jkk_bin1->Draw();
  TCanv_jkk_3bin->cd(2);
  h_v2_jkk_bin2->Draw();
  TCanv_jkk_3bin->cd(3);
  h_v2_jkk_bin3->Draw();

  // ========================== (9) v1 3 bin jackknife 7.2 GeV ========================
  TCanvas *TCanv_jkk_rap = new TCanvas("TCanv_jkk_rap","TCanv_jkk_rap",200,10,1024,768);
  TCanv_jkk_rap->Divide(3,1);
  Double_t jkk0_rap[10] = {-0.0153919,0.0103276,-0.00969343,-0.00250536,-0.000626034,
    0.00819924,0.000214182,0.00193367,0.00878884,-0.00066745};
  Double_t jkk1_rap[10] =  {0.0272558,0.0178734,0.0269522,0.0329685,0.0317695,
    0.035321,0.0295679,0.0378912,0.0301462,0.0254561};
  Double_t jkk2_rap[10] =   {0.0716877,0.111469,0.0604407,0.0561723,0.0998641,
    0.109151,0.0839622,0.109282,0.057333,0.0915568};

  TH1D *h_v1_jkk_bin1 = new TH1D("h_v1_jkk_bin1","h_v1_jkk_bin1",8, -0.05,0.05);
  TH1D *h_v1_jkk_bin2 = new TH1D("h_v1_jkk_bin2","h_v1_jkk_bin2",8, -0.05,0.05);
  TH1D *h_v1_jkk_bin3 = new TH1D("h_v1_jkk_bin3","h_v1_jkk_bin3",8, 0.0,0.15);
  for(int i=0;i<10;i++){
    h_v1_jkk_bin1->Fill(jkk0_rap[i]);
    h_v1_jkk_bin2->Fill(jkk1_rap[i]);
    h_v1_jkk_bin3->Fill(jkk2_rap[i]);
  }
  TCanv_jkk_rap->cd(1);
  h_v1_jkk_bin1->Draw();
  TCanv_jkk_rap->cd(2);
  h_v1_jkk_bin2->Draw();
  TCanv_jkk_rap->cd(3);
  h_v1_jkk_bin3->Draw();
}

Double_t proportion(Double_t *x, Double_t *p)
{
  return ( p[0] * x[0]);
}

Double_t Pol3plus1(Double_t *x, Double_t *p)
{
  return ( p[0] * x[0] + p[1] * x[0] * x[0] * x[0]);
}

Double_t d_v2nq(Double_t d_v2, Int_t ncq ){
  Double_t v2nq = d_v2/(Double_t)ncq;
  return v2nq;
}
Double_t d_mTm0nq(Double_t d_pT, Int_t ncq, Double_t mass ){ // pT to (mT - m0)/nq
  Double_t mTm0nq =  (sqrt(mass * mass + d_pT * d_pT) - mass) / (Double_t)ncq;
  return mTm0nq;
}
