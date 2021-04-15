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

void v2_bkg(){
  gStyle->SetOptStat(0);
  gStyle->SetOptDate(0);
  gStyle->SetEndErrorSize(6);
  // gStyle->SetOptTitle(0);
  TCanvas *can_v2_bkg = new TCanvas("can_v2_bkg","v_{2} vs. pT",200,10,1024,768);
  can_v2_bkg->DrawFrame(2., -0.06, 250., 0.2);
  // 7.2 GeV v2 vs. pT 10-40%
  double x_ptSetA[2]    = { 0.9,1.8};
  double xErr_ptSetA[2] = {0.3,0.6};

  // double y_v2_10_40[2]      = {-0.00966106, -0.0205263}; // w bin
  // double yErr_stat_v2_10_40[2] = {0.014085270558323343, 0.009257925680356804};

  double y_v2_10_40[2]      = {0.0199758, -0.0232837}; // w bin
  double yErr_stat_v2_10_40[2] = {0.0277097, 0.0285158};

  // double y_v2_bkg_10_40[2]      = {0.0113, 0.0106}; // w bin
  // double yErr_stat_v2_bkg_10_40[2] = {0.0028, 0.0027};
  double y_v2_bkg_10_40[2]      = {0.0046, 0.0034}; // w bin
  double yErr_stat_v2_bkg_10_40[2] = {0.0027, 0.0027};

  // ---------------- Particle species: KM ----------------
  Double_t pt_bin_center_KM[7] = {0.315,0.495,0.705,0.885,1.155,1.545,1.935};
  Double_t v2_values_KM[7] = {0.00187643,0.027042,0.0401026,0.0520291,0.0636185,0.0922199,0.0851067};
  Double_t v2_stat_error_KM[7] = {0.00416041,0.00282332,0.00306204,0.00390723,0.0043618,0.00900285,0.0219624};
  Double_t v2_syst_low_error_KM[7] = {3.75182e-05,2.14134e-05,7.42303e-05,0.000240313,0.000298115,0.00046376,0.00201082};
  Double_t v2_syst_high_error_KM[7] = {1.87695e-05,1.2589e-05,0.000140415,0.000293047,0.000149896,0.00028154,0.0022204};

  // ---------------- Particle species: KP ----------------
  Double_t pt_bin_center_KP[7] = {0.315,0.495,0.705,0.885,1.155,1.545,1.965};
  Double_t v2_values_KP[7] = {0.0123349,0.0284567,0.0493235,0.0674997,0.080886,0.100798,0.0878592};
  Double_t v2_stat_error_KP[7] = {0.00259688,0.00169075,0.00177886,0.00220443,0.00239161,0.00478715,0.0106033};
  Double_t v2_syst_low_error_KP[7] = {0.000213675,1.68026e-05,8.83694e-05,9.47739e-05,5.02883e-05,0.000454596,0.00393434};
  Double_t v2_syst_high_error_KP[7] = {0.000181583,3.12802e-05,5.80217e-05,0.000169471,8.27577e-05,0.000264862,0.00378213};
  // ---------------- Particle species: Phi ----------------
  Double_t pt_bin_center_10_40[1] = {0.9945};
  Double_t pt_bin_center_errCen_10_40[1] = {0};
  Double_t v2_values_10_40[1] = {0.0399507};
  Double_t v2_stat_error_10_40[1] = {0.0193815};
  TH2D * temp = new TH2D("temp","10-40% Au+Au",1000,0,2.5,1000,-0.04,0.12);
  temp->GetYaxis()->SetTitle("v_{2}");
  temp->GetYaxis()->SetTitleOffset(1.2);
  temp->GetXaxis()->SetTitle("pT (GeV/c)");
  temp->Draw();
  TGraphErrors* gr = new TGraphErrors(2, x_ptSetA, y_v2_bkg_10_40, xErr_ptSetA, yErr_stat_v2_bkg_10_40);
  gr->SetMarkerStyle(23);
  gr->SetMarkerColor(kBlack);
  gr->SetLineColor(kBlack);
  gr->SetMarkerSize(2);
  gr->Draw("P");
  TGraphErrors* gr1 = new TGraphErrors(7, pt_bin_center_KM, v2_values_KM, 0, v2_stat_error_KM);
  gr1->SetMarkerStyle(20);
  gr1->SetMarkerColor(kRed+1);
  gr1->SetLineColor(kRed+1);
  gr1->SetMarkerSize(2);
  gr1->Draw("P");
  TGraphErrors* gr2 = new TGraphErrors(7, pt_bin_center_KP, v2_values_KP, 0, v2_stat_error_KP);
  gr2->SetMarkerStyle(20);
  gr2->SetMarkerColor(kBlue+1);
  gr2->SetLineColor(kBlue+1);
  gr2->SetMarkerSize(2);
  gr2->Draw("P");
  TGraphErrors* gr3 = new TGraphErrors(2, x_ptSetA, y_v2_10_40, xErr_ptSetA, yErr_stat_v2_10_40);
  gr3->SetMarkerStyle(30);
  gr3->SetMarkerColor(kBlue);
  gr3->SetLineColor(kBlue);
  gr3->SetMarkerSize(2);
  gr3->Draw("P");
  TGraphErrors* gr7p7Cen_10_40 = new TGraphErrors(1, pt_bin_center_10_40, v2_values_10_40, pt_bin_center_errCen_10_40, v2_stat_error_10_40);
  gr7p7Cen_10_40->SetMarkerStyle(23);
  gr7p7Cen_10_40->SetMarkerColor(kGreen+2);
  gr7p7Cen_10_40->SetLineColor(kGreen+2);
  gr7p7Cen_10_40->SetMarkerSize(2);
  gr7p7Cen_10_40->Draw("P");
  TLegend *leg = new TLegend(0.1,0.65,0.5,0.9);
  leg->AddEntry(gr3,"7.2 GeV #phi","p");
  leg->AddEntry(gr,"7.2 GeV K^{+}K^{-} background","p");
  leg->AddEntry(gr7p7Cen_10_40,"7.7 GeV #phi","p");
  leg->AddEntry(gr1,"7.7 GeV K^{-}","p");
  leg->AddEntry(gr2,"7.7 GeV K^{+}","p");
  leg->Draw("same");
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  TLine *line1 = new TLine(0, 0, 2.5, 0);
  line1->SetLineStyle(7);
  line1->Draw("same");
  can_v2_bkg->SaveAs("v2_bkg.pdf");
}
