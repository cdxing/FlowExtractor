// From Guannan
// Date: Oct 14, 2020
TLatex* drawLatex(Double_t x, Double_t y, char* text, Int_t textFont, Double_t textSize, Int_t colorIndex)
{
   TLatex *latex = new TLatex(x, y, text);
   latex->SetNDC();
   latex->SetTextFont(textFont);
   latex->SetTextSize(textSize);
   latex->SetTextColor(colorIndex);
   latex->Draw("same");
   return latex;
}

void plot7p7Slop()
{
   //fit and plot
   TCanvas *c1 = new TCanvas("c1", "c1", 10, 10, 800, 800);
   gStyle->SetEndErrorSize(0.01);

   gStyle->SetOptStat(0);
   gStyle->SetOptFit(1);

   double x1[1] = {0};
   double y1[1] = {0 };
   double x2[1] = {1 };
   double y2[1] = {1 };
   double left[1] = {0.2  };
   double right[1] = { 0.02 };
   double top[1] = {0.06   };
   double bottom[1] = {0.16 };

   for (int i = 0; i < 1; i++)
   {
      c1->SetFillStyle(4000);
      c1->SetFillColor(10);
      c1->SetBorderMode(0);
      c1->SetBorderSize(0);
      c1->SetFrameFillColor(0);
      c1->SetFrameBorderMode(0);
      c1->SetFrameBorderSize(0);
      c1->SetFrameLineWidth(1);
      c1->SetLeftMargin(left[i]);
      c1->SetRightMargin(right[i]);
      c1->SetTopMargin(top[i]);
      c1->SetBottomMargin(bottom[i]);
      c1->SetTickx();
      c1->SetTicky();
      c1->Draw();
   }
   // # xy scan the data point
   // # Format: x y -dx +dx -dy +dy
   // -0.459459 -0.0173913  0.128378  0.128378  0.0275362 0.0275362
   // -0.148649 0.00217391  0.121622  0.121622  0.0210145 0.0210145
   // 0.148649  -0.00797101 0.121622  0.121622  0.0195652 0.0195652
   // 0.445946  0.0057971 0.121622  0.121622  0.0271739 0.0271739


   double yy[4] = { -0.45, -0.15, 0.15, 0.45};
   double yyerr[4] = {0.15, 0.15, 0.15, 0.15};
   double v11[4] = { -0.0173913, 0.00217391, -0.00797101, 0.0057971};
   double v11err[4] = {0.0275362, 0.0210145, 0.0195652 , 0.0271739};

   TGraphErrors* gr7p7Cen_v11_slop;
   // gr7p7Cen_v11_slop = new TGraphErrors(4, yy, v11, yyerr, v11err);
   gr7p7Cen_v11_slop = new TGraphErrors(4, yy, v11, 0, v11err);


   TH1D* h0 = new TH1D("h1", "", 100, -10, 10);
   h0->GetYaxis()->SetRangeUser(-0.06, 0.06);
   h0->SetTitle("");
   h0->GetXaxis()->SetTitle("y");
   h0->GetYaxis()->SetTitle("v_{1}");
   h0->GetXaxis()->SetRangeUser(-1., +1.);
   h0->GetXaxis()->CenterTitle();
   h0->GetYaxis()->CenterTitle();
   h0->GetXaxis()->SetTitleSize(0.06);
   h0->GetXaxis()->SetLabelOffset(0.005);
   h0->GetXaxis()->SetLabelSize(0.05);
   h0->GetXaxis()->SetLabelFont(42);
   h0->GetXaxis()->SetTitleFont(42);
   h0->GetXaxis()->SetTitleOffset(0.9);
   h0->GetYaxis()->SetTitleOffset(1.2);
   h0->GetYaxis()->SetTitleSize(0.06);
   h0->GetYaxis()->SetLabelSize(0.05);
   h0->GetYaxis()->SetLabelOffset(0.018);
   h0->GetYaxis()->SetLabelFont(42);
   h0->GetYaxis()->SetTitleFont(42);
   h0->GetXaxis()->SetNdivisions(505);
   h0->GetYaxis()->SetNdivisions(505);
   h0->Draw();


   gr7p7Cen_v11_slop->SetMarkerStyle(20);
   gr7p7Cen_v11_slop->SetMarkerSize(0.8);
   gr7p7Cen_v11_slop->SetMarkerColor(2);
   gr7p7Cen_v11_slop->SetLineWidth(2);
   gr7p7Cen_v11_slop->SetLineColor(2);
   gr7p7Cen_v11_slop->Draw("Psame");

   TF1* f1 = new TF1("f1","[0]+[1]*x",-1,1);
   f1->SetParNames("p0","p1");
   TF1* f11 = new TF1("f11", "[0]*x", -1, 1);
   f11->SetParNames("p0");
   gr7p7Cen_v11_slop->Fit("f1", "0+", "", -0.6, 0.6);
   gr7p7Cen_v11_slop->Fit("f11", "0+", "", -0.6, 0.6);
   f1->Draw("same");
   f11->Draw("same");

   drawLatex(0.22, 0.88, "STAR Au+Au #sqrt{s_{NN}} 7.7GeV", 42, 0.05, 1);
   drawLatex(0.32, 0.80, Form("#phi meson"), 42, 0.042, 4);
   drawLatex(0.32, 0.65, Form("p0 = %3.3f#pm%3.3f", f1->GetParameter(0), f1->GetParError(0)), 132, 0.042, 2);


   TLegend* legend2 = new TLegend(0.73, .65, .93, .80, Form("  v_{1}    "));
   legend2->SetTextFont(42);
   legend2->SetTextSize(0.042);
   legend2->SetBorderSize(0);
   legend2->SetFillStyle(0);
   legend2->AddEntry(gr7p7Cen_v11_slop, " #phi", "pl");
   legend2->Draw("same");


   c1->SaveAs(Form("7p7_phi_v1Slop11.png"));


}
