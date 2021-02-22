#include <TStyle.h>


void flowPlotter_4p5GeV()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptDate(0);
  // gStyle->SetEndErrorSize(6);
  gStyle->SetOptTitle(0);


    TCanvas *c1 = new TCanvas("c1","c1",200,10,1024,768);
    c1->DrawFrame(3., -0.05, 40., 0.06);
    // TH2D * histTemp = new TH2D("histTemp","histTemp",1000,3.,250,1000,-0.05,0.05);
    TH2D * histTemp = new TH2D("histTemp","histTemp",38,3.,41,1000,-0.05,0.06);
    histTemp->GetYaxis()->SetTitle("Directed Flow Slope dv_{1} /dy|_{y=0}");
    histTemp->GetYaxis()->SetTitleOffset(1.4);
    histTemp->GetXaxis()->SetTitleOffset(1.3);
    // histTemp->GetXaxis()->SetNdivisions(438);
    histTemp->GetXaxis()->SetBinLabel(1,"label");
    histTemp->GetXaxis()->SetTitle("Collision Energy #sqrt{s_{NN}} (GeV)");
    histTemp->Draw();
    TPaveText * ptxt_3_40 = new TPaveText(0.15,0.7,0.38,0.85,"NDCARC");
    ptxt_3_40->SetFillColor(0);
    ptxt_3_40 -> AddText("Au+Au 4.5, 7.2 GeV FXT");
    ptxt_3_40 -> AddText("STAR preliminary");
    ptxt_3_40->Draw("same");
    // TGaxis *axis1 = new TGaxis(3,-0.05,41,-0.05,3,40,3,"G");
    // axis1->SetName("axis1");
    // axis1->Draw();
    // TPaveText * ptxt_4_xlabel = new TPaveText(0.178,0.055,0.198,0.095,"NDCARC");
    // ptxt_4_xlabel->SetFillColor(0);
    // ptxt_4_xlabel -> AddText("4");
    // ptxt_4_xlabel->Draw("same");
    // TPaveText * ptxt_40_xlabel = new TPaveText(0.878,0.055,0.908,0.095,"NDCARC");
    // ptxt_40_xlabel->SetFillColor(0);
    // ptxt_40_xlabel -> AddText("40");
    // ptxt_40_xlabel->SetAllWith("40", "font", 5);
    // ptxt_40_xlabel->Draw("same");
    //4.5 GeV phi
    double x0[1]    = {4.5};
    double zero0[1] = {0};
    double py0[1]      = {-0.01301};
    double ey_stat0[1] = {0.01566};
    ey_stat0[0] *= sqrt(2);
    double ey_sys0[1]  = {0.013012};
    ey_sys0[0] *= sqrt(2);
    std::cout<< "ey_sys0[0] = "<< ey_sys0[0] << std::endl;
    //0.013012
    // 7.2 GeV
    double x_phi_7p2[1]    = {7.2};
    double xErr_phi_7p2[1] = {0};
    double y_phi_7p2[1]      = {0.0138904};
    double yErr_stat_phi_7p2[1] = {0.00916586};
    double yErr_sys_phi_7p2[1] = {0.00647123};

    // 7.2 7.7 GeV together
    double x_phi_7p2_7p7[1]    = {7.6};
    double xErr_phi_7p2_7p7[1] = {0};
    double y_phi_7p2_7p7[1]      = {0.0139926};//{0.0159611};
    double yErr_stat_phi_7p2_7p7[1] = {0.00873405};//{0.0403535};

    // 7.7 GeV BES-I data, credit Guannan
    double x_phi_7p7[1]    = {8.2};
    double xErr_phi_7p7[1] = {0};
    double y_phi_7p7[1]      = {0.0153485};//{0.0159611};
    double yErr_stat_phi_7p7[1] = {0.0391923};//{0.0403535};

    // BES-I phi
    double x1[6]    = { 11.5, 14.5, 19.6, 27, 39, 200};
    double zero1[6] = {0};
    double py1[6]      = {0.0199, -0.0311, -0.0304, -0.02,   -0.013,   -0.00366944};
    double ey_stat1[6] = {0.0183, 0.00938, 0.00763, 0.00644, 0.006142, 0.00203012};
    double ey_sys1[6]  = {0.03,   0.021,   0.014,   0.015,   0.012,    0.000533};

    // BES-I Proton
    double x2[8]    = {7.7, 11.5, 14.5, 19.6, 27, 39, 62.4, 200};
    double zero2[8] = {0};
    double py2[8]      = {0.022534   , -0.000486069, -0.00684062, -0.00690107, -0.00561899, -0.00415636, -0.00127671, 0.000601299};
    double ey_stat2[8] = {0.000863638, 0.000618966 , 0.000843574, 0.000853584, 0.000836661, 0.000565664, 0.00168219 , 0.000738392};
    double ey_sys2[8]  = {0.00447883 , 0.00225644  , 0.0021     , 0.000933001, 0.0008     , 0.0001     , 0.0031039  , 0.00011};

    // BES-I KPlus
    double x3[8]       = {7.7, 11.5, 14.5, 19.6, 27, 39, 62.4, 200};
    double zero3[8]    = {0};
    double py3[8]      = {-0.0137338, -0.0135527,  -0.0127973, -0.0124038,  -0.0104937,  -0.0082035,  -0.00771108, -0.00270561};
    double ey_stat3[8] = {0.00166579, 0.00100918,  0.001172,   0.00103033,  0.000927541, 0.000581183, 0.00185789,  0.000290576};
    double ey_sys3[8]  = {0.00076095, 0.000472041, 0.00104204, 0.000836145, 0.00021,     0.0002,      0.0029149,   0.000137};

    // BES-I KMinus
    double x4[8]       = {7.7, 11.5, 14.5, 19.6, 27, 39, 62.4, 200};
    double zero4[8]    = {0};
    double py4[8]      = {-0.00756657, -0.0163133,  -0.0160707,  -0.0174652,  -0.0147908, -0.011649,   -0.010862,  -0.00336346};
    double ey_stat4[8] = {0.00290359,  0.00150192,  0.00162578,  0.00132712,  0.00112134, 0.000672674, 0.00205018, 0.000302512};
    double ey_sys4[8]  = {0.000950781, 0.000974057, 0.000974057, 0.000458308, 0.0002,     0.0001,      0.0030239,  0.000194};

    double d_offset = 0.1;
    for(int i = 0;i<6;i++){
      x1[i] += 5*d_offset;
      // x2[i] -= d_offset;
      // x3[i] -= d_offset * 2;
      // x4[i] -= d_offset * 3;
    }
    // for(int i = 0;i<8;i++){
    //   x2[i] += d_offset;
    //   x3[i] += d_offset * 2;
    //   x4[i] += d_offset * 3;
    // }

    // data set (1) with stat and sys errors

    //

    // Now draw data set (1)
    // We first have to draw it only with the stat errors
    gStyle->SetOptDate(0);
    gStyle->SetEndErrorSize(6);
    c1->SetGridx(0);
    c1->SetGridy(0);
    c1->SetLogx();

    //4.5 GeV phi
    TGraphErrors *graph0 = new TGraphErrors(1, x0, py0, zero0, ey_stat0);
    graph0->SetMarkerStyle(34);
    graph0->SetMarkerColor(kBlack);
    graph0->SetLineColor(kBlack);
    graph0->SetMarkerSize(2);
    graph0->Draw("P");

    // 7.2 GeV
    TGraphErrors *graph_7p2 = new TGraphErrors(1, x_phi_7p2, y_phi_7p2, xErr_phi_7p2, yErr_stat_phi_7p2);
    graph_7p2->SetMarkerStyle(29);
    graph_7p2->SetMarkerColor(kBlue);
    graph_7p2->SetLineColor(kBlue);
    graph_7p2->SetMarkerSize(2.5);
    graph_7p2->Draw("P");
    TGraphErrors *graph_7p2_sys = new TGraphErrors(1, x_phi_7p2, y_phi_7p2, xErr_phi_7p2, yErr_sys_phi_7p2);
    graph_7p2_sys->SetMarkerColor(kBlue);
    graph_7p2_sys->SetLineColor(kBlue);
    graph_7p2_sys->Draw("[]");

    // 7.2 GeV
    TGraphErrors *graph_7p2_7p7 = new TGraphErrors(1, x_phi_7p2_7p7, y_phi_7p2_7p7, xErr_phi_7p2_7p7, yErr_stat_phi_7p2_7p7);
    graph_7p2_7p7->SetMarkerStyle(29);
    graph_7p2_7p7->SetMarkerColor(kRed);
    graph_7p2_7p7->SetLineColor(kRed);
    graph_7p2_7p7->SetMarkerSize(2.5);
    graph_7p2_7p7->Draw("P");

    // BES-I phi
    TGraphErrors *graph1 = new TGraphErrors(6, x1, py1, zero1, ey_stat1);
    graph1->SetMarkerStyle(28);
    graph1->SetMarkerColor(kBlack);
    graph1->SetLineColor(kBlack);
    graph1->SetMarkerSize(2);
    graph1->Draw("P");

    // BES-I Proton
    TGraphErrors *graph2 = new TGraphErrors(8, x2, py2, zero2, ey_stat2);
    graph2->SetMarkerStyle(25);
    graph2->SetMarkerColor(kBlack);
    graph2->SetLineColor(kBlack);
    graph2->SetMarkerSize(1);
    graph2->Draw("P");

    // BES-I KPlus
    TGraphErrors *graph3 = new TGraphErrors(8, x3, py3, zero3, ey_stat3);
    graph3->SetMarkerStyle(26);
    graph3->SetMarkerColor(kBlack);
    graph3->SetLineColor(kBlack);
    graph3->SetMarkerSize(1);
    graph3->Draw("P");

    // BES-I KMinus
    TGraphErrors *graph4 = new TGraphErrors(8, x4, py4, zero4, ey_stat4);
    graph4->SetMarkerStyle(32);
    graph4->SetMarkerColor(kBlack);
    graph4->SetLineColor(kBlack);
    graph4->SetMarkerSize(1);
    graph4->Draw("P");

    TLegend *legend = new TLegend(0.54,0.6,0.88,0.88);
    legend->SetBorderSize(0);
    legend->AddEntry(graph0,"#phi  4.5 GeV FXT 0-30%","p");
    legend->AddEntry(graph_7p2,"#phi  7.2 GeV FXT 10-40%","p");
    legend->AddEntry(graph_7p2_7p7,"#phi  7.2 FXT 7.7 COL combined 10-40%","p");
    legend->AddEntry(graph1,"#phi  BES-I 10-40%","p");
    legend->AddEntry(graph2,"p  BES-I 10-40%","p");
    legend->AddEntry(graph3,"K^{+}  BES-I 10-40%","p");
    legend->AddEntry(graph4,"K^{-}  BES-I 10-40%","p");
    //, offset -0.2 GeV  , 7.7 GeV offset 0.2 GeV
    legend->Draw("same");
    // Now we have to somehow depict the sys errors

    //4.5 GeV phi
    TGraphErrors *graph0_sys = new TGraphErrors(1, x0, py0, zero0, ey_sys0);
    graph0_sys->SetLineColor(kBlack);
    graph0_sys->Draw("[]");

    cout << "py0 = "     << py0[0] << endl;
    cout << "ey_stat0 = "<< ey_stat0[0] << endl;
    cout << "ey_sys0 = " << ey_sys0[0] << endl;
    // 7.7 GeV BES-I data, credit Guannan
    TGraphErrors *graph_7p7 = new TGraphErrors(1, x_phi_7p7, y_phi_7p7, xErr_phi_7p7, yErr_stat_phi_7p7);
    graph_7p7->SetMarkerStyle(28);
    graph_7p7->SetMarkerColor(kBlack);
    graph_7p7->SetLineColor(kBlack);
    graph_7p7->SetMarkerSize(2);
    graph_7p7->Draw("P");
    // BES-I phi
    TGraphErrors *graph1_sys = new TGraphErrors(6, x1, py1, zero1, ey_sys1);
    graph1_sys->SetLineColor(kBlack);
    graph1_sys->Draw("[]");

    // BES-I Proton
    TGraphErrors *graph2_sys = new TGraphErrors(8, x2, py2, zero2, ey_sys2);
    graph2_sys->SetLineColor(kBlack);
    graph2_sys->Draw("[]");

    // BES-I KPlus
    TGraphErrors *graph3_sys = new TGraphErrors(8, x3, py3, zero3, ey_sys3);
    graph3_sys->SetLineColor(kBlack);
    graph3_sys->Draw("[]");

    // BES-I KMinus
    TGraphErrors *graph4_sys = new TGraphErrors(8, x4, py4, zero4, ey_sys4);
    graph4_sys->SetLineColor(kBlack);
    graph4_sys->Draw("[]");

    TLine *line = new TLine(3, 0, 41, 0);
    line->SetLineStyle(7);
    line->Draw();
}
