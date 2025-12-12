#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TStyle.h>
#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <cmath>
#include <sstream>
#include <algorithm>
#include "ePIC_style.C"

void plot_multiplot_sigma_red_vs_x() { 
  // gStyle->SetPadLeftMargin(0.15);
  // gROOT->ForceStyle();
  gROOT->ProcessLine("set_ePIC_style()");

  std::ifstream infile("../FinalGrids/NNPDF31_nnlo_as_0118_proton_StandardBinning_ele_Ee10_Ep130_lumi1.0fb-1_pessimistic_grid.csv"); 
  std::string line;

  std::getline(infile, line); // Skip header

  struct DataRow {
    double Q2, x, y, purity, stability, sigma_red, sigma_red_pdf_unc, sigma_stat, sigma_p2p, sigma_norm;
  };
  std::vector<DataRow> rows;
  std::vector<double> Q2_values;

  while (std::getline(infile, line)) {
    std::stringstream ss(line);
    std::string token;
    std::vector<std::string> tokens;

    while (std::getline(ss, token, ',')) {
      tokens.push_back(token);
    }

    double Q2 = std::stod(tokens[0]);
    double x = std::stod(tokens[1]);
    double y = std::stod(tokens[2]);
    double purity = std::stod(tokens[3]);
    double stability = std::stod(tokens[4]);
    double sigma_red = std::stod(tokens[5]);
    double sigma_red_pdf_unc = std::stod(tokens[6]);
    double sigma_stat = std::stod(tokens[7]);
    double sigma_p2p = std::stod(tokens[8]);
    double sigma_norm = std::stod(tokens[9]);

    if (y < 0.01 || y > 0.95) continue;
    // if (purity < 0.2 || stability < 0.2) continue;

    // rows.push_back({Q2, x, sigma_red, sigma_stat, sigma_red_pdf, sigma_stat_pdf});
    rows.push_back({Q2, x, y, purity, stability, sigma_red, sigma_red_pdf_unc, sigma_stat, sigma_p2p, sigma_norm});
    Q2_values.push_back(Q2);
  }

  std::sort(Q2_values.begin(), Q2_values.end());
  Q2_values.erase(std::unique(Q2_values.begin(), Q2_values.end()), Q2_values.end());

  std::map<double, int> Q2_index_map;
  for (size_t i = 0; i < Q2_values.size(); ++i) {
    Q2_index_map[Q2_values[i]] = i;
  }

  // Gradient colors
  std::vector<int> colors;
  int startColor = 51;
  int endColor = 99;
  int nBins = Q2_values.size();
  for (int i = 0; i < nBins; ++i) {
    int color = startColor + i * (endColor - startColor) / std::max(1, nBins - 1);
    colors.push_back(color);
  }

  std::vector<TGraphErrors*> graphs;
  std::vector<TGraphErrors*> graphs_pdf;

  for (size_t i = 0; i < Q2_values.size(); ++i) {
    double Q2 = Q2_values[i];
    std::vector<double> x_vals, sigma_vals, sigma_errs, sigma_vals_pdf, sigma_errs_pdf;

    for (const auto& row : rows) {
      if (std::abs(row.Q2 - Q2) < 1e-6) {
        double scale = 1;
        x_vals.push_back(row.x);
        // sigma_vals.push_back(row.sigma_red * scale);
        // sigma_errs.push_back(0.01 * row.sigma_red * std::sqrt(row.sigma_stat*row.sigma_stat+1.9*1.9+3.4*3.4) * scale);
        // sigma_vals_pdf.push_back(row.sigma_red_pdf * scale);
        sigma_vals.push_back(row.sigma_red * scale);
        // sigma_errs.push_back(0.01 * row.sigma_red * row.sigma_stat * scale); // stat only
        sigma_errs.push_back(0.01 * row.sigma_red * sqrt(row.sigma_stat*row.sigma_stat+row.sigma_p2p*row.sigma_p2p+row.sigma_norm*row.sigma_norm) * scale); // stat + sys
        sigma_vals_pdf.push_back(row.sigma_red * scale);
        sigma_errs_pdf.push_back(row.sigma_red_pdf_unc * scale);
      }
    }

    TGraphErrors* graph = new TGraphErrors(x_vals.size(), x_vals.data(), sigma_vals.data(), nullptr, sigma_errs.data());
    TGraphErrors* graph_pdf = new TGraphErrors(x_vals.size(), x_vals.data(), sigma_vals_pdf.data(), nullptr, sigma_errs_pdf.data());
    graph->SetMarkerStyle(20);
    graph->SetMarkerSize(1);
    graph->SetMarkerColor(colors[i]);
    graph->SetLineColor(colors[i]);
    graphs.push_back(graph);
    graph_pdf->SetLineWidth(2);
    graph_pdf->SetFillColor(kBlack);
    graph_pdf->SetFillStyle(3002);
    graphs_pdf.push_back(graph_pdf);
  }

  TCanvas* c1 = new TCanvas("c1", "sigma_red_vs_x", 1600, 1000);
  c1->Divide(5, 4);

  double x_min = 1e-5;
  double x_max = 1;
  double y_min = 0;
  double y_max = 2;

  for (size_t i = 0; i < graphs.size(); ++i) {
    c1->cd(i + 1);
    gPad->SetLogx();
    gPad->SetTopMargin(0.05);
    gPad->SetRightMargin(0.05);
    gPad->SetBottomMargin(0.15);
    gPad->SetLeftMargin(0.15);

    graphs_pdf[i]->SetTitle("");
    graphs_pdf[i]->GetXaxis()->SetTitle("x");
    graphs_pdf[i]->GetYaxis()->SetTitle("#sigma_{red}");
    graphs_pdf[i]->GetXaxis()->SetLimits(x_min, x_max);
    graphs_pdf[i]->SetMinimum(y_min);
    graphs_pdf[i]->SetMaximum(2.2 * graphs_pdf[i]->GetHistogram()->GetMaximum());

    graphs_pdf[i]->GetXaxis()->SetTitleSize(0.06);
    graphs_pdf[i]->GetXaxis()->SetLabelSize(0.05);
    graphs_pdf[i]->GetYaxis()->SetTitleSize(0.06);
    graphs_pdf[i]->GetYaxis()->SetLabelSize(0.05);

    graphs_pdf[i]->Draw("3AL");
    graphs[i]->Draw("P SAME");

    TLegend* leg = new TLegend(0.48, 0.7, 0.93, 0.93);
    leg->SetTextSize(0.045);
    leg->SetBorderSize(0);
    leg->AddEntry(graphs[i], "ePIC Simulation", "ep");
    leg->AddEntry(graphs_pdf[i], "NNPDF 3.1 Prediction", "fl");
    leg->Draw();

    TLatex latex;
    latex.SetTextSize(0.06);
    latex.SetTextFont(62); // bold font
    latex.SetNDC(true);
    latex.DrawLatex(0.18, 0.85, Form("Q^{2} = %.1f GeV^{2}", Q2_values[i]));
  }

 c1->cd(graphs.size() + 1);
 
 TLatex Text_com;
 Text_com.SetTextAlign(13);  //align at top
 Text_com.SetTextSize(0.12);
 // Manually change c.o.m. energy
 Text_com.DrawLatexNDC(.1,.81,"e+p, #sqrt{s} = 72 GeV");
 Text_com.DrawLatexNDC(.1,.66,"L_{proj} = 1 fb^{-1}");
 
 
 TLatex Text_ePIC;
 Text_ePIC.SetTextSize(0.12);
 Text_ePIC.SetTextFont(62);
 Text_ePIC.DrawLatexNDC(.1,.9,"ePIC Performance");  // performance plot
 //Text_ePIC.DrawLatexNDC(.15,.88,"ePIC Internal");  // for internal use only
 //Text_ePIC.DrawLatexNDC(.15,.88,"ePIC Preliminary"); // preliminary released version 
 //Text_ePIC.DrawLatexNDC(.15,.88,"ePIC Work in Progress"); // work in progress to be shown outside
 //Text_ePIC.DrawLatexNDC(.15,.88,"ePIC"); // final published version

 // Add dates: needed for performance plots
 TLatex Text_date;
 Text_date.SetTextSize(0.1);
 Text_date.SetTextFont(52);
 Text_date.DrawLatexNDC(.1,.1,"Simu campaign: 10/2025");  // performance plot

  c1->Update();
}
