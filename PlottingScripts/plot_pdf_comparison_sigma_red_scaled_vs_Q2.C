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

void plot_pdf_comparison_sigma_red_scaled_vs_Q2() {
  // gStyle->SetPadLeftMargin(0.15);
  // gStyle->SetPalette(kBird);
  // gROOT->ForceStyle();
  gROOT->ProcessLine("set_ePIC_style()");

  std::ifstream infile("../FinalGrids/NNPDF31_nnlo_as_0118_proton_StandardBinning_ele_Ee10_Ep130_lumi1.0fb-1_pessimistic_y_degraded_below_0.1_grid.csv"); 
  std::string line;

  std::getline(infile, line); // Skip header

  // struct DataRow {
  // double Q2, x, sigma_red, sigma_stat;
  // };
  struct DataRow {
    double Q2, x, y, purity, stability, sigma_red, sigma_red_pdf_unc, sigma_stat, sigma_p2p, sigma_norm;
  };
  std::vector<DataRow> rows;
  std::vector<double> x_values;

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

    // Drawing cuts
    if (y < 0.01 || y > 0.95) continue;
    // if (purity < 0.3 || stability < 0.3) continue;

    // rows.push_back({Q2, x, sigma_red, sigma_stat});
    rows.push_back({Q2, x, y, purity, stability, sigma_red, sigma_red_pdf_unc, sigma_stat, sigma_p2p, sigma_norm});
    x_values.push_back(x);
  }

  std::sort(x_values.begin(), x_values.end());
  x_values.erase(std::unique(x_values.begin(), x_values.end()), x_values.end());
  std::reverse(x_values.begin(), x_values.end());

  std::map<double, int> x_index_map;
  for (size_t i = 0; i < x_values.size(); ++i) {
    x_index_map[x_values[i]] = i;
  }

  // Gradient colors: from blue to red using default ROOT palette (start at 51)
  std::vector<int> colors;
  int startColor = 51;
  int endColor = 99;
  int nBins = x_values.size();
  for (int i = 0; i < nBins; ++i) {
    int color = startColor + i * (endColor - startColor) / std::max(1, nBins - 1);
    colors.push_back(color);
  }
  // Uncomment to use palette
  // int numColors = x_values.size();
  // for (int i = 0; i < numColors; ++i) {
  // int colorIndex = i * gStyle->GetNumberOfColors() / std::max(1, numColors - 1);
  // int color = TColor::GetColorPalette(colorIndex);
  // colors.push_back(color);
  // }

  std::vector<TGraphErrors*> graphs;
  std::vector<TGraphErrors*> graphs_pdf;

  for (size_t i = 0; i < x_values.size(); ++i) {
    double x = x_values[i];
    std::vector<double> Q2_vals, sigma_vals, sigma_errs, sigma_vals_pdf, sigma_errs_pdf;

    for (const auto& row : rows) {
      if (std::abs(row.x - x) < 1e-10) {
        double scale = std::pow(2, i);
        Q2_vals.push_back(row.Q2);
        sigma_vals.push_back(row.sigma_red * scale);
        // sigma_errs.push_back(0.01 * row.sigma_red * row.sigma_stat * scale); // stat only
        sigma_errs.push_back(0.01 * row.sigma_red * sqrt(row.sigma_stat*row.sigma_stat+row.sigma_p2p*row.sigma_p2p+row.sigma_norm*row.sigma_norm) * scale); // stat + sys
        sigma_vals_pdf.push_back(row.sigma_red * scale);
        sigma_errs_pdf.push_back(row.sigma_red_pdf_unc * scale);
        // sigma_errs.push_back(0.01 * row.sigma_red * std::sqrt(row.sigma_stat*row.sigma_stat+1*1) * scale);
        // sigma_errs.push_back(0.01 * row.sigma_red * std::sqrt(row.sigma_stat*row.sigma_stat+1.9*1.9+3.4*3.4) * scale);
      }
    }

    TGraphErrors* graph = new TGraphErrors(Q2_vals.size(), Q2_vals.data(), sigma_vals.data(), nullptr, sigma_errs.data());
    TGraphErrors* graph_pdf = new TGraphErrors(Q2_vals.size(), Q2_vals.data(), sigma_vals_pdf.data(), nullptr, sigma_errs_pdf.data());
    graph->SetMarkerStyle(20);
    graph->SetMarkerSize(1);
    // graph->SetMarkerSize(0.5);
    graph->SetMarkerColor(colors[i]);
    graph->SetLineColor(colors[i]);
    graphs.push_back(graph);
    graph_pdf->SetLineWidth(2);
    graph_pdf->SetFillColor(kBlack);
    graph_pdf->SetFillStyle(3002);
    graphs_pdf.push_back(graph_pdf);

  }

  // TCanvas* c1 = new TCanvas("c1", "sigma_red * 2^i vs Q2", 800, 700);
  TCanvas* c1 = new TCanvas("c1", "sigma_red * 2^i vs Q2", 1200, 700);
  c1->SetLogx();
  c1->SetLogy();

  double x_min = 0.9;
  double x_max = 1e5;
  double y_min = 1e-3;
  double y_max = 1e7;

  // Draw first to set up axes
  graphs_pdf[0]->SetTitle("");
  graphs_pdf[0]->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");
  graphs_pdf[0]->GetYaxis()->SetTitle("#sigma_{red} #times 2^{i}");
  graphs_pdf[0]->GetXaxis()->SetLimits(x_min, x_max);
  graphs_pdf[0]->SetMinimum(y_min);
  graphs_pdf[0]->SetMaximum(y_max);
  graphs_pdf[0]->Draw("3AL");
  graphs[0]->Draw("P SAME");

  // Draw remaining graphs
  for (size_t i = 1; i < graphs.size(); ++i) {
    graphs_pdf[i]->Draw("3L SAME");
    graphs[i]->Draw("P SAME");
  }

  // Add legend
  // TLegend* legend = new TLegend(0.45, 0.75, 0.98, 0.88);
  TLegend* legend = new TLegend(0.55, 0.8, 0.96, 0.93);
  legend->SetFillStyle(0); // transparent background
  legend->SetTextSize(0.035);
  legend->AddEntry(graphs[0], "ePIC Simulation", "ep");
  legend->AddEntry(graphs_pdf[0], "NNPDF 3.1 Prediction", "fl");
  legend->SetBorderSize(0);
  legend->Draw();

  // Add annotation
  // TLatex latex;
  // latex.SetNDC();
  // latex.SetTextSize(0.035);
  // latex.DrawLatex(0.16, 0.92, "ePIC ep 10#times130 GeV^{2} (10.0 fb^{-1})");

  // Labels per line
  // latex.SetTextSize(0.025);
  // latex.SetTextFont(42);
  // latex.SetTextAlign(12);
  // latex.SetNDC(false); // Use axis coords

  for (size_t i = 0; i < graphs.size(); ++i) {
    TGraphErrors* g = graphs[i];
    int n = g->GetN();
    if (n == 0) continue;

    double q2_avg = 0, y_avg = 0;
    for (int j = 0; j < n; ++j) {
      double x, y;
      g->GetPoint(j, x, y);
      q2_avg = x;
      y_avg = y;
    }
    // q2_avg /= n;
    // y_avg /= n;

    TLatex Text_com;
    Text_com.SetTextAlign(13);  //align at top
    Text_com.SetTextSize(0.035);
    Text_com.DrawLatexNDC(.15,.85,"e+p, #sqrt{s} = 100 GeV, L_{proj} = 1 fb^{-1}");
    
    TLatex Text_ePIC;
    Text_ePIC.SetTextSize(0.035);
    Text_ePIC.SetTextFont(62);
    Text_ePIC.DrawLatexNDC(.15,.88,"ePIC Performance");  // performance plot
    
    // Add dates: needed for performance plots
    TLatex Text_date;
    Text_date.SetTextSize(0.035);
    Text_date.SetTextFont(52);
    Text_date.DrawLatexNDC(.65,.96,"Simu campaign: 10/2025");  // performance plot

    TLatex latex;
    latex.SetTextColor(colors[i]);
    latex.SetTextSize(0.025);
    latex.SetTextFont(42);
    latex.SetTextAlign(12);
    latex.SetNDC(false); // Use axis coords
    latex.DrawLatex(q2_avg * 2, y_avg, Form("x=%.1e, i=%zu", x_values[i], i));
  }


  c1->Update();
}
