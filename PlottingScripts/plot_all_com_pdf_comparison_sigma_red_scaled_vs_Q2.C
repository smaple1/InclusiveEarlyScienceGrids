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
// Style file
#include "ePIC_style.C"


void make_10x250_canvas(std::vector<int> colors) {
  gStyle->SetPadLeftMargin(0.15);
  // gStyle->SetPalette(kBird);
  gROOT->ForceStyle();

  std::ifstream infile("../FinalGrids/NNPDF31_nnlo_as_0118_proton_StandardBinning_ele_Ee10_Ep250_lumi1.0fb-1_pessimistic_grid.csv"); 
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
    Q2 = Q2 / 0.9;
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
  // std::vector<int> colors;
  // int startColor = 51;
  // int endColor = 99;
  // int nBins = x_values.size();
  // for (int i = 0; i < nBins; ++i) {
  // int color = startColor + i * (endColor - startColor) / std::max(1, nBins - 1);
  // colors.push_back(color);
  // }
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
    graph->SetMarkerStyle(24);
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


  double x_min = 0.9;
  double x_max = 1e5;
  double y_min = 1e-4;
  double y_max = 1e6;

  // Draw first to set up axes
  graphs_pdf[0]->SetTitle("");
  graphs_pdf[0]->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");
  graphs_pdf[0]->GetYaxis()->SetTitle("#sigma_{red} #times 2^{i}");
  graphs_pdf[0]->GetXaxis()->SetLimits(x_min, x_max);
  graphs_pdf[0]->SetMinimum(y_min);
  graphs_pdf[0]->SetMaximum(y_max);
  graphs_pdf[0]->Draw("3L SAME");
  graphs[0]->Draw("P SAME");

  // Draw remaining graphs
  for (size_t i = 1; i < graphs.size(); ++i) {
    graphs_pdf[i]->Draw("3L SAME");
    graphs[i]->Draw("P SAME");
  }

  TMarker* marker1 = new TMarker(0, 0, 20);  // style 20 = full circle
  marker1->SetMarkerStyle(20);
  
  TMarker* marker2 = new TMarker(0, 0, 21);  // style 21 = full square
  marker2->SetMarkerStyle(21);

  TMarker* marker3 = new TMarker(0, 0, 24);  // style 21 = full square
  marker2->SetMarkerStyle(21);

  // Add legend
  // TLegend* legend = new TLegend(0.5, 0.75, 0.91, 0.88);
  TLegend* legend = new TLegend(0.6, 0.70, 0.91, 0.93);// use with ePIC Simulation on left
  // TLegend* legend = new TLegend(0.6, 0.55, 0.91, 0.78);// use with ePIC Simulation on right
  legend->SetFillStyle(0); // transparent background
  // legend->AddEntry(graphs[0], "ePIC Full Simulation", "p");
  // legend->AddEntry(marker1, "e+p 5#times41 GeV^{2}", "p");
  // legend->AddEntry(marker2, "e+p 10#times100 GeV^{2}", "p");
  // legend->AddEntry(marker3, "e+p 18#times275 GeV^{2}", "p");
  // legend->AddEntry((TGraph*)0, "L_{proj} = 10 fb^{-1}", "");
  legend->AddEntry(marker1, "e+p, #sqrt{s} = 51 GeV", "p");
  legend->AddEntry(marker2, "e+p, #sqrt{s} = 72 GeV", "p");
  legend->AddEntry(marker3, "e+p, #sqrt{s} = 100 GeV", "p");
  // legend->AddEntry(graphs[0], "Statistical Error (ePIC Full Simulation)", "e");
  legend->AddEntry(graphs_pdf[0], "NNPDF 3.1 Prediction", "FL");
  legend->SetBorderSize(0);
  legend->Draw();

  // Add annotation
  TLatex latex;
  latex.SetNDC();
  latex.SetTextSize(0.035);
  // latex.DrawLatex(0.16, 0.92, "ePIC ep 18#times275 GeV^{2} (10.0 fb^{-1})");

  // Labels per line
  latex.SetTextSize(0.025);
  latex.SetTextFont(42);
  latex.SetTextAlign(12);
  latex.SetNDC(false); // Use axis coords

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

    latex.SetTextColor(colors[i]);
    latex.DrawLatex(q2_avg * 2, y_avg, Form("x=%.1e, i=%zu", x_values[i], i));
  }

}

void make_10x130_canvas(std::vector<int> colors) {
  gStyle->SetPadLeftMargin(0.15);
  // gStyle->SetPalette(kBird);
  gROOT->ForceStyle();

  std::ifstream infile("../FinalGrids/NNPDF31_nnlo_as_0118_proton_StandardBinning_ele_Ee10_Ep130_lumi1.0fb-1_pessimistic_grid.csv"); 
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
  // std::vector<int> colors;
  // int startColor = 51;
  // int endColor = 99;
  // int nBins = x_values.size();
  // for (int i = 0; i < nBins; ++i) {
  // int color = startColor + i * (endColor - startColor) / std::max(1, nBins - 1);
  // colors.push_back(color);
  // }
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
      }
    }

    TGraphErrors* graph = new TGraphErrors(Q2_vals.size(), Q2_vals.data(), sigma_vals.data(), nullptr, sigma_errs.data());
    TGraphErrors* graph_pdf = new TGraphErrors(Q2_vals.size(), Q2_vals.data(), sigma_vals_pdf.data(), nullptr, sigma_errs_pdf.data());
    graph->SetMarkerStyle(21);
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

  // c1->SetLogx();
  // c1->SetLogy();

  double x_min = 0.9;
  double x_max = 1e5;
  double y_min = 1e-4;
  double y_max = 1e6;

  // Draw first to set up axes
  graphs_pdf[0]->SetTitle("");
  graphs_pdf[0]->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");
  graphs_pdf[0]->GetYaxis()->SetTitle("#sigma_{red} #times 2^{i}");
  graphs_pdf[0]->GetXaxis()->SetLimits(x_min, x_max);
  graphs_pdf[0]->SetMinimum(y_min);
  graphs_pdf[0]->SetMaximum(y_max);
  graphs_pdf[0]->Draw("3L SAME");
  graphs[0]->Draw("P SAME");

  // Draw remaining graphs
  for (size_t i = 1; i < graphs.size(); ++i) {
    graphs_pdf[i]->Draw("3L SAME");
    graphs[i]->Draw("P SAME");
  }

  // c1->Update();
}

void make_5x100_canvas(std::vector<int> colors) {
  gStyle->SetPadLeftMargin(0.15);
  // gStyle->SetPalette(kBird);
  gROOT->ForceStyle();

  std::ifstream infile("../FinalGrids/NNPDF31_nnlo_as_0118_proton_StandardBinning_ele_Ee5_Ep130_lumi1.0fb-1_pessimistic_grid.csv"); 
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
    Q2 = Q2 * 0.9;
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
  // std::vector<int> colors;
  // int startColor = 51;
  // int endColor = 99;
  // int nBins = x_values.size();
  // for (int i = 0; i < nBins; ++i) {
  // int color = startColor + i * (endColor - startColor) / std::max(1, nBins - 1);
  // colors.push_back(color);
  // }
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

  // c1->SetLogx();
  // c1->SetLogy();

  double x_min = 0.9;
  double x_max = 5e4;
  double y_min = 5e-4;
  double y_max = 1e8;

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

  // c1->Update();
}

void plot_all_com_pdf_comparison_sigma_red_scaled_vs_Q2() {
  gROOT->ProcessLine("set_ePIC_style()");
  // colour palette
  std::vector<int> colors;
  int startColor = 51;
  int endColor = 99;
  int nBins = 20;
  for (int i = 0; i < nBins; ++i) {
    int color = startColor + i * (endColor - startColor) / std::max(1, nBins - 1);
    colors.push_back(color);
  }
  // TCanvas* c1 = new TCanvas("c1", "sigma_red * 2^i vs Q2", 600, 800);
  TCanvas* c1 = new TCanvas("c1", "sigma_red * 2^i vs Q2", 1200, 700);
  c1->SetLogx();
  c1->SetLogy();
  make_5x100_canvas(colors);
  c1->Update();
  make_10x130_canvas(colors);
  c1->Update();
  make_10x250_canvas(colors);
  c1->Update();
  // ePIC Simulation on left
  TLatex Text_ePIC;
  Text_ePIC.SetTextSize(0.05);
  Text_ePIC.SetTextFont(62);
  Text_ePIC.DrawLatexNDC(.15,.88,"ePIC Performance");
  TLatex Text_com;
  Text_com.SetTextAlign(13);  //align at top
  Text_com.DrawLatexNDC(.15,.85,"L_{proj} = 1 fb^{-1} per beam setting");
  // ePIC internal on right
  // TLatex Text_ePIC;
  // Text_ePIC.SetTextSize(0.05);
  // Text_ePIC.SetTextFont(62);
  // Text_ePIC.DrawLatexNDC(.61,.88,"ePIC Simulation");
  // TLatex Text_com;
  // Text_com.SetTextAlign(13);  //align at top
  // Text_com.DrawLatexNDC(.61,.85,"L_{proj} = 10 fb^{-1}");
  // Add dates: needed for performance plots
  // TLatex Text_date;
  // Text_date.SetTextSize(0.035);
  // Text_date.SetTextFont(52);
  // Text_date.DrawLatexNDC(.72,.96,"Released: MM/YYYY");
}
