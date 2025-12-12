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

void plot_errors_vs_Q2() {
  // gStyle->SetPadLeftMargin(0.15);
  gROOT->ProcessLine("set_ePIC_style()");

  std::ifstream infile("../FinalGrids/NNPDF31_nnlo_as_0118_proton_StandardBinning_ele_Ee10_Ep250_lumi1.0fb-1_all_columns_grid.csv"); 
  std::string line;

  std::getline(infile, line); // Skip header

  struct DataRow {
    double Q2, x, y, purity, stability, sigma_red, sigma_red_pdf_unc, sigma_stat, sigma_p2p_opt, sigma_p2p_pess, sigma_norm_opt, sigma_norm_pess;
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
    double purity = std::stod(tokens[6]);
    double stability = std::stod(tokens[7]);
    double sigma_red = std::stod(tokens[8]);
    double sigma_red_pdf_unc = std::stod(tokens[9]);
    double sigma_stat = std::stod(tokens[10]);
    double sigma_p2p_opt = std::stod(tokens[11]);
    double sigma_p2p_pess = std::stod(tokens[12]);
    double sigma_norm_opt = std::stod(tokens[13]);
    double sigma_norm_pess= std::stod(tokens[14]);

    // Drawing cuts
    if (y < 0.01 || y > 0.95) continue;
    // if (purity < 0.3 || stability < 0.3) continue;

    rows.push_back({Q2, x, y, purity, stability, sigma_red, sigma_red_pdf_unc, sigma_stat, sigma_p2p_opt, sigma_p2p_pess, sigma_norm_opt, sigma_norm_pess});
    x_values.push_back(x);
  }

  std::sort(x_values.begin(), x_values.end());
  x_values.erase(std::unique(x_values.begin(), x_values.end()), x_values.end());
  std::reverse(x_values.begin(), x_values.end());

  std::map<double, int> x_index_map;
  for (size_t i = 0; i < x_values.size(); ++i) {
    x_index_map[x_values[i]] = i;
  }

  for (size_t i = 0; i < x_values.size(); ++i) {
    double x = x_values[i];
    std::vector<double> Q2_vals, sigma_vals, sigma_errs, sigma_errs_opt, sigma_errs_pess;

    for (const auto& row : rows) {
      if (std::abs(row.x - x) < 1e-10) {
        Q2_vals.push_back(row.Q2);
        sigma_vals.push_back(0);
        sigma_errs.push_back(row.sigma_stat);
        sigma_errs_opt.push_back(sqrt(row.sigma_stat*row.sigma_stat+row.sigma_p2p_opt*row.sigma_p2p_opt+row.sigma_norm_opt*row.sigma_norm_opt));
        sigma_errs_pess.push_back(sqrt(row.sigma_stat*row.sigma_stat+row.sigma_p2p_pess*row.sigma_p2p_pess+row.sigma_norm_pess*row.sigma_norm_pess));
      }
    }

    TGraphErrors* graph = new TGraphErrors(Q2_vals.size(), Q2_vals.data(), sigma_vals.data(), nullptr, sigma_errs.data());
    graph->SetMarkerStyle(20);
    // graph->SetLineColor(kBlack);
    graph->SetMarkerSize(1);

    TGraphErrors* graph_opt = new TGraphErrors(Q2_vals.size(), Q2_vals.data(), sigma_vals.data(), nullptr, sigma_errs_opt.data());
    graph_opt->SetMarkerStyle(20);
    // graph_opt->SetLineColor(kGreen);
    graph_opt->SetMarkerSize(1);

    TGraphErrors* graph_pess= new TGraphErrors(Q2_vals.size(), Q2_vals.data(), sigma_vals.data(), nullptr, sigma_errs_pess.data());
    graph_pess->SetMarkerStyle(20);
    // graph_pess->SetLineColor(kBlue);
    graph_pess->SetMarkerSize(1);

    // TCanvas* c1 = new TCanvas("c1", "sigma_red * 2^i vs Q2", 600, 800);
    TCanvas* c1 = new TCanvas();
    c1->SetLogx();

    graph->SetFillColorAlpha(kGreen, 0.3);
    graph_opt->SetFillColorAlpha(5, 0.3);
    graph_pess->SetFillColorAlpha(kRed, 0.3);

    graph_pess->SetTitle("");
    graph_pess->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");
    graph_pess->GetXaxis()->SetMoreLogLabels();
    graph_pess->GetYaxis()->SetTitle("#delta [%]");
    double max_val = *std::max_element(sigma_errs_pess.begin(), sigma_errs_pess.end());
    // graph_pess->GetYaxis()->SetRangeUser(-6.8,6.8);
    graph_pess->GetYaxis()->SetRangeUser(-max_val/0.6,max_val/0.6);
    graph_pess->Draw("A3");
    graph_opt->Draw("3SAME");
    graph->Draw("3SAME");

    // Add annotation (old - pre ePIC style)
    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.05);
    latex.SetTextFont(62);
    latex.DrawLatex(0.16, 0.88, Form("x=%.1e", x_values[i]));

    TLatex Text_com;
    Text_com.SetTextAlign(13);  //align at top
    Text_com.DrawLatexNDC(.5,.85,"e+p, #sqrt{s} = 100 GeV, L_{proj} = 1 fb^{-1}");
    // Text_com.DrawLatexNDC(.65,.8,"L_{proj} = 10 fb^{-1}");
    
    
    TLatex Text_ePIC;
    Text_ePIC.SetTextSize(0.05);
    Text_ePIC.SetTextFont(62);
    Text_ePIC.DrawLatexNDC(.5,.88,"ePIC Performance");  // performance plot
    
    // Add dates: needed for performance plots
    TLatex Text_date;
    Text_date.SetTextSize(0.035);
    Text_date.SetTextFont(52);
    Text_date.DrawLatexNDC(.65,.96,"Simu campaign: 10/2025");  // performance plot

    // Add legend
    // TLegend* legend = new TLegend(0.17, 0.75, 0.85, 0.88);// top
    TLegend* legend = new TLegend(0.17, 0.12, 0.85, 0.25);// bottom
    legend->SetNColumns(3);
    legend->SetFillStyle(0); // transparent background
    legend->AddEntry(graph, "#delta_{stat} ", "f");
    legend->AddEntry(graph_opt, "#delta_{stat} #oplus #delta_{sys,optimistic}  ", "f");
    legend->AddEntry(graph_pess, "#delta_{stat} #oplus #delta_{sys,pessimistic}", "f");
    legend->SetBorderSize(0);
    legend->Draw();

    c1->Update();
    c1->Print(Form("error_plots/Ee10_Ep250/Error_vs_Q2_x%.1e.png",x_values[i]));
  }


}
