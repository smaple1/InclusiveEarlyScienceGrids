#ifndef BENCHMARKS_HEADER_ONLY_H
#define BENCHMARKS_HEADER_ONLY_H

#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TString.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h"

// Global histograms
#define NEW1D(name,title) new TH1F(name,title,100,-1,1)
#define NEW2D(name,title) new TH2F(name,title,30,0.001,1,30,-1,1)

static TH1F *hResoX_electron;
static TH1F *hResoX_jb;
static TH1F *hResoX_da;
static TH1F *hResoX_sigma;
static TH1F *hResoX_esigma;

static TH1F *hResoY_electron;
static TH1F *hResoY_jb;
static TH1F *hResoY_da;
static TH1F *hResoY_sigma;
static TH1F *hResoY_esigma;

static TH1F *hResoQ2_electron;
static TH1F *hResoQ2_jb;
static TH1F *hResoQ2_da;
static TH1F *hResoQ2_sigma;
static TH1F *hResoQ2_esigma;

static TH2F *hResoX_2D_electron;
static TH2F *hResoX_2D_jb;
static TH2F *hResoX_2D_da;
static TH2F *hResoX_2D_sigma;
static TH2F *hResoX_2D_esigma;

static TH2F *hResoY_2D_electron;
static TH2F *hResoY_2D_jb;
static TH2F *hResoY_2D_da;
static TH2F *hResoY_2D_sigma;
static TH2F *hResoY_2D_esigma;

static TH2F *hResoQ2_2D_electron;
static TH2F *hResoQ2_2D_jb;
static TH2F *hResoQ2_2D_da;
static TH2F *hResoQ2_2D_sigma;
static TH2F *hResoQ2_2D_esigma;

static TH1F *hReco_E_electron;
static TH1F *hReco_theta_electron;
static TH1F *hReco_Sigma_HFS;
static TH1F *hReco_pt_HFS;
static TH1F *hReco_Sigma_total;

static TH1F *hReso_Sigma_HFS;
static TH1F *hReso_pt_HFS;

// Initialise histograms
inline void init_benchmarks()
{
    hResoX_electron = NEW1D("hResoX_electron","Electron method;#Delta x/x;Counts");
    hResoX_jb       = NEW1D("hResoX_jb","JB method;#Delta x/x;Counts");
    hResoX_da       = NEW1D("hResoX_da","Double Angle method;#Delta x/x;Counts");
    hResoX_sigma    = NEW1D("hResoX_sigma","#Sigma method;#Delta x/x;Counts");
    hResoX_esigma   = NEW1D("hResoX_esigma","e-#Sigma method;#Delta x/x;Counts");

    hResoY_electron = NEW1D("hResoY_electron","Electron method;#Delta y/y;Counts");
    hResoY_jb       = NEW1D("hResoY_jb","JB method;#Delta y/y;Counts");
    hResoY_da       = NEW1D("hResoY_da","Double Angle method;#Delta y/y;Counts");
    hResoY_sigma    = NEW1D("hResoY_sigma","#Sigma method;#Delta y/y;Counts");
    hResoY_esigma   = NEW1D("hResoY_esigma","e-#Sigma method;#Delta y/y;Counts");

    hResoQ2_electron = NEW1D("hResoQ2_electron","Electron method;#Delta Q^{2}/Q^{2};Counts");
    hResoQ2_jb       = NEW1D("hResoQ2_jb","JB method;#Delta Q^{2}/Q^{2};Counts");
    hResoQ2_da       = NEW1D("hResoQ2_da","Double Angle method;#Delta Q^{2}/Q^{2};Counts");
    hResoQ2_sigma    = NEW1D("hResoQ2_sigma","#Sigma method;#Delta Q^{2}/Q^{2};Counts");
    hResoQ2_esigma   = NEW1D("hResoQ2_esigma","e-#Sigma method;#Delta Q^{2}/Q^{2};Counts");

    hResoX_2D_electron = NEW2D("hResoX_2D_electron","Electron method;y;#Delta x/x");
    hResoX_2D_jb       = NEW2D("hResoX_2D_jb","JB method;y;#Delta x/x");
    hResoX_2D_da       = NEW2D("hResoX_2D_da","Double Angle method;y;#Delta x/x");
    hResoX_2D_sigma    = NEW2D("hResoX_2D_sigma","#Sigma method;y;#Delta x/x");
    hResoX_2D_esigma   = NEW2D("hResoX_2D_esigma","e-#Sigma method;y;#Delta x/x");

    hResoY_2D_electron = NEW2D("hResoY_2D_electron","Electron method;y;#Delta y/y");
    hResoY_2D_jb       = NEW2D("hResoY_2D_jb","JB method;y;#Delta y/y");
    hResoY_2D_da       = NEW2D("hResoY_2D_da","Double Angle method;y;#Delta y/y");
    hResoY_2D_sigma    = NEW2D("hResoY_2D_sigma","#Sigma method;y;#Delta y/y");
    hResoY_2D_esigma   = NEW2D("hResoY_2D_esigma","e-#Sigma method;y;#Delta y/y");

    hResoQ2_2D_electron = NEW2D("hResoQ2_2D_electron","Electron method;y;#Delta Q^{2}/Q^{2}");
    hResoQ2_2D_jb       = NEW2D("hResoQ2_2D_jb","JB method;y;#Delta Q^{2}/Q^{2}");
    hResoQ2_2D_da       = NEW2D("hResoQ2_2D_da","Double Angle method;y;#Delta Q^{2}/Q^{2}");
    hResoQ2_2D_sigma    = NEW2D("hResoQ2_2D_sigma","#Sigma method;y;#Delta Q^{2}/Q^{2}");
    hResoQ2_2D_esigma   = NEW2D("hResoQ2_2D_esigma","e-#Sigma method;y;#Delta Q^{2}/Q^{2}");

    hReco_E_electron = new TH1F("hReco_E_electron","Electron energy;Energy [GeV];",100,0,25);
    hReco_theta_electron = new TH1F("hReco_theta_electron","Electron polar angle;Polar angle [rad];",100,0,3.5);
    hReco_Sigma_HFS = new TH1F("hReco_Sigma_HFS","HFS E-p_{z};E-p_{z} [GeV];",100,0,50);
    hReco_pt_HFS = new TH1F("hReco_pt_HFS","HFS p_{T};p_{T} [GeV];",100,0,100);
    hReco_Sigma_total = new TH1F("hReco_Sigma_total","Total E-p_{z};E-p_{z} [GeV];",100,0,50);

    hReso_Sigma_HFS = NEW1D("hReso_Sigma_HFS","HFS E-p_{z} Resolution;#Delta#Sigma_{h}/#Sigma_{h};");
    hReso_pt_HFS = NEW1D("hReso_pt_HFS","HFS p_{T} Resolution;#Delta p_{T,h}/p_{T,h};");
}


// Fill function
inline float rel(float reco, float truth)
{
    return (truth != 0 ? (reco - truth) / truth : 0);
}

inline void fill_benchmarks(
    float x_truth, float x_ele, float x_jb, float x_da, float x_sig, float x_esig,
    float y_truth, float y_ele, float y_jb, float y_da, float y_sig, float y_esig,
    float Q2_truth, float Q2_ele, float Q2_jb, float Q2_da, float Q2_sig, float Q2_esig,
    float E, float theta, float sigma_h, float pt_had,
    float sigma_h_true, float pt_had_true,
    float weight = 1.0)
{
    hResoX_electron->Fill(rel(x_ele,x_truth),weight);
    hResoX_jb->Fill(rel(x_jb,x_truth),weight);
    hResoX_da->Fill(rel(x_da,x_truth),weight);
    hResoX_sigma->Fill(rel(x_sig,x_truth),weight);
    hResoX_esigma->Fill(rel(x_esig,x_truth),weight);

    hResoY_electron->Fill(rel(y_ele,y_truth),weight);
    hResoY_jb->Fill(rel(y_jb,y_truth),weight);
    hResoY_da->Fill(rel(y_da,y_truth),weight);
    hResoY_sigma->Fill(rel(y_sig,y_truth),weight);
    hResoY_esigma->Fill(rel(y_esig,y_truth),weight);

    hResoQ2_electron->Fill(rel(Q2_ele,Q2_truth),weight);
    hResoQ2_jb->Fill(rel(Q2_jb,Q2_truth),weight);
    hResoQ2_da->Fill(rel(Q2_da,Q2_truth),weight);
    hResoQ2_sigma->Fill(rel(Q2_sig,Q2_truth),weight);
    hResoQ2_esigma->Fill(rel(Q2_esig,Q2_truth),weight);

    hResoX_2D_electron->Fill(y_truth,rel(x_ele,x_truth),weight);
    hResoX_2D_jb->Fill(y_truth,rel(x_jb,x_truth),weight);
    hResoX_2D_da->Fill(y_truth,rel(x_da,x_truth),weight);
    hResoX_2D_sigma->Fill(y_truth,rel(x_sig,x_truth),weight);
    hResoX_2D_esigma->Fill(y_truth,rel(x_esig,x_truth),weight);

    hResoY_2D_electron->Fill(y_truth,rel(y_ele,y_truth),weight);
    hResoY_2D_jb->Fill(y_truth,rel(y_jb,y_truth),weight);
    hResoY_2D_da->Fill(y_truth,rel(y_da,y_truth),weight);
    hResoY_2D_sigma->Fill(y_truth,rel(y_sig,y_truth),weight);
    hResoY_2D_esigma->Fill(y_truth,rel(y_esig,y_truth),weight);

    hResoQ2_2D_electron->Fill(y_truth,rel(Q2_ele,Q2_truth),weight);
    hResoQ2_2D_jb->Fill(y_truth,rel(Q2_jb,Q2_truth),weight);
    hResoQ2_2D_da->Fill(y_truth,rel(Q2_da,Q2_truth),weight);
    hResoQ2_2D_sigma->Fill(y_truth,rel(Q2_sig,Q2_truth),weight);
    hResoQ2_2D_esigma->Fill(y_truth,rel(Q2_esig,Q2_truth),weight);

    hReco_E_electron->Fill(E,weight);
    hReco_theta_electron->Fill(theta,weight);
    hReco_Sigma_HFS->Fill(sigma_h,weight);
    hReco_pt_HFS->Fill(pt_had,weight);
    hReco_Sigma_total->Fill(sigma_h+E*(1-cos(theta)),weight);

    hReso_Sigma_HFS->Fill(rel(sigma_h,sigma_h_true),weight);
    hReso_pt_HFS->Fill(rel(pt_had,pt_had_true),weight);
}


// Write out PDF file 
inline void write_benchmarks_pdf(const TString &pdfname, bool logx=false)
{
    TCanvas *c = new TCanvas("c","Benchmark Results",1200,900);

    c->Print(pdfname + "[");

    auto draw_1D = [&](std::initializer_list<TH1*> h) {
        c->Clear();
        c->Divide(3,2);
        int i=1;
        for (auto *hist : h) {
            c->cd(i++);
            hist->Draw("hist");
        }
        c->Print(pdfname);
    };

    auto draw_2D = [&](std::initializer_list<TH2*> h) {
        c->Clear();
        c->Divide(3,2);
        int i=1;
        for (auto *hist : h) {
            c->cd(i++);
            if (logx) gPad->SetLogx();
            hist->Draw("colz");
        }
        c->Print(pdfname);
    };

    // 1D kinematics pages
    draw_1D({hResoX_electron,hResoX_jb,hResoX_da,hResoX_sigma,hResoX_esigma});
    draw_1D({hResoY_electron,hResoY_jb,hResoY_da,hResoY_sigma,hResoY_esigma});
    draw_1D({hResoQ2_electron,hResoQ2_jb,hResoQ2_da,hResoQ2_sigma,hResoQ2_esigma});

    // 2D kinematics pages
    draw_2D({hResoX_2D_electron,hResoX_2D_jb,hResoX_2D_da,hResoX_2D_sigma,hResoX_2D_esigma});
    draw_2D({hResoY_2D_electron,hResoY_2D_jb,hResoY_2D_da,hResoY_2D_sigma,hResoY_2D_esigma});
    draw_2D({hResoQ2_2D_electron,hResoQ2_2D_jb,hResoQ2_2D_da,hResoQ2_2D_sigma,hResoQ2_2D_esigma});

    // input variables page
    draw_1D({hReco_Sigma_HFS,hReco_pt_HFS,hReco_Sigma_total,hReso_Sigma_HFS,hReso_pt_HFS});
    draw_1D({hReco_E_electron,hReco_theta_electron});

    c->Print(pdfname + "]");
}

#endif

