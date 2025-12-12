#include "Boost.h"
#include "Beam.h"

// PODIO
#include "podio/Frame.h"
#include "podio/ROOTReader.h"

// DATA MODEL
#include "edm4eic/InclusiveKinematicsCollection.h"
#include "edm4hep/MCParticleCollection.h"
#include "edm4eic/ReconstructedParticleCollection.h"
#include "edm4eic/MCRecoClusterParticleAssociationCollection.h"
#include "edm4eic/MCRecoParticleAssociationCollection.h"
#include "edm4eic/HadronicFinalStateCollection.h"
#include "edm4eic/ClusterCollection.h"
#include "edm4hep/Vector3f.h"

// GENERAL
#include <filesystem>
#include <vector>

// BENCHMARKS
#include "benchmarks.h"

using ROOT::Math::PxPyPzEVector;

std::vector<float> calc_elec_method(float E, float theta, float pt_had, float sigma_h, float E_ebeam, float E_pbeam);
std::vector<float> calc_jb_method(float E, float theta, float pt_had, float sigma_h, float E_ebeam, float E_pbeam);
std::vector<float> calc_da_method(float E, float theta, float pt_had, float sigma_h, float E_ebeam, float E_pbeam);
std::vector<float> calc_sig_method(float E, float theta, float pt_had, float sigma_h, float E_ebeam, float E_pbeam);
std::vector<float> calc_esig_method(float E, float theta, float pt_had, float sigma_h, float E_ebeam, float E_pbeam);

TLorentzVector TransformLabToHeadOnFrame(TLorentzVector eBeam, TLorentzVector pBeam, TLorentzVector Lvec);
std::vector<std::string> get_filevector(std::string filename);

void Podio_ExtractYield_HeraBinning_5x130() {
  // Settings
  Float_t E_ebeam = 5;
  Float_t E_pbeam = 130;
  Float_t target_lumi = 1.;
  Float_t m_e = 0.000511;
  Float_t m_p = 0.938;
  double xAngle = 25e-3;
  Float_t min_Empz = 6;
  Float_t max_Empz = 12;
  PxPyPzEVector pni, ei;
  ei.SetPxPyPzE(0, 0, -E_ebeam, sqrt(E_ebeam*E_ebeam+m_e*m_e));
  pni.SetPxPyPzE(-1*E_pbeam*TMath::Sin(xAngle), 0, E_pbeam*TMath::Cos(xAngle), sqrt(E_pbeam*E_pbeam+m_p*m_p));
  
  // paths to file lists
  std::vector<TString> filenames = {
    "../../RECO/ep_5x130_Q2_1_10/filelist.txt",
    "../../RECO/ep_5x130_Q2_10_100/filelist.txt",
    "../../RECO/ep_5x130_Q2_100_1000/filelist.txt",
    "../../RECO/ep_5x130_Q2_1000_10000/filelist.txt"
  };
  
  std::vector<double> gen_xsecs = {0.54441532506348078,3.5926395427926944e-2,8.8862696458699047e-4,8.2239104878847066e-7}; // ub
  
  
  double Q2Bins[] = {1.721, 2.324, 3.074, 3.969, 5.408, 7.433, 9.219,
		     10.954, 13.412, 16.432, 19.899, 24.372, 30.741, 39.686, 51.962, 64.807, 79.373,
		     103.923, 134.164, 173.205, 223.607, 273.861, 346.41, 447.214, 570.088, 721.11, 894.443,
		     1095.45, 1341.64, 1677.05, 2236.07, 2738.61, 3464.1, 4472.14, 5700.88, 7211.1, 8944.43};
  
  double xBins[] = {0.0001, 0.0001585, 0.0002512, 0.000398, 0.000631,
		    0.001, 0.001585, 0.002512, 0.00398, 0.00631,
		    0.01, 0.01585, 0.02512, 0.0398, 0.0631,
		    0.1, 0.145, 0.209, 0.316, 0.501, 1};
  const int nBinsX = (sizeof(xBins)/sizeof(*xBins)) - 1;
  const int nBinsQ2 = (sizeof(Q2Bins)/sizeof(*Q2Bins)) - 1;
  
  // Create 2D yield histograms
  TString histTitle = TString::Format("Luminosity = %.3f fb-1;x; Q^{2}",target_lumi);
  TH2D* hQ2vsX_yield_gen = new TH2D("hQ2vsX_yield_gen", histTitle, nBinsX, xBins, nBinsQ2, Q2Bins);
  
  TH2D* hQ2vsX_yield_rec_ele = new TH2D("hQ2vsX_yield_rec_ele", histTitle, nBinsX, xBins, nBinsQ2, Q2Bins);
  TH2D* hQ2vsX_yield_rec_jb= new TH2D("hQ2vsX_yield_rec_jb", histTitle, nBinsX, xBins, nBinsQ2, Q2Bins);
  TH2D* hQ2vsX_yield_rec_da  = new TH2D("hQ2vsX_yield_rec_da", histTitle, nBinsX, xBins, nBinsQ2, Q2Bins);
  TH2D* hQ2vsX_yield_rec_esig  = new TH2D("hQ2vsX_yield_rec_esig", histTitle, nBinsX, xBins, nBinsQ2, Q2Bins);
  TH2D* hQ2vsX_yield_rec_sig  = new TH2D("hQ2vsX_yield_rec_sig", histTitle, nBinsX, xBins, nBinsQ2, Q2Bins);
  
  TH2D* hQ2vsX_yield_gen_rec_ele = new TH2D("hQ2vsX_yield_gen_rec_ele", histTitle, nBinsX, xBins, nBinsQ2, Q2Bins);
  TH2D* hQ2vsX_yield_gen_rec_jb = new TH2D("hQ2vsX_yield_gen_rec_jb", histTitle, nBinsX, xBins, nBinsQ2, Q2Bins);
  TH2D* hQ2vsX_yield_gen_rec_da  = new TH2D("hQ2vsX_yield_gen_rec_da", histTitle, nBinsX, xBins, nBinsQ2, Q2Bins);
  TH2D* hQ2vsX_yield_gen_rec_esig  = new TH2D("hQ2vsX_yield_gen_rec_esig", histTitle, nBinsX, xBins, nBinsQ2, Q2Bins);
  TH2D* hQ2vsX_yield_gen_rec_sig  = new TH2D("hQ2vsX_yield_gen_rec_sig", histTitle, nBinsX, xBins, nBinsQ2, Q2Bins);
  
  Float_t x_truth, x_ele, x_jb, x_da, x_sig, x_esig;
  Float_t y_truth, y_ele, y_jb, y_da, y_sig, y_esig;
  Float_t Q2_truth, Q2_ele, Q2_jb, Q2_da, Q2_sig, Q2_esig;
  Float_t E, theta, sigma_h, pt_had;
  Float_t sigma_h_true, pt_had_true;
  
  
  init_benchmarks();

  int file_index = 0;
  for (auto filename : filenames){
    auto inFiles = get_filevector(filenames[file_index].Data());
    
    auto reader = podio::ROOTReader();
    reader.openFiles(inFiles);
    double nEntries = reader.getEntries("events");
    cout << nEntries << " events found" << endl;
    
    double weight = target_lumi/((nEntries)/(gen_xsecs[file_index]*1e9));
    cout << "Weight for file set " << file_index << ": " << weight << endl;
    
    
    for (size_t i = 0; i < reader.getEntries("events"); i++) {// begin event loop
      const auto event = podio::Frame(reader.readNextEntry("events"));
      if (i%100==0) cout << i << " events processed" << endl;
      
      // Retrieve Inclusive Kinematics Collections
      auto& kin_truth = event.get<edm4eic::InclusiveKinematicsCollection>("InclusiveKinematicsTruth");
      auto& kin_electron = event.get<edm4eic::InclusiveKinematicsCollection>("InclusiveKinematicsElectron");
      auto& kin_jb = event.get<edm4eic::InclusiveKinematicsCollection>("InclusiveKinematicsJB");
      
      // Retrieve Scattered electron and HFS
      // replace "Truth" with "EMinusPz" below for basic eicrecon electron-finder implementation
      auto& eleCollection = event.get<edm4eic::ReconstructedParticleCollection>("ScatteredElectronsTruth");
      // auto& hfsCollection = event.get<edm4eic::HadronicFinalStateCollection>("HadronicFinalState");
      auto& ecalClusters = event.get<edm4eic::ClusterCollection>("EcalClusters");
      
      auto& rcparts = event.get<edm4eic::ReconstructedParticleCollection>("ReconstructedParticles");
      
      if (kin_truth.empty()) continue;
      
      // Store kinematics from InclusiveKinematics branches
      x_truth = kin_truth.x()[0];
      y_truth = kin_truth.y()[0];
      Q2_truth = kin_truth.Q2()[0];
      
      // Fill truth histograms before any cuts/selection criteria
      hQ2vsX_yield_gen->Fill(x_truth, Q2_truth, weight);
      
      auto boost = eicrecon::determine_boost(ei, pni);

      // Only consider events where it's possible to reconstruct both electron and HFS
      if (kin_electron.empty() || kin_jb.empty()) continue;
      
      PxPyPzEVector scat_ele;
      E = eleCollection[0].getEnergy();
      auto& ele_momentum = eleCollection[0].getMomentum();
      scat_ele.SetPxPyPzE(ele_momentum.x, ele_momentum.y, ele_momentum.z, E);
      theta = scat_ele.Theta();
      
      // Comment these out since we're getting the HFS manually
      // sigma_h = hfsCollection[0].getSigma();
      // pt_had = hfsCollection[0].getPT();
      
      PxPyPzEVector ebeam_boosted = boost(ei);
      PxPyPzEVector pbeam_boosted = boost(pni);
      
      double pxsum = 0;
      double pysum = 0;
      double pzsum = 0;
      double Esum  = 0;
      const auto ef_rc_id{eleCollection[0].getObjectID().index};
      for (const auto p : rcparts) {
	bool isHadron = true;
	// Check if it's the scattered electron
      if (p.getObjectID().index != ef_rc_id) {
	// Lorentz vector in lab frame
	PxPyPzEVector hf_lab(p.getMomentum().x, p.getMomentum().y, p.getMomentum().z, p.getEnergy());
	// Boost to collinear frame
	PxPyPzEVector hf_boosted = boost(hf_lab);
	// PxPyPzEVector hf_boosted = hf_lab;// HFS no boost

	pxsum += hf_boosted.Px();
	pysum += hf_boosted.Py();
	pzsum += hf_boosted.Pz();
	Esum += hf_boosted.E();
      }
    }
    sigma_h = Esum - pzsum;
    pt_had = sqrt(pxsum*pxsum + pysum*pysum);
    sigma_h_true = y_truth*2*E_ebeam;// stored for benchmarking only
    pt_had_true = sqrt(Q2_truth*(1 - y_truth));// stored for benchmarking only
   
    // Calculate kinematics manually
    std::vector<float> elec_reco = calc_elec_method(E, theta, pt_had, sigma_h, E_ebeam, E_pbeam);
    std::vector<float> jb_reco = calc_jb_method(E, theta, pt_had, sigma_h, E_ebeam, E_pbeam);
    std::vector<float> da_reco = calc_da_method(E, theta, pt_had, sigma_h, E_ebeam, E_pbeam);
    std::vector<float> sigma_reco = calc_sig_method(E, theta, pt_had, sigma_h, E_ebeam, E_pbeam);
    std::vector<float> esigma_reco = calc_esig_method(E, theta, pt_had, sigma_h, E_ebeam, E_pbeam);

    x_ele = elec_reco[0];
    x_jb = jb_reco[0];
    x_da = da_reco[0];
    x_sig = sigma_reco[0];
    x_esig = esigma_reco[0];
   
    y_ele = elec_reco[1];
    y_jb = jb_reco[1];
    y_da = da_reco[1];
    y_sig = sigma_reco[1];
    y_esig = esigma_reco[1];
   
    Q2_ele = elec_reco[2];
    Q2_jb = jb_reco[2];
    Q2_da = da_reco[2];
    Q2_sig = sigma_reco[2];
    Q2_esig = esigma_reco[2];
   

    // Some minimal cuts
    bool cuts = true;
    // will cut on e.g. y_ele later
    // cuts = cuts && (y_ele < 0.95);
    // cuts = cuts && (y_ele > 0.001);
    // cuts = cuts && (Q2_ele > 1);
    cuts = cuts && (sigma_h + E*(1-cos(theta)) > min_Empz);
    cuts = cuts && (sigma_h + E*(1-cos(theta)) < max_Empz);
    
    if (!cuts) continue;

    fill_benchmarks(x_truth, x_ele, x_jb, x_da, x_sig, x_esig,
		    y_truth, y_ele, y_jb, y_da, y_sig, y_esig,
		    Q2_truth, Q2_ele, Q2_jb, Q2_da, Q2_sig, Q2_esig,
		    E, theta, sigma_h, pt_had, sigma_h_true, pt_had_true,
		    weight);

    hQ2vsX_yield_rec_ele->Fill(x_ele, Q2_ele, weight);
    hQ2vsX_yield_rec_jb->Fill(x_jb, Q2_jb, weight);
    hQ2vsX_yield_rec_da->Fill(x_da, Q2_da, weight);
    hQ2vsX_yield_rec_esig->Fill(x_esig, Q2_esig, weight);
    hQ2vsX_yield_rec_sig->Fill(x_sig, Q2_sig, weight);

    int gen_bin = hQ2vsX_yield_gen->FindBin(x_truth, Q2_truth);
    if (hQ2vsX_yield_rec_ele->FindBin(x_ele,Q2_ele) == gen_bin) hQ2vsX_yield_gen_rec_ele->Fill(x_ele,Q2_ele,weight);
    if (hQ2vsX_yield_rec_jb->FindBin(x_jb,Q2_jb) == gen_bin) hQ2vsX_yield_gen_rec_jb->Fill(x_jb,Q2_jb,weight);
    if (hQ2vsX_yield_rec_da->FindBin(x_da,Q2_da) == gen_bin) hQ2vsX_yield_gen_rec_da->Fill(x_da,Q2_da,weight);
    if (hQ2vsX_yield_rec_esig->FindBin(x_esig,Q2_esig) == gen_bin) hQ2vsX_yield_gen_rec_esig->Fill(x_esig,Q2_esig,weight);
    if (hQ2vsX_yield_rec_sig->FindBin(x_sig,Q2_sig) == gen_bin) hQ2vsX_yield_gen_rec_sig->Fill(x_sig,Q2_sig,weight);
   
  }// end event loop
  file_index++;
  }
  // Drawing the histograms
  TFile *histFile = new TFile(Form("yield_histograms_podio_HeraBinning_%gx%g.root",E_ebeam,E_pbeam),"recreate");
  hQ2vsX_yield_gen->Write();
  hQ2vsX_yield_rec_ele->Write();
  hQ2vsX_yield_rec_jb->Write();
  hQ2vsX_yield_rec_da->Write();
  hQ2vsX_yield_rec_esig->Write();
  hQ2vsX_yield_rec_sig->Write();
  hQ2vsX_yield_gen_rec_ele->Write();
  hQ2vsX_yield_gen_rec_jb->Write();
  hQ2vsX_yield_gen_rec_da->Write();
  hQ2vsX_yield_gen_rec_esig->Write();
  hQ2vsX_yield_gen_rec_sig->Write();

  // close the output file
  histFile->Close();

  
  // write_benchmarks_pdf("benchmarks.pdf", /*logx=*/true);
  write_benchmarks_pdf(Form("benchmarks_HeraBinning_%gx%g.pdf",E_ebeam, E_pbeam), /*logx=*/false);
 
  cout << "Done!" << endl;
}

// electron method
std::vector<float> calc_elec_method(float E, float theta, float pt_had, float sigma_h, float E_ebeam, float E_pbeam) {
  float Q2  = 2.*E_ebeam*E*(1+TMath::Cos(theta));
  float y = 1. - (E/E_ebeam)*TMath::Sin(theta/2)*TMath::Sin(theta/2);
  float x = Q2/(4*E_ebeam*E_pbeam*y);
  return {x, y, Q2};
}

// jb method
std::vector<float> calc_jb_method(float E, float theta, float pt_had, float sigma_h, float E_ebeam, float E_pbeam) {
  float y = sigma_h/(2*E_ebeam);
  float Q2 = pt_had*pt_had / (1-y);
  float x = Q2/(4*E_ebeam*E_pbeam*y);
  return {x, y, Q2};
}

// float angle method
std::vector<float> calc_da_method(float E, float theta, float pt_had, float sigma_h, float E_ebeam, float E_pbeam) {
  float alpha_h = sigma_h/pt_had;
  float alpha_e = TMath::Tan(theta/2);
  float y = alpha_h / (alpha_e + alpha_h);
  float Q2 = 4*E_ebeam*E_ebeam / (alpha_e * (alpha_h + alpha_e));
  float x = Q2/(4*E_ebeam*E_pbeam*y);
  return {x, y, Q2};
}

// sigma method
std::vector<float> calc_sig_method(float E, float theta, float pt_had, float sigma_h, float E_ebeam, float E_pbeam) {
  float y = sigma_h/(sigma_h + E*(1 - TMath::Cos(theta))); 
  float Q2 = E*E*TMath::Sin(theta)*TMath::Sin(theta) / (1-y);
  float x = Q2/(4*E_ebeam*E_pbeam*y);
  return {x, y, Q2};
}

// e-sigma method
std::vector<float> calc_esig_method(float E, float theta, float pt_had, float sigma_h, float E_ebeam, float E_pbeam) {
  float Q2  = 2.*E_ebeam*E*(1+TMath::Cos(theta));
  float x = calc_sig_method(E,theta,pt_had,sigma_h,E_ebeam,E_pbeam)[0];
  float y = Q2/(4*E_ebeam*E_pbeam*x);
  return {x, y, Q2};
}

std::vector<std::string> get_filevector(std::string filename){
  std::vector<std::string> fileVector;
  std::ifstream in(filename);
  std::string file("");
  while (in >> file) fileVector.push_back(file.data());
  return fileVector;
}
