#include <vector>
#include <iostream>
#include <cmath>
#include <TMath.h>

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TStyle.h"

using namespace std;

void parton_qa(const char* inputFileName = "/eos/cms/store/group/phys_heavyions/xiaoyul/wenbin/sample/sample.root", 
               const char* outputDir = "/eos/cms/store/group/phys_heavyions/xiaoyul/wenbin/anaOutput/qa/") {
    
    // Open input file
    TFile* inFile = TFile::Open(inputFileName, "READ");
    if (!inFile || inFile->IsZombie()) {
        cout << "Error: Cannot open input file " << inputFileName << endl;
        return;
    }
    
    // Get the tree
    TTree* inTree = (TTree*)inFile->Get("trackTree");
    if (!inTree) {
        cout << "Error: Cannot find trackTree in input file" << endl;
        return;
    }
    
    // Set up branch addresses
    vector<int>* b_par_pdgid = nullptr;
    vector<float>* b_par_px = nullptr;
    vector<float>* b_par_py = nullptr;
    vector<float>* b_par_pz = nullptr;
    vector<float>* b_par_e = nullptr;
    vector<float>* b_par_x = nullptr;
    vector<float>* b_par_y = nullptr;
    vector<float>* b_par_z = nullptr;
    vector<float>* b_par_t = nullptr;
    vector<int>* b_par_color1 = nullptr;
    vector<int>* b_par_color2 = nullptr;
    
    vector<int>* b_par_pdgid_after_zpc = nullptr;
    vector<float>* b_par_px_after_zpc = nullptr;
    vector<float>* b_par_py_after_zpc = nullptr;
    vector<float>* b_par_pz_after_zpc = nullptr;
    vector<float>* b_par_e_after_zpc = nullptr;
    vector<float>* b_par_x_after_zpc = nullptr;
    vector<float>* b_par_y_after_zpc = nullptr;
    vector<float>* b_par_z_after_zpc = nullptr;
    vector<float>* b_par_t_after_zpc = nullptr;
    vector<int>* b_par_color1_after_zpc = nullptr;
    vector<int>* b_par_color2_after_zpc = nullptr;
    
    int b_total_collisions = 0;
    
    vector<float>* b_px = nullptr;
    vector<float>* b_py = nullptr;
    vector<float>* b_pz = nullptr;
    vector<float>* b_m = nullptr;
    vector<int>* b_pid = nullptr;
    vector<int>* b_chg = nullptr;
    
    vector<float>* b_genJetEta = nullptr;
    vector<float>* b_genJetPt = nullptr;
    vector<float>* b_genJetPhi = nullptr;
    vector<int>* b_genJetChargedMultiplicity = nullptr;
    
    // Set branch addresses
    inTree->SetBranchAddress("par_pdgid", &b_par_pdgid);
    inTree->SetBranchAddress("par_px", &b_par_px);
    inTree->SetBranchAddress("par_py", &b_par_py);
    inTree->SetBranchAddress("par_pz", &b_par_pz);
    inTree->SetBranchAddress("par_e", &b_par_e);
    inTree->SetBranchAddress("par_x", &b_par_x);
    inTree->SetBranchAddress("par_y", &b_par_y);
    inTree->SetBranchAddress("par_z", &b_par_z);
    inTree->SetBranchAddress("par_t", &b_par_t);
    inTree->SetBranchAddress("par_color1", &b_par_color1);
    inTree->SetBranchAddress("par_color2", &b_par_color2);
    
    inTree->SetBranchAddress("par_pdgid_after_zpc", &b_par_pdgid_after_zpc);
    inTree->SetBranchAddress("par_px_after_zpc", &b_par_px_after_zpc);
    inTree->SetBranchAddress("par_py_after_zpc", &b_par_py_after_zpc);
    inTree->SetBranchAddress("par_pz_after_zpc", &b_par_pz_after_zpc);
    inTree->SetBranchAddress("par_e_after_zpc", &b_par_e_after_zpc);
    inTree->SetBranchAddress("par_x_after_zpc", &b_par_x_after_zpc);
    inTree->SetBranchAddress("par_y_after_zpc", &b_par_y_after_zpc);
    inTree->SetBranchAddress("par_z_after_zpc", &b_par_z_after_zpc);
    inTree->SetBranchAddress("par_t_after_zpc", &b_par_t_after_zpc);
    inTree->SetBranchAddress("par_color1_after_zpc", &b_par_color1_after_zpc);
    inTree->SetBranchAddress("par_color2_after_zpc", &b_par_color2_after_zpc);
    
    inTree->SetBranchAddress("total_collisions", &b_total_collisions);
    
    inTree->SetBranchAddress("px", &b_px);
    inTree->SetBranchAddress("py", &b_py);
    inTree->SetBranchAddress("pz", &b_pz);
    inTree->SetBranchAddress("m", &b_m);
    inTree->SetBranchAddress("pid", &b_pid);
    inTree->SetBranchAddress("chg", &b_chg);
    
    inTree->SetBranchAddress("genJetEta", &b_genJetEta);
    inTree->SetBranchAddress("genJetPt", &b_genJetPt);
    inTree->SetBranchAddress("genJetPhi", &b_genJetPhi);
    inTree->SetBranchAddress("genJetChargedMultiplicity", &b_genJetChargedMultiplicity);
    
    // Create histograms
    TH2D* hNchjVsNcollision = new TH2D("hNchjVsNcollision", "Leading Jet Charged Multiplicity vs Total Collisions; Total Collisions; N_{ch}^{jet}", 
                                       20, 0, 20, 100, 0, 100);//x-axis is total collisions, y-axis is leading jet charged multiplicity
    
    TH2D* hNparVsNparticles = new TH2D("hNparVsNparticles", "Number of Partons vs Number of Particles; N_{particles}; N_{partons}", 
                                       600, 0, 600, 200, 0, 200);//x-axis is number of particles, y-axis is number of partons, in an event.
    
    // genJet distributions
    TH1D* hGenJetPt = new TH1D("hGenJetPt", "GenJet p_{T} Distribution; p_{T} (GeV/c); Counts", 100, 0, 1000);
    TH1D* hGenJetEta = new TH1D("hGenJetEta", "GenJet #eta Distribution; #eta; Counts", 100, -3, 3);
    TH1D* hGenJetPhi = new TH1D("hGenJetPhi", "GenJet #phi Distribution; #phi; Counts", 100, -TMath::Pi(), TMath::Pi());
    TH1D* hGenJetChargedMultiplicity = new TH1D("hGenJetChargedMultiplicity", "GenJet Charged Multiplicity Distribution; N_{ch}; Counts", 100, 0, 100);
    
    // Total collisions distribution
    TH1D* hTotalCollisions = new TH1D("hTotalCollisions", "Total Collisions Distribution; N_{collisions}; Counts", 50, 0, 50);
    
    // Np vs leading jet multiplicity
    TH2D* hNpVsLeadingJetMult = new TH2D("hNpVsLeadingJetMult", "Number of Partons vs Leading Jet Multiplicity; Leading Jet Multiplicity; N_{partons}", 100, 0, 100, 200, 0, 200);
    
    // For events with 0 collisions
    TH2D* hPxBeforeVsAfter_0coll = new TH2D("hPxBeforeVsAfter_0coll", "Parton p_{x} Before vs After ZPC (0 collisions); p_{x} Before ZPC (GeV/c); p_{x} After ZPC (GeV/c)", 
                                            50, -10, 10, 50, -10, 10);
    TH2D* hPyBeforeVsAfter_0coll = new TH2D("hPyBeforeVsAfter_0coll", "Parton p_{y} Before vs After ZPC (0 collisions); p_{y} Before ZPC (GeV/c); p_{y} After ZPC (GeV/c)", 
                                            50, -10, 10, 50, -10, 10);
    TH2D* hPzBeforeVsAfter_0coll = new TH2D("hPzBeforeVsAfter_0coll", "Parton p_{z} Before vs After ZPC (0 collisions); p_{z} Before ZPC (GeV/c); p_{z} After ZPC (GeV/c)", 
                                            50, -50, 50, 50, -50, 50);
    
    // For events with nonzero collisions
    TH2D* hPxBeforeVsAfter_ncoll = new TH2D("hPxBeforeVsAfter_ncoll", "Parton p_{x} Before vs After ZPC (nonzero collisions); p_{x} Before ZPC (GeV/c); p_{x} After ZPC (GeV/c)", 
                                            50, -10, 10, 50, -10, 10);
    TH2D* hPyBeforeVsAfter_ncoll = new TH2D("hPyBeforeVsAfter_ncoll", "Parton p_{y} Before vs After ZPC (nonzero collisions); p_{y} Before ZPC (GeV/c); p_{y} After ZPC (GeV/c)", 
                                            50, -10, 10, 50, -10, 10);
    TH2D* hPzBeforeVsAfter_ncoll = new TH2D("hPzBeforeVsAfter_ncoll", "Parton p_{z} Before vs After ZPC (nonzero collisions); p_{z} Before ZPC (GeV/c); p_{z} After ZPC (GeV/c)", 
                                            50, -50, 50, 50, -50, 50);
    
    // Event loop
    Long64_t nEntries = inTree->GetEntries();
    cout << "Processing " << nEntries << " events..." << endl;
    
    for (Long64_t ientry = 0; ientry < nEntries; ientry++) {
        inTree->GetEntry(ientry);
        
        if (ientry % 1000 == 0) {
            cout << "Processing event " << ientry << "/" << nEntries << endl;
        }
        
        // Check if we have jets
        // Check if first jet (leading jet) meets criteria: |eta|<1.6, pt>550
        if (b_genJetEta->size() == 0) continue;
        
        if (fabs((*b_genJetEta)[0]) >= 1.6 || (*b_genJetPt)[0] <= 550) continue;
        
        // Calculate Nchj (jet charged multiplicity) for leading jet
        int Nchj = 0;
        for (size_t i = 0; i < b_px->size(); i++) {
            // Check if particle is charged
            if ((*b_chg)[i] != 0) {
                float pt = sqrt((*b_px)[i] * (*b_px)[i] + (*b_py)[i] * (*b_py)[i]);
                float eta = -log(tan(atan2(pt, (*b_pz)[i]) / 2.0));
                
                // Check pt and eta cuts only
                if (pt > 0.3 && fabs(eta) < 2.4) {
                    Nchj++;
                }
            }
        }
        
        // Fill Nchj vs Ncollision
        hNchjVsNcollision->Fill(b_total_collisions, Nchj);
        
        // Fill Npar vs Nparticles
        int Npar = b_par_pdgid->size();
        int Nparticles = b_px->size();
        hNparVsNparticles->Fill(Nparticles, Npar);
        
        // Fill Np vs leading jet multiplicity (first jet entry is leading jet)
        hNpVsLeadingJetMult->Fill(Nchj, Npar);
        
        // Fill genJet distributions
        for (size_t i = 0; i < b_genJetEta->size(); i++) {
            hGenJetPt->Fill((*b_genJetPt)[i]);
            hGenJetEta->Fill((*b_genJetEta)[i]);
            hGenJetPhi->Fill((*b_genJetPhi)[i]);
            hGenJetChargedMultiplicity->Fill((*b_genJetChargedMultiplicity)[i]);
        }
        
        // Fill total collisions distribution
        hTotalCollisions->Fill(b_total_collisions);
        
        // Fill momentum correlation plots per parton based on collision count
        if (b_total_collisions == 0) {
            // Fill per parton for events with 0 collisions
            for (size_t i = 0; i < b_par_pdgid->size(); i++) {
                hPxBeforeVsAfter_0coll->Fill((*b_par_px)[i], (*b_par_px_after_zpc)[i]);
                hPyBeforeVsAfter_0coll->Fill((*b_par_py)[i], (*b_par_py_after_zpc)[i]);
                hPzBeforeVsAfter_0coll->Fill((*b_par_pz)[i], (*b_par_pz_after_zpc)[i]);
            }
        } else {
            // Fill per parton for events with nonzero collisions
            for (size_t i = 0; i < b_par_pdgid->size(); i++) {
                hPxBeforeVsAfter_ncoll->Fill((*b_par_px)[i], (*b_par_px_after_zpc)[i]);
                hPyBeforeVsAfter_ncoll->Fill((*b_par_py)[i], (*b_par_py_after_zpc)[i]);
                hPzBeforeVsAfter_ncoll->Fill((*b_par_pz)[i], (*b_par_pz_after_zpc)[i]);
            }
        }
    }
    
    // Create output file based on input filename (same as parton_v2.C)
    TString baseFileName = TString(inputFileName);
    baseFileName.ReplaceAll(".root", "");
    baseFileName.ReplaceAll("/", "_");
    baseFileName.ReplaceAll("eos_cms_store_group_phys_heavyions_xiaoyul_wenbin_sample_", "");
    
    TString outputFileName = Form("%sparton_qa_output_%s.root", outputDir, baseFileName.Data());
    TFile* outFile = new TFile(outputFileName, "RECREATE");
    
    // Write histograms
    hNchjVsNcollision->Write();
    hNparVsNparticles->Write();
    
    // Write Np vs leading jet multiplicity
    hNpVsLeadingJetMult->Write();
    
    // Write genJet distributions
    hGenJetPt->Write();
    hGenJetEta->Write();
    hGenJetPhi->Write();
    hGenJetChargedMultiplicity->Write();
    
    // Write total collisions distribution
    hTotalCollisions->Write();
    
    hPxBeforeVsAfter_0coll->Write();
    hPyBeforeVsAfter_0coll->Write();
    hPzBeforeVsAfter_0coll->Write();
    hPxBeforeVsAfter_ncoll->Write();
    hPyBeforeVsAfter_ncoll->Write();
    hPzBeforeVsAfter_ncoll->Write();
    
    outFile->Close();
    inFile->Close();
    
    cout << "QA analysis completed. Output saved to: " << outputFileName << endl;
}
