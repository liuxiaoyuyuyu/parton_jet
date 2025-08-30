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
                                       50, 0, 100, 50, 0, 100);
    
    TH2D* hNparVsNparticles = new TH2D("hNparVsNparticles", "Number of Partons vs Number of Particles; N_{particles}; N_{partons}", 
                                       50, 0, 200, 50, 0, 200);
    
    // For events with 0 collisions
    TH2D* hPxBeforeVsAfter_0coll = new TH2D("hPxBeforeVsAfter_0coll", "Parton <p_{x}> Before vs After ZPC (0 collisions); <p_{x}> Before ZPC (GeV/c); <p_{x}> After ZPC (GeV/c)", 
                                            50, -10, 10, 50, -10, 10);
    TH2D* hPyBeforeVsAfter_0coll = new TH2D("hPyBeforeVsAfter_0coll", "Parton <p_{y}> Before vs After ZPC (0 collisions); <p_{y}> Before ZPC (GeV/c); <p_{y}> After ZPC (GeV/c)", 
                                            50, -10, 10, 50, -10, 10);
    TH2D* hPzBeforeVsAfter_0coll = new TH2D("hPzBeforeVsAfter_0coll", "Parton <p_{z}> Before vs After ZPC (0 collisions); <p_{z}> Before ZPC (GeV/c); <p_{z}> After ZPC (GeV/c)", 
                                            50, -50, 50, 50, -50, 50);
    
    // For events with nonzero collisions
    TH2D* hPxBeforeVsAfter_ncoll = new TH2D("hPxBeforeVsAfter_ncoll", "Parton <p_{x}> Before vs After ZPC (nonzero collisions); <p_{x}> Before ZPC (GeV/c); <p_{x}> After ZPC (GeV/c)", 
                                            50, -10, 10, 50, -10, 10);
    TH2D* hPyBeforeVsAfter_ncoll = new TH2D("hPyBeforeVsAfter_ncoll", "Parton <p_{y}> Before vs After ZPC (nonzero collisions); <p_{y}> Before ZPC (GeV/c); <p_{y}> After ZPC (GeV/c)", 
                                            50, -10, 10, 50, -10, 10);
    TH2D* hPzBeforeVsAfter_ncoll = new TH2D("hPzBeforeVsAfter_ncoll", "Parton <p_{z}> Before vs After ZPC (nonzero collisions); <p_{z}> Before ZPC (GeV/c); <p_{z}> After ZPC (GeV/c)", 
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
        if (b_genJetEta->size() == 0) continue;
        
        // Find leading jet that meets criteria: |eta|<1.6, pt>550
        int leadingJetIdx = -1;
        float maxPt = 0;
        
        for (size_t i = 0; i < b_genJetEta->size(); i++) {
            if (fabs((*b_genJetEta)[i]) < 1.6 && (*b_genJetPt)[i] > 550) {
                if ((*b_genJetPt)[i] > maxPt) {
                    maxPt = (*b_genJetPt)[i];
                    leadingJetIdx = i;
                }
            }
        }
        
        // Skip if no qualifying jet found
        if (leadingJetIdx == -1) continue;
        
        // Calculate Nchj (jet charged multiplicity) for leading jet
        int Nchj = 0;
        for (size_t i = 0; i < b_px->size(); i++) {
            // Check if particle is charged
            if ((*b_chg)[i] != 0) {
                float pt = sqrt((*b_px)[i] * (*b_px)[i] + (*b_py)[i] * (*b_py)[i]);
                float eta = -log(tan(atan2(pt, (*b_pz)[i]) / 2.0));
                
                // Check pt and eta cuts
                if (pt > 0.3 && fabs(eta) < 2.4) {
                    // Check if particle belongs to leading jet (simple deltaR cut)
                    float jetEta = (*b_genJetEta)[leadingJetIdx];
                    float jetPhi = (*b_genJetPhi)[leadingJetIdx];
                    float partPhi = atan2((*b_py)[i], (*b_px)[i]);
                    
                    float deltaEta = eta - jetEta;
                    float deltaPhi = partPhi - jetPhi;
                    if (deltaPhi > TMath::Pi()) deltaPhi -= 2 * TMath::Pi();
                    if (deltaPhi < -TMath::Pi()) deltaPhi += 2 * TMath::Pi();
                    
                    float deltaR = sqrt(deltaEta * deltaEta + deltaPhi * deltaPhi);
                    if (deltaR < 0.4) { // R=0.4 jet cone
                        Nchj++;
                    }
                }
            }
        }
        
        // Fill Nchj vs Ncollision
        hNchjVsNcollision->Fill(b_total_collisions, Nchj);
        
        // Fill Npar vs Nparticles
        int Npar = b_par_pdgid->size();
        int Nparticles = b_px->size();
        hNparVsNparticles->Fill(Nparticles, Npar);
        
        // Calculate average parton momenta before and after ZPC
        double avgPxBefore = 0, avgPyBefore = 0, avgPzBefore = 0;
        double avgPxAfter = 0, avgPyAfter = 0, avgPzAfter = 0;
        
        if (Npar > 0) {
            for (size_t i = 0; i < b_par_pdgid->size(); i++) {
                avgPxBefore += (*b_par_px)[i];
                avgPyBefore += (*b_par_py)[i];
                avgPzBefore += (*b_par_pz)[i];
            }
            avgPxBefore /= Npar;
            avgPyBefore /= Npar;
            avgPzBefore /= Npar;
        }
        
        if (b_par_pdgid_after_zpc->size() > 0) {
            for (size_t i = 0; i < b_par_pdgid_after_zpc->size(); i++) {
                avgPxAfter += (*b_par_px_after_zpc)[i];
                avgPyAfter += (*b_par_py_after_zpc)[i];
                avgPzAfter += (*b_par_pz_after_zpc)[i];
            }
            avgPxAfter /= b_par_pdgid_after_zpc->size();
            avgPyAfter /= b_par_pdgid_after_zpc->size();
            avgPzAfter /= b_par_pdgid_after_zpc->size();
        }
        
        // Fill momentum correlation plots based on collision count
        if (b_total_collisions == 0) {
            hPxBeforeVsAfter_0coll->Fill(avgPxBefore, avgPxAfter);
            hPyBeforeVsAfter_0coll->Fill(avgPyBefore, avgPyAfter);
            hPzBeforeVsAfter_0coll->Fill(avgPzBefore, avgPzAfter);
        } else {
            hPxBeforeVsAfter_ncoll->Fill(avgPxBefore, avgPxAfter);
            hPyBeforeVsAfter_ncoll->Fill(avgPyBefore, avgPyAfter);
            hPzBeforeVsAfter_ncoll->Fill(avgPzBefore, avgPzAfter);
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
