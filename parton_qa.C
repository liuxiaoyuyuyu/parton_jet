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
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TTreeReaderValue.h"

using namespace std;

void parton_qa(const char* inputFileName = "/eos/cms/store/group/phys_heavyions/xiaoyul/wenbin/sample/batch2/pp_parton_cascade_19607.root", 
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
    
    // Set up TTreeReader
    TTreeReader reader(inTree);
    
    // Set up TTreeReaderArrays for all branches
    TTreeReaderArray<int> b_par_pdgid(reader, "par_pdgid");
    TTreeReaderArray<float> b_par_px(reader, "par_px");
    TTreeReaderArray<float> b_par_py(reader, "par_py");
    TTreeReaderArray<float> b_par_pz(reader, "par_pz");
    TTreeReaderArray<float> b_par_e(reader, "par_e");
    TTreeReaderArray<float> b_par_x(reader, "par_x");
    TTreeReaderArray<float> b_par_y(reader, "par_y");
    TTreeReaderArray<float> b_par_z(reader, "par_z");
    TTreeReaderArray<float> b_par_t(reader, "par_t");
    TTreeReaderArray<int> b_par_color1(reader, "par_color1");
    TTreeReaderArray<int> b_par_color2(reader, "par_color2");
    
    TTreeReaderArray<int> b_par_pdgid_after_zpc(reader, "par_pdgid_after_zpc");
    TTreeReaderArray<float> b_par_px_after_zpc(reader, "par_px_after_zpc");
    TTreeReaderArray<float> b_par_py_after_zpc(reader, "par_py_after_zpc");
    TTreeReaderArray<float> b_par_pz_after_zpc(reader, "par_pz_after_zpc");
    TTreeReaderArray<float> b_par_e_after_zpc(reader, "par_e_after_zpc");
    TTreeReaderArray<float> b_par_x_after_zpc(reader, "par_x_after_zpc");
    TTreeReaderArray<float> b_par_y_after_zpc(reader, "par_y_after_zpc");
    TTreeReaderArray<float> b_par_z_after_zpc(reader, "par_z_after_zpc");
    TTreeReaderArray<float> b_par_t_after_zpc(reader, "par_t_after_zpc");
    TTreeReaderArray<int> b_par_color1_after_zpc(reader, "par_color1_after_zpc");
    TTreeReaderArray<int> b_par_color2_after_zpc(reader, "par_color2_after_zpc");
    
    TTreeReaderValue<int> b_total_collisions(reader, "total_collisions");
    
    TTreeReaderArray<float> b_px(reader, "px");
    TTreeReaderArray<float> b_py(reader, "py");
    TTreeReaderArray<float> b_pz(reader, "pz");
    TTreeReaderArray<float> b_m(reader, "m");
    TTreeReaderArray<int> b_pid(reader, "pid");
    TTreeReaderArray<int> b_chg(reader, "chg");
    
    TTreeReaderArray<float> b_genJetEta(reader, "genJetEta");
    TTreeReaderArray<float> b_genJetPt(reader, "genJetPt");
    TTreeReaderArray<float> b_genJetPhi(reader, "genJetPhi");
    TTreeReaderArray<int> b_genJetChargedMultiplicity(reader, "genJetChargedMultiplicity");
    
    // genDau branches for jet constituents
    TTreeReaderArray<vector<int> > b_genDau_chg(reader, "genDau_chg");
    TTreeReaderArray<vector<int> > b_genDau_pid(reader, "genDau_pid");
    TTreeReaderArray<vector<float> > b_genDau_pt(reader, "genDau_pt");
    TTreeReaderArray<vector<float> > b_genDau_eta(reader, "genDau_eta");
    TTreeReaderArray<vector<float> > b_genDau_phi(reader, "genDau_phi");
    
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
    
    // Event multiplicity vs leading jet Nchj
    TH2D* hEventMultVsLeadingJetNchj = new TH2D("hEventMultVsLeadingJetNchj", "Event Multiplicity vs Leading Jet N_{ch}^{jet}; Leading Jet N_{ch}^{jet}; Event Multiplicity", 100, 0, 100, 600, 0, 600);
    
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
    
    // Event loop using TTreeReader
    Long64_t nEntries = inTree->GetEntries();
    cout << "Processing " << nEntries << " events..." << endl;
    
    Long64_t ientry = 0;
    while (reader.Next()) {
        if (ientry % 1000 == 0) {
            cout << "Processing event " << ientry << "/" << nEntries << endl;
        }
        
        // Check if we have jets
        if (b_genJetEta.GetSize() == 0) {
            ientry++;
            continue;
        }
        
        // Check if first jet (leading jet) meets criteria: |eta|<1.6, pt>550
        if (fabs(b_genJetEta[0]) >= 1.6 || b_genJetPt[0] <= 550) {
            ientry++;
            continue;
        }
        
        // Calculate Nchj (jet charged multiplicity) for leading jet using genDau branches
        int Nchj = 0;
        if (b_genDau_chg.GetSize() > 0) {
            for (size_t i = 0; i < b_genDau_chg[0].size(); i++) {
                // Check if particle is charged
                if (b_genDau_chg[0][i] != 0) {
                    float pt = b_genDau_pt[0][i];
                    float eta = b_genDau_eta[0][i];
                    
                    // Check pt and eta cuts
                    if (pt > 0.3 && fabs(eta) < 2.4) {
                        Nchj++;
                    }
                }
            }
        }
        
        // Debug output for first few events
        if (ientry < 5) {
            cout << "Event " << ientry << ": genDau_chg size = " << b_genDau_chg.GetSize() << endl;
            if (b_genDau_chg.GetSize() > 0) {
                cout << "  Leading jet constituents: " << b_genDau_chg[0].size() << endl;
                cout << "  Leading jet genJetChargedMultiplicity: " << b_genJetChargedMultiplicity[0] << endl;
                cout << "  Calculated Nchj: " << Nchj << endl;
            }
        }
        
        // Fill Nchj vs Ncollision
        hNchjVsNcollision->Fill(*b_total_collisions, Nchj);
        
        // Fill Npar vs Nparticles
        int Npar = b_par_pdgid.GetSize();
        int Nparticles = b_px.GetSize();
        hNparVsNparticles->Fill(Nparticles, Npar);
        
        // Fill Np vs leading jet multiplicity (first jet entry is leading jet)
        hNpVsLeadingJetMult->Fill(Nchj, Npar);
        
        // Fill event multiplicity vs leading jet Nchj
        hEventMultVsLeadingJetNchj->Fill(Nchj, Nparticles);
        
        // Fill genJet distributions
        for (size_t i = 0; i < b_genJetEta.GetSize(); i++) {
            hGenJetPt->Fill(b_genJetPt[i]);
            hGenJetEta->Fill(b_genJetEta[i]);
            hGenJetPhi->Fill(b_genJetPhi[i]);
            hGenJetChargedMultiplicity->Fill(b_genJetChargedMultiplicity[i]);
        }
        
        // Fill total collisions distribution
        hTotalCollisions->Fill(*b_total_collisions);
        
        // Fill momentum correlation plots per parton based on collision count
        if (*b_total_collisions == 0) {
            // Fill per parton for events with 0 collisions
            for (size_t i = 0; i < b_par_pdgid.GetSize(); i++) {
                hPxBeforeVsAfter_0coll->Fill(b_par_px[i], b_par_px_after_zpc[i]);
                hPyBeforeVsAfter_0coll->Fill(b_par_py[i], b_par_py_after_zpc[i]);
                hPzBeforeVsAfter_0coll->Fill(b_par_pz[i], b_par_pz_after_zpc[i]);
            }
        } else {
            // Fill per parton for events with nonzero collisions
            for (size_t i = 0; i < b_par_pdgid.GetSize(); i++) {
                hPxBeforeVsAfter_ncoll->Fill(b_par_px[i], b_par_px_after_zpc[i]);
                hPyBeforeVsAfter_ncoll->Fill(b_par_py[i], b_par_py_after_zpc[i]);
                hPzBeforeVsAfter_ncoll->Fill(b_par_pz[i], b_par_pz_after_zpc[i]);
            }
        }
        
        ientry++;
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
    
    // Write event multiplicity vs leading jet Nchj
    hEventMultVsLeadingJetNchj->Write();
    
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
