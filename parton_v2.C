#include <vector>
#include <unordered_map>
#include <iostream>
#include <cmath>
#include <TMath.h>

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TVector2.h"
#include "coordinateTools.h"
#include "binning.h"

using namespace std;

const double maxTime = 9;
const int nTimeBins = 300;
const double tauForv2vsEcc = 3.0; // tau value for v2 vs eccentricity correlation plot

// Global histograms for binned analysis
TProfile* hv2VsTau[trackbin][etasbin];       // Average v2 vs tau with uncertainties
TProfile* hEccVsTau[trackbin][etasbin];      // Average eccentricity vs tau with uncertainties

// TH1D* hv2RMSvsTau[trackbin][etasbin];        // RMS v2 vs tau - COMMENTED OUT
// TH1D* hEccRMSvsTau[trackbin][etasbin];       // RMS eccentricity vs tau - COMMENTED OUT
TH2D* hv2vsEcc[trackbin][etasbin];           // v2 vs eccentricity correlation at tau=3fm/c

// QA plots
TH2D* hNpartonVsJetMult;  // Nparton vs jet charged multiplicity
TProfile* hNpartonVsTau;  // Nparton vs tau (changed to TProfile)

// Binned QA plots
TH2D* hNpartonVsJetMult_binned[trackbin][etasbin];  // Nparton vs jet charged multiplicity per bin
TProfile* hNpartonVsTau_binned[trackbin][etasbin];  // Nparton vs tau per bin (changed to TProfile)
TH2D* hPsi2_binned[trackbin][etasbin];              // Psi2 distribution vs tau per bin

// Initialize histograms for each bin
void initializeHistograms() {
    // QA plots
    hNpartonVsJetMult = new TH2D("hNpartonVsJetMult", "Nparton vs Jet Charged Multiplicity; Jet Charged Multiplicity; Nparton", 100, 0, 100, 100, 0, 100);
    hNpartonVsTau = new TProfile("hNpartonVsTau", "Nparton vs Tau; #tau (fm/c); <Nparton>", nTimeBins, 0, maxTime);
    
    for (int i = 0; i < trackbin; i++) {
        for (int j = 0; j < etasbin; j++) {
            // TProfile for average values with uncertainties
            TString name = Form("hv2VsTau_Nch%d_etas%d", i, j);
            TString title = Form("Average v2 vs tau (%d < N_{ch} < %d, %.1f < #eta_{s} < %.1f); #tau (fm/c); <v2>", 
                                trackbinbounds_MC[i], trackbinboundsUpper_MC[i], etasbinbounds[j], etasbinboundsUpper[j]);
            hv2VsTau[i][j] = new TProfile(name, title, nTimeBins, 0, maxTime, 0, 1);
            
            name = Form("hEccVsTau_Nch%d_etas%d", i, j);
            title = Form("Average Eccentricity vs tau (%d < N_{ch} < %d, %.1f < #eta_{s} < %.1f); #tau (fm/c); <#varepsilon_{2}>", 
                        trackbinbounds_MC[i], trackbinboundsUpper_MC[i], etasbinbounds[j], etasbinboundsUpper[j]);
            hEccVsTau[i][j] = new TProfile(name, title, nTimeBins, 0, maxTime, 0, 1);
            
            // RMS histograms - COMMENTED OUT
            // name = Form("hv2RMSvsTau_Nch%d_etas%d", i, j);
            // title = Form("RMS v2 vs tau (%d < N_{ch} < %d, %.1f < #eta_{s} < %.1f); #tau (fm/c); v2^{RMS}", 
            //             trackbinbounds_MC[i], trackbinboundsUpper_MC[i], etasbinbounds[j], etasbinboundsUpper[j]);
            // hv2RMSvsTau[i][j] = new TH1D(name, title, nTimeBins, 0, maxTime);
            
            // name = Form("hEccRMSvsTau_Nch%d_etas%d", i, j);
            // title = Form("RMS Eccentricity vs tau (%d < N_{ch} < %d, %.1f < #eta_{s} < %.1f); #tau (fm/c); #varepsilon_{2}^{RMS}", 
            //             trackbinbounds_MC[i], trackbinboundsUpper_MC[i], etasbinbounds[j], etasbinboundsUpper[j]);
            // hEccRMSvsTau[i][j] = new TH1D(name, title, nTimeBins, 0, maxTime);
            
            // Binned QA plots
            name = Form("hNpartonVsJetMult_Nch%d_etas%d", i, j);
            title = Form("Nparton vs Jet Charged Multiplicity (%d < N_{ch} < %d, %.1f < #eta_{s} < %.1f); Jet Charged Multiplicity; Nparton", 
                        trackbinbounds_MC[i], trackbinboundsUpper_MC[i], etasbinbounds[j], etasbinboundsUpper[j]);
            hNpartonVsJetMult_binned[i][j] = new TH2D(name, title, 100, 0, 100, 100, 0, 100);
            
            name = Form("hNpartonVsTau_Nch%d_etas%d", i, j);
            title = Form("Nparton vs Tau (%d < N_{ch} < %d, %.1f < #eta_{s} < %.1f); #tau (fm/c); <Nparton>", 
                        trackbinbounds_MC[i], trackbinboundsUpper_MC[i], etasbinbounds[j], etasbinboundsUpper[j]);
            hNpartonVsTau_binned[i][j] = new TProfile(name, title, nTimeBins, 0, maxTime);
            
            name = Form("hPsi2_Nch%d_etas%d", i, j);
            title = Form("Psi2 Distribution vs Tau (%d < N_{ch} < %d, %.1f < #eta_{s} < %.1f); #tau (fm/c); #Psi_{2}; Counts", 
                        trackbinbounds_MC[i], trackbinboundsUpper_MC[i], etasbinbounds[j], etasbinboundsUpper[j]);
            hPsi2_binned[i][j] = new TH2D(name, title, nTimeBins, 0, maxTime, 100, -TMath::Pi()/2, TMath::Pi()/2);
            
            // v2 vs eccentricity correlation plot at specific tau
            name = Form("hv2vsEcc_Nch%d_etas%d", i, j);
            title = Form("v2 vs Eccentricity at #tau<%.1f fm/c (%d < N_{ch} < %d, %.1f < #eta_{s} < %.1f); #varepsilon_{2}; v2", 
                        tauForv2vsEcc, trackbinbounds_MC[i], trackbinboundsUpper_MC[i], etasbinbounds[j], etasbinboundsUpper[j]);
            hv2vsEcc[i][j] = new TH2D(name, title, 50, 0, 1, 50, -1, 1);
        }
    }
}

struct Parton {
    int   pdgid;
    double px, py, pz, e;
    double x,  y,  z,  t;
    double pt, eta, phi;
    double eta_s, tau; // defined in jet frame
    int   jetID;
    //defined in jet frame
    double jet_par_x;
    double jet_par_y;
    double jet_par_z;
    double jet_par_px;
    double jet_par_py;
    double jet_par_pz;
};

struct Jet {
    double Pt;
    double Eta;
    double Phi;
    int   genJetChargedMultiplicity;
};



// Match Partons to Jets, match parton to the first jet that is within 0.8 of the parton's eta and phi.
void matchPartonsToJets(vector<Parton>& partons, const vector<Jet>& jets) {
    for (auto& p : partons) {
        for (int j = 0; j < (int)jets.size(); ++j) {
            double dphi = TVector2::Phi_mpi_pi(p.phi - jets[j].Phi);
            double dEta = p.eta - jets[j].Eta;
            double dR   = sqrt(dEta * dEta + dphi * dphi);
            if (dR < 0.8) {
                p.jetID = j;
                break;
            }
        }
    }
}

// Coordinate Transformation from Lab to Jet Frame
void transformToJetFrame(vector<Parton>& partons, const vector<Jet>& jets) {
    for (auto& p : partons) {
        if (p.jetID < 0) continue;
        
        TVector3 jv;
        jv.SetPtEtaPhi(jets[p.jetID].Pt, jets[p.jetID].Eta, jets[p.jetID].Phi);
        
        // rotate momentum into the jet frame
        TVector3 mom(p.px, p.py, p.pz);
        //the naming is a bit confusing, ptJ_m is the momentum of the parton in the jet frame.
        double ptJ_m = ptWRTJet(jv, mom);
        double thJ_m = thetaWRTJet(jv, mom);
        double phJ_m = phiWRTJet(jv, mom);
        
        p.jet_par_px = ptJ_m * cos(phJ_m);
        p.jet_par_py = ptJ_m * sin(phJ_m);
        p.jet_par_pz = ptJ_m / tan(thJ_m);
        
        // Use small tolerance for double comparison
        const double tolerance = 1e-10;
        if ((fabs(p.x) < tolerance) && (fabs(p.y) < tolerance) && (fabs(p.z) < tolerance)) {
            p.jet_par_x = 0;
            p.jet_par_y = 0;
            p.jet_par_z = 0;
            continue;
        }
        
        TVector3 pos(p.x, p.y, p.z);
        double ptJ = ptWRTJet(jv, pos);
        double thJ = thetaWRTJet(jv, pos);
        double phJ = phiWRTJet(jv, pos);
        
        p.jet_par_x = ptJ * cos(phJ);
        p.jet_par_y = ptJ * sin(phJ);
        p.jet_par_z = ptJ / tan(thJ);
    }
}

// Group by JetID and compute observables for each bin
void fillBinnedObservables(std::map<int, std::vector<int>>& partonsByJet, std::vector<Parton>& partons, const vector<Jet>& jets) {
    for (const auto& [jetID, idxs] : partonsByJet) {
        if (idxs.size() < 2) continue;
        
        // Use jet charged multiplicity (not number of partons)
        int jetChargedMult = jets[jetID].genJetChargedMultiplicity;
        int npartons = idxs.size();
        
        // Fill QA plots
        hNpartonVsJetMult->Fill(jetChargedMult, npartons);
        
        // Determine multiplicity bin based on jet charged multiplicity
        int multBin = -1;
        for (int i = 0; i < trackbin; i++) {
            if (jetChargedMult >= trackbinbounds_MC[i] && jetChargedMult < trackbinboundsUpper_MC[i]) {
                multBin = i;
                break;
            }
        }
        if (multBin == -1) continue; // Skip if outside multiplicity range
        
        // Calculate eta_s for each parton and group them by eta_s bins
        std::vector<std::vector<int>> partonsByEtaBin(etasbin);
        
        for (int idx : idxs) {
            const auto& P = partons[idx];
            // Use the pre-calculated eta_s value
            double eta_s = P.eta_s;
            
            // Assign parton to appropriate eta_s bin
            for (int j = 0; j < etasbin; j++) {
                if (eta_s >= etasbinbounds[j] && eta_s < etasbinboundsUpper[j]) {
                    partonsByEtaBin[j].push_back(idx);
                    break;
                }
            }
        }
        
        // Fill binned QA plots that don't depend on tau (once per jet)
        for (int etaBin = 0; etaBin < etasbin; etaBin++) {
            const auto& etaPartons = partonsByEtaBin[etaBin];
            if (etaPartons.size() >= 2) {
                hNpartonVsJetMult_binned[multBin][etaBin]->Fill(jetChargedMult, etaPartons.size());
            }
        }
        
        // Fill histograms for each tau bin
        for (int it = 1; it <= nTimeBins; ++it) {
            double tauTarget = it * maxTime / nTimeBins; // Upper edge of the bin
            
            // Count partons whose formation time is smaller than tauTarget
            int nPartonsAtTau = 0;
            for (int idx : idxs) {
                const auto& P = partons[idx];
                if (P.tau < tauTarget) nPartonsAtTau++;
            }
            // For TProfile, we want to average the number of partons per jet at each tau
            // So we fill with the count for this specific jet
            hNpartonVsTau->Fill(tauTarget, nPartonsAtTau);
            

            

            

            
            // Calculate observables for each eta_s bin separately
            for (int etaBin = 0; etaBin < etasbin; etaBin++) {
                const auto& etaPartons = partonsByEtaBin[etaBin];
                if (etaPartons.size() < 2) continue; // Need at least 2 partons
                

                
                // Count partons in this eta_s bin whose formation time is smaller than tauTarget
                int nPartonsInEtaBin = 0;
                for (int idx : etaPartons) {
                    const auto& P = partons[idx];
                    if (P.tau < tauTarget) nPartonsInEtaBin++;
                }
                

                
                if (nPartonsInEtaBin < 2) continue; // Need at least 2 partons
                
                // Compute spatial centroid and eccentricity in one loop
                double sumE = 0, sumx = 0, sumy = 0;
                double Re = 0, Im = 0, Wtot = 0;
                
                for (int idx : etaPartons) {
                    const auto& P = partons[idx];
                    if (P.tau >= tauTarget) continue;
                    
                    // Use lab time for free streaming (from evolution_jet_bins.cc)
                    double t_upperEdge = tauTarget * cosh(P.eta_s);
                    double x_jetzt = P.jet_par_x + (t_upperEdge - P.t) * P.jet_par_px / P.e;
                    double y_jetzt = P.jet_par_y + (t_upperEdge - P.t) * P.jet_par_py / P.e;
                    
                    sumE += P.e;
                    sumx += P.e * x_jetzt;
                    sumy += P.e * y_jetzt;
                }
                
                if (sumE < 1e-10) continue; // Use small tolerance instead of exact zero comparison
                
                double xm = sumx / sumE;
                double ym = sumy / sumE;
                
                // Compute eccentricity in the same loop as centroid calculation
                for (int idx : etaPartons) {
                    const auto& P = partons[idx];
                    if (P.tau >= tauTarget) continue;
                    
                    // Use lab time for free streaming
                    double t_upperEdge = tauTarget * cosh(P.eta_s);
                    double x_jetzt = P.jet_par_x + (t_upperEdge - P.t) * P.jet_par_px / P.e;
                    double y_jetzt = P.jet_par_y + (t_upperEdge - P.t) * P.jet_par_py / P.e;
                    
                    double dx = x_jetzt - xm;
                    double dy = y_jetzt - ym;
                    double r2 = dx * dx + dy * dy;
                    double phi = atan2(dy, dx);
                    double w = P.e * r2;
                    Re += w * cos(2 * phi);
                    Im += w * sin(2 * phi);
                    Wtot += w;
                }
                
                if (Wtot < 1e-10) continue; // Use small tolerance instead of exact zero comparison
                
                double ecc = sqrt(Re * Re + Im * Im) / Wtot;
                double psi2 = 0.5 * atan2(Im, Re);
                hEccVsTau[multBin][etaBin]->Fill(tauTarget, ecc); 
                
                // Calculate v2 by averaging over partons in this jet first, then average over jets
                double sumv2 = 0;
                int nPartonsForv2 = 0;
                for (int idx : etaPartons) {
                    const auto& P = partons[idx];
                    if (P.tau >= tauTarget) continue;
                    
                    double phi_mom = atan2(P.jet_par_py, P.jet_par_px);
                    double v2 = cos(2 * (phi_mom - psi2));
                    sumv2 += v2;
                    nPartonsForv2++;
                }
                
                if (nPartonsForv2 > 0) { // This comparison is fine since it's an integer
                    double avgv2 = sumv2 / nPartonsForv2; // Average over partons in this jet
                    hv2VsTau[multBin][etaBin]->Fill(tauTarget, avgv2); // Average over jets
                    
                    // Fill v2 vs eccentricity correlation at specific tau
                    if (fabs(tauTarget - tauForv2vsEcc) < maxTime/nTimeBins/2) { // Check if we're at the target tau bin
                        hv2vsEcc[multBin][etaBin]->Fill(ecc, avgv2);
                    }
                }
                
                // Fill binned QA plots
                hNpartonVsTau_binned[multBin][etaBin]->Fill(tauTarget, nPartonsInEtaBin); // Number of partons with tau < tauTarget
                hPsi2_binned[multBin][etaBin]->Fill(tauTarget, psi2);
            }//loop over eta_s
        }//loop over tau
    }//loop over jet
}

// Compute RMS correctly: RMS = sqrt(<x^2> - <x>^2)
void computeRMSvsTau(TProfile* hInput, TH1D* hRMS) {
    for (int i = 1; i <= hInput->GetNbinsX(); ++i) {
        double tau = hInput->GetXaxis()->GetBinCenter(i);
        int entries = hInput->GetBinEntries(i);
        
        if (entries > 0) {
            // For TProfile, we can get RMS directly from the error
            // TProfile error = RMS/sqrt(N), so RMS = error * sqrt(N)
            double error = hInput->GetBinError(i);
            double rms = error * sqrt(entries);
            hRMS->SetBinContent(i, rms);
        } else {
            hRMS->SetBinContent(i, 0);
        }
    }
}

void parton_v2(const char* inputFileName = "/eos/cms/store/group/phys_heavyions/huangxi/PC/pp_parton_cascade_0.root", 
               const char* outputDir = "/eos/cms/store/group/phys_heavyions/xiaoyul/wenbin/anaOutput/") {
    // Initialize histograms
    initializeHistograms();
    
    // Open input file
    TFile* inFile = TFile::Open(inputFileName, "READ");
    if (!inFile || inFile->IsZombie()) {
        cout << "Error: Cannot open file " << inputFileName << endl;
        return;
    }
    
    TTree *inTree = (TTree*)inFile->Get("trackTree");
    if (!inTree) {
        cout << "Error: Cannot find TTree 'trackTree' in file " << inputFileName << endl;
        inFile->Close();
        return;
    }
    

    
    int nentries = inTree->GetEntries();
    
    //----------------------------------------------------------------------
    // Set up branch addresses (like evolution_jet_bins.cc)
    //----------------------------------------------------------------------
    
    // Parton branches
    vector<int>*    b_par_pdgid = nullptr;
    vector<double>* b_par_px    = nullptr;
    vector<double>* b_par_py    = nullptr;
    vector<double>* b_par_pz    = nullptr;
    vector<double>* b_par_e     = nullptr;
    vector<double>* b_par_x     = nullptr;
    vector<double>* b_par_y     = nullptr;
    vector<double>* b_par_z     = nullptr;
    vector<double>* b_par_t     = nullptr;

    inTree->SetBranchAddress("par_pdgid", &b_par_pdgid);
    inTree->SetBranchAddress("par_px",    &b_par_px);
    inTree->SetBranchAddress("par_py",    &b_par_py);
    inTree->SetBranchAddress("par_pz",    &b_par_pz);
    inTree->SetBranchAddress("par_e",     &b_par_e);
    inTree->SetBranchAddress("par_x",     &b_par_x);
    inTree->SetBranchAddress("par_y",     &b_par_y);
    inTree->SetBranchAddress("par_z",     &b_par_z);
    inTree->SetBranchAddress("par_t",     &b_par_t);

    // Jet branches
    vector<double>* b_genJetPt  = nullptr;
    vector<double>* b_genJetEta = nullptr;
    vector<double>* b_genJetPhi = nullptr;
    vector<int>* b_genJetChargedMult = nullptr;

    inTree->SetBranchAddress("genJetPt",  &b_genJetPt);
    inTree->SetBranchAddress("genJetEta", &b_genJetEta);
    inTree->SetBranchAddress("genJetPhi", &b_genJetPhi);
    inTree->SetBranchAddress("genJetChargedMultiplicity", &b_genJetChargedMult);
    
    // Prepare containers
    vector<Parton> partons;
    vector<Jet> jets;
    std::map<int, vector<int>> partonsByJet;
    
    // Event loop
    for (int ientry = 0; ientry < nentries; ientry++) {
    //for (int ientry = 0; ientry < 10; ientry++) {
        // Get the current entry (like evolution_jet_bins.cc)
        inTree->GetEntry(ientry);
        
        // Read jets (like evolution_jet_bins.cc)
        jets.clear();
        for (int j = 0; j < (int)b_genJetPt->size(); ++j) {
            Jet jet;
            jet.Pt = b_genJetPt->at(j);
            jet.Eta = b_genJetEta->at(j);
            jet.Phi = b_genJetPhi->at(j);
            jet.genJetChargedMultiplicity = b_genJetChargedMult->at(j);
            jets.push_back(jet);
        }
        
        // Read partons (like evolution_jet_bins.cc)
        partons.clear();
        for (int j = 0; j < (int)b_par_pdgid->size(); ++j) {
            Parton p;
            p.pdgid = b_par_pdgid->at(j);
            p.px    = b_par_px->at(j);
            p.py    = b_par_py->at(j);
            p.pz    = b_par_pz->at(j);
            p.e     = b_par_e->at(j);
            p.x     = b_par_x->at(j);
            p.y     = b_par_y->at(j);
            p.z     = b_par_z->at(j);
            p.t     = b_par_t->at(j);
            
            // Calculate parton kinematics
            p.pt   = sqrt(p.px * p.px + p.py * p.py);
            p.phi  = atan2(p.py, p.px);
            double theta = atan2(p.pt, p.pz);
            p.eta  = -log(tan(theta / 2.0));
            
            p.jetID       = -1;
            p.jet_par_x   = 0.0;
            p.jet_par_y   = 0.0;
            p.jet_par_z   = 0.0;
            p.jet_par_px  = 0.0;
            p.jet_par_py  = 0.0;
            p.jet_par_pz  = 0.0;
            
            partons.push_back(p);
        }
        
        matchPartonsToJets(partons, jets);
        transformToJetFrame(partons, jets);
        
        // Calculate proper time and eta_s for each parton
        for (auto& p : partons) {
            double z2 = p.jet_par_z * p.jet_par_z;
            p.tau = (p.t * p.t > z2) ? sqrt(p.t * p.t - z2) : 0.0;
            
            // Spacetime rapidity (same as evolution_jet_bins.cc)
            if (p.t > p.jet_par_z) {
                p.eta_s = 0.5 * log((p.t + p.jet_par_z) / (p.t - p.jet_par_z));
            } else {
                double vz = p.jet_par_pz / p.e;
                p.eta_s = 0.5 * log((1 + vz) / (1 - vz));
            }
        }
        

        // Group partons by jet
        partonsByJet.clear();
        for (int ip = 0; ip < partons.size(); ++ip) {
            auto& p = partons[ip];
            if (p.jetID >= 0) {
                partonsByJet[p.jetID].push_back(ip);
            }
        }
        // Fill binned observables
        fillBinnedObservables(partonsByJet, partons, jets);
        
        // Debug: Print some info for first few events
        if (ientry < 3) {
            cout << "Event " << ientry << ": " << jets.size() << " jets, " << partons.size() << " partons" << endl;
            for (int i = 0; i < min(3, (int)partons.size()); i++) {
                cout << "  Parton " << i << ": tau=" << partons[i].tau << ", eta_s=" << partons[i].eta_s << endl;
            }
        }
    }
    
    // Compute RMS for each bin - COMMENTED OUT
    // for (int i = 0; i < trackbin; i++) {
    //     for (int j = 0; j < etasbin; j++) {
    //         computeRMSvsTau(hv2VsTau[i][j], hv2RMSvsTau[i][j]);
    //         computeRMSvsTau(hEccVsTau[i][j], hEccRMSvsTau[i][j]);
    //     }
    // }
    
    // Create output file based on input filename
    TString baseFileName = TString(inputFileName);
    baseFileName.ReplaceAll(".root", "");
    baseFileName.ReplaceAll("/", "_");
    baseFileName.ReplaceAll("eos_cms_store_group_phys_heavyions_huangxi_PC_", "");
    
    TString outputFileName = Form("%sparton_v2_output_%s.root", outputDir, baseFileName.Data());
    TFile* outFile = TFile::Open(outputFileName, "RECREATE");
    
    // Write histograms
    for (int i = 0; i < trackbin; i++) {
        for (int j = 0; j < etasbin; j++) {
            hv2VsTau[i][j]->Write();
            hEccVsTau[i][j]->Write();
            // hv2RMSvsTau[i][j]->Write();  // COMMENTED OUT
            // hEccRMSvsTau[i][j]->Write(); // COMMENTED OUT
        }
    }
    
    // Write QA plots
    hNpartonVsJetMult->Write();
    hNpartonVsTau->Write();
    
    // Write binned QA plots
    for (int i = 0; i < trackbin; i++) {
        for (int j = 0; j < etasbin; j++) {
            hNpartonVsJetMult_binned[i][j]->Write();
            hNpartonVsTau_binned[i][j]->Write();
            hPsi2_binned[i][j]->Write();
            hv2vsEcc[i][j]->Write();
        }
    }
    
    outFile->Close();
    inFile->Close();
}