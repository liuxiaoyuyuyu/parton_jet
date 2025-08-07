#include "trackTree.C"
#include "coordinateTools.h"
#include "binning.h"

const double maxTime = 10;
const int nTauBins = 40;

// Global histograms for binned analysis
TProfile* hV2VsTau[trackbin][etasbin];       // Average v2 vs tau with uncertainties
TProfile* hEccVsTau[trackbin][etasbin];      // Average eccentricity vs tau with uncertainties
TH1D* hV2RMSvsTau[trackbin][etasbin];        // RMS v2 vs tau
TH1D* hEccRMSvsTau[trackbin][etasbin];       // RMS eccentricity vs tau

// QA plots
TH2D* hNpartonVsJetMult;  // Nparton vs jet charged multiplicity
TH1D* hNpartonVsTau;      // Nparton vs tau

// Binned QA plots
TH2D* hNpartonVsJetMult_binned[trackbin][etasbin];  // Nparton vs jet charged multiplicity per bin
TH1D* hNpartonVsTau_binned[trackbin][etasbin];      // Nparton vs tau per bin
TH1D* hPsi2_binned[trackbin][etasbin];              // Psi2 distribution per bin

// Initialize histograms for each bin
void initializeHistograms() {
    // QA plots
    hNpartonVsJetMult = new TH2D("hNpartonVsJetMult", "Nparton vs Jet Charged Multiplicity; Jet Charged Multiplicity; Nparton", 50, 0, 100, 50, 0, 100);
    hNpartonVsTau = new TH1D("hNpartonVsTau", "Nparton vs Tau; #tau (fm/c); Nparton", 40, 0, maxTime);
    
    for (int i = 0; i < trackbin; i++) {
        for (int j = 0; j < etasbin; j++) {
            // TProfile for average values with uncertainties
            TString name = Form("hV2VsTau_%d_%d", i, j);
            TString title = Form("Average v2 vs tau (Nch bin %d, eta_s bin %d); #tau (fm/c); <v2>", i, j);
            hV2VsTau[i][j] = new TProfile(name, title, 40, 0, maxTime, 0, 1);
            
            name = Form("hEccVsTau_%d_%d", i, j);
            title = Form("Average Eccentricity vs tau (Nch bin %d, eta_s bin %d); #tau (fm/c); <#varepsilon_{2}>", i, j);
            hEccVsTau[i][j] = new TProfile(name, title, 40, 0, maxTime, 0, 1);
            
            // RMS histograms
            name = Form("hV2RMSvsTau_%d_%d", i, j);
            title = Form("RMS v2 vs tau (Nch bin %d, eta_s bin %d); #tau (fm/c); v2^{RMS}", i, j);
            hV2RMSvsTau[i][j] = new TH1D(name, title, 40, 0, maxTime);
            
            name = Form("hEccRMSvsTau_%d_%d", i, j);
            title = Form("RMS Eccentricity vs tau (Nch bin %d, eta_s bin %d); #tau (fm/c); #varepsilon_{2}^{RMS}", i, j);
            hEccRMSvsTau[i][j] = new TH1D(name, title, 40, 0, maxTime);
            
            // Binned QA plots
            name = Form("hNpartonVsJetMult_%d_%d", i, j);
            title = Form("Nparton vs Jet Charged Multiplicity (Nch bin %d, eta_s bin %d); Jet Charged Multiplicity; Nparton", i, j);
            hNpartonVsJetMult_binned[i][j] = new TH2D(name, title, 100, 0, 100, 100, 0, 100);
            
            name = Form("hNpartonVsTau_%d_%d", i, j);
            title = Form("Nparton vs Tau (Nch bin %d, eta_s bin %d); #tau (fm/c); Nparton", i, j);
            hNpartonVsTau_binned[i][j] = new TH1D(name, title, 40, 0, maxTime);
            
            name = Form("hPsi2_%d_%d", i, j);
            title = Form("Psi2 Distribution (Nch bin %d, eta_s bin %d); #Psi_{2}; Counts", i, j);
            hPsi2_binned[i][j] = new TH1D(name, title, 50, -TMath::Pi()/2, TMath::Pi()/2);
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

// Read and store jets
void readJets(trackTree* tree, vector<Jet>& jets) {
    jets.clear();
    for (int j = 0; j < (int)tree->genJetPt->size(); ++j) {
        Jet jet;
        jet.Pt = tree->genJetPt->at(j);
        jet.Eta = tree->genJetEta->at(j);
        jet.Phi = tree->genJetPhi->at(j);
        jet.genJetChargedMultiplicity = tree->genJetChargedMultiplicity->at(j);
        jets.push_back(jet);
    }
}

// Read and store partons
void readPartons(trackTree* tree, vector<Parton>& partons) {
    partons.clear();
    for (int j = 0; j < (int)tree->par_pdgid->size(); ++j) {
        Parton p;
        p.pdgid = (double)tree->par_pdgid->at(j);
        p.px    = (double)tree->par_px->at(j);
        p.py    = (double)tree->par_py->at(j);
        p.pz    = (double)tree->par_pz->at(j);
        p.e     = (double)tree->par_e->at(j);
        p.x     = (double)tree->par_x->at(j);
        p.y     = (double)tree->par_y->at(j);
        p.z     = (double)tree->par_z->at(j);
        p.t     = (double)tree->par_t->at(j);
        
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
}

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
            // Calculate eta_s in jet frame (similar to tau calculation)
            double z2 = P.jet_par_z * P.jet_par_z;
            double eta_s = (P.t * P.t > z2) ? 0.5 * log((P.t + P.jet_par_z) / (P.t - P.jet_par_z)) : 0.0;
            
            // Assign parton to appropriate eta_s bin
            for (int j = 0; j < etasbin; j++) {
                if (eta_s >= etasbinbounds[j] && eta_s < etasbinboundsUpper[j]) {
                    partonsByEtaBin[j].push_back(idx);
                    break;
                }
            }
        }
        
        // Fill histograms for each tau bin
        for (int it = 1; it <= nTauBins; ++it) {
            double tauTarget = it * maxTime / nTauBins;
            
            // Fill QA plot for Nparton vs tau
            hNpartonVsTau->Fill(tauTarget, npartons);
            
            // Calculate observables for each eta_s bin separately
            for (int etaBin = 0; etaBin < etasbin; etaBin++) {
                const auto& etaPartons = partonsByEtaBin[etaBin];
                if (etaPartons.size() < 2) continue; // Need at least 2 partons
                
                // Compute spatial centroid at tauTarget for this eta_s bin
                double sumE = 0, sumx = 0, sumy = 0;
                for (int idx : etaPartons) {
                    const auto& P = partons[idx];
                    if (P.tau > tauTarget) continue;
                    double dt = tauTarget - P.tau;
                    double x = P.jet_par_x + dt * P.jet_par_px / P.e;
                    double y = P.jet_par_y + dt * P.jet_par_py / P.e;
                    sumE += P.e;
                    sumx += P.e * x;
                    sumy += P.e * y;
                }
                if (sumE == 0) continue;
                
                double xm = sumx / sumE;
                double ym = sumy / sumE;
                
                // Compute eccentricity for this eta_s bin
                double Re = 0, Im = 0, Wtot = 0;
                for (int idx : etaPartons) {
                    const auto& P = partons[idx];
                    if (P.tau > tauTarget) continue;
                    double dt = tauTarget - P.tau;
                    double x = P.jet_par_x + dt * P.jet_par_px / P.e;
                    double y = P.jet_par_y + dt * P.jet_par_py / P.e;
                    double dx = x - xm;
                    double dy = y - ym;
                    double r2 = dx * dx + dy * dy;
                    double phi = atan2(dy, dx);
                    double w = P.e * r2;
                    Re += w * cos(2 * phi);
                    Im += w * sin(2 * phi);
                    Wtot += w;
                }
                if (Wtot == 0) continue;
                
                double eps2 = sqrt(Re * Re + Im * Im) / Wtot;
                double psi2 = 0.5 * atan2(Im, Re);
                
                // Compute v2 for this eta_s bin
                double sumV2 = 0;
                int n = 0;
                for (int idx : etaPartons) {
                    const auto& P = partons[idx];
                    if (P.tau > tauTarget) continue;
                    double phi_mom = atan2(P.jet_par_py, P.jet_par_px);
                    sumV2 += cos(2 * (phi_mom - psi2));
                    ++n;
                }
                
                if (n > 0) {
                    double v2 = sumV2 / n;
                    // Use TProfile which automatically handles uncertainties
                    hV2VsTau[multBin][etaBin]->Fill(tauTarget, v2);
                    hEccVsTau[multBin][etaBin]->Fill(tauTarget, eps2);
                    
                    // Fill binned QA plots
                    hNpartonVsJetMult_binned[multBin][etaBin]->Fill(jetChargedMult, etaPartons.size());
                    hNpartonVsTau_binned[multBin][etaBin]->Fill(tauTarget, etaPartons.size());
                    hPsi2_binned[multBin][etaBin]->Fill(psi2);
                }
            }
        }
    }
}

// Compute RMS correctly: RMS = sqrt(<x^2> - <x>^2)
void computeRMSvsTau(TProfile* hInput, TH1D* hRMS) {
    for (int i = 1; i <= hInput->GetNbinsX(); ++i) {
        double tau = hInput->GetXaxis()->GetBinCenter(i);
        double mean = hInput->GetBinContent(i);
        double meanSq = hInput->GetBinSumw2()->At(i) / hInput->GetBinEntries(i);
        double rms = sqrt(meanSq - mean * mean);
        hRMS->SetBinContent(i, rms);
    }
}

void parton_v2(const char* inputFileName = "pp_parton_cascade_1.root") {
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
    
    trackTree *tree = new trackTree(inTree);
    int nentries = inTree->GetEntries();
    
    // Prepare containers
    vector<Parton> partons;
    vector<Jet> jets;
    std::map<int, vector<int>> partonsByJet;
    
    // Event loop
    for (int ientry = 0; ientry < nentries; ientry++) {
        tree->GetEntry(ientry);
        
        readJets(tree, jets);
        readPartons(tree, partons);
        
        matchPartonsToJets(partons, jets);
        transformToJetFrame(partons, jets);
        
        // Calculate proper time for each parton
        for (auto& p : partons) {
            double z2 = p.jet_par_z * p.jet_par_z;
            p.tau = (p.t * p.t > z2) ? sqrt(p.t * p.t - z2) : 0.0;
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
    }
    
    // Compute RMS for each bin (TProfile already handles averages and uncertainties)
    for (int i = 0; i < trackbin; i++) {
        for (int j = 0; j < etasbin; j++) {
            computeRMSvsTau(hV2VsTau[i][j], hV2RMSvsTau[i][j]);
            computeRMSvsTau(hEccVsTau[i][j], hEccRMSvsTau[i][j]);
        }
    }
    
    // Create output file
    TString baseFileName = TString(inputFileName);
    baseFileName.ReplaceAll(".root", "");
    baseFileName.ReplaceAll("/", "_");
    baseFileName.ReplaceAll("eos_cms_store_group_phys_heavyions_huangxi_PC_", "");
    
    TString outputFileName = Form("parton_v2_output_%s.root", baseFileName.Data());
    TFile* outFile = TFile::Open(outputFileName, "RECREATE");
    
    // Write histograms
    for (int i = 0; i < trackbin; i++) {
        for (int j = 0; j < etasbin; j++) {
            hV2VsTau[i][j]->Write();
            hEccVsTau[i][j]->Write();
            hV2RMSvsTau[i][j]->Write();
            hEccRMSvsTau[i][j]->Write();
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
        }
    }
    
    outFile->Close();
    inFile->Close();
    delete tree;
}