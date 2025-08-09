#include <vector>
#include <unordered_map>
#include <iostream>
#include <cmath>
#include <TMath.h>

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TVector2.h"
#include "coordinateTools.h"
#include "binning.h"

using namespace std;

//----------------------------------------------------------------------
// 数据结构定义
//----------------------------------------------------------------------

struct Parton {
    int   pdgid;

    double px, py, pz, e;
    double x,  y,  z,  t;
    double pt, eta, phi;
    double eta_s, tau; // defined in jet frame

    int   jetID;

    double jet_par_px, jet_par_py, jet_par_pz;
    double jet_par_x, jet_par_y, jet_par_z;
};

struct Jet {
    double Pt;
    double Eta;
    double Phi;
    int   genJetChargedMultiplicity;
    int multBin;
};

//----------------------------------------------------------------------
// 主程序
//----------------------------------------------------------------------

int evolution_jet_bins(int idx) 
{
    // 打开 ROOT 文件
    TString inputName = Form("/eos/cms/store/group/phys_heavyions/huangxi/PC/pp_parton_cascade_%d.root", idx * 200);
    TFile *file = TFile::Open(inputName);
    if (!file || file->IsZombie())
    {
        cerr << "Error: cannot open file pp_parton_cascade.root" << endl;
        return 1;
    }

    // 获取 TTree
    TTree* tree = dynamic_cast<TTree*>(file->Get("trackTree"));
    if (!tree) {
        cerr << "Error: cannot find TTree 'trackTree' in file" << endl;
        file->Close();
        return 1;
    }

    //----------------------------------------------------------------------
    // 设置分支地址
    //----------------------------------------------------------------------
    
    vector<int>*    b_par_pdgid = nullptr;
    vector<double>* b_par_px    = nullptr;
    vector<double>* b_par_py    = nullptr;
    vector<double>* b_par_pz    = nullptr;
    vector<double>* b_par_e     = nullptr;
    vector<double>* b_par_x     = nullptr;
    vector<double>* b_par_y     = nullptr;
    vector<double>* b_par_z     = nullptr;
    vector<double>* b_par_t     = nullptr;

    tree->SetBranchAddress("par_pdgid", &b_par_pdgid);
    tree->SetBranchAddress("par_px",    &b_par_px);
    tree->SetBranchAddress("par_py",    &b_par_py);
    tree->SetBranchAddress("par_pz",    &b_par_pz);
    tree->SetBranchAddress("par_e",     &b_par_e);
    tree->SetBranchAddress("par_x",     &b_par_x);
    tree->SetBranchAddress("par_y",     &b_par_y);
    tree->SetBranchAddress("par_z",     &b_par_z);
    tree->SetBranchAddress("par_t",     &b_par_t);

    vector<double>* b_genJetPt  = nullptr;
    vector<double>* b_genJetEta = nullptr;
    vector<double> *b_genJetPhi = nullptr;
    vector<int> *b_genJetChargedMult = nullptr;

    tree->SetBranchAddress("genJetPt",  &b_genJetPt);
    tree->SetBranchAddress("genJetEta", &b_genJetEta);
    tree->SetBranchAddress("genJetPhi", &b_genJetPhi);
    tree->SetBranchAddress("genJetChargedMultiplicity", &b_genJetChargedMult);

    //----------------------------------------------------------------------
    // 创建直方图
    //----------------------------------------------------------------------

    int nTimeBins = 100;
    double maxTime = 3;
    double timeBinWidth = maxTime / nTimeBins;

    TH1D *hTime = new TH1D("hTime", "Generation proper time distribution; #tau (fm/c); Entries", 100, 0, 100);
    TH1D *hMult = new TH1D("hMult", "Jet multiplicity distribution; N_{ch}^{j}; Entries", 12, 0, 120);
    TH2D *hJetMultParton = new TH2D("hJetMultParton", ";N_{ch}^{j};N_{p}^{j}", 50, 0, 100, 25, 0, 50);

    TH2D *hMultTime[trackbin];
    TH2D *hEccTime[trackbin];
    TH2D *hDimTime[trackbin];
    TH3D *hDenTime[trackbin];
    TH3D *hNpTime[trackbin];
    TH3D *hAreaTime[trackbin];
    TH2D *hATest[trackbin];

    for (int ibin = 0; ibin < trackbin; ibin++)
    {
        hMultTime[ibin] = new TH2D(Form("hMultTime_%d", ibin), Form("Jet parton multiplicity vs time, %d < N_{ch}^{j} < %d; t (fm/c); N_{parton}", trackbinbounds_MC[ibin], trackbinboundsUpper_MC[ibin]), nTimeBins, 0, maxTime, 100, 0, 100);
        hEccTime[ibin] = new TH2D(Form("hEccTime_%d", ibin), Form("Eccentricity vs time, %d < N_{ch}^{j} < %d; t (fm/c); #epsilon_{2}", trackbinbounds_MC[ibin], trackbinboundsUpper_MC[ibin]), nTimeBins, 0, maxTime, 50, 0, 1);
        hDimTime[ibin] = new TH2D(Form("hDimTime_%d", ibin), Form("Dimension vs time, %d < N_{ch}^{j} < %d; t (fm/c); r (fm)", trackbinbounds_MC[ibin], trackbinboundsUpper_MC[ibin]), nTimeBins, 0, maxTime, 50, 0, 2);
        hDenTime[ibin] = new TH3D(Form("hDenTime_%d", ibin), Form("Density vs time vs #eta_{s}, %d < N_{ch}^{j} < %d; t (fm/c); #eta_{s}; density (fm^{-2})", trackbinbounds_MC[ibin], trackbinboundsUpper_MC[ibin]), nTimeBins, 0, maxTime, etasbin, etasbinbounds[0], etasbinboundsUpper[etasbin - 1], 400, 0, 100);
        hNpTime[ibin] = new TH3D(Form("hNpTime_%d", ibin), Form("Parton number vs time vs #eta_{s}, %d < N_{ch}^{j} < %d; t (fm/c); #eta_{s}; # parton", trackbinbounds_MC[ibin], trackbinboundsUpper_MC[ibin]), nTimeBins, 0, maxTime, etasbin, etasbinbounds[0], etasbinboundsUpper[etasbin - 1], 50, 0, 50);
        hAreaTime[ibin] = new TH3D(Form("hAreaTime_%d", ibin), Form("Area vs time vs #eta_{s}, %d < N_{ch}^{j} < %d; t (fm/c); #eta_{s}; S (fm^{2})", trackbinbounds_MC[ibin], trackbinboundsUpper_MC[ibin]), nTimeBins, 0, maxTime, etasbin, etasbinbounds[0], etasbinboundsUpper[etasbin - 1], 100, 0, 20);
        hATest[ibin] = new TH2D(Form("hATest_%d", ibin), Form("Area vs time, %d < N_{ch}^{j} < %d, 1 < #eta_{s} < 2; t (fm/c); S (fm^{2})", trackbinbounds_MC[ibin], trackbinboundsUpper_MC[ibin]), nTimeBins, 0, maxTime, 100, 0, 20);
    }

    // 风格设置
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kBird);

    //----------------------------------------------------------------------
    // 准备容器
    //----------------------------------------------------------------------

    vector<Parton>           partons;
    vector<Jet>              jets;

    //----------------------------------------------------------------------
    // 事件循环
    //----------------------------------------------------------------------

    Long64_t nentries = tree->GetEntries();

    for (Long64_t ie = 0; ie < nentries; ++ie) {
        tree->GetEntry(ie);

        // 读取并存储 jets
        jets.clear();
        for (int j = 0; j < (int)b_genJetPt->size(); ++j) {
            Jet jet;
            jet.Pt  = b_genJetPt->at(j);
            jet.Eta = b_genJetEta->at(j);
            jet.Phi = b_genJetPhi->at(j);
            jet.genJetChargedMultiplicity = b_genJetChargedMult->at(j);
            jet.multBin = -1;
            for (int ibin = 0; ibin < trackbin; ibin++)
            {
                if (jet.genJetChargedMultiplicity >= trackbinbounds_MC[ibin] && jet.genJetChargedMultiplicity < trackbinboundsUpper_MC[ibin]) 
                {
                    jet.multBin = ibin;
                    break;
                }
            }
            hMult->Fill(jet.genJetChargedMultiplicity);

            jets.push_back(jet);
        }

        // 读取并存储 partons
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

            // 计算动量相关量
            p.pt   = sqrt(p.px*p.px + p.py*p.py);
            p.phi  = atan2(p.py, p.px);
            double theta = atan2(p.pt, p.pz);
            p.eta  = -log(tan(theta/2.0f));

            p.jetID       = -1;

            partons.push_back(p);
        }

        //----------------------------------------------------------------------
        // 关联 parton 与 jet
        //----------------------------------------------------------------------

        for (auto& p : partons) {
            for (int j = 0; j < (int)jets.size(); ++j) {
                double dphi = TVector2::Phi_mpi_pi(p.phi - jets[j].Phi);
                double dEta = p.eta - jets[j].Eta;
                double dR   = sqrt(dEta*dEta + dphi*dphi);

                if (dR < 0.8f) {
                    p.jetID = j;
                    break;
                }
            }
        }
        
        // std::cout << "Event " << ie
        //           << ":  partons=" << partons.size()
        //           << ", jets="    << jets.size() <<"\n";
            //      << ", first jet mult="    << jets[0].genJetChargedMultiplicity
            //      <<", second jet mult="    << jets[1].genJetChargedMultiplicity<<"\n";

        //----------------------------------------------------------------------
        // 坐标转换：实验室系 -> jet 参考系
        //----------------------------------------------------------------------

        for (auto& p : partons) 
        {
            // Filter partons not belonging to any jet
            if ((p.px == 0) && (p.py == 0) && (p.pz == 0)) p.jetID = -1;
            if (p.jetID < 0) continue;

            // Rotate parton momenta to jet frame
            TVector3 jv;
            jv.SetPtEtaPhi(jets[p.jetID].Pt, jets[p.jetID].Eta, jets[p.jetID].Phi);

            TVector3 mom(p.px, p.py, p.pz);

            double ptJ_mom = ptWRTJet  (jv, mom);
            double thJ_mom = thetaWRTJet(jv, mom);
            double phJ_mom = phiWRTJet  (jv, mom);

            p.jet_par_px = ptJ_mom * cos(phJ_mom);
            p.jet_par_py = ptJ_mom * sin(phJ_mom);
            p.jet_par_pz = ptJ_mom / tan(thJ_mom);

            // Jet frame coordinates and proper time
            if ((p.x == 0) && (p.y == 0) && (p.z == 0))
            {
                p.jet_par_x = 0;
                p.jet_par_y = 0;
                p.jet_par_z = 0;
                p.tau = 0;
            }
            else 
            {
                TVector3 pos(p.x, p.y, p.z);

                double ptJ = ptWRTJet(jv, pos);
                double thJ = thetaWRTJet(jv, pos);
                double phJ = phiWRTJet(jv, pos);

                p.jet_par_x = ptJ * cos(phJ);
                p.jet_par_y = ptJ * sin(phJ);
                p.jet_par_z = ptJ / tan(thJ);
                p.tau = sqrt(p.t * p.t - p.jet_par_z * p.jet_par_z);
            }

            // Spacetime rapidity
            if (p.t > p.jet_par_z) 
            {
                p.eta_s = 0.5 * log((p.t + p.jet_par_z) / (p.t - p.jet_par_z));
            }
            else 
            {
                double vz = p.jet_par_pz / p.e;
                p.eta_s = 0.5 * log((1 + vz) / (1 - vz));
            }
        }

        //----------------------------------------------------------------------
        // 按 jetID 分组，并计算 eccentricity
        //----------------------------------------------------------------------

        unordered_map<int, vector<int>> partonsByJet;

        for (int ip = 0; ip < (int)partons.size(); ++ip) 
        {
            auto& p = partons[ip];

            if (p.jetID < 0) continue;

            partonsByJet[p.jetID].push_back(ip);
            hTime->Fill(p.tau);
        }

        for (auto &kv : partonsByJet) 
        {
            int jID = kv.first;
            hJetMultParton->Fill(jets[jID].genJetChargedMultiplicity, kv.second.size());
        }

        for (int itau = 1; itau <= nTimeBins; itau++)
        {
            double tau_lowEdge = (itau - 1) * timeBinWidth;
            double tau_center = (itau - 0.5) * timeBinWidth;

            for (auto &kv : partonsByJet)
            {
                int jID = kv.first;
                int jBin = jets[jID].multBin;

                auto &idxs = kv.second;
                if (idxs.empty()) continue;

                int nPartons = 0;
                double sumE = 0.0;
                double sumx = 0.0;
                double sumy = 0.0;
                double sumE_test = 0, sumx_test = 0, sumy_test = 0;

                for (int idx : idxs)
                {
                    auto &P = partons[idx];
                    if (P.tau < tau_lowEdge) 
                    {
                        nPartons++;

                        double t_lowEdge = tau_lowEdge * cosh(P.eta_s);
                        double x_jetzt = P.jet_par_x + (t_lowEdge - P.t) * P.jet_par_px / P.e;
                        double y_jetzt = P.jet_par_y + (t_lowEdge - P.t) * P.jet_par_py / P.e;

                        double w = P.e;
                        sumE += w;
                        sumx += w * x_jetzt;
                        sumy += w * y_jetzt;

                        if (P.eta_s >= 1 && P.eta_s < 2) 
                        {
                            sumE_test += w;
                            sumx_test += w * x_jetzt;
                            sumy_test += w * y_jetzt;
                        }
                    }
                }
                hMultTime[jBin]->Fill(tau_center, nPartons);

                // eccentricity calculation, only proceed if it is well-defined
                if (nPartons > 1) 
                {
                    double xm = sumx / sumE;
                    double ym = sumy / sumE;

                    double xm_test = sumx_test / sumE_test;
                    double ym_test = sumy_test / sumE_test;

                    double real = 0.0;
                    double imag = 0.0;
                    double norm = 0.0;
                    double x2sum_test = 0, y2sum_test = 0;

                    for (int idx : idxs)
                    {
                        auto &P = partons[idx];
                        if (P.tau < tau_lowEdge)
                        {
                            double t_lowEdge = tau_lowEdge * cosh(P.eta_s);
                            double x_jetzt = P.jet_par_x + (t_lowEdge - P.t) * P.jet_par_px / P.e;
                            double y_jetzt = P.jet_par_y + (t_lowEdge - P.t) * P.jet_par_py / P.e;

                            double dx = x_jetzt - xm;
                            double dy = y_jetzt - ym;
                            double r2 = dx * dx + dy * dy;
                            double phi_loc = atan2(dy, dx);
                            double w = P.e * r2;

                            real += w * cos(2 * phi_loc);
                            imag += w * sin(2 * phi_loc);
                            norm += w;

                            if (P.eta_s >= 1 && P.eta_s < 2)
                            {
                                double dx_test = x_jetzt - xm_test;
                                double dy_test = y_jetzt - ym_test;
                                x2sum_test += P.e * pow(dx_test, 2);
                                y2sum_test += P.e * pow(dy_test, 2);
                            }
                        }
                    }
                    double S_test = TMath::Pi() * sqrt(x2sum_test * y2sum_test) / sumE_test;
                    //hATest[jBin]->Fill(t, S_test);

                    double ecc = (norm > 0.0) ? sqrt(real * real + imag * imag) / norm : 0.0;
                    double dimen = sqrt(norm / sumE);

                    hEccTime[jBin]->Fill(tau_center, ecc);
                    hDimTime[jBin]->Fill(tau_center, dimen);
                    if (ecc > 1.0000001)
                        cout << "WARNING: ecc = " << ecc << endl;
                }

                for (int ieta = 0; ieta < etasbin; ieta++)
                {
                    vector<int> partonsIDbyEta;
                    for (int idx : idxs)
                    {
                        auto &P = partons[idx];
                        if (P.tau < tau_lowEdge && P.eta_s >= etasbinbounds[ieta] && P.eta_s < etasbinboundsUpper[ieta])
                        {
                            partonsIDbyEta.push_back(idx);
                        }
                    }
                    int nsubpartons = partonsIDbyEta.size();
                    hNpTime[jBin]->Fill(tau_center, (etasbinbounds[ieta] + etasbinboundsUpper[ieta]) / 2, nsubpartons);

                    if (nsubpartons < 2) continue;

                    double subsume = 0, subsumx = 0, subsumy = 0;

                    for (int idx : partonsIDbyEta)
                    {
                        auto &P = partons[idx];

                        double t_lowEdge = tau_lowEdge * cosh(P.eta_s);
                        double x_jetzt = P.jet_par_x + (t_lowEdge - P.t) * P.jet_par_px / P.e;
                        double y_jetzt = P.jet_par_y + (t_lowEdge - P.t) * P.jet_par_py / P.e;

                        double w = P.e;
                        subsume += w;
                        subsumx += w * x_jetzt;
                        subsumy += w * y_jetzt;
                    }

                    double subxm = subsumx / subsume;
                    double subym = subsumy / subsume;

                    double x2sum = 0, y2sum = 0;

                    for (int idx : partonsIDbyEta) 
                    {
                        auto &P = partons[idx];

                        double t_lowEdge = tau_lowEdge * cosh(P.eta_s);
                        double x_jetzt = P.jet_par_x + (t_lowEdge - P.t) * P.jet_par_px / P.e;
                        double y_jetzt = P.jet_par_y + (t_lowEdge - P.t) * P.jet_par_py / P.e;

                        double dx = x_jetzt - subxm;
                        double dy = y_jetzt - subym;

                        x2sum += P.e * pow(dx, 2);
                        y2sum += P.e * pow(dy, 2);
                    }

                    double area = TMath::Pi() * sqrt(x2sum * y2sum) / subsume;
                    if (ieta == 1) hATest[jBin]->Fill(tau_center, area);
                    //if (area < 1e-6) continue;
                    double dens = nsubpartons / (etasbinboundsUpper[ieta] - etasbinbounds[ieta]) / area;
                    /*
                    if (0&&area < 1e-3) 
                    {
                        cout << "area = " << area << "\tdensity = " << dens << endl;
                        for (int idx : partonsIDbyEta)
                        {
                            auto &P = partons[idx];

                            double x_jetzt = P.jet_par_x + (t - P.t) * P.jet_par_px / P.e;
                            double y_jetzt = P.jet_par_y + (t - P.t) * P.jet_par_py / P.e;

                            double dx = x_jetzt - subxm;
                            double dy = y_jetzt - subym;
                            std::cout << "orig_pos=(" << P.jet_par_x << "," << P.jet_par_y << ") curr_pos=(" << x_jetzt << "," << y_jetzt << ") mom=(" << P.jet_par_px << "," << P.jet_par_py << ") e=" << P.e << " dx=" << dx << " dy=" << dy << std::endl;
                        }
                        cout << endl;
                    }*/
                    //if (jBin == 0 && dens > 1e3)
                    //    cout << "Look: dens = " << dens << "\tnsubpartons = " << nsubpartons << "\tarea = " << area << endl;

                    hDenTime[jBin]->Fill(tau_center, (etasbinbounds[ieta] + etasbinboundsUpper[ieta]) / 2, dens);
                    hAreaTime[jBin]->Fill(tau_center, (etasbinbounds[ieta] + etasbinboundsUpper[ieta]) / 2, area);
                }
            }
        }
    }

    //----------------------------------------------------------------------
    // 绘图并保存
    //----------------------------------------------------------------------

    TString outputName = Form("/eos/cms/store/group/phys_heavyions/huangxi/parton_jet_bins/parton_jet_bins_%d.root", idx);
    TFile *outFile = TFile::Open(outputName, "RECREATE");
    if (!outFile || outFile->IsZombie()) {
        cerr << "Cannot create output file " << outputName << endl;
        return 1;
    }

    hTime->Write();
    hMult->Write();
    hJetMultParton->Write();
    for (int ibin = 0; ibin < trackbin; ibin++)
    {
        hMultTime[ibin]->Write();
        hEccTime[ibin]->Write();
        hDimTime[ibin]->Write();
        hDenTime[ibin]->Write();
        hNpTime[ibin]->Write();
        hAreaTime[ibin]->Write();
        hATest[ibin]->Write();
    }

    outFile->Close();
    file->Close();
    return 0;
}
