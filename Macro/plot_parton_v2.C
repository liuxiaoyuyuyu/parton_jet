#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TProfile.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TString.h>
#include <TMath.h>
#include <iostream>
#include <vector>
#include "../binning.h"

using namespace std;

void plot_parton_v2(const char* inputFileName = "/Users/xl155/Documents/wenbin/root/parton_v2_output_round2.root") {
    
    // Set ROOT style
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kBird);
    gStyle->SetPadGridX(0);
    gStyle->SetPadGridY(0);
    gStyle->SetTitleSize(0.06, "t"); // Title size
    gStyle->SetLabelSize(0.05, "x"); // X-axis label size
    gStyle->SetLabelSize(0.05, "y"); // Y-axis label size
    gStyle->SetTitleSize(0.06, "x"); // X-axis title size
    gStyle->SetTitleSize(0.06, "y"); // Y-axis title size
    gStyle->SetPadLeftMargin(0.12);
    gStyle->SetPadRightMargin(0.05);
    gStyle->SetPadTopMargin(0.08);
    gStyle->SetPadBottomMargin(0.12);
    
    // Rebinning factor for tau axis (increase for coarser bins)
    const int tauRebinFactor = 4; // Change this to rebin more or less
    
    // Open input file
    TFile* inFile = TFile::Open(inputFileName, "READ");
    if (!inFile || inFile->IsZombie()) {
        cout << "Error: Cannot open file " << inputFileName << endl;
        return;
    }
    
    // Create output PDF file path
    TString pdfName = "/Users/xl155/Documents/wenbin/pdf/parton_v2_plots.pdf";
    
    // Page 1: Global QA plots
    TCanvas *c1 = new TCanvas("c1", "Global QA Plots", 1200, 800);
    c1->Divide(2, 2);
    
    // Nparton vs Jet Charged Multiplicity
    c1->cd(1);
    TH2D* hNpartonVsJetMult = (TH2D*)inFile->Get("hNpartonVsJetMult");
    if (hNpartonVsJetMult) {
        hNpartonVsJetMult->Draw("COLZ");
        hNpartonVsJetMult->SetTitle("Nparton vs Jet Charged Multiplicity");
    }
    
    // Nparton vs Tau
    c1->cd(2);
    TH1D* hNpartonVsTau = (TH1D*)inFile->Get("hNpartonVsTau");
    if (hNpartonVsTau) {
        // Rebin the tau axis
        hNpartonVsTau->Rebin(tauRebinFactor);
        hNpartonVsTau->Draw();
        hNpartonVsTau->SetTitle("Nparton vs Tau");
    }
    
    c1->SaveAs(pdfName + "(");
    
    // Page 2: v2 vs tau for all Nch and eta_s bins (2x4 layout for eta_s bins)
    TCanvas *c2 = new TCanvas("c2", "v2 vs Tau - All Eta_s Bins", 1200, 800);
    c2->Divide(2, 4);
    
    for (int j = 0; j < etasbin; j++) {
        c2->cd(j + 1);
        TString name = Form("hv2VsTau_Nch0_etas%d", j);
        TProfile* hv2 = (TProfile*)inFile->Get(name);
        if (hv2) {
            // Rebin the tau axis
            hv2->RebinX(tauRebinFactor);
            hv2->Draw();
            hv2->SetTitle(Form("v2 vs tau (%d < N_{ch} < %d, %.1f < #eta_{s} < %.1f)", 
                               trackbinbounds_MC[0], trackbinboundsUpper_MC[0], 
                               etasbinbounds[j], etasbinboundsUpper[j]));
            hv2->GetXaxis()->SetTitle("#tau (fm/c)");
            hv2->GetYaxis()->SetTitle("<v2>");
        }
    }
    
    c2->SaveAs(pdfName);
    
    // Page 3: Eccentricity vs tau for all eta_s bins (multiplicity bin 0)
    TCanvas *c3 = new TCanvas("c3", "Eccentricity vs Tau - All Eta_s Bins", 1200, 800);
    c3->Divide(2, 4);
    
    for (int j = 0; j < etasbin; j++) {
        c3->cd(j + 1);
        TString name = Form("hEccVsTau_Nch0_etas%d", j);
        TProfile* hEcc = (TProfile*)inFile->Get(name);
        if (hEcc) {
            // Rebin the tau axis
            hEcc->RebinX(tauRebinFactor);
            hEcc->Draw();
            hEcc->SetTitle(Form("Eccentricity vs tau (%d < N_{ch} < %d, %.1f < #eta_{s} < %.1f)", 
                                trackbinbounds_MC[0], trackbinboundsUpper_MC[0], 
                                etasbinbounds[j], etasbinboundsUpper[j]));
            hEcc->GetXaxis()->SetTitle("#tau (fm/c)");
            hEcc->GetYaxis()->SetTitle("<#varepsilon_{2}>");
        }
    }
    
    c3->SaveAs(pdfName);
    
    // Page 4: v2 vs tau for all multiplicity bins (eta_s bin 0)
    TCanvas *c4 = new TCanvas("c4", "v2 vs Tau - All Multiplicity Bins", 1200, 800);
    c4->Divide(3, 3);
    
    for (int i = 0; i < trackbin; i++) {
        c4->cd(i + 1);
        TString name = Form("hv2VsTau_Nch%d_etas0", i);
        TProfile* hv2 = (TProfile*)inFile->Get(name);
        if (hv2) {
            // Rebin the tau axis
            hv2->RebinX(tauRebinFactor);
            hv2->Draw();
            hv2->SetTitle(Form("v2 vs tau (%d < N_{ch} < %d, %.1f < #eta_{s} < %.1f)", 
                               trackbinbounds_MC[i], trackbinboundsUpper_MC[i], 
                               etasbinbounds[0], etasbinboundsUpper[0]));
            hv2->GetXaxis()->SetTitle("#tau (fm/c)");
            hv2->GetYaxis()->SetTitle("<v2>");
        }
    }
    
    c4->SaveAs(pdfName);
    
    // Page 5: Eccentricity vs tau for all multiplicity bins (eta_s bin 0)
    TCanvas *c5 = new TCanvas("c5", "Eccentricity vs Tau - All Multiplicity Bins", 1200, 800);
    c5->Divide(3, 3);
    
    for (int i = 0; i < trackbin; i++) {
        c5->cd(i + 1);
        TString name = Form("hEccVsTau_Nch%d_etas0", i);
        TProfile* hEcc = (TProfile*)inFile->Get(name);
        if (hEcc) {
            // Rebin the tau axis
            hEcc->RebinX(tauRebinFactor);
            hEcc->Draw();
            hEcc->SetTitle(Form("Eccentricity vs tau (%d < N_{ch} < %d, %.1f < #eta_{s} < %.1f)", 
                                trackbinbounds_MC[i], trackbinboundsUpper_MC[i], 
                                etasbinbounds[0], etasbinboundsUpper[0]));
            hEcc->GetXaxis()->SetTitle("#tau (fm/c)");
            hEcc->GetYaxis()->SetTitle("<#varepsilon_{2}>");
        }
    }
    
    c5->SaveAs(pdfName);
    
    // Page 6: v2 vs tau for all Nch bins (eta_s bin 1)
    TCanvas *c6a = new TCanvas("c6a", "v2 vs Tau - All Nch Bins (eta_s bin 1)", 1200, 800);
    c6a->Divide(3, 3);
    
    for (int i = 0; i < trackbin; i++) {
        c6a->cd(i + 1);
        TString name = Form("hv2VsTau_Nch%d_etas1", i);
        TProfile* hv2 = (TProfile*)inFile->Get(name);
        if (hv2) {
            // Rebin the tau axis
            hv2->RebinX(tauRebinFactor);
            hv2->Draw();
            hv2->SetTitle(Form("v2 vs tau (%d < N_{ch} < %d, %.1f < #eta_{s} < %.1f)", 
                               trackbinbounds_MC[i], trackbinboundsUpper_MC[i], 
                               etasbinbounds[1], etasbinboundsUpper[1]));
            hv2->GetXaxis()->SetTitle("#tau (fm/c)");
            hv2->GetYaxis()->SetTitle("<v2>");
        }
    }
    
    c6a->SaveAs(pdfName);
    
    // Page 7: v2 vs tau for all Nch bins (eta_s bin 2)
    TCanvas *c7a = new TCanvas("c7a", "v2 vs Tau - All Nch Bins (eta_s bin 2)", 1200, 800);
    c7a->Divide(3, 3);
    
    for (int i = 0; i < trackbin; i++) {
        c7a->cd(i + 1);
        TString name = Form("hv2VsTau_Nch%d_etas2", i);
        TProfile* hv2 = (TProfile*)inFile->Get(name);
        if (hv2) {
            // Rebin the tau axis
            hv2->RebinX(tauRebinFactor);
            hv2->Draw();
            hv2->SetTitle(Form("v2 vs tau (%d < N_{ch} < %d, %.1f < #eta_{s} < %.1f)", 
                               trackbinbounds_MC[i], trackbinboundsUpper_MC[i], 
                               etasbinbounds[2], etasbinboundsUpper[2]));
            hv2->GetXaxis()->SetTitle("#tau (fm/c)");
            hv2->GetYaxis()->SetTitle("<v2>");
        }
    }
    
    c7a->SaveAs(pdfName);
    
    // Page 8: v2 vs tau for all Nch bins (eta_s bin 3)
    TCanvas *c8a = new TCanvas("c8a", "v2 vs Tau - All Nch Bins (eta_s bin 3)", 1200, 800);
    c8a->Divide(3, 3);
    
    for (int i = 0; i < trackbin; i++) {
        c8a->cd(i + 1);
        TString name = Form("hv2VsTau_Nch%d_etas3", i);
        TProfile* hv2 = (TProfile*)inFile->Get(name);
        if (hv2) {
            // Rebin the tau axis
            hv2->RebinX(tauRebinFactor);
            hv2->Draw();
            hv2->SetTitle(Form("v2 vs tau (%d < N_{ch} < %d, %.1f < #eta_{s} < %.1f)", 
                               trackbinbounds_MC[i], trackbinboundsUpper_MC[i], 
                               etasbinbounds[3], etasbinboundsUpper[3]));
            hv2->GetXaxis()->SetTitle("#tau (fm/c)");
            hv2->GetYaxis()->SetTitle("<v2>");
        }
    }
    
    c8a->SaveAs(pdfName);
    
    // Page 9: v2 vs tau for all Nch bins (eta_s bin 4)
    TCanvas *c9a = new TCanvas("c9a", "v2 vs Tau - All Nch Bins (eta_s bin 4)", 1200, 800);
    c9a->Divide(3, 3);
    
    for (int i = 0; i < trackbin; i++) {
        c9a->cd(i + 1);
        TString name = Form("hv2VsTau_Nch%d_etas4", i);
        TProfile* hv2 = (TProfile*)inFile->Get(name);
        if (hv2) {
            // Rebin the tau axis
            hv2->RebinX(tauRebinFactor);
            hv2->Draw();
            hv2->SetTitle(Form("v2 vs tau (%d < N_{ch} < %d, %.1f < #eta_{s} < %.1f)", 
                               trackbinbounds_MC[i], trackbinboundsUpper_MC[i], 
                               etasbinbounds[4], etasbinboundsUpper[4]));
            hv2->GetXaxis()->SetTitle("#tau (fm/c)");
            hv2->GetYaxis()->SetTitle("<v2>");
        }
    }
    
    c9a->SaveAs(pdfName);
    
    // Page 10: v2 vs tau for all Nch bins (eta_s bin 5)
    TCanvas *c10a = new TCanvas("c10a", "v2 vs Tau - All Nch Bins (eta_s bin 5)", 1200, 800);
    c10a->Divide(3, 3);
    
    for (int i = 0; i < trackbin; i++) {
        c10a->cd(i + 1);
        TString name = Form("hv2VsTau_Nch%d_etas5", i);
        TProfile* hv2 = (TProfile*)inFile->Get(name);
        if (hv2) {
            // Rebin the tau axis
            hv2->RebinX(tauRebinFactor);
            hv2->Draw();
            hv2->SetTitle(Form("v2 vs tau (%d < N_{ch} < %d, %.1f < #eta_{s} < %.1f)", 
                               trackbinbounds_MC[i], trackbinboundsUpper_MC[i], 
                               etasbinbounds[5], etasbinboundsUpper[5]));
            hv2->GetXaxis()->SetTitle("#tau (fm/c)");
            hv2->GetYaxis()->SetTitle("<v2>");
        }
    }
    
    c10a->SaveAs(pdfName);
    
    // Page 11: v2 vs tau for all Nch bins (eta_s bin 6)
    TCanvas *c11a = new TCanvas("c11a", "v2 vs Tau - All Nch Bins (eta_s bin 6)", 1200, 800);
    c11a->Divide(3, 3);
    
    for (int i = 0; i < trackbin; i++) {
        c11a->cd(i + 1);
        TString name = Form("hv2VsTau_Nch%d_etas6", i);
        TProfile* hv2 = (TProfile*)inFile->Get(name);
        if (hv2) {
            // Rebin the tau axis
            hv2->RebinX(tauRebinFactor);
            hv2->Draw();
            hv2->SetTitle(Form("v2 vs tau (%d < N_{ch} < %d, %.1f < #eta_{s} < %.1f)", 
                               trackbinbounds_MC[i], trackbinboundsUpper_MC[i], 
                               etasbinbounds[6], etasbinboundsUpper[6]));
            hv2->GetXaxis()->SetTitle("#tau (fm/c)");
            hv2->GetYaxis()->SetTitle("<v2>");
        }
    }
    
    c11a->SaveAs(pdfName);
    
    // Page 12: v2 vs tau for all Nch bins (eta_s bin 7)
    TCanvas *c12a = new TCanvas("c12a", "v2 vs Tau - All Nch Bins (eta_s bin 7)", 1200, 800);
    c12a->Divide(3, 3);
    
    for (int i = 0; i < trackbin; i++) {
        c12a->cd(i + 1);
        TString name = Form("hv2VsTau_Nch%d_etas7", i);
        TProfile* hv2 = (TProfile*)inFile->Get(name);
        if (hv2) {
            // Rebin the tau axis
            hv2->RebinX(tauRebinFactor);
            hv2->Draw();
            hv2->SetTitle(Form("v2 vs tau (%d < N_{ch} < %d, %.1f < #eta_{s} < %.1f)", 
                               trackbinbounds_MC[i], trackbinboundsUpper_MC[i], 
                               etasbinbounds[7], etasbinboundsUpper[7]));
            hv2->GetXaxis()->SetTitle("#tau (fm/c)");
            hv2->GetYaxis()->SetTitle("<v2>");
        }
    }
    
    c12a->SaveAs(pdfName);
    
    // Page 13: Eccentricity vs tau for all Nch bins (eta_s bin 1)
    TCanvas *c13a = new TCanvas("c13a", "Eccentricity vs Tau - All Nch Bins (eta_s bin 1)", 1200, 800);
    c13a->Divide(3, 3);
    
    for (int i = 0; i < trackbin; i++) {
        c13a->cd(i + 1);
        TString name = Form("hEccVsTau_Nch%d_etas1", i);
        TProfile* hEcc = (TProfile*)inFile->Get(name);
        if (hEcc) {
            // Rebin the tau axis
            hEcc->RebinX(tauRebinFactor);
            hEcc->Draw();
            hEcc->SetTitle(Form("Eccentricity vs tau (%d < N_{ch} < %d, %.1f < #eta_{s} < %.1f)", 
                                trackbinbounds_MC[i], trackbinboundsUpper_MC[i], 
                                etasbinbounds[1], etasbinboundsUpper[1]));
            hEcc->GetXaxis()->SetTitle("#tau (fm/c)");
            hEcc->GetYaxis()->SetTitle("<#varepsilon_{2}>");
        }
    }
    
    c13a->SaveAs(pdfName);
    
    // Page 14: Eccentricity vs tau for all Nch bins (eta_s bin 2)
    TCanvas *c14a = new TCanvas("c14a", "Eccentricity vs Tau - All Nch Bins (eta_s bin 2)", 1200, 800);
    c14a->Divide(3, 3);
    
    for (int i = 0; i < trackbin; i++) {
        c14a->cd(i + 1);
        TString name = Form("hEccVsTau_Nch%d_etas2", i);
        TProfile* hEcc = (TProfile*)inFile->Get(name);
        if (hEcc) {
            // Rebin the tau axis
            hEcc->RebinX(tauRebinFactor);
            hEcc->Draw();
            hEcc->SetTitle(Form("Eccentricity vs tau (%d < N_{ch} < %d, %.1f < #eta_{s} < %.1f)", 
                                trackbinbounds_MC[i], trackbinboundsUpper_MC[i], 
                                etasbinbounds[2], etasbinboundsUpper[2]));
            hEcc->GetXaxis()->SetTitle("#tau (fm/c)");
            hEcc->GetYaxis()->SetTitle("<#varepsilon_{2}>");
        }
    }
    
    c14a->SaveAs(pdfName);
    
    // Page 15: Eccentricity vs tau for all Nch bins (eta_s bin 3)
    TCanvas *c15a = new TCanvas("c15a", "Eccentricity vs Tau - All Nch Bins (eta_s bin 3)", 1200, 800);
    c15a->Divide(3, 3);
    
    for (int i = 0; i < trackbin; i++) {
        c15a->cd(i + 1);
        TString name = Form("hEccVsTau_Nch%d_etas3", i);
        TProfile* hEcc = (TProfile*)inFile->Get(name);
        if (hEcc) {
            // Rebin the tau axis
            hEcc->RebinX(tauRebinFactor);
            hEcc->Draw();
            hEcc->SetTitle(Form("Eccentricity vs tau (%d < N_{ch} < %d, %.1f < #eta_{s} < %.1f)", 
                                trackbinbounds_MC[i], trackbinboundsUpper_MC[i], 
                                etasbinbounds[3], etasbinboundsUpper[3]));
            hEcc->GetXaxis()->SetTitle("#tau (fm/c)");
            hEcc->GetYaxis()->SetTitle("<#varepsilon_{2}>");
        }
    }
    
    c15a->SaveAs(pdfName);
    
    // Page 16: Eccentricity vs tau for all Nch bins (eta_s bin 4)
    TCanvas *c16a = new TCanvas("c16a", "Eccentricity vs Tau - All Nch Bins (eta_s bin 4)", 1200, 800);
    c16a->Divide(3, 3);
    
    for (int i = 0; i < trackbin; i++) {
        c16a->cd(i + 1);
        TString name = Form("hEccVsTau_Nch%d_etas4", i);
        TProfile* hEcc = (TProfile*)inFile->Get(name);
        if (hEcc) {
            // Rebin the tau axis
            hEcc->RebinX(tauRebinFactor);
            hEcc->Draw();
            hEcc->SetTitle(Form("Eccentricity vs tau (%d < N_{ch} < %d, %.1f < #eta_{s} < %.1f)", 
                                trackbinbounds_MC[i], trackbinboundsUpper_MC[i], 
                                etasbinbounds[4], etasbinboundsUpper[4]));
            hEcc->GetXaxis()->SetTitle("#tau (fm/c)");
            hEcc->GetYaxis()->SetTitle("<#varepsilon_{2}>");
        }
    }
    
    c16a->SaveAs(pdfName);
    
    // Page 17: Eccentricity vs tau for all Nch bins (eta_s bin 5)
    TCanvas *c17a = new TCanvas("c17a", "Eccentricity vs Tau - All Nch Bins (eta_s bin 5)", 1200, 800);
    c17a->Divide(3, 3);
    
    for (int i = 0; i < trackbin; i++) {
        c17a->cd(i + 1);
        TString name = Form("hEccVsTau_Nch%d_etas5", i);
        TProfile* hEcc = (TProfile*)inFile->Get(name);
        if (hEcc) {
            // Rebin the tau axis
            hEcc->RebinX(tauRebinFactor);
            hEcc->Draw();
            hEcc->SetTitle(Form("Eccentricity vs tau (%d < N_{ch} < %d, %.1f < #eta_{s} < %.1f)", 
                                trackbinbounds_MC[i], trackbinboundsUpper_MC[i], 
                                etasbinbounds[5], etasbinboundsUpper[5]));
            hEcc->GetXaxis()->SetTitle("#tau (fm/c)");
            hEcc->GetYaxis()->SetTitle("<#varepsilon_{2}>");
        }
    }
    
    c17a->SaveAs(pdfName);
    
    // Page 18: Eccentricity vs tau for all Nch bins (eta_s bin 6)
    TCanvas *c18a = new TCanvas("c18a", "Eccentricity vs Tau - All Nch Bins (eta_s bin 6)", 1200, 800);
    c18a->Divide(3, 3);
    
    for (int i = 0; i < trackbin; i++) {
        c18a->cd(i + 1);
        TString name = Form("hEccVsTau_Nch%d_etas6", i);
        TProfile* hEcc = (TProfile*)inFile->Get(name);
        if (hEcc) {
            // Rebin the tau axis
            hEcc->RebinX(tauRebinFactor);
            hEcc->Draw();
            hEcc->SetTitle(Form("Eccentricity vs tau (%d < N_{ch} < %d, %.1f < #eta_{s} < %.1f)", 
                                trackbinbounds_MC[i], trackbinboundsUpper_MC[i], 
                                etasbinbounds[6], etasbinboundsUpper[6]));
            hEcc->GetXaxis()->SetTitle("#tau (fm/c)");
            hEcc->GetYaxis()->SetTitle("<#varepsilon_{2}>");
        }
    }
    
    c18a->SaveAs(pdfName);
    
    // Page 19: Eccentricity vs tau for all Nch bins (eta_s bin 7)
    TCanvas *c19a = new TCanvas("c19a", "Eccentricity vs Tau - All Nch Bins (eta_s bin 7)", 1200, 800);
    c19a->Divide(3, 3);
    
    for (int i = 0; i < trackbin; i++) {
        c19a->cd(i + 1);
        TString name = Form("hEccVsTau_Nch%d_etas7", i);
        TProfile* hEcc = (TProfile*)inFile->Get(name);
        if (hEcc) {
            // Rebin the tau axis
            hEcc->RebinX(tauRebinFactor);
            hEcc->Draw();
            hEcc->SetTitle(Form("Eccentricity vs tau (%d < N_{ch} < %d, %.1f < #eta_{s} < %.1f)", 
                                trackbinbounds_MC[i], trackbinboundsUpper_MC[i], 
                                etasbinbounds[7], etasbinboundsUpper[7]));
            hEcc->GetXaxis()->SetTitle("#tau (fm/c)");
            hEcc->GetYaxis()->SetTitle("<#varepsilon_{2}>");
        }
    }
    
    c19a->SaveAs(pdfName);
    
    // Page 20: v2 RMS vs tau for different multiplicity bins (eta_s bin 0)
    TCanvas *c6 = new TCanvas("c6", "v2 RMS vs Tau - Multiplicity Bins", 1200, 800);
    c6->Divide(3, 3);
    
    for (int i = 0; i < trackbin; i++) {
        c6->cd(i + 1);
        TString name = Form("hv2RMSvsTau_Nch%d_etas0", i);
        TH1D* hv2RMS = (TH1D*)inFile->Get(name);
        if (hv2RMS) {
            // Rebin the tau axis
            hv2RMS->Rebin(tauRebinFactor);
            hv2RMS->Draw();
            hv2RMS->SetTitle(Form("v2 RMS vs tau (%d < N_{ch} < %d, %.1f < #eta_{s} < %.1f)", 
                                  trackbinbounds_MC[i], trackbinboundsUpper_MC[i], 
                                  etasbinbounds[0], etasbinboundsUpper[0]));
            hv2RMS->GetXaxis()->SetTitle("#tau (fm/c)");
            hv2RMS->GetYaxis()->SetTitle("v2^{RMS}");
        }
    }
    
    c6->SaveAs(pdfName);
    
    // Page 7: Eccentricity RMS vs tau for different multiplicity bins (eta_s bin 0)
    TCanvas *c7 = new TCanvas("c7", "Eccentricity RMS vs Tau - Multiplicity Bins", 1200, 800);
    c7->Divide(3, 3);
    
    for (int i = 0; i < trackbin; i++) {
        c7->cd(i + 1);
        TString name = Form("hEccRMSvsTau_Nch%d_etas0", i);
        TH1D* hEccRMS = (TH1D*)inFile->Get(name);
        if (hEccRMS) {
            // Rebin the tau axis
            hEccRMS->Rebin(tauRebinFactor);
            hEccRMS->Draw();
            hEccRMS->SetTitle(Form("Eccentricity RMS vs tau (%d < N_{ch} < %d, %.1f < #eta_{s} < %.1f)", 
                                   trackbinbounds_MC[i], trackbinboundsUpper_MC[i], 
                                   etasbinbounds[0], etasbinboundsUpper[0]));
            hEccRMS->GetXaxis()->SetTitle("#tau (fm/c)");
            hEccRMS->GetYaxis()->SetTitle("#varepsilon_{2}^{RMS}");
        }
    }
    
    c7->SaveAs(pdfName);
    
    // Page 8: v2 vs Eccentricity correlation plots
    TCanvas *c8 = new TCanvas("c8", "v2 vs Eccentricity Correlation", 1200, 800);
    c8->Divide(3, 3);
    
    for (int i = 0; i < trackbin; i++) {
        c8->cd(i + 1);
        TString name = Form("hv2vsEcc_Nch%d_etas0", i);
        TH2D* hv2vsEcc = (TH2D*)inFile->Get(name);
        if (hv2vsEcc) {
            hv2vsEcc->Draw("COLZ");
            hv2vsEcc->SetTitle(Form("v2 vs Eccentricity (%d < N_{ch} < %d, %.1f < #eta_{s} < %.1f)", 
                                    trackbinbounds_MC[i], trackbinboundsUpper_MC[i], 
                                    etasbinbounds[0], etasbinboundsUpper[0]));
            hv2vsEcc->GetXaxis()->SetTitle("#varepsilon_{2}");
            hv2vsEcc->GetYaxis()->SetTitle("v2");
        }
    }
    
    c8->SaveAs(pdfName);
    
    // Page 9: Psi2 distribution vs tau (multiplicity bin 0, all eta_s bins)
    TCanvas *c9 = new TCanvas("c9", "Psi2 Distribution vs Tau", 1200, 800);
    c9->Divide(2, 4);
    
    for (int j = 0; j < etasbin; j++) {
        c9->cd(j + 1);
        TString name = Form("hPsi2_Nch0_etas%d", j);
        TH2D* hPsi2 = (TH2D*)inFile->Get(name);
        if (hPsi2) {
            // Rebin the tau axis (X-axis for TH2D)
            hPsi2->RebinX(tauRebinFactor);
            hPsi2->Draw("COLZ");
            hPsi2->SetTitle(Form("Psi2 Distribution vs Tau (%d < N_{ch} < %d, %.1f < #eta_{s} < %.1f)", 
                                 trackbinbounds_MC[0], trackbinboundsUpper_MC[0], 
                                 etasbinbounds[j], etasbinboundsUpper[j]));
            hPsi2->GetXaxis()->SetTitle("#tau (fm/c)");
            hPsi2->GetYaxis()->SetTitle("#Psi_{2}");
        }
    }
    
    c9->SaveAs(pdfName);
    
    // Page 10: Binned QA plots - Nparton vs Jet Mult
    TCanvas *c11 = new TCanvas("c11", "Binned QA - Nparton vs Jet Mult", 1200, 800);
    c11->Divide(3, 3);
    
    for (int i = 0; i < trackbin; i++) {
        c11->cd(i + 1);
        TString name = Form("hNpartonVsJetMult_Nch%d_etas0", i);
        TH2D* hQA = (TH2D*)inFile->Get(name);
        if (hQA) {
            hQA->Draw("COLZ");
            hQA->SetTitle(Form("Nparton vs Jet Mult (%d < N_{ch} < %d, %.1f < #eta_{s} < %.1f)", 
                               trackbinbounds_MC[i], trackbinboundsUpper_MC[i], 
                               etasbinbounds[0], etasbinboundsUpper[0]));
            hQA->GetXaxis()->SetTitle("Jet Charged Multiplicity");
            hQA->GetYaxis()->SetTitle("Nparton");
        }
    }
    
    c11->SaveAs(pdfName);
    
    // Page 11: Binned QA plots - Nparton vs Tau
    TCanvas *c12 = new TCanvas("c12", "Binned QA - Nparton vs Tau", 1200, 800);
    c12->Divide(3, 3);
    
    for (int i = 0; i < trackbin; i++) {
        c12->cd(i + 1);
        TString name = Form("hNpartonVsTau_Nch%d_etas0", i);
        TH1D* hQA = (TH1D*)inFile->Get(name);
        if (hQA) {
            // Rebin the tau axis
            hQA->Rebin(tauRebinFactor);
            hQA->Draw();
            hQA->SetTitle(Form("Nparton vs Tau (%d < N_{ch} < %d, %.1f < #eta_{s} < %.1f)", 
                               trackbinbounds_MC[i], trackbinboundsUpper_MC[i], 
                               etasbinbounds[0], etasbinboundsUpper[0]));
            hQA->GetXaxis()->SetTitle("#tau (fm/c)");
            hQA->GetYaxis()->SetTitle("Nparton");
        }
    }
    
    c12->SaveAs(pdfName + ")");
    
    // Clean up
    inFile->Close();
    delete c1;
    delete c2;
    delete c3;
    delete c4;
    delete c5;
    delete c6a;
    delete c7a;
    delete c8a;
    delete c9a;
    delete c10a;
    delete c11a;
    delete c12a;
    delete c13a;
    delete c14a;
    delete c15a;
    delete c16a;
    delete c17a;
    delete c18a;
    delete c19a;
    delete c6;
    delete c7;
    delete c8;
    delete c9;
    delete c11;
    delete c12;
    
    cout << "Plots saved to: " << pdfName << endl;
}

// Function to plot multiple files
void plot_multiple_files(const char* fileList = "filelist.txt") {
    ifstream file(fileList);
    string line;
    
    while (getline(file, line)) {
        if (!line.empty() && line[0] != '#') {
            cout << "Plotting: " << line << endl;
            plot_parton_v2(line.c_str());
        }
    }
    
    file.close();
} 