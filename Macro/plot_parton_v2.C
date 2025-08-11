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

// Generic function to plot histograms with Nch and eta_s binning
void plotBinnedHistograms(TFile* inFile, const char* histPattern, const char* plotTitle, 
                          const char* xTitle, const char* yTitle, int tauRebinFactor, 
                          const char* pdfName, bool isTProfile = false, bool isTH2D = false, 
                          bool needsTauRebin = false, double yMin = 0, double yMax = 0, 
                          const vector<int>& etaBinsToPlot = {}) {
    
    for (int etaBin = 0; etaBin < etasbin; etaBin++) {
        // Skip if specific eta bins are requested and this one is not in the list
        if (!etaBinsToPlot.empty() && find(etaBinsToPlot.begin(), etaBinsToPlot.end(), etaBin) == etaBinsToPlot.end()) {
            continue;
        }
        
        TCanvas *c = new TCanvas(Form("c_eta%d", etaBin), Form("%s - Eta_s bin %d", plotTitle, etaBin), 1200, 800);
        c->Divide(3, 3);
        
        for (int multBin = 0; multBin < trackbin; multBin++) {
            c->cd(multBin + 1);
            TString name = Form(histPattern, multBin, etaBin);
            
            if (isTProfile) {
                TProfile* h = (TProfile*)inFile->Get(name);
                if (h) {
                    if (needsTauRebin) h->RebinX(tauRebinFactor);
                    h->Draw();
                    h->SetTitle(Form("%s (%d < N_{ch} < %d, %.1f < #eta_{s} < %.1f)", 
                                    plotTitle, trackbinbounds_MC[multBin], trackbinboundsUpper_MC[multBin], 
                                    etasbinbounds[etaBin], etasbinboundsUpper[etaBin]));
                    h->GetXaxis()->SetTitle(xTitle);
                    h->GetYaxis()->SetTitle(yTitle);
                    // Set Y-axis range for v2 and eccentricity plots
                    if (yMin != 0 || yMax != 0) {
                        h->SetMinimum(yMin);
                        h->SetMaximum(yMax);
                    }
                }
            } else if (isTH2D) {
                TH2D* h = (TH2D*)inFile->Get(name);
                if (h) {
                    if (needsTauRebin) h->RebinX(tauRebinFactor);
                    h->Draw("COLZ");
                    h->SetTitle(Form("%s (%d < N_{ch} < %d, %.1f < #eta_{s} < %.1f)", 
                                    plotTitle, trackbinbounds_MC[multBin], trackbinboundsUpper_MC[multBin], 
                                    etasbinbounds[etaBin], etasbinboundsUpper[etaBin]));
                    h->GetXaxis()->SetTitle(xTitle);
                    h->GetYaxis()->SetTitle(yTitle);
                }
            } else {
                TH1D* h = (TH1D*)inFile->Get(name);
                if (h) {
                    if (needsTauRebin) h->Rebin(tauRebinFactor);
                    h->Draw();
                    h->SetTitle(Form("%s (%d < N_{ch} < %d, %.1f < #eta_{s} < %.1f)", 
                                    plotTitle, trackbinbounds_MC[multBin], trackbinboundsUpper_MC[multBin], 
                                    etasbinbounds[etaBin], etasbinboundsUpper[etaBin]));
                    h->GetXaxis()->SetTitle(xTitle);
                    h->GetYaxis()->SetTitle(yTitle);
                }
            }
        }
        
        c->SaveAs(pdfName);
        delete c;
    }
}

void plot_parton_v2(const char* inputFileName = "/Users/xl155/Documents/wenbin/root/parton_v2_output_round4.root", 
                   const vector<int>& etaBinsToPlot = {}) {
    
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
    const int tauRebinFactor = 1; // Change this to rebin more or less
    
    // Open input file
    TFile* inFile = TFile::Open(inputFileName, "READ");
    if (!inFile || inFile->IsZombie()) {
        cout << "Error: Cannot open file " << inputFileName << endl;
        return;
    }
    
    // Create output PDF file path
    TString pdfName = "/Users/xl155/Documents/wenbin/pdf/parton_v2_plots_round4.pdf";
    
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
    TProfile* hNpartonVsTau = (TProfile*)inFile->Get("hNpartonVsTau");
    if (hNpartonVsTau) {
        // Rebin the tau axis
        hNpartonVsTau->RebinX(tauRebinFactor);
        hNpartonVsTau->Draw();
        hNpartonVsTau->SetTitle("Nparton vs Tau");
    }
    
    c1->SaveAs(pdfName + "(");
    
    // Pages 2-9: v2 vs tau
    plotBinnedHistograms(inFile, "hv2VsTau_Nch%d_etas%d", "v2 vs tau", "#tau (fm/c)", "<v2>", 
                        tauRebinFactor, pdfName, true, false, true, 0.0, 1.0, etaBinsToPlot);
    
    // Pages 10-17: Eccentricity vs tau
    plotBinnedHistograms(inFile, "hEccVsTau_Nch%d_etas%d", "Eccentricity vs tau", "#tau (fm/c)", "<#varepsilon_{2}>", 
                        tauRebinFactor, pdfName, true, false, true, 0.0, 1.0, etaBinsToPlot);
    
    // Pages 18-25: v2 RMS vs tau - COMMENTED OUT
    // plotBinnedHistograms(inFile, "hv2RMSvsTau_Nch%d_etas%d", "v2 RMS vs tau", "#tau (fm/c)", "v2^{RMS}", 
    //                     tauRebinFactor, pdfName, false, false, true);
    
    // Pages 26-33: Eccentricity RMS vs tau - COMMENTED OUT
    // plotBinnedHistograms(inFile, "hEccRMSvsTau_Nch%d_etas%d", "Eccentricity RMS vs tau", "#tau (fm/c)", "#varepsilon_{2}^{RMS}", 
    //                     tauRebinFactor, pdfName, false, false, true);
    
    // Pages 18-25: Nparton vs Jet Mult
    plotBinnedHistograms(inFile, "hNpartonVsJetMult_Nch%d_etas%d", "Nparton vs Jet Mult", "Jet Charged Multiplicity", "Nparton", 
                        tauRebinFactor, pdfName, false, true, false, 0.0, 0.0, etaBinsToPlot);
    
    // Pages 26-33: Nparton vs Tau
    plotBinnedHistograms(inFile, "hNpartonVsTau_Nch%d_etas%d", "Nparton vs Tau", "#tau (fm/c)", "<Nparton>", 
                        tauRebinFactor, pdfName, true, false, true, 0.0, 0.0, etaBinsToPlot);
    
    // Pages 34-41: v2 vs Eccentricity correlation
    plotBinnedHistograms(inFile, "hv2vsEcc_Nch%d_etas%d", "v2 vs Eccentricity", "#varepsilon_{2}", "v2", 
                        tauRebinFactor, pdfName, false, true, false, 0.0, 0.0, etaBinsToPlot);
    
    // Pages 42-49: Psi2 Distribution vs Tau
    plotBinnedHistograms(inFile, "hPsi2_Nch%d_etas%d", "Psi2 Distribution vs Tau", "#tau (fm/c)", "#Psi_{2}", 
                        tauRebinFactor, pdfName, false, true, true, 0.0, 0.0, etaBinsToPlot);
    
    // Close the PDF file properly
    TCanvas *cClose = new TCanvas("cClose", "Close", 1, 1);
    cClose->SaveAs(pdfName + ")");
    delete cClose;
    
    // Clean up
    inFile->Close();
    delete c1;
    
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

// Helper functions for plotting specific eta_s bins
void plot_parton_v2_eta0(const char* inputFileName = "/Users/xl155/Documents/wenbin/root/parton_v2_output_round2.root") {
    vector<int> etaBins = {0};
    plot_parton_v2(inputFileName, etaBins);
}

void plot_parton_v2_eta1(const char* inputFileName = "/Users/xl155/Documents/wenbin/root/parton_v2_output_round2.root") {
    vector<int> etaBins = {1};
    plot_parton_v2(inputFileName, etaBins);
}

void plot_parton_v2_eta2(const char* inputFileName = "/Users/xl155/Documents/wenbin/root/parton_v2_output_round2.root") {
    vector<int> etaBins = {2};
    plot_parton_v2(inputFileName, etaBins);
}

void plot_parton_v2_eta3(const char* inputFileName = "/Users/xl155/Documents/wenbin/root/parton_v2_output_round2.root") {
    vector<int> etaBins = {3};
    plot_parton_v2(inputFileName, etaBins);
}

void plot_parton_v2_eta4(const char* inputFileName = "/Users/xl155/Documents/wenbin/root/parton_v2_output_round2.root") {
    vector<int> etaBins = {4};
    plot_parton_v2(inputFileName, etaBins);
}

void plot_parton_v2_eta5(const char* inputFileName = "/Users/xl155/Documents/wenbin/root/parton_v2_output_round2.root") {
    vector<int> etaBins = {5};
    plot_parton_v2(inputFileName, etaBins);
}

void plot_parton_v2_eta6(const char* inputFileName = "/Users/xl155/Documents/wenbin/root/parton_v2_output_round2.root") {
    vector<int> etaBins = {6};
    plot_parton_v2(inputFileName, etaBins);
}

void plot_parton_v2_eta7(const char* inputFileName = "/Users/xl155/Documents/wenbin/root/parton_v2_output_round2.root") {
    vector<int> etaBins = {7};
    plot_parton_v2(inputFileName, etaBins);
}

// Function to plot multiple eta_s bins
void plot_parton_v2_multiple_etas(const char* inputFileName = "/Users/xl155/Documents/wenbin/root/parton_v2_output_round2.root", 
                                 const vector<int>& etaBins = {0, 1, 2, 3}) {
    plot_parton_v2(inputFileName, etaBins);
} 