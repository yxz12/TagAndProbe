#include <iostream>
#include <cmath>
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH1.h"
#include "TF1.h"
#include "TMinuit.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TStyle.h"
#include "TAttAxis.h"
#include "TFitResult.h"
#include "TLatex.h"
#include <fstream>
#include <TSystem.h>
#include "TLegend.h"
#include "TGaxis.h"

int tagandprobe() {
    
    double DeltaR(double eta_1, double phi_1, double eta_2, double phi_2) {
        double deltaR=sqrt((eta_1-eta_2)*(eta_1-eta_2) + (phi_1-phi_2)*(phi_1-phi_2));
        return deltaR;
    }
    
    //--------------------------------------------------------------------------------------
    //                      1. Extract the root file and the tree
    //--------------------------------------------------------------------------------------
    TFile *file = new TFile("../Trees111115/DYJetsToLL_tpzee_2015.root");
    TTree *tree = dynamic_cast<TTree*>(file->Get("ntuple"));
    
    //--------------------------------------------------------------------------------------
    //                      1.1 Extract the properties of the particles
    //--------------------------------------------------------------------------------------
    double id_1, id_2, iso_1, iso_2;
    int q_1, q_2;
    tree->SetBranchAddress("id_1",&id_1);
    tree->SetBranchAddress("id_2",&id_2);
    tree->SetBranchAddress("iso_1",&iso_1);
    tree->SetBranchAddress("iso_2",&iso_2);
    tree->SetBranchAddress("q_1",&q_1);
    tree->SetBranchAddress("q_2",&q_2);
    double pt_1,pt_2,E_1,E_2,eta_1,eta_2,phi_1,phi_2,m_vis;
    tree->SetBranchAddress("pt_1", &pt_1);
    tree->SetBranchAddress("pt_2", &pt_2);
    tree->SetBranchAddress("E_1", &E_1);
    tree->SetBranchAddress("E_2", &E_2);
    tree->SetBranchAddress("eta_1", &eta_1);
    tree->SetBranchAddress("eta_2", &eta_2);
    tree->SetBranchAddress("phi_1", &phi_1);
    tree->SetBranchAddress("phi_2", &phi_2);
    tree->SetBranchAddress("m_vis", &m_vis);
    bool os;
    tree->SetBranchAddress("os", &os);
    double pt_tag,pt_probe,E_tag,E_probe,eta_tag,eta_probe,phi_tag,phi_probe;
    int q_tag, q_probe;

    int nPtBins = 8;
    double ptBins[9] = {10,13,16,20,25,30,40,60,10000};
    int nEtaBins = 3;
    float etaBins[4] = {0,0.9,1.2,2.4};
    TString PtBins[8] = {"Pt10to13",
                      "Pt13to16",
                      "Pt16to20",
                      "Pt20to25",
                      "Pt25to30",
                      "Pt30to40",
                      "Pt40to60",
                      "PtGt60"};

    TString EtaBins[3] = {"EtaLt0p9",
                       "Eta0p9to1p2",
                       "EtaGt1p2"};
      
 
    //--------------------------------------------------------------------------------------
    //                      1.1 Create the histograms for storing the properties
    //--------------------------------------------------------------------------------------
    
    TH1F *id_pass=new TH1F("id_pass","",500,0,1000);
    TH1F *id_fail=new TH1F("id_fail","",500,0,1000);
    TH1F *iso_pass=new TH1F("iso_pass","",500,0,1000);
    TH1F *iso_fail=new TH1F("iso_fail","",500,0,1000);
    TH1F *id_iso_pass=new TH1F("id_iso_pass","",800,0,1000);
    TH1F *id_iso_fail=new TH1F("id_iso_fail","",800,0,1000);
    
    TH1F *id_iso_pass_1=new TH1F("id_iso_pass_1","",800,0,1000);
    TH1F *id_iso_fail_1=new TH1F("id_iso_fail_1","",800,0,1000);
    TH1F *id_iso_pass_2=new TH1F("id_iso_pass_2","",800,0,1000);
    TH1F *id_iso_fail_2=new TH1F("id_iso_fail_2","",800,0,1000);
    TH1F *id_iso_pass_3=new TH1F("id_iso_pass_3","",800,0,1000);
    TH1F *id_iso_fail_3=new TH1F("id_iso_fail_3","",800,0,1000);
    TH1F *id_iso_pass_4=new TH1F("id_iso_pass_4","",800,0,1000);
    TH1F *id_iso_fail_4=new TH1F("id_iso_fail_4","",800,0,1000);

//    TH1F *eta_tag=new TH1F("eta_tag","",100,-2.5,2.5);
//    TH1F *eta_probe=new TH1F("eta_probe","",100,-2.5,2.5);
    
    // keep track of the statistic errors on the histogram, before filling
    //m_vis_pass->Sumw2();
    //m_vis_fail->Sumw2();
    TH1F* ZMassEtaPtPass[3][9];
    TH1F* ZMassEtaPtFail[3][9];
    for (int iEta=0; iEta<nEtaBins; ++iEta) {
        for (int iPt=0; iPt<nPtBins; ++iPt) {
            ZMassEtaPtPass[iEta][iPt] = new TH1F("ZMass"+EtaBins[iEta]+PtBins[iPt]+"Pass","",80,50,130);
            ZMassEtaPtFail[iEta][iPt] = new TH1F("ZMass"+EtaBins[iEta]+PtBins[iPt]+"Fail","",80,50,130);
        }
    }
    


    //--------------------------------------------------------------------------------------
    //                      1.3 Applying selection criteria - Tag/probe
    //--------------------------------------------------------------------------------------
    for (int i = 0; i < tree->GetEntries(); ++i) {
        tree->GetEntry(i);
        if(id_1>0.5 && iso_1<0.15 && pt_1>22 && fabs(eta_1)<2.4) {
        //Here we apply the tag condition. This could be just id and iso, although Alexei applies some other things like high pt
        //and eta cuts, possible dxy and dz cuts and trigger matching
            // Candidate 1 is a valid tag
            if(id_2>0.5 && iso_2<0.15 && pt_2>22 && fabs(eta_2)<2.4) {
            // Candidate 2 is also a valid tag
                //Here we check the pair is opposite sign. We could also check they are well separated, using Delta R(commented out) 
                if(os /*&& DeltaR(eta_1,phi_1,eta_2,phi_2)>0.5*/){
                    if(q_1 == 1) {
                    // Candidate 1 is the tag
                        if(id_2>0.5 && iso_2<0.15) {
                            //Passing probe events fill the passing histogram
                            id_iso_pass->Fill(m_vis);
                            //Separate events using the probe conditions
                            for (int iEta=0; iEta<nEtaBins; ++iEta) {
                                for (int iPt=0; iPt<nPtBins; ++iPt) {
                                    if( fabs(eta_2) > etaBins[iEta] && fabs(eta_2) < etaBins[iEta+1] 
                                        && pt_2 > ptBins[iPt] && pt_2 < ptBins[iPt+1]) {
                                            ZMassEtaPtPass[iEta][iPt]->Fill(m_vis);
                                    }
                                }
                            }
                        } else {
                            //Failing probe events fill the failing histogram
                            id_iso_fail->Fill(m_vis);
                            for (int iEta=0; iEta<nEtaBins-1; ++iEta) {
                                for (int iPt=0; iPt<nPtBins-1; ++iPt) {
                                    if( fabs(eta_2) > etaBins[iEta] && fabs(eta_2) < etaBins[iEta+1] 
                                        && pt_2 > ptBins[iPt] && pt_2 < ptBins[iPt+1]) {
                                            ZMassEtaPtFail[iEta][iPt]->Fill(m_vis);
                                    }
                                }
                            }
                        }
                    } else if (q_1 == -1) {
                    // Candidate 2 is the tag
                        if(id_1>0.5 && iso_1<0.15) {
                            id_iso_pass->Fill(m_vis);
                            for (int iEta=0; iEta<nEtaBins-1; ++iEta) {
                                for (int iPt=0; iPt<nPtBins-1; ++iPt) {
                                    if( fabs(eta_1) > etaBins[iEta] && fabs(eta_1) < etaBins[iEta+1] 
                                        && pt_1 > ptBins[iPt] && pt_1 < ptBins[iPt+1]) {
                                            ZMassEtaPtPass[iEta][iPt]->Fill(m_vis);
                                    }
                                }
                            }
                        } else {
                            id_iso_fail->Fill(m_vis);
                            for (int iEta=0; iEta<nEtaBins-1; ++iEta) {
                                for (int iPt=0; iPt<nPtBins-1; ++iPt) {
                                    if( fabs(eta_1) > etaBins[iEta] && fabs(eta_1) < etaBins[iEta+1] 
                                        && pt_1 > ptBins[iPt] && pt_1 < ptBins[iPt+1]) {
                                            ZMassEtaPtFail[iEta][iPt]->Fill(m_vis);
                                    }
                                }
                            }
                        }
                    }
                }
            } else {
                //Candidate 1 is the tag, check conditions on candidate 2
                if(os /*&& DeltaR(eta_1,phi_1,eta_2,phi_2)>0.5*/){
                    //Only use events where charge == 1 for the tag, even in events with only one tag
                    if(q_1==1){    
                        if(id_2>0.5 && iso_2<0.15) {
                            id_iso_pass->Fill(m_vis);
                            for (int iEta=0; iEta<nEtaBins-1; ++iEta) {
                                for (int iPt=0; iPt<nPtBins-1; ++iPt) {
                                    if( fabs(eta_2) > etaBins[iEta] && fabs(eta_2) < etaBins[iEta+1] 
                                        && pt_2 > ptBins[iPt] && pt_2 < ptBins[iPt+1]) {
                                            ZMassEtaPtPass[iEta][iPt]->Fill(m_vis);
                                    }
                                }
                            }
                        } else {
                            id_iso_fail->Fill(m_vis);
                            for (int iEta=0; iEta<nEtaBins-1; ++iEta) {
                                for (int iPt=0; iPt<nPtBins-1; ++iPt) {
                                    if( fabs(eta_2) > etaBins[iEta] && fabs(eta_2) < etaBins[iEta+1] 
                                        && pt_2 > ptBins[iPt] && pt_2 < ptBins[iPt+1]) {
                                            ZMassEtaPtFail[iEta][iPt]->Fill(m_vis);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    TFile* file0=new TFile("electronIDIso_alexei.root","RECREATE");
    
    for (int iEta=0; iEta<nEtaBins; ++iEta) {
        for (int iPt=0; iPt<nPtBins; ++iPt) {
            ZMassEtaPtPass[iEta][iPt]->Write();
            ZMassEtaPtFail[iEta][iPt]->Write();
        }
    }
    return 0;
}
