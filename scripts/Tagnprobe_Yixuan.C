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

int cut() {
    double DeltaR(double eta_1, double phi_1, double eta_2, double phi_2) {
        double deltaR=sqrt((eta_1-eta_2)*(eta_1-eta_2) + (phi_1-phi_2)*(phi_1-phi_2));
        return deltaR;
    }
    
    //--------------------------------------------------------------------------------------
    //                      1. Extract the root file and the tree
    //--------------------------------------------------------------------------------------
    TFile *file = new TFile("DYJetsToLL_tpzmm_2015.root");
    TTree *tree = dynamic_cast<TTree*>(file->Get("ntuple"));
    
    //--------------------------------------------------------------------------------------
    //                      1.1 Extract the properties of the particles
    //--------------------------------------------------------------------------------------
    double id_1, id_2, iso_1, iso_2; int q_1, q_2; bool trigger_match_1, trigger_match_2;
    tree->SetBranchAddress("id_1",&id_1);
    tree->SetBranchAddress("id_2",&id_2);
    tree->SetBranchAddress("iso_1",&iso_1);
    tree->SetBranchAddress("iso_2",&iso_2);
    tree->SetBranchAddress("q_1",&q_1);
    tree->SetBranchAddress("q_2",&q_2);
    tree->SetBranchAddress("trigger_match_1",&trigger_match_1);
    tree->SetBranchAddress("trigger_match_2",&trigger_match_2);
    
    double pt_1,pt_2,E_1,E_2,eta_1,eta_2,phi_1,phi_2,m_vis,dxy_1,dxy_2,dz_1,dz_2;
    tree->SetBranchAddress("pt_1", &pt_1);
    tree->SetBranchAddress("pt_2", &pt_2);
    tree->SetBranchAddress("E_1", &E_1);
    tree->SetBranchAddress("E_2", &E_2);
    tree->SetBranchAddress("eta_1", &eta_1);
    tree->SetBranchAddress("eta_2", &eta_2);
    tree->SetBranchAddress("phi_1", &phi_1);
    tree->SetBranchAddress("phi_2", &phi_2);
    tree->SetBranchAddress("m_vis", &m_vis);
    tree->SetBranchAddress("dxy_1", &dxy_1);
    tree->SetBranchAddress("dxy_2", &dxy_2);
    tree->SetBranchAddress("dz_1", &dz_1);
    tree->SetBranchAddress("dz_2", &dz_2);
    
    //--------------------------------------------------------------------------------------
    //                      1.1 Create the histograms for storing the properties
    //--------------------------------------------------------------------------------------
    
    TH1F *id_pass=new TH1F("id_pass","",500,0,1000);
    TH1F *id_fail=new TH1F("id_fail","",500,0,1000);
    TH1F *iso_pass=new TH1F("iso_pass","",500,0,1000);
    TH1F *iso_fail=new TH1F("iso_fail","",500,0,1000);
    TH1F *id_iso_pass=new TH1F("id_iso_pass","",800,0,1000);
    TH1F *id_iso_fail=new TH1F("id_iso_fail","",800,0,1000);
    
    TH1F *id_iso_pass_pt1013=new TH1F("id_iso_pass_pt1013","",800,0,1000);
    TH1F *id_iso_fail_pt1013=new TH1F("id_iso_fail_pt1013","",800,0,1000);
    TH1F *id_iso_pass_pt1316=new TH1F("id_iso_pass_pt1316","",800,0,1000);
    TH1F *id_iso_fail_pt1316=new TH1F("id_iso_fail_pt1316","",800,0,1000);
    TH1F *id_iso_pass_pt1620=new TH1F("id_iso_pass_pt1620","",800,0,1000);
    TH1F *id_iso_fail_pt1620=new TH1F("id_iso_fail_pt1620","",800,0,1000);
    TH1F *id_iso_pass_pt2025=new TH1F("id_iso_pass_pt2025","",800,0,1000);
    TH1F *id_iso_fail_pt2025=new TH1F("id_iso_fail_pt2025","",800,0,1000);
    TH1F *id_iso_pass_pt2530=new TH1F("id_iso_pass_pt2530","",800,0,1000);
    TH1F *id_iso_fail_pt2530=new TH1F("id_iso_fail_pt2530","",800,0,1000);
    TH1F *id_iso_pass_pt3040=new TH1F("id_iso_pass_pt3040","",800,0,1000);
    TH1F *id_iso_fail_pt3040=new TH1F("id_iso_fail_pt3040","",800,0,1000);
    TH1F *id_iso_pass_pt4060=new TH1F("id_iso_pass_pt4060","",800,0,1000);
    TH1F *id_iso_fail_pt4060=new TH1F("id_iso_fail_pt4060","",800,0,1000);
    
    // keep track of the statistic errors on the histogram, before filling
    //m_vis_pass->Sumw2();
    //m_vis_fail->Sumw2();
    
    //--------------------------------------------------------------------------------------
    //                      1.3 Applying selection criteria - Tag/probe
    //--------------------------------------------------------------------------------------
    // fill the histograms with the entries from the leaves, inside the .Fill() is the equation
    /*for (int i = 0; i < tree->GetEntries(); ++i) {
        tree->GetEntry(i);
                if (id_1>0.5 && iso_1<0.1 && q_1==1 && q_2==-1 && pt_1>10. && pt_1<15.){
                    if (id_2>0.5){
                        id_pass->Fill(sqrt((E_1+E_2)*(E_1+E_2)-(pt_1*cos(phi_1)+pt_2*cos(phi_2))*(pt_1*cos(phi_1)+pt_2*cos(phi_2))-(pt_1*sin(phi_1)+pt_2*sin(phi_2))*(pt_1*sin(phi_1)+pt_2*sin(phi_2))-(pt_1*sinh(eta_1)+pt_2*sinh(eta_2))*(pt_1*sinh(eta_1)+pt_2*sinh(eta_1))));
                    }
                    else {
                        id_fail->Fill(sqrt((E_1+E_2)*(E_1+E_2)-(pt_1*cos(phi_1)+pt_2*cos(phi_2))*(pt_1*cos(phi_1)+pt_2*cos(phi_2))-(pt_1*sin(phi_1)+pt_2*sin(phi_2))*(pt_1*sin(phi_1)+pt_2*sin(phi_2))-(pt_1*sinh(eta_1)+pt_2*sinh(eta_2))*(pt_1*sinh(eta_1)+pt_2*sinh(eta_1))));
                    }
                }
                else if (id_2>0.5 && iso_2<0.1 && q_2==1 && q_1==-1 && pt_2>10. && pt_2<15.){
           
                    if (id_1>0.5){
                        id_pass->Fill(sqrt((E_1+E_2)*(E_1+E_2)-(pt_1*cos(phi_1)+pt_2*cos(phi_2))*(pt_1*cos(phi_1)+pt_2*cos(phi_2))-(pt_1*sin(phi_1)+pt_2*sin(phi_2))*(pt_1*sin(phi_1)+pt_2*sin(phi_2))-(pt_1*sinh(eta_1)+pt_2*sinh(eta_2))*(pt_1*sinh(eta_1)+pt_2*sinh(eta_1))));
                    }
                    else {
                        id_fail->Fill(sqrt((E_1+E_2)*(E_1+E_2)-(pt_1*cos(phi_1)+pt_2*cos(phi_2))*(pt_1*cos(phi_1)+pt_2*cos(phi_2))-(pt_1*sin(phi_1)+pt_2*sin(phi_2))*(pt_1*sin(phi_1)+pt_2*sin(phi_2))-(pt_1*sinh(eta_1)+pt_2*sinh(eta_2))*(pt_1*sinh(eta_1)+pt_2*sinh(eta_1))));
                    }
                }
    }
    
    // fill the histograms with the entries from the leaves, inside the .Fill() is the equation
    for (int i = 0; i < tree->GetEntries(); ++i) {
        tree->GetEntry(i);
        if (pt_1>15. && pt_1<20. && pt_2>15. && pt_2<20.){
            if (id_1>0.5 && iso_1<0.1 && q_1==1 && q_2==-1){
                if (iso_2<0.1){
                    iso_pass->Fill(sqrt((E_1+E_2)*(E_1+E_2)-(pt_1*cos(phi_1)+pt_2*cos(phi_2))*(pt_1*cos(phi_1)+pt_2*cos(phi_2))-(pt_1*sin(phi_1)+pt_2*sin(phi_2))*(pt_1*sin(phi_1)+pt_2*sin(phi_2))-(pt_1*sinh(eta_1)+pt_2*sinh(eta_2))*(pt_1*sinh(eta_1)+pt_2*sinh(eta_1))));
                }
                else {
                    iso_fail->Fill(sqrt((E_1+E_2)*(E_1+E_2)-(pt_1*cos(phi_1)+pt_2*cos(phi_2))*(pt_1*cos(phi_1)+pt_2*cos(phi_2))-(pt_1*sin(phi_1)+pt_2*sin(phi_2))*(pt_1*sin(phi_1)+pt_2*sin(phi_2))-(pt_1*sinh(eta_1)+pt_2*sinh(eta_2))*(pt_1*sinh(eta_1)+pt_2*sinh(eta_1))));
                }
            }
            else if (id_2>0.5 && iso_2<0.1 && q_2==1 && q_1==-1){

                if (iso_1<0.1){
                    iso_pass->Fill(sqrt((E_1+E_2)*(E_1+E_2)-(pt_1*cos(phi_1)+pt_2*cos(phi_2))*(pt_1*cos(phi_1)+pt_2*cos(phi_2))-(pt_1*sin(phi_1)+pt_2*sin(phi_2))*(pt_1*sin(phi_1)+pt_2*sin(phi_2))-(pt_1*sinh(eta_1)+pt_2*sinh(eta_2))*(pt_1*sinh(eta_1)+pt_2*sinh(eta_1))));
                }
                else {
                    iso_fail->Fill(sqrt((E_1+E_2)*(E_1+E_2)-(pt_1*cos(phi_1)+pt_2*cos(phi_2))*(pt_1*cos(phi_1)+pt_2*cos(phi_2))-(pt_1*sin(phi_1)+pt_2*sin(phi_2))*(pt_1*sin(phi_1)+pt_2*sin(phi_2))-(pt_1*sinh(eta_1)+pt_2*sinh(eta_2))*(pt_1*sinh(eta_1)+pt_2*sinh(eta_1))));
                }
            }
        }
    }*/
    for (int i = 0; i < tree->GetEntries(); ++i) {
        tree->GetEntry(i);
        // each i is a pair of muon and contain information like pt E eta phi...
        //if(os && DeltaR(eta_1,phi_1,eta_2,phi_2)>0.5){
        if (id_1>0.5 && iso_1<0.15 && q_1==1 && q_2==-1 && pt_1>22 && fabs(eta_1)<2.1 && trigger_match_1 && fabs(dxy_1) < 0.045 && fabs(dz_1) < 0.2){
        //first check the tag in particle 1 and requiring particle 2 to have a minus charge
            if (fabs(dxy_2) < 0.2 && fabs(dz_2) < 0.5){
            //now particle 1 is tag, 2 is probe, then select the pt range and eta range on the probe only.
                if (id_2>0.5 && iso_2<0.15){
            //check if probe passes the id+iso the fill the m_vis
                    if (pt_2>10. && pt_2<13. && fabs(eta_2)>1.2 && fabs(eta_2)<2.4){
                        id_iso_pass_pt1013->Fill(m_vis);
                    }
                    if (pt_2>13. && pt_2<16. && fabs(eta_2)>1.2 && fabs(eta_2)<2.4){
                        id_iso_pass_pt1316->Fill(m_vis);
                    }
                    if (pt_2>16. && pt_2<20. && fabs(eta_2)>1.2 && fabs(eta_2)<2.4){
                        id_iso_pass_pt1620->Fill(m_vis);
                    }
                    if (pt_2>20. && pt_2<25. && fabs(eta_2)>1.2 && fabs(eta_2)<2.4){
                        id_iso_pass_pt2025->Fill(m_vis);
                    }
                    if (pt_2>25. && pt_2<30. && fabs(eta_2)>1.2 && fabs(eta_2)<2.4){
                        id_iso_pass_pt2530->Fill(m_vis);
                    }
                    if (pt_2>30. && pt_2<40. && fabs(eta_2)>1.2 && fabs(eta_2)<2.4){
                        id_iso_pass_pt3040->Fill(m_vis);
                    }
                    if (pt_2>40. && pt_2<60. && fabs(eta_2)>1.2 && fabs(eta_2)<2.4){
                        id_iso_pass_pt4060->Fill(m_vis);
                    }
                }
                else {
                    if (pt_2>10. && pt_2<13. && fabs(eta_2)>1.2 && fabs(eta_2)<2.4){
                        id_iso_fail_pt1013->Fill(m_vis);
                    }
                    if (pt_2>13. && pt_2<16. && fabs(eta_2)>1.2 && fabs(eta_2)<2.4){
                        id_iso_fail_pt1316->Fill(m_vis);
                    }
                    if (pt_2>16. && pt_2<20. && fabs(eta_2)>1.2 && fabs(eta_2)<2.4){
                        id_iso_fail_pt1620->Fill(m_vis);
                    }
                    if (pt_2>20. && pt_2<25. && fabs(eta_2)>1.2 && fabs(eta_2)<2.4){
                        id_iso_fail_pt2025->Fill(m_vis);
                    }
                    if (pt_2>25. && pt_2<30. && fabs(eta_2)>1.2 && fabs(eta_2)<2.4){
                        id_iso_fail_pt2530->Fill(m_vis);
                    }
                    if (pt_2>30. && pt_2<40. && fabs(eta_2)>1.2 && fabs(eta_2)<2.4){
                        id_iso_fail_pt3040->Fill(m_vis);
                    }
                    if (pt_2>40. && pt_2<60. && fabs(eta_2)>1.2 && fabs(eta_2)<2.4){
                        id_iso_fail_pt4060->Fill(m_vis);
                    }
                }
            }
        }
        else if (id_2>0.5 && iso_2<0.15 && q_2==1 && q_1==-1 && pt_2>22 && fabs(eta_2)<2.1 && trigger_match_2 && fabs(dxy_2) < 0.045 && fabs(dz_2) < 0.2){
        // similarly letting particle 2 to be the tag and 1 to be the probe
            if (fabs(dxy_2) < 0.2 && fabs(dz_2) < 0.5){
                if (id_1>0.5 && iso_1<0.15){
                //check if probe passes the id+iso the fill the m_vis
                    if (pt_1>10. && pt_1<13. && fabs(eta_1)>1.2 && fabs(eta_1)<2.4){
                        id_iso_pass_pt1013->Fill(m_vis);
                    }
                    if (pt_1>13. && pt_1<16. && fabs(eta_1)>1.2 && fabs(eta_1)<2.4){
                        id_iso_pass_pt1316->Fill(m_vis);
                    }
                    if (pt_1>16. && pt_1<20. && fabs(eta_1)>1.2 && fabs(eta_1)<2.4){
                        id_iso_pass_pt1620->Fill(m_vis);
                    }
                    if (pt_1>20. && pt_1<25. && fabs(eta_1)>1.2 && fabs(eta_1)<2.4){
                        id_iso_pass_pt2025->Fill(m_vis);
                    }
                    if (pt_1>25. && pt_1<30. && fabs(eta_1)>1.2 && fabs(eta_1)<2.4){
                        id_iso_pass_pt2530->Fill(m_vis);
                    }
                    if (pt_1>30. && pt_1<40. && fabs(eta_1)>1.2 && fabs(eta_1)<2.4){
                        id_iso_pass_pt3040->Fill(m_vis);
                    }
                    if (pt_1>40. && pt_1<60. && fabs(eta_1)>1.2 && fabs(eta_1)<2.4){
                        id_iso_pass_pt4060->Fill(m_vis);
                    }
                }
                else {
                    if (pt_1>10. && pt_1<13. && fabs(eta_1)>1.2 && fabs(eta_1)<2.4){
                        id_iso_fail_pt1013->Fill(m_vis);
                    }
                    if (pt_1>13. && pt_1<16. && fabs(eta_1)>1.2 && fabs(eta_1)<2.4){
                        id_iso_fail_pt1316->Fill(m_vis);
                    }
                    if (pt_1>16. && pt_1<20. && fabs(eta_1)>1.2 && fabs(eta_1)<2.4){
                        id_iso_fail_pt1620->Fill(m_vis);
                    }
                    if (pt_1>20. && pt_1<25. && fabs(eta_1)>1.2 && fabs(eta_1)<2.4){
                        id_iso_fail_pt2025->Fill(m_vis);
                    }
                    if (pt_1>25. && pt_1<30. && fabs(eta_1)>1.2 && fabs(eta_1)<2.4){
                        id_iso_fail_pt2530->Fill(m_vis);
                    }
                    if (pt_1>30. && pt_1<40. && fabs(eta_1)>1.2 && fabs(eta_1)<2.4){
                        id_iso_fail_pt3040->Fill(m_vis);
                    }
                    if (pt_1>40. && pt_1<60. && fabs(eta_1)>1.2 && fabs(eta_1)<2.4){
                        id_iso_fail_pt4060->Fill(m_vis);
                    }
                }
            }
        }
    }
    
    TFile* file0=new TFile("ZMassEta1p2to2p4Ptcut.root","RECREATE");
    
    /*id_pass->Write();
    id_fail->Write();
    iso_pass->Write();
    iso_fail->Write();*/
    id_iso_pass_pt1013->Write();
    id_iso_fail_pt1013->Write();
    id_iso_pass_pt1316->Write();
    id_iso_fail_pt1316->Write();
    id_iso_pass_pt1620->Write();
    id_iso_fail_pt1620->Write();
    id_iso_pass_pt2025->Write();
    id_iso_fail_pt2025->Write();
    id_iso_pass_pt2530->Write();
    id_iso_fail_pt2530->Write();
    id_iso_pass_pt3040->Write();
    id_iso_fail_pt3040->Write();
    id_iso_pass_pt4060->Write();
    id_iso_fail_pt4060->Write();
    
    return 0;
}
