#include "KinematicReconstruction.h"


void testing(){


    TLorentzVector lep_1;
    TLorentzVector lep_2;
    TLorentzVector vecJet1;
    TLorentzVector vecJet2;
    TLorentzVector vecJet3;
    TLorentzVector vecJet4;
    TLorentzVector met;
    vector<TLorentzVector> vecJets;

    int test_var = 1;
    std::cout << "teste " << test_var << std::endl;


    if ( test_var == 0 )
    {
        lep_1.SetPtEtaPhiM( 61.979, 0.61853, -2.63232, 0.105658);
        lep_2.SetPtEtaPhiM( 60.1229, 0.185944, -0.762817, 0.000510999);
        vecJet1.SetPtEtaPhiM( 150.884, -0.316589, 2.74658, -16.0837);
        vecJet2.SetPtEtaPhiM( 41.95, -0.547974, 0.543579, 7.84399);
        double MET_pt = 102.042; double MET_phi = -0.190522; 
        met.SetPxPyPzE ( MET_pt*cos(MET_phi) , MET_pt*sin(MET_phi) , 0., MET_pt);
        vecJets.push_back(vecJet1);
        vecJets.push_back(vecJet2);
        // ttbar_reco: 1
    }
    if ( test_var == 1 )
    {
        lep_1.SetPtEtaPhiM( 109.308, 0.161591, -0.274414, 0.000511);
        lep_2.SetPtEtaPhiM( 31.2923, -0.709839, 1.73462, 0.105658);
        vecJet1.SetPtEtaPhiM( 141.798, 1.42139, 1.30029, 19.8736);
        vecJet2.SetPtEtaPhiM( 119.502, -0.067276, -2.27637, -11.1585);
        vecJet3.SetPtEtaPhiM( 29.545, -0.671753, 0.882935,6.43838);
        vecJet4.SetPtEtaPhiM( 26.8387, -1.9751,-0.34137, 5.77417);
        double MET_pt = 132.666; double MET_phi = -2.41464; 
        met.SetPxPyPzE ( MET_pt*cos(MET_phi) , MET_pt*sin(MET_phi) , 0., MET_pt);
        vecJets.push_back(vecJet1);
        vecJets.push_back(vecJet2);
        vecJets.push_back(vecJet3);
        vecJets.push_back(vecJet4);
        // ttbar_reco: 0
    }
    if ( test_var == 2 )
    {
        lep_1.SetPtEtaPhiM( 77.0179, -0.681641, 0.358215, 0.105658);
        lep_2.SetPtEtaPhiM( 30.2312, -1.34985, -1.52686, 0.000510999);
        vecJet1.SetPtEtaPhiM( 96.799, -1.90283, 1.93408, -11.1372);
        vecJet2.SetPtEtaPhiM( 72.1742, -1.89502, -0.858765, -10.6295);
        double MET_pt = 91.1357; double MET_phi = -2.96994; 
        met.SetPxPyPzE ( MET_pt*cos(MET_phi) , MET_pt*sin(MET_phi) , 0., MET_pt);
        vecJets.push_back(vecJet1);
        vecJets.push_back(vecJet2);
        // ttbar_reco: 1
    }
    if ( test_var == 3 )
    {
        lep_1.SetPtEtaPhiM( 101.426, -0.337646, 2.95215, 0.000510999);
        lep_2.SetPtEtaPhiM( 28.868, -0.634766, 0.333191, 0.000510999);
        vecJet1.SetPtEtaPhiM( 83.6766, -1.30469, -3.08984, -11.833);
        vecJet2.SetPtEtaPhiM( 64.6089, 0.081131, 2.13574, 9.17354);
        double MET_pt = 117.4; double MET_phi = -0.460152; 
        met.SetPxPyPzE ( MET_pt*cos(MET_phi) , MET_pt*sin(MET_phi) , 0., MET_pt);
        vecJets.push_back(vecJet1);
        vecJets.push_back(vecJet2);
        // ttbar_reco: 0
    }
    if ( test_var == 4 )
    {
        lep_1.SetPtEtaPhiM( 55.9799, -0.593994, 0.011961, 0.105658);
        lep_2.SetPtEtaPhiM( 28.5938, -0.984375, -0.0715027, 0.000510999);
        vecJet1.SetPtEtaPhiM( 116.389, -1.57788, -1.92725, -24.6782);
        vecJet2.SetPtEtaPhiM( 38.8176, -0.158997, 1.38745, 5.81478);
        double MET_pt = 101.086; double MET_phi = 2.20604; 
        met.SetPxPyPzE ( MET_pt*cos(MET_phi) , MET_pt*sin(MET_phi) , 0., MET_pt);
        vecJets.push_back(vecJet1);
        vecJets.push_back(vecJet2);
        // ttbar_reco: 0
    }


    //std::cout << " lep_1.Pt() " << lep_1.Pt() << " lep_1.Eta() " << lep_1.Eta() << " lep_1.Phi() " << lep_1.Phi() << " lep_1.M() " << lep_1.M() <<   std::endl;
    //std::cout << " lep_2.Pt() " << lep_2.Pt() << " lep_2.Eta() " << lep_2.Eta() << " lep_2.Phi() " << lep_2.Phi() <<  " lep_2.M() " << lep_2.M() << std::endl;
    //std::cout << " vecJet1.Pt() " << vecJet1.Pt() << " vecJet1.Eta() " << vecJet1.Eta() << " vecJet1.Phi() " << vecJet1.Phi() << " vecJet1.M() " << vecJet1.M() <<  std::endl;
    //std::cout << " vecJet2.Pt() " << vecJet2.Pt() << " vecJet2.Eta() " << vecJet2.Eta() << " vecJet2.Phi() " << vecJet2.Phi() << " vecJet2.M() " << vecJet2.M() <<  std::endl;
    //std::cout << " met.Pt() " << met.Pt() << " met.Eta() " << met.Eta() << " met.Phi() " << met.Phi() << " met.M() " << met.M() <<  std::endl;

    std::cout << " ============================================== " <<  std::endl;

    // loop over 1st jet
    for(std::vector<TLorentzVector>::const_iterator jet1 = vecJets.begin(); jet1 != vecJets.end(); jet1++)
    {
        std::cout << " --------------------------------------------------- " <<  std::endl;
        // loop over 2nd jet
        for(std::vector<TLorentzVector>::const_iterator jet2 = vecJets.begin(); jet2 != vecJets.end(); jet2++)
        {
            // skip same jets
            if(jet1 == jet2) continue;
            
            TLorentzVector jetB, jetBbar;
            if(jet1->M() < 0)
            {
                TLorentzVector jet;
                jet.SetPtEtaPhiM(jet1->Pt(), jet1->Eta(), jet1->Phi(), -1 * jet1->M());
                jetB = jet;
            }
            else jetB = *jet1;

            if(jet2->M() < 0)
            {
                TLorentzVector jet;
                jet.SetPtEtaPhiM(jet2->Pt(), jet2->Eta(), jet2->Phi(), -1 * jet2->M());
                jetBbar = jet;
            }
            else jetBbar = *jet2;

            //std::cout << " lep_1.Pt() " << lep_1.Pt() << " lep_1.Eta() " << lep_1.Eta() << " lep_1.Phi() " << lep_1.Phi() << " lep_1.M() " << lep_1.M() <<   std::endl;
            //std::cout << " lep_2.Pt() " << lep_2.Pt() << " lep_2.Eta() " << lep_2.Eta() << " lep_2.Phi() " << lep_2.Phi() <<  " lep_2.M() " << lep_2.M() << std::endl;
            //std::cout << " jetB.Pt() " << jetB.Pt() << " jetB.Eta() " << jetB.Eta() << " jetB.Phi() " << jetB.Phi() << " jetB.M() " << jetB.M() <<  std::endl;
            //std::cout << " jetBbar.Pt() " << jetBbar.Pt() << " jetBbar.Eta() " << jetBbar.Eta() << " jetBbar.Phi() " << jetBbar.Phi() << " jetBbar.M() " << jetBbar.M() <<  std::endl;
            //std::cout << " met.Pt() " << met.Pt() << " met.Eta() " << met.Eta() << " met.Phi() " << met.Phi() << " met.M() " << met.M() <<  std::endl;

            // get solution
            KinematicReconstruction myTTbarObject;
            myTTbarObject.setConstraints( lep_1, lep_2, jetB, jetBbar, met) ;
            bool hasSolution = myTTbarObject.solutionSmearing( ) ;
            std::cout << "It has Solution " << hasSolution << std::endl;
           
        }
    }
    

    // std::cout << "-------------------------- " << ttbar_reco << std::endl;
    std::cout << "END PROGRAM" << std::endl;
}