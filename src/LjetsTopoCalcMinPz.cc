/*
 Calculator for the old D0-style LJets topo variables
 
 Author: Gena Kukartsev, 2012
 */



#include <iostream>

#include "DataFormats/Math/interface/LorentzVector.h"
#include "LJMet/Com/interface/BaseCalc.h"
#include "LJMet/Com/interface/LJetsTopoVarsNew.h" // needs work on compatibility and refactoring
#include "LJMet/Com/interface/LjmetEventContent.h"
#include "LJMet/Com/interface/LjmetFactory.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

class LjmetFactory;

class LjetsTopoCalcMinPz : public BaseCalc{
    
public:
    
    LjetsTopoCalcMinPz();
    virtual ~LjetsTopoCalcMinPz(){}
    
    virtual int BeginJob(){
        
        if (mPset.exists("debug"))      debug_ = mPset.getParameter<bool>("debug");
        else                            debug_ = false;
        
        return 0;
    }
    virtual int AnalyzeEvent(edm::Event const & event, BaseEventSelector * selector);
    virtual int EndJob(){return 0;}
    
    
private:
    
    bool debug_;
    
    int FillLjetsBranches( std::vector<edm::Ptr<pat::Muon> > const & vTightMuons,
                          std::vector<edm::Ptr<pat::Electron> > const & vTightElectrons,
                          std::vector<std::pair<TLorentzVector,bool> >  const & vCorrBtagJets,
                          //edm::Ptr<pat::MET> const & pMet,
                          TLorentzVector const & corrMET,
                          bool isMuon,
                          bool bestTop);
    
};



//static int reg = LjmetFactory::GetInstance()->Register(new LjetsTopoCalcMinPz(), "LjetsTopoCalcMinPz");



LjetsTopoCalcMinPz::LjetsTopoCalcMinPz(){
    mLegend = "[LjetsTopoCalcMinPz]: ";
}


int LjetsTopoCalcMinPz::AnalyzeEvent(edm::Event const & event,
                                     BaseEventSelector * selector){
    //
    // compute event variables here
    //
    
    //
    // _____ Get objects from the selector _____________________
    //
    std::vector<edm::Ptr<pat::Jet> >      const & vSelJets = selector->GetSelectedJets();
    std::vector<edm::Ptr<pat::Jet> >      const & vSelBtagJets = selector->GetSelectedBtagJets();
    std::vector<std::pair<TLorentzVector,bool>> const & vCorrBtagJets = selector->GetCorrJetsWithBTags();
    std::vector<edm::Ptr<pat::Jet> >      const & vAllJets = selector->GetAllJets();
    std::vector<edm::Ptr<pat::Muon> >     const & vSelMuons = selector->GetSelectedMuons();
    std::vector<edm::Ptr<pat::Electron> > const & vSelElectrons = selector->GetSelectedElectrons();
    edm::Ptr<pat::MET>                    const & pMet = selector->GetMet();
    TLorentzVector                        const & corrMET = selector->GetCorrectedMet();
    
    //
    //_____ Old D0-style LJetsTopoVarsNew: mu ______
    //
    // (safety checks inside)
    
    bool muonchannel = false;
    bool electronchannel = false;
    if (vSelMuons.size()>0) muonchannel = true;
    if (vSelElectrons.size()>0) electronchannel = true;
    if (muonchannel && electronchannel && debug_)
        std::cout << mLegend <<"WARNING: Two Leptons in the event"<<std::endl;
    
    if (vSelMuons.size()==0 and vSelElectrons.size()==0) {
        if (debug_)
            std::cout << mLegend << "No Lepton in event! "<<std::endl;
        return 0;
    }
    
    if (vCorrBtagJets.size()==0) {
        if (debug_) std::cout << mLegend <<"No Corrected jets in event! "<<std::endl;
        return 0;
    }
    
    if (corrMET.Pt()==0) {
        if (debug_) std::cout << mLegend <<"No Corrected MET in event! "<<std::endl;
        return 0;
    }
    
    
    FillLjetsBranches(vSelMuons,
                      vSelElectrons,
                      vCorrBtagJets,
                      corrMET,
                      muonchannel, // isMuon
                      false); //bestTop for neutrino pz
    
    
    return 0;
}



int LjetsTopoCalcMinPz::FillLjetsBranches( std::vector<edm::Ptr<pat::Muon> > const & vSelMuons,
                                          std::vector<edm::Ptr<pat::Electron> > const & vSelElectrons,
                                          std::vector<std::pair<TLorentzVector,bool>> const & vCorrBtagJets,
                                          //edm::Ptr<pat::MET> const & pMet,
                                          TLorentzVector const & corrMET,
                                          bool isMuon,
                                          bool bestTop){
    //
    // Fills old D0-style LJetsTopoVarsNew
    // Returns 0 if success,
    //        -1 if failed (no needed objects in the event)
    //
    
    int result = -1;
    
    while(1){
        
        TLorentzVector tlv_lepton;
        
        if ( vSelMuons.size() == 0 && vSelElectrons.size() == 0) break;
        if ( vSelMuons.size() > 0 && vSelElectrons.size() > 0) {
            std::cout<<mLegend<<"Two Leptons, using the Muon...."<<std::endl;
            tlv_lepton.SetPxPyPzE( vSelMuons[0]->px(),
                                  vSelMuons[0]->py(),
                                  vSelMuons[0]->pz(),
                                  vSelMuons[0]->energy() );
        }
        if ( vSelMuons.size() > 0 ) {
            tlv_lepton.SetPxPyPzE( vSelMuons[0]->px(),
                                  vSelMuons[0]->py(),
                                  vSelMuons[0]->pz(),
                                  vSelMuons[0]->energy() );
        }
        if ( vSelElectrons.size() > 0 ) {
            tlv_lepton.SetPxPyPzE( vSelElectrons[0]->px(),
                                  vSelElectrons[0]->py(),
                                  vSelElectrons[0]->pz(),
                                  vSelElectrons[0]->energy() );
        }
        
        // make a vector of jets
        /*
         std::vector<TLorentzVector> tlv_jets;
         for ( std::vector<std::pair<TLorentzVector,bool>>::const_iterator iJet = vCorrBtagJets.begin();
         iJet != vCorrBtagJets.end();
         ++iJet ){
         //
         
         //
         TLorentzVector tlv_jet( (*iJet)->first.px(),
         (*iJet)->first.py(),
         (*iJet)->first.pz(),
         (*iJet)->first.energy() );
         
         tlv_jets.push_back(tlv_jet);
         }
         */
        
        // check if MET exists
        //if ( pMet.isNonnull() && pMet.isAvailable() ){ }
        //else break;
        
        // check if corrected MET exists
        if ( corrMET.Pt() > 0 ){ }
        else break;
        
        TLorentzVector tlv_met( corrMET.Px(),
                               corrMET.Py(),
                               corrMET.Pz(),
                               corrMET.Energy() );
        
        // topovars calculator
        //LJetsTopoVarsNew topovars(tlv_jets, tlv_muon, tlv_met, isMuon, bestTop);
        LJetsTopoVarsNew topovars(vCorrBtagJets, tlv_lepton, tlv_met, isMuon,  bestTop);
        
        // compute branches
        SetValue("Jet1Jet2W_M", topovars.Jet1Jet2W_M());
        SetValue("Jet1Jet2W_Pt", topovars.Jet1Jet2W_Pt());
        SetValue("W_MT", topovars.W_MT());
        SetValue("W_M", topovars.W_M());
        SetValue("W_Pt", topovars.W_Pt());
        SetValue("Jet1TagJet2TagW_M", topovars.Jet1TagJet2TagW_M());
        SetValue("AllJetsW_M", topovars.AllJetsW_M());//sqrt_shat
        SetValue("PzNu", topovars.pznu());//sqrt_shat
        SetValue("PzOtherNu", topovars.pzothernu());//sqrt_shat
        
        result = 0;
        break;
        
    } // end of while(1)
    
    return result;
}
