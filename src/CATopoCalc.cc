/*
 Calculator for topo variables needing Fat Jets
 
 Author: Joshua Swanson, 2014
 */

#include <iostream>
#include <limits>   // std::numeric_limits

#include "DataFormats/Math/interface/LorentzVector.h"
#include "LJMet/Com/interface/BaseCalc.h"
#include "LJMet/Com/interface/LjmetEventContent.h"
#include "LJMet/Com/interface/LjmetFactory.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

class LjmetFactory;

class CATopoCalc : public BaseCalc {
    
public:
    CATopoCalc();
    virtual ~CATopoCalc();
    virtual int BeginJob();
    virtual int AnalyzeEvent(edm::Event const & event, BaseEventSelector * selector);
    virtual int EndJob();
    
private:
    bool debug_;
    edm::InputTag AK8slimmedJetColl_it;
    
    int FillBranches( std::vector<edm::Ptr<pat::Muon> > const & vTightMuons,
                     std::vector<edm::Ptr<pat::Electron> > const & vTightElectrons,
                     std::vector<std::pair<TLorentzVector,bool> > const & vCorrBtagJets,
                     TLorentzVector const & corrMET,
                     std::vector<TLorentzVector> const & vCAWJets,
                     bool isMuon );
};

static int reg = LjmetFactory::GetInstance()->Register(new CATopoCalc(), "CATopoCalc");

CATopoCalc::CATopoCalc()
{
    mLegend = "[CATopoCalc]: ";
}

CATopoCalc::~CATopoCalc()
{
}

int CATopoCalc::BeginJob()
{
    if (mPset.exists("debug")) debug_ = mPset.getParameter<bool>("debug");
    else debug_ = false;
    
    if (mPset.exists("AK8slimmedJetColl")) AK8slimmedJetColl_it = mPset.getParameter<edm::InputTag>("AK8slimmedJetColl");
    else AK8slimmedJetColl_it = edm::InputTag("slimmedJetsAK8");
    
    return 0;
}

int CATopoCalc::AnalyzeEvent(edm::Event const & event, BaseEventSelector * selector)
{
    //
    // compute event variables here
    //
    
    //
    // _____ Get objects from the selector _____________________
    //
    std::vector<std::pair<TLorentzVector,bool>> const & vCorrBtagJets = selector->GetCorrJetsWithBTags();
    std::vector<edm::Ptr<pat::Muon> > const & vSelMuons = selector->GetSelectedMuons();
    std::vector<edm::Ptr<pat::Electron> > const & vSelElectrons = selector->GetSelectedElectrons();
    TLorentzVector const & corrMET = selector->GetCorrectedMet();
    
    
    bool muonchannel = false;
    bool electronchannel = false;
    if (vSelMuons.size()>0) muonchannel = true;
    if (vSelElectrons.size()>0) electronchannel = true;
    if (muonchannel && electronchannel && debug_)
        std::cout << mLegend << "WARNING: Two Leptons in the event" << std::endl;
    
    if (vSelMuons.size()==0 and vSelElectrons.size()==0) {
        if (debug_)
            std::cout << mLegend << "No Lepton in event! " << std::endl;
        return 0;
    }
    
    if (vCorrBtagJets.size()==0) {
        if (debug_) std::cout << mLegend << "No Corrected jets in event! " << std::endl;
        return 0;
    }
    
    if (corrMET.Pt()==0) {
        if (debug_) std::cout << mLegend << "No Corrected MET in event! " << std::endl;
        return 0;
    }
    
    edm::Handle<std::vector<pat::Jet> > CAWJets;
    event.getByLabel(AK8slimmedJetColl_it, CAWJets);
    std::vector<TLorentzVector> CAWP4;
    TLorentzVector CAJet;
    
    for (std::vector<pat::Jet>::const_iterator ijet = CAWJets->begin(); ijet != CAWJets->end(); ijet++) {
        CAJet.SetPxPyPzE(ijet->px(),
                         ijet->py(),
                         ijet->pz(),
                         ijet->energy());
        CAWP4.push_back(CAJet);
    }
    
    FillBranches(vSelMuons,
                 vSelElectrons,
                 vCorrBtagJets,
                 corrMET,
                 CAWP4,
                 muonchannel // isMuon
                 );
    return 0;
}

int CATopoCalc::FillBranches( std::vector<edm::Ptr<pat::Muon> > const & vSelMuons,
                             std::vector<edm::Ptr<pat::Electron> > const & vSelElectrons,
                             std::vector<std::pair<TLorentzVector,bool>> const & vCorrBtagJets,
                             TLorentzVector const & corrMET,
                             std::vector<TLorentzVector> const & vCAWJets,
                             bool isMuon
                             )
{
    while(1) {
        //Create lepton and met four vectors
        TLorentzVector tlv_lepton;
        
        if ( vSelMuons.size() == 0 && vSelElectrons.size() == 0) break;
        if ( vSelMuons.size() > 0 && vSelElectrons.size() > 0) {
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
        
        if ( corrMET.Pt() > 0 ){ }
        else break;
        
        TLorentzVector tlv_met( corrMET.Px(),
                               corrMET.Py(),
                               corrMET.Pz(),
                               corrMET.Energy() );
        
        std::vector<TLorentzVector> jets;
        std::vector<TLorentzVector> bjets;
        int nJets = 0;
        int nBJets = 0;
        std::vector <double> bJetPt;
        std::vector <double> bJetEta;
        std::vector <double> bJetPhi;
        
        //Remove jets/bjets overlapping with leading CA Jet
        for (vector<std::pair<TLorentzVector,bool>>::const_iterator jet = vCorrBtagJets.begin(); jet != vCorrBtagJets.end(); ++jet){
            
            if( vCAWJets.size() > 0 ){
                double CAtoAKJetDR = vCAWJets[0].DeltaR((*jet).first);
                if( CAtoAKJetDR > 0.65 ){
                    if((*jet).second){
                        bjets.push_back((*jet).first);
                        bJetPt.push_back((*jet).first.Pt());
                        bJetEta.push_back((*jet).first.Eta());
                        bJetPhi.push_back((*jet).first.Phi());
                        
                        ++nBJets;
                    }
                    else{
                        jets.push_back((*jet).first);
                        ++nJets;
                    }
                }
            }
        }
        
        // Calculate topological variables with resulting collections
        double tPrimeMass = -std::numeric_limits<double>::max();
        double minDRCAtoB = std::numeric_limits<double>::max();
        double CAMindrBMass = -std::numeric_limits<double>::max();
        double dR = std::numeric_limits<double>::max();
        TLorentzVector bestTop;
        double topMass = -std::numeric_limits<double>::max();
        double massDiff = std::numeric_limits<double>::max();
        double tPrimeMassBestTop = -std::numeric_limits<double>::max();
        double bestTopMass = -std::numeric_limits<double>::max();
        if( bjets.size() > 0 && vCAWJets.size() > 0) {
            tPrimeMass = double(( tlv_met + tlv_lepton + vCAWJets[0] + bjets[0] ).M());
            
            for (unsigned int i = 0; i < bjets.size(); ++i){
                
                //Find the best l-nu-b top candidate
                topMass = double((tlv_met + tlv_lepton + bjets[i]).M());
                // Orduna: 192.2 into a parameter in the config file instead of a hardcoded value?
                if( fabs(topMass - 192.2) < massDiff ) {
                    massDiff = fabs(topMass - 192.2);
                    bestTop = tlv_met + tlv_lepton + bjets[i];
                    tPrimeMassBestTop = double((bestTop + vCAWJets[0]).M());
                    bestTopMass = topMass;
                }
                
                //Find the bjet nearest to the CA jet but not overlapping
                dR = vCAWJets[0].DeltaR(bjets[i]);
                if( dR < minDRCAtoB ){
                    minDRCAtoB = dR;
                    CAMindrBMass = double((vCAWJets[0] + bjets[i]).M());
                }
            }
        }
        
        // Create branches
        SetValue("tPrimeMass",        tPrimeMass);
        SetValue("tPrimeMassBestTop", tPrimeMassBestTop );
        SetValue("bestTopMasslnub",   bestTopMass);
        SetValue("minDRCAtoB",        minDRCAtoB);
        SetValue("CAMindrBMass",      CAMindrBMass);
        SetValue("nJets",             nJets);
        SetValue("nBJets",            nBJets);
        SetValue("bJetPt" ,           bJetPt);
        SetValue("bJetEta" ,          bJetEta);
        SetValue("bJetPhi" ,          bJetPhi);
        break;
    }
    
    return 0;
}

int CATopoCalc::EndJob()
{
    return 0;
}
