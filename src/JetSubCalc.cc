/*
 Calculator for substructure variables
 
 Author: Joshua Swanson, 2014
 */

#include <iostream>
#include <limits>   // std::numeric_limits
#include <vector>
#include <string>

#include "LJMet/Com/interface/BaseCalc.h"
#include "LJMet/Com/interface/LjmetFactory.h"
#include "LJMet/Com/interface/LjmetEventContent.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/PATObject.h"
#include "DataFormats/JetReco/interface/CATopJetTagInfo.h"

class LjmetFactory;

class JetSubCalc : public BaseCalc {
    
public:
    JetSubCalc();
    virtual ~JetSubCalc();
    virtual int BeginJob();
    virtual int AnalyzeEvent(edm::Event const & event, BaseEventSelector * selector);
    virtual int EndJob();
    
private:
    edm::InputTag slimmedJetColl_it;
    edm::InputTag slimmedJetsAK8Coll_it;
    std::string bDiscriminant;
    std::string tagInfo;
};

static int reg = LjmetFactory::GetInstance()->Register(new JetSubCalc(), "JetSubCalc");

JetSubCalc::JetSubCalc()
{
}

JetSubCalc::~JetSubCalc()
{
}

int JetSubCalc::BeginJob()
{
    if (mPset.exists("slimmedJetColl")) slimmedJetColl_it = mPset.getParameter<edm::InputTag>("slimmedJetColl");
    else slimmedJetColl_it = edm::InputTag("slimmedJets");
    
    if (mPset.exists("slimmedJetsAK8Coll")) slimmedJetsAK8Coll_it = mPset.getParameter<edm::InputTag>("slimmedJetsAK8Coll");
    else slimmedJetsAK8Coll_it = edm::InputTag("slimmedJetsAK8");
    
    if (mPset.exists("bDiscriminant")) bDiscriminant = mPset.getParameter<std::string>("bDiscriminant");
    else bDiscriminant = "pfCombinedSecondaryVertexBJetTags";
    
    if (mPset.exists("tagInfo")) tagInfo = mPset.getParameter<std::string>("tagInfo");
    else tagInfo = "caTop";
    
    return 0;
}

int JetSubCalc::AnalyzeEvent(edm::Event const & event, BaseEventSelector * selector)
{
    float subjetCSV;
    int CSVL, CSVM, CSVT;

    // I think these are AK4
    edm::Handle<std::vector<pat::Jet> > theJets;
    event.getByLabel(slimmedJetColl_it, theJets);
    
    // Available variables
    std::vector<double> theJetPt;
    std::vector<double> theJetEta;
    std::vector<double> theJetPhi;
    std::vector<double> theJetEnergy;
    std::vector<double> theJetCSV;
    
    // Additional variables related to the associated secondary vertex if there is one
    // Mass of the vertex
    std::vector<double> theJetVtxMass;
    // Number of tracks
    std::vector<double> theJetVtxNtracks;
    // Decay length value and significance
    std::vector<double> theJetVtx3DVal;
    std::vector<double> theJetVtx3DSig;
    
    // Discriminator for the MVA PileUp id.
    // NOTE: Training used is for ak5PFJetsCHS in CMSSW 5.3.X and Run 1 pileup
    std::vector<double> theJetPileupJetId;
    
    //Identity
    std::vector<int> theJetIndex;
    std::vector<int> theJetnDaughters;
    
    //Daughter four vector and index
    std::vector<double> theJetDaughterPt;
    std::vector<double> theJetDaughterEta;
    std::vector<double> theJetDaughterPhi;
    std::vector<double> theJetDaughterEnergy;
    
    std::vector<int> theJetDaughterMotherIndex;
    
    std::vector<int> theJetCSVLSubJets;
    std::vector<int> theJetCSVMSubJets;
    std::vector<int> theJetCSVTSubJets;
    
    double theVtxMass, theVtxNtracks, theVtx3DVal, theVtx3DSig, thePileupJetId;
    
    for (std::vector<pat::Jet>::const_iterator ijet = theJets->begin(); ijet != theJets->end(); ijet++) {
        int index = (int)(ijet-theJets->begin());
        
        theVtxNtracks  = -std::numeric_limits<double>::max();
        theVtxNtracks  = (double)ijet->userFloat("vtxNtracks");
        
        if (theVtxNtracks > 0) {
            theVtxMass     = -std::numeric_limits<double>::max();
            theVtx3DVal    = -std::numeric_limits<double>::max();
            theVtx3DSig    = -std::numeric_limits<double>::max();

            theVtxMass     = (double)ijet->userFloat("vtxMass");
            theVtx3DVal    = (double)ijet->userFloat("vtx3DVal");
            theVtx3DSig    = (double)ijet->userFloat("vtx3DSig");
            
            theJetVtxNtracks.push_back(theVtxNtracks);
            theJetVtxMass.push_back(theVtxMass);
            theJetVtx3DVal.push_back(theVtx3DVal);
            theJetVtx3DSig.push_back(theVtx3DSig);
        }
        
        thePileupJetId = -std::numeric_limits<double>::max();
        thePileupJetId = (double)ijet->userFloat("pileupJetId:fullDiscriminant");
        theJetPileupJetId.push_back(thePileupJetId);
        
        theJetPt    .push_back(ijet->pt());
        theJetEta   .push_back(ijet->eta());
        theJetPhi   .push_back(ijet->phi());
        theJetEnergy.push_back(ijet->energy());
        theJetCSV.push_back(ijet->bDiscriminator( bDiscriminant ));
        
        theJetIndex.push_back(index);
        theJetnDaughters.push_back((int)ijet->numberOfDaughters());
        
        CSVL = 0;
        CSVM = 0;
        CSVT = 0;
        subjetCSV = -std::numeric_limits<float>::max();
        
        for (size_t ui = 0; ui < ijet->numberOfDaughters(); ui++) {
            pat::PackedCandidate const * theDaughter = dynamic_cast<pat::PackedCandidate const *>(ijet->daughter(ui));
            
            theJetDaughterPt    .push_back(theDaughter->pt());
            theJetDaughterEta   .push_back(theDaughter->eta());
            theJetDaughterPhi   .push_back(theDaughter->phi());
            theJetDaughterEnergy.push_back(theDaughter->energy());
            
            theJetDaughterMotherIndex.push_back(index);
            
            
//            pat::Jet const * subjet = dynamic_cast<pat::Jet const *>(ijet->daughter(ui));
//            subjetCSV = subjet->bDiscriminator(bDiscriminant);
//            if (subjetCSV > 0.244 && subjet->pt() > 20.) {
//                CSVL++;
//            }
//            if (subjetCSV > 0.679 && subjet->pt() > 20.) {
//                CSVM++;
//            }
//            if (subjetCSV > 0.898 && subjet->pt() > 20.) {
//                CSVT++;
//            }
        }
        theJetCSVLSubJets.push_back(CSVL);
        theJetCSVMSubJets.push_back(CSVM);
        theJetCSVTSubJets.push_back(CSVT);
    }
    
    SetValue("theJetPt",     theJetPt);
    SetValue("theJetEta",    theJetEta);
    SetValue("theJetPhi",    theJetPhi);
    SetValue("theJetEnergy", theJetEnergy);
    SetValue("theJetCSV",    theJetCSV);
    
    SetValue("theJetVtxMass",     theJetVtxMass);
    SetValue("theJetVtxNtracks",  theJetVtxNtracks);
    SetValue("theJetVtx3DVal",    theJetVtx3DVal);
    SetValue("theJetVtx3DSig",    theJetVtx3DSig);
    SetValue("theJetPileupJetId", theJetPileupJetId);
    
    SetValue("theJetIndex",      theJetIndex);
    SetValue("theJetnDaughters", theJetnDaughters);
    
    SetValue("theJetDaughterPt",     theJetDaughterPt);
    SetValue("theJetDaughterEta",    theJetDaughterEta);
    SetValue("theJetDaughterPhi",    theJetDaughterPhi);
    SetValue("theJetDaughterEnergy", theJetDaughterEnergy);
    
    SetValue("theJetDaughterMotherIndex", theJetDaughterMotherIndex);
    
    SetValue("theJetCSVLSubJets", theJetCSVLSubJets);
    SetValue("theJetCSVMSubJets", theJetCSVMSubJets);
    SetValue("theJetCSVTSubJets", theJetCSVTSubJets);
    
    // I think these are AK8 jets so topMass, minMass and nSubJets make sense
    edm::Handle<std::vector<pat::Jet> > theAK8Jets;
    event.getByLabel(slimmedJetsAK8Coll_it, theAK8Jets);
    
    // Four vector
    std::vector<double> theJetAK8Pt;
    std::vector<double> theJetAK8Eta;
    std::vector<double> theJetAK8Phi;
    std::vector<double> theJetAK8Energy;
    std::vector<double> theJetAK8CSV;
    
    // Pruned, trimmed and filtered masses available
    std::vector<double> theJetAK8PrunedMass;
    std::vector<double> theJetAK8TrimmedMass;
    std::vector<double> theJetAK8FilteredMass;
    
    // n-subjettiness variables tau1, tau2, and tau3 available
    std::vector<double> theJetAK8NjettinessTau1;
    std::vector<double> theJetAK8NjettinessTau2;
    std::vector<double> theJetAK8NjettinessTau3;
    
    std::vector<double> theJetAK8caTopTopMass;
    std::vector<double> theJetAK8caTopMinMass;
    std::vector<int> theJetAK8caTopnSubJets;

    std::vector<double> theJetAK8Mass;
    std::vector<int>    theJetAK8Index;
    std::vector<int>    theJetAK8nDaughters;
    
    // Daughter four vector and index
    std::vector<double> theJetAK8DaughterPt;
    std::vector<double> theJetAK8DaughterEta;
    std::vector<double> theJetAK8DaughterPhi;
    std::vector<double> theJetAK8DaughterEnergy;
    
    std::vector<int> theJetAK8DaughterMotherIndex;
    
    std::vector<int> theJetAK8CSVLSubJets;
    std::vector<int> theJetAK8CSVMSubJets;
    std::vector<int> theJetAK8CSVTSubJets;
    
    double topMass, minMass;
    int nSubJets;
    double thePrunedMass, theTrimmedMass, theFilteredMass;
    double theNjettinessTau1, theNjettinessTau2, theNjettinessTau3;
    
    for (std::vector<pat::Jet>::const_iterator ijet = theAK8Jets->begin(); ijet != theAK8Jets->end(); ijet++) {
        int index = (int)(ijet-theAK8Jets->begin());
        
        thePrunedMass   = -std::numeric_limits<double>::max();
        theTrimmedMass  = -std::numeric_limits<double>::max();
        theFilteredMass = -std::numeric_limits<double>::max();
        thePrunedMass   = (double)ijet->userFloat("ak8PFJetsCHSPrunedLinks");
        theTrimmedMass  = (double)ijet->userFloat("ak8PFJetsCHSTrimmedLinks");
        theFilteredMass = (double)ijet->userFloat("ak8PFJetsCHSFilteredLinks");
        
        theNjettinessTau1 = -std::numeric_limits<double>::max();
        theNjettinessTau2 = -std::numeric_limits<double>::max();
        theNjettinessTau3 = -std::numeric_limits<double>::max();
        theNjettinessTau1 = (double)ijet->userFloat("NjettinessAK8:tau1");
        theNjettinessTau2 = (double)ijet->userFloat("NjettinessAK8:tau2");
        theNjettinessTau3 = (double)ijet->userFloat("NjettinessAK8:tau3");

        theJetAK8Pt    .push_back(ijet->pt());
        theJetAK8Eta   .push_back(ijet->eta());
        theJetAK8Phi   .push_back(ijet->phi());
        theJetAK8Energy.push_back(ijet->energy());
        theJetAK8CSV.push_back(ijet->bDiscriminator( bDiscriminant ));
        
        theJetAK8PrunedMass  .push_back(thePrunedMass);
        theJetAK8TrimmedMass .push_back(theTrimmedMass);
        theJetAK8FilteredMass.push_back(theFilteredMass);
        
        theJetAK8NjettinessTau1.push_back(theNjettinessTau1);
        theJetAK8NjettinessTau2.push_back(theNjettinessTau2);
        theJetAK8NjettinessTau3.push_back(theNjettinessTau3);
        
        theJetAK8Mass.push_back(ijet->mass());
        theJetAK8nDaughters.push_back((int)ijet->numberOfDaughters());
        
        theJetAK8Index.push_back(index);
        
        reco::CATopJetTagInfo const * jetInfo = dynamic_cast<reco::CATopJetTagInfo const *>( ijet->tagInfo( tagInfo ));

        if ( jetInfo != 0 ) {
            topMass   = -std::numeric_limits<double>::max();
            minMass   = -std::numeric_limits<double>::max();
            nSubJets  = std::numeric_limits<int>::min();
            
            topMass = jetInfo->properties().topMass;
            minMass = jetInfo->properties().minMass;
            nSubJets = jetInfo->properties().nSubJets;
            
            theJetAK8caTopTopMass.push_back(topMass);
            theJetAK8caTopMinMass.push_back(minMass);
            theJetAK8caTopnSubJets.push_back(nSubJets);
        }
        
        subjetCSV = -std::numeric_limits<float>::max();
        CSVL = 0;
        CSVM = 0;
        CSVT = 0;
        
        for (size_t ui = 0; ui < ijet->numberOfDaughters(); ui++) {
            pat::PackedCandidate const * theDaughter = dynamic_cast<pat::PackedCandidate const *>(ijet->daughter(ui));
            theJetAK8DaughterPt    .push_back(theDaughter->pt());
            theJetAK8DaughterEta   .push_back(theDaughter->eta());
            theJetAK8DaughterPhi   .push_back(theDaughter->phi());
            theJetAK8DaughterEnergy.push_back(theDaughter->energy());
            
            theJetAK8DaughterMotherIndex.push_back(index);
            
//            subjetCSV = subjet->bDiscriminator(bDiscriminant);
//            if (subjetCSV > 0.244 && ijet->daughter(ui)->pt() > 20.) {
//                CSVL++;
//            }
//            if (subjetCSV > 0.679 && ijet->daughter(ui)->pt() > 20.){
//                CSVM++;
//            }
//            if (subjetCSV > 0.898 && ijet->daughter(ui)->pt() > 20.){
//                CSVT++;
//            }
        }
        theJetAK8CSVLSubJets.push_back(CSVL);
        theJetAK8CSVMSubJets.push_back(CSVM);
        theJetAK8CSVTSubJets.push_back(CSVT);
    }
    
    SetValue("theJetAK8Pt",     theJetAK8Pt);
    SetValue("theJetAK8Eta",    theJetAK8Eta);
    SetValue("theJetAK8Phi",    theJetAK8Phi);
    SetValue("theJetAK8Energy", theJetAK8Energy);
    SetValue("theJetAK8CSV",    theJetAK8CSV);
    
    SetValue("theJetAK8PrunedMass",   theJetAK8PrunedMass);
    SetValue("theJetAK8TrimmedMass",  theJetAK8TrimmedMass);
    SetValue("theJetAK8FilteredMass", theJetAK8FilteredMass);
    
    SetValue("theJetAK8NjettinessTau1", theJetAK8NjettinessTau1);
    SetValue("theJetAK8NjettinessTau2", theJetAK8NjettinessTau2);
    SetValue("theJetAK8NjettinessTau3", theJetAK8NjettinessTau3);
    
    SetValue("theJetAK8Mass",   theJetAK8Mass);
    
    SetValue("theJetAK8Index",      theJetAK8Index);
    SetValue("theJetAK8nDaughters", theJetAK8nDaughters);
    
    SetValue("theJetAK8caTopTopMass", theJetAK8caTopTopMass);
    SetValue("theJetAK8caTopMinMass", theJetAK8caTopMinMass);
    SetValue("theJetAK8caTopnSubJets", theJetAK8caTopnSubJets);
    
    SetValue("theJetAK8DaughterPt",     theJetAK8DaughterPt);
    SetValue("theJetAK8DaughterEta",    theJetAK8DaughterEta);
    SetValue("theJetAK8DaughterPhi",    theJetAK8DaughterPhi);
    SetValue("theJetAK8DaughterEnergy", theJetAK8DaughterEnergy);
    
    SetValue("theJetAK8DaughterMotherIndex", theJetAK8DaughterMotherIndex);
    
    SetValue("theJetAK8CSVLSubJets", theJetAK8CSVLSubJets);
    SetValue("theJetAK8CSVMSubJets", theJetAK8CSVMSubJets);
    SetValue("theJetAK8CSVTSubJets", theJetAK8CSVTSubJets);
    
    return 0;
}

int JetSubCalc::EndJob()
{
    return 0;
}
