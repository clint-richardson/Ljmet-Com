#include "LJMet/Com/interface/LjmetEventContent.h"

LjmetEventContent::LjmetEventContent():
mName("LjmetEventContent"),
mLegend("[LjmetEventContent]: "),
mpTree(0),
mFirstEntry(true),
mVerbosity(0)
{
}

LjmetEventContent::LjmetEventContent(std::map<std::string, edm::ParameterSet const> mPar):
mName("LjmetEventContent"),
mLegend("[LjmetEventContent]: "),
mpTree(0),
mFirstEntry(true),
mVerbosity(0)
{
    if (mPar.find("ljmet") != mPar.end()) {
        if (mPar["ljmet"].exists("verbosity")) {
            mVerbosity = mPar["ljmet"].getParameter<int>("verbosity");
        }
    }
}

LjmetEventContent::~LjmetEventContent()
{
}

void LjmetEventContent::SetTree(TTree * tree)
{
    mpTree = tree;
}

void LjmetEventContent::SetHistogram(std::string modname, std::string histname, int nbins, double low, double high)
{
    // Create histogram entry in event content, so it is created by the LjmetFactory
    mLegend = "[" + mName + "]: " ;
    
    //mDoubleHist[modname][histname]
    
    std::map<std::string, std::map<std::string, HistMetadata>>::iterator iMod;
    std::map<std::string, HistMetadata>::iterator iHist;
    
    iMod = mDoubleHist.lower_bound(modname);
    
    if (iMod != mDoubleHist.end() && !(mDoubleHist.key_comp()(modname, iMod->first))) {
        // Key already exists. Update iMod->second
        iHist = iMod->second.lower_bound(histname);
        if(iHist != iMod->second.end() && !(iMod->second.key_comp()(histname, iHist->first))) {
            // Key already exists. Update iHist->second - not needed
            std::cout << mLegend << "Histogram " << modname << "/" << histname << " is already set" << std::endl;
        } else {
            // The key does not exist in the map. Add it.
            iMod->second.insert(iHist, std::pair<std::string, HistMetadata>(histname, HistMetadata(histname, nbins, low, high)));
        }
    } else {
        // The key does not exist in the map. Add it.
        std::map<std::string,HistMetadata> _map;
        _map.insert(std::pair<std::string, HistMetadata>(histname, HistMetadata(histname, nbins, low, high)));
        mDoubleHist.insert(iMod, std::pair<std::string, std::map<std::string, HistMetadata>>(modname, _map));
        // Use iMod as a hint to insert, so it can avoid another lookup
    }
}

void LjmetEventContent::SetValue(std::string key, bool value)
{
    mBoolBranch[key] = value;
}

void LjmetEventContent::SetValue(std::string key, int value)
{
    mIntBranch[key] = value;
}

void LjmetEventContent::SetValue(std::string key, long long value)
{
  mLongIntBranch[key] = value;
}
 
void LjmetEventContent::SetValue(std::string key, double value)
{
    mDoubleBranch[key] = value;
}

void LjmetEventContent::SetValue(std::string key, std::vector<bool> value)
{
    mVectorBoolBranch[key] = value;
}

void LjmetEventContent::SetValue(std::string key, std::vector<int> value)
{
    mVectorIntBranch[key] = value;
}

void LjmetEventContent::SetValue(std::string key, std::vector<double> value)
{
    mVectorDoubleBranch[key] = value;
}


void LjmetEventContent::SetValue(std::string key,std::vector<std::string> value){
    mVectorStringBranch[key] = value;
}


void LjmetEventContent::SetHistValue(std::string modname, std::string histname, double value)
{
    // Assign current hist value to hist metadata collection
    
    mLegend = "[" + mName + "]: ";
    
    std::map<std::string, std::map<std::string, HistMetadata>>::iterator iMod;
    std::map<std::string, HistMetadata>::iterator iHist;
    
    iMod = mDoubleHist.lower_bound(modname);
    
    if(iMod != mDoubleHist.end() && !(mDoubleHist.key_comp()(modname, iMod->first))) {
        // Key already exists. Update iMod->second
        iHist = iMod->second.lower_bound(histname);
        if(iHist != iMod->second.end() && !(iMod->second.key_comp()(histname, iHist->first))) {
            // Key already exists. Update iHist->second - not needed
            iHist->second.SetValue(value);
        } else {
            // The key does not exist in the map. Add it.
            std::cout << mLegend << "Cannot set value, histogram " << histname << " in module " << modname << " does not exist" << std::endl;
        }
    } else {
        // The key does not exist in the map. Add it to the map
        std::cout << mLegend << "Cannot set value, histogram " << modname << "/" << histname << " does not exist" << std::endl;
    }
}

void LjmetEventContent::Fill()
{
    if (mFirstEntry) {
        createBranches();
        mFirstEntry = false;
    }
    mpTree->Fill();
    
    // fill histograms
    // create histograms
    std::map<std::string, std::map<std::string, LjmetEventContent::HistMetadata>>::iterator iMod;
    std::map<std::string, LjmetEventContent::HistMetadata>::iterator iHist;
    for (iMod = mDoubleHist.begin(); iMod != mDoubleHist.end(); ++iMod) {
        for (iHist = iMod->second.begin(); iHist != iMod->second.end(); ++iHist) {
            TH1 * _hist = iHist->second.GetHist();
            if (_hist) _hist->Fill(iHist->second.GetValue());
        }
    }
}

int LjmetEventContent::createBranches()
{
    // Create branches in the tree according to maps
    
    mLegend = "[" + mName + "]: " ;
    std::string name_type;
    
    std::cout << mLegend << "Creating branches in output tree" << std::endl;
    
    // Boolean branches
    for (std::map<std::string, bool>::iterator br = mBoolBranch.begin(); br != mBoolBranch.end(); ++br) {
        name_type = br->first + "/O";
        mpTree->Branch(br->first.c_str(), &(br->second), name_type.c_str());
        
        if (mVerbosity > 0) {
            std::cout << mLegend << "Branch " << name_type << " created" << std::endl;
        }
    }
    std::cout << mLegend << "bool branches created: " << mBoolBranch.size() << std::endl;
    
    // Integer branches
    for (std::map<std::string, int>::iterator br = mIntBranch.begin(); br != mIntBranch.end(); ++br) {
        name_type = br->first + "/I";
        mpTree->Branch(br->first.c_str(), &(br->second), name_type.c_str());
        
        if (mVerbosity > 0) {
            std::cout << mLegend << "Branch " << name_type << " created" << std::endl;
        }
    }
    std::cout << mLegend << "integer branches created: " << mIntBranch.size() << std::endl;
    

    // Long Integer branches
    for (std::map<std::string, long long>::iterator br = mLongIntBranch.begin(); br != mLongIntBranch.end(); ++br) {
      name_type = br->first + "/L";
      mpTree->Branch(br->first.c_str(), &(br->second), name_type.c_str());
      
      if (mVerbosity > 0) {
	std::cout << mLegend << "Branch " << name_type << " created" << std::endl;
      }
    }
    std::cout << mLegend << "integer branches created: " << mLongIntBranch.size() << std::endl;
    
 

    // Double branches
    for (std::map<std::string, double>::iterator br = mDoubleBranch.begin(); br != mDoubleBranch.end(); ++br) {
        name_type = br->first + "/D";
        mpTree->Branch(br->first.c_str(), &(br->second), name_type.c_str());
        
        if (mVerbosity > 0) {
            std::cout << mLegend << "Branch " << name_type << " created" << std::endl;
        }
    }
    std::cout << mLegend << "double branches created: " << mDoubleBranch.size() << std::endl;
    
    // Vector-of-bool branches
    for (std::map<std::string, std::vector<bool>>::iterator br = mVectorBoolBranch.begin(); br != mVectorBoolBranch.end(); ++br) {
        mpTree->Branch(br->first.c_str(), &(br->second));
        
        if (mVerbosity > 0) {
            std::cout << mLegend << "Branch " << br->first << " std::vector<bool> created" << std::endl;
        }
    }
    std::cout << mLegend << "std::vector<bool> branches created: " << mVectorBoolBranch.size() << std::endl;
    
    // Vector-of-int branches
    for (std::map<std::string, std::vector<int>>::iterator br = mVectorIntBranch.begin(); br != mVectorIntBranch.end(); ++br) {
        mpTree-> Branch(br->first.c_str(), &(br->second));
        
        if (mVerbosity > 0) {
            std::cout << mLegend << "Branch " << br->first << " std::vector<int> created" << std::endl;
        }
    }
    std::cout << mLegend << "std::vector<int> branches created: " << mVectorIntBranch.size() << std::endl;
    
    // Vector-of-double branches
    for (std::map<std::string, std::vector<double>>::iterator br = mVectorDoubleBranch.begin(); br != mVectorDoubleBranch.end(); ++br) {
        mpTree->Branch(br->first.c_str(), &(br->second));
        
        if (mVerbosity > 0) {
            std::cout << mLegend << "Branch " << br->first << " std::vector<double> created" << std::endl;
        }
    }

    std::cout << mLegend << "std::vector<double> branches created: "
    << mVectorDoubleBranch.size() << std::endl;

    // std::vector std::vector std::string branches
    for(std::map<std::string,std::vector<std::string> >::iterator br = mVectorStringBranch.begin();
	br != mVectorStringBranch.end();
	++br){
        std::string name_type = br->first+" std::vector<std::string>";
      //      std::string type = "VVString";
      mpTree -> Branch(br->first.c_str(),
		       &(br->second));

    if (mVerbosity>0){
            std::cout << mLegend << "Branch " << name_type
            << " created" << std::endl;
        }
    }

    return 0;
}
