#ifndef LJMet_Com_interface_LjmetEventContent_h
#define LJMet_Com_interface_LjmetEventContent_h

/*
 Container class for passing internal event content
 
 Author: Gena Kukartsev, 2012
 */

#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <limits>
#include "TH1.h"
#include "TTree.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

class LjmetEventContent {
public:
    /// Container for histogram basic info and current value to be held in event content and filled into hist
    class HistMetadata {
    public:
        HistMetadata(std::string name, int nbins, double xmin, double xmax):
        mName(name),
        mNBins(nbins),
        mXMin(xmin),
        mXMax(xmax),
        mpHist(0),
        mValue(std::numeric_limits<double>::max()) { }
        
        ~HistMetadata() { }
        std::string GetName() { return mName; }
        int GetNBins() { return mNBins; }
        double GetXMin() { return mXMin; }
        double GetXMax() { return mXMax; }
        TH1 * GetHist() { return mpHist; }
        double GetValue() { return mValue; }
        void SetHist(TH1 * pHist){ mpHist = pHist; }
        void SetValue(double value) { mValue = value; }
        
    private:
        HistMetadata() { }
        std::string mName;
        int mNBins;
        double mXMin;
        double mXMax;
        TH1 * mpHist;
        double mValue;
    };
    
    LjmetEventContent();
    LjmetEventContent(std::map<std::string, edm::ParameterSet const> mPar);
    virtual ~LjmetEventContent();
    LjmetEventContent(const LjmetEventContent &); // stop default
    
    void SetTree(TTree * tree);
    
    /// Create histogram entry in event content, so it is created by the LjmetFactory
    void SetHistogram(std::string modname, std::string histname, int nbins, double low, double high);
    
    void SetValue(std::string key, bool value);
    void SetValue(std::string key, int value);
    void SetValue(std::string key, long long value);
    void SetValue(std::string key, double value);
    void SetValue(std::string key, std::vector<bool> value);
    void SetValue(std::string key, std::vector<int> value);
    void SetValue(std::string key, std::vector<double> value);
    void SetValue(std::string key, std::vector<std::string> value);
    // histograms: mDoubleHist[module][histname]
    // actual histograms get created by TFileService in the main application
    // based on info in this container
    std::map<std::string, std::map<std::string, HistMetadata>> & GetHistMap() { return mDoubleHist; }
    
    /// Assign current hist value to hist metadata collection
    void SetHistValue(std::string modname, std::string histname, double value);
    void Fill();
    
private:
    /// Create branches in the tree according to maps
    int createBranches();
    std::string mName;
    std::string mLegend;
    TTree * mpTree;

    
    std::map<std::string,bool> mBoolBranch;
    std::map<std::string,int> mIntBranch;
    std::map<std::string,long long> mLongIntBranch;
    std::map<std::string,double> mDoubleBranch;
    std::map<std::string,std::vector<bool> > mVectorBoolBranch;
    std::map<std::string,std::vector<int> > mVectorIntBranch;
    std::map<std::string,std::vector<double> > mVectorDoubleBranch;
    std::map<std::string,std::vector<std::string> > mVectorStringBranch;
    // mDoubleHist[module][histname]=value
    std::map<std::string,std::map<std::string,HistMetadata> > mDoubleHist;
    bool mFirstEntry;
    int mVerbosity;
};

#endif
