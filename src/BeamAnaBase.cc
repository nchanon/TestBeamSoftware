/*!
        \file                BeamAnaBase.cc
        \brief               Base Analysis class, reads input files, sets the hit vectors etc. Provides
                             access to hits, clusters, stubs, condition, telescope data. All individual applications
                             should inherit from this class.
        \author              Suvankar Roy Chowdhury
        \date                05/07/16
        Support :            mail to : suvankar.roy.chowdhury@cern.ch
*/

#include "BeamAnaBase.h"
#include "Utility.h"
#include "TSystem.h"
#include "TChain.h"
#include<algorithm>
#include <fstream>

BeamAnaBase::BeamAnaBase() :
  fin_(nullptr),
  analysisTree_(nullptr),
  dutEv_(new tbeam::dutEvent()),
  condEv_(new tbeam::condEvent()),
  telEv_(new  tbeam::TelescopeEvent()),
  periodcictyF_(false),
  isGood_(false),
  hasTelescope_(false),
  doTelMatching_(false),
  doChannelMasking_(false),  
  sw_(-1),
  offset1_(-1),
  offset2_(-1),
  cwd_(-1),
  dut0_chtempC0_(new std::vector<int>()),
  dut0_chtempC1_(new std::vector<int>()),
  dut1_chtempC0_(new std::vector<int>()),
  dut1_chtempC1_(new std::vector<int>()),
  dutRecoClmap_(new std::map<std::string,std::vector<tbeam::cluster>>),
  dutRecoStubmap_(new std::map<std::string,std::vector<tbeam::stub>>),
  recostubChipids_(new std::map<std::string,std::vector<unsigned int>>()),
  cbcstubChipids_(new std::map<std::string,std::vector<unsigned int>>()),
  dut_maskedChannels_(new std::map<std::string,std::vector<int>>()),
  nStubsrecoSword_(0),
  nStubscbcSword_(0)
{
  dutRecoClmap_->insert({("det0C0"),std::vector<tbeam::cluster>()});
  dutRecoClmap_->insert({("det0C1"),std::vector<tbeam::cluster>()});
  dutRecoClmap_->insert({("det1C0"),std::vector<tbeam::cluster>()});
  dutRecoClmap_->insert({("det1C1"),std::vector<tbeam::cluster>()});
  dutRecoStubmap_->insert({("C0"),std::vector<tbeam::stub>()});
  dutRecoStubmap_->insert({("C1"),std::vector<tbeam::stub>()});
  recostubChipids_->insert({("C0"),std::vector<unsigned int>()});
  recostubChipids_->insert({("C1"),std::vector<unsigned int>()});
  cbcstubChipids_->insert({("C0"),std::vector<unsigned int>()});
  cbcstubChipids_->insert({("C1"),std::vector<unsigned int>()});
}

bool BeamAnaBase::setInputFile(const std::string& fname) {
  fin_ = TFile::Open(fname.c_str());
  if(!fin_)    {
    std::cout <<  "File " << fname << " could not be opened!!" << std::endl;
    return false;
  }
  analysisTree_ = dynamic_cast<TTree*>(fin_->Get("analysisTree"));
  if(analysisTree_)    return true;
  return false; 
}

void BeamAnaBase::setTelMatching(const bool mtel) {
  doTelMatching_ = mtel;
}

void BeamAnaBase::setChannelMasking(const bool mch, const std::string cFile) {
  doChannelMasking_ = mch;
  if(doChannelMasking_)   readChannelMaskData(cFile);
}
bool BeamAnaBase::branchFound(const string& b)
{
  TBranch* branch = analysisTree_->GetBranch(b.c_str());
  if (!branch) {
    cout << ">>> SetBranchAddress: <" << b << "> not found!" << endl;
    return false;
  }
  cout << ">>> SetBranchAddress: <" << b << "> found!" << endl;
  hasTelescope_ = true;
  return true;
}


void BeamAnaBase::setAddresses() {
  //set the address of the DUT tree
  if(branchFound("DUT"))    analysisTree_->SetBranchAddress("DUT", &dutEv_);
  if(branchFound("Condition"))    analysisTree_->SetBranchAddress("Condition", &condEv_);
  if(branchFound("TelescopeEvent"))    analysisTree_->SetBranchAddress("TelescopeEvent",&telEv_);
  if(branchFound("periodicityFlag"))    analysisTree_->SetBranchAddress("periodicityFlag",&periodcictyF_);
  if(branchFound("goodEventFlag"))    analysisTree_->SetBranchAddress("goodEventFlag",&isGood_);
  analysisTree_->SetBranchStatus("*",1);
}

void BeamAnaBase::setDetChannelVectors() {
  if(doChannelMasking_) {
    if( dutEv_->dut_channel.find("det0") != dutEv_->dut_channel.end() )
      Utility::getChannelMaskedHits(dutEv_->dut_channel.at("det0"), dut_maskedChannels_->at("det0")); 
    if( dutEv_->dut_channel.find("det1") != dutEv_->dut_channel.end() )
      Utility::getChannelMaskedHits(dutEv_->dut_channel.at("det1"), dut_maskedChannels_->at("det1")); 
    if( dutEv_->clusters.find("det1") != dutEv_->clusters.end() )  
      Utility::getChannelMaskedClusters(dutEv_->clusters.at("det0"), dut_maskedChannels_->at("det0"));
    if( dutEv_->clusters.find("det1") != dutEv_->clusters.end() )
      Utility::getChannelMaskedClusters(dutEv_->clusters.at("det1"), dut_maskedChannels_->at("det1"));
    //stub seeding layer os det1
    Utility::getChannelMaskedStubs(dutEv_->stubs,dut_maskedChannels_->at("det1"));
  }
  //std::cout << "setP1" << std::endl;
  if( dutEv_->dut_channel.find("det0") != dutEv_->dut_channel.end() ) {
      for( unsigned int j = 0; j<(dutEv_->dut_channel.at("det0")).size(); j++ ) {
        int ch = (dutEv_->dut_channel.at("det0")).at(j);
	if( ch <= 1015 )  dut0_chtempC0_->push_back(ch);
	else dut0_chtempC1_->push_back(ch-1016);
      }
  }
  if( dutEv_->dut_channel.find("det1") != dutEv_->dut_channel.end() ) {
      for( unsigned int j = 0; j<(dutEv_->dut_channel.at("det1")).size(); j++ ) {
        int ch = (dutEv_->dut_channel.at("det1")).at(j);
        if( ch <= 1015 )  dut1_chtempC0_->push_back(ch);
        else  dut1_chtempC1_->push_back(ch-1016);
      }
  }
  //std::cout << "setP2" << std::endl;
  for(auto& cl : (dutEv_->clusters)){
    std::string ckey = cl.first;//keys are det0 and det1
    for(auto& c : cl.second)  {
      if(c->x <= 1015)  dutRecoClmap_->at(ckey +"C0").push_back(*c);
      else {
        auto ctemp = *c;
        ctemp.x -= 1016;//even for column 1 we fill histograms between 0 and 1015 
        dutRecoClmap_->at( ckey + "C1").push_back(ctemp);
      }
    }    
  }
  //std::cout << "setP3" << std::endl;
  for(auto& s : dutEv_->stubs) {
    tbeam::stub st = *s;
    if(st.x <= 1015)   dutRecoStubmap_->at("C0").push_back(st);
    else dutRecoStubmap_->at("C1").push_back(st);
  }
  
  //for(auto& t:*recostubChipids_) std::cout << t.first << ",";
  //for(auto& t:*cbcstubChipids_) std::cout << t.first << ",";
  nStubsrecoSword_ = Utility::readStubWord(*recostubChipids_,dutEv_->stubWordReco);
  nStubscbcSword_ = Utility::readStubWord(*cbcstubChipids_,dutEv_->stubWord);
  //std::cout << "Leaving set" << std::endl;
}


void BeamAnaBase::getCbcConfig(uint32_t cwdWord, uint32_t windowWord){
  sw_ = windowWord >>4;
  offset1_ = (cwdWord)%4;
  if ((cwdWord>>2)%2) offset1_ = -offset1_;
  offset2_ = (cwdWord>>3)%4;
  if ((cwdWord>>5)%2) offset2_ = -offset2_;
  cwd_ = (cwdWord>>6)%4;
}

void BeamAnaBase::getExtrapolatedTracks(std::vector<double>& xTkdut0, std::vector<double>& xTkdut1) {
  for(unsigned int itrk = 0; itrk<telEv_->nTrackParams;itrk++) {
    double XTkatDUT0_itrk = (z_DUT0-z_FEI4)*telEv_->dxdz->at(itrk) + telEv_->xPos->at(itrk);
    double XTkatDUT1_itrk = (z_DUT1-z_FEI4)*telEv_->dxdz->at(itrk) + telEv_->xPos->at(itrk);
    //check for duplicate tracks
    if(std::find(xTkdut0.begin(), xTkdut0.end(), XTkatDUT0_itrk) == xTkdut0.end() 
       && std::find(xTkdut1.begin(), xTkdut1.end(), XTkatDUT1_itrk) == xTkdut1.end() ) {
      xTkdut0.push_back(XTkatDUT0_itrk);
      xTkdut1.push_back(XTkatDUT1_itrk);
    }
  }
}

void BeamAnaBase::readChannelMaskData(const std::string cmaskF) {
  std::ifstream fin(cmaskF.c_str(),std::ios::in);
  if(!fin) {
    std::cout << "Channel Mask File could not be opened!!" << std::endl;
    return;
  }
  while(fin) {
    std::string line;
    std::getline(fin,line);
    //std::cout << "Line=" << line << ">>" << fin << std::endl;
    if(fin) {  
      if (line.substr(0,1) == "#" || line.substr(0,2) == "//") continue;
      std::vector<std::string> tokens;
      //first split against :
      Utility::tokenize(line,tokens,":");
      std::vector<std::string> maskedCh;
      //first split against , to get masked hits
      Utility::tokenize(tokens[1],maskedCh,",");
      int cbcid = std::atoi(tokens[0].c_str());
      for(auto& ch : maskedCh) {
        cbcMaskedChannelsMap_[cbcid].push_back(std::atoi(ch.c_str()));
      }
    }
  }
  fin.close();  
   
  std::cout << "Masked Channels List" << std::endl;
  dut_maskedChannels_->insert({("det0"),std::vector<int>()});
  dut_maskedChannels_->insert({("det1"),std::vector<int>()});
  for(auto& cbc : cbcMaskedChannelsMap_) {
    int cbcId = cbc.first;
    std::cout << "CBCid=" << cbcId << "  MaskedCh>>";
    int hitposX = -1;
    for(auto& ch : cbc.second) {
      std::cout << ch << ",";
      int ichan = ch / 2;
      if(cbcId <= 7) {
	hitposX = 127*cbcId + ichan;
      } else {
	  hitposX = 2032 - (127*cbcId + ichan);  
      }
      if( ch % 2 == 0 ) 
        for(int ic = hitposX-2; ic <= hitposX+2; ic++)
          dut_maskedChannels_->at("det1").push_back(ic);
      else 
        for(int ic = hitposX-2; ic <= hitposX+2; ic++)
          dut_maskedChannels_->at("det0").push_back(ic);
    }
    std::cout << std::endl;
  }
  
  std::cout << "Masked Channels Unfolded>>" << std::endl;
  for( auto& d : *dut_maskedChannels_) {
    std::cout << "DET=" << d.first << "  Masked Channels>>";
    for(auto& ch : d.second) 
      std::cout << ch << ",";
    std::cout << std::endl;
  } 
}
void BeamAnaBase::endJob() {
  
}
void BeamAnaBase::clearEvent() {
  dut0_chtempC0_->clear();
  dut0_chtempC1_->clear();
  dut1_chtempC0_->clear();
  dut1_chtempC1_->clear();
  for(auto& c: *dutRecoClmap_)
    c.second.clear();
  for(auto& s: *dutRecoStubmap_)
    s.second.clear();
  for(auto& rs : *recostubChipids_)
    rs.second.clear();
  for(auto& rs : *cbcstubChipids_)
    rs.second.clear();
  nStubsrecoSword_ = 0;
  nStubscbcSword_ = 0;
}

BeamAnaBase::~BeamAnaBase() {
}
