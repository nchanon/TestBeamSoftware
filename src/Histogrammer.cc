/*!
        \file                Histogrammer.cc
        \brief               Class provides user to book, fill histograms
        \author              Suvankar Roy Chowdhury
        \date                05/07/16
        Support :            mail to : suvankar.roy.chowdhury@cern.ch
*/
#include "Histogrammer.h"
#include <climits>
#include <cmath>
#include <cassert>
#include <iomanip>
#include <sstream>
#include <algorithm>
#include <iterator>
#include "TLorentzVector.h"
#include "TFile.h"

using std::cout;
using std::cerr;
using std::endl;
using std::ios;
using std::setw;

Histogrammer::Histogrammer(std::string& outFile) {
  fout_ = new TFile(TString(outFile),"RECREATE");
  isFileopen_ = true;
}

void Histogrammer::bookEventHistograms() {
  fout_->cd();
  fout_->mkdir("EventInfo");
  fout_->cd("EventInfo");
  new TH1I("nevents", "#Events", 10000001, -0.5, 10000000.5);
  new TH1I("dutAngle", "DUT Angle;DUTAngle;#Events", 3100, -0.5, 3099.5);  
  new TH1I("hvSettings", "High Voltage settings;HV;#Events", 1000,-0.5,999.5);
  new TH1I("vcth", "Vcth value;vcth;#Events", 200,-0.5,199.5);
  new TH1I("offset", ";offset;#Events", 200,-0.5,199.5);
  new TH1I("window", ";window;#Events", 200,-0.5,199.5);
  new TH1I("tilt", ";tilt;#Events", 200,-0.5,199.5);
  new TH1I("condData", ";condData;#Events", 20,-0.5,19.5);
  new TH1I("tdcPhase",";tdc;#Events",17,-0.5,16.5);
  new TH1I("isPeriodic",";isPeriodic;#Events",2,-0.5,1.5);
  new TH1I("isGoodFlag",";isGood;#Events",2,-0.5,1.5);
}

void Histogrammer::bookDUTHistograms(std::string det) {
  TString d(det);
  fout_->cd();
  fout_->mkdir(d);
  fout_->cd(d);
  new TH2I("hitmapfull",d + " hitmap;strip no.;#Events",1016,-0.5,1015.5,2,-0.5,1.5);
  bookDUTHistoForColumn(d,"C0");
  //bookDUTHistoForColumn(d,"C1");
}
void Histogrammer::bookDUTHistoForColumn(TString& d, TString c) {
  fout_->cd(d);    
  new TH1I("chsize" + c ,"dut0 channel occupancy per event" + c + ";#Channels;#Events",51,-0.5,50.5);
  new TH1I("hitmap" + c,"dut0 hitmap " + c + ";strip no.;#Events",1016,-0.5,1015.5);
  new TH1D("ncluster" + c,"#cluster dut0 " + c + ";#Clusters;#Events",51,-0.5,50.5);
  new TH1I("clusterWidth" + c,"dut0 clusterWidth " + c +";#ClusterWidth;#Events",128,-0.5,127.5);
  new TH1D("clusterPos" + c,"dut0 clusterPos " + c + ";Strip Number;#Events",1016,-0.5,1015.5);
  new TProfile("clusterWidthVsPosProf" + c,"dut0 clusterWidthVsPos " + c + ";Strip Number;Cluster Width",1016,-0.5,1015.5);
  new TH2D("clusterWidthVsPos2D" + c ,"dut0 clusterWidthVsPos " + c + ";Strip Number;Cluster Width",1016,-0.5,1015.5, 20,-0.5,19.5);
   new TH2D("nhitvsnclus" + c,"#Clusters vs #Hits;#Hits;#Clusters",50,-0.5,49.5, 50,-0.5,49.5);
   new TH2D("nhitvsHitClusPosDiff" + c,"Cluster-Hit MinPosDiff vs #Hits;#Hits;#Pos Diff",50,-0.5,49.5, 256,-0.5,255.5);
   new TH2D("propertyVsTDC2D" + c, "Hit Property vs TDC " + d + c + ";TDC;",17,-0.5,16.5, 10, 0.5, 10.5);
}

void Histogrammer::bookStubHistograms() {
  fout_->cd();
  fout_->mkdir("StubInfo");
  fout_->cd("StubInfo");

  new TH1I("cbcStubWord","Stub Bit from CBC",16,-0.5,15.5);
  new TH1I("recoStubWord","Stub Bit from offline CBC logic emulation",16,-0.5,15.5);
  new TH1I("nstubsFromCBCSword","Total number of stubs from CBC stub word",20,-.5,19.5);
  new TH1I("nstubsFromRecoSword","Total number of stubs from Reco Stub word",20,-.5,19.5); 
  new TH1I("nstubsFromReco","Total number of stubs from Reconstruction",20,-.5,19.5); 
  new TH1I("stubMatch","Matching between CBC Stub and Reco Stub",4,0.5,4.5);
  new TH1I("nstubsdiffSword","#StubsRecoStubword - #StubsfromStubWord",20,-0.5,19.5);
  new TH1I("nstubsdiff","#StubsReco - #StubsfromStubWord",20,-0.5,19.5);
  bookStubHistoForColumn("C0");
  //bookStubHistoForColumn("C1");
}

void Histogrammer::bookStubHistoForColumn(TString c) {
  fout_->cd("StubInfo");
  new TH1I("nstubReco" + c,"Number of stubs for " + c + " from offline reconstruction;#stubs;Events",20,-0.5,19.5);
}

void Histogrammer::bookCorrelationHistograms() {
  fout_->mkdir("Correlation");
  fout_->cd("Correlation");
  bookCorrelationHistoForColumn("C0");
  //bookCorrelationHistoForColumn("C1");
}
    
void Histogrammer::bookCorrelationHistoForColumn(TString c) {
  new TH1D("cor_hit" + c,"Sensor Hit Correlation " + c,4,0.5,4.5);
  new TH1I("nclusterdiff" + c,"Difference in #clusters between dut0 and dut1() for " + c + ";#cluster_{det0} - #cluster_  {det1_};Events",20,-0.5,19.5);
  //new TH1D("nhitdiff" + c,"#Hits det0 - #Hits det1 " + c,100,-.5,99.5);
}


void Histogrammer::bookTrackMatchHistograms() {
  fout_->mkdir("TrackMatch");
  fout_->cd("TrackMatch");
  new TH1D("nTrackParams","#Tracks Telescope;#tracks;#events",30,-0.5,29.5);
  new TH1D("nTrackParamsNodupl","#Tracks Telescope after duplicate removal;#tracks;#events",30,-0.5,29.5);
  new TH1D("hposxTkDUT0","Xpos of extrapolated track at DUT0 plane after alignment(#fid trk=1)",100,-20.,20.);
  new TH1D("hposxTkDUT1","Xpos of extrapolated track at DUT0 plane after alignment(#fid trk=1)",100,-20.,20.);
  new TH1D("hminposClsDUT0","Xpos of the cluster hit at DUT0 plane with min dist from track(#fid trk=1)",100,-20.,20.);
  new TH1D("hminposClsDUT1","Xpos of the cluster hit at DUT1 plane with min dist from track(#fid trk=1)",100,-20.,20.);
  new TH1D("hminposStub","Xpos of the stub hit with min dist from trk(#fid trk=1)",100,-20.,20.);
  new TH2D("minclsTrkPoscorrD0","Cluster xTrk Pos Correlation;trk;cluster",255,-0.5,254.5, 255,-0.5,254.5);
  new TH2D("minclsTrkPoscorrD1","Cluster xTrk Pos Correlation;trk;cluster",255,-0.5,254.5, 255,-0.5,254.5);
  TH2D* h2d = dynamic_cast<TH2D*>(Utility::getHist2D("minclsTrkPoscorrD0"));
  h2d->SetOption("colz");
  h2d = dynamic_cast<TH2D*>(Utility::getHist2D("minclsTrkPoscorrD1"));
  h2d->SetOption("colz");


  //new TH1D("residualDUT0multitrkfidNodupl","ClusterResidual at DUT0 plane(fiducial)(#trk>1, no duplicate tracks)",400,-20.,20.);
  //new TH1D("residualDUT1multitrkfidNodupl","ClusterResidual at DUT1 plane(fiducial)(#trk>1, no duplicate tracks)",400,-20.,20.);
  new TH1D("minresidualDUT0_1trkfid","ClusterResidual at DUT0 plane(fiducial)(#trk=1)",1000,-10.,10.);
  new TH1D("minresidualDUT1_1trkfid","ClusterResidual at DUT1 plane(fiducial)(#trk=1)",1000,-10.,10.);
  new TH1D("clswidthDUT0_1trkfid","ClusterWidth(cluster matched to track) at DUT0 plane(fiducial)(#trk=1)",50,-0.5,49.5);
  new TH1D("clswidthDUT1_1trkfid","ClusterWidth(cluster matched to track) at DUT1 plane(fiducial)(#trk=1)",50,-0.5,49.5);

  //new TH1D("sresidualC0multitrkfidNodupl","Stub Residual at DUT0 plane(fiducial)(#trk>1, no duplicate tracks)",800,-20.,20.);
  new TH2D("minstubTrkPoscorrD1","Cluster xTrk Pos Correlation;trk;stub",255,-0.5,254.5, 255,-0.5,254.5);
  h2d = dynamic_cast<TH2D*>(Utility::getHist2D("minstubTrkPoscorrD1"));
  h2d->SetOption("colz");

  new TH1D("sminresidualC0_1trkfid","Stub Residual at DUT1 plane(fiducial)(#trk=1)",1000,-10.,10.);

  new TH1I("trkcluseff","",9,-0.5,8.5);
  TH1I* h = dynamic_cast<TH1I*>(Utility::getHist1D("trkcluseff"));
  h->GetXaxis()->SetBinLabel(1,"xtkNodupl(=1)");
  h->GetXaxis()->SetBinLabel(2,"xtkFidD0");
  h->GetXaxis()->SetBinLabel(3,"xtkFidD1");
  h->GetXaxis()->SetBinLabel(4,"xtkFidD0 && xtkFidD1");
  h->GetXaxis()->SetBinLabel(5,"xtkClsMatchD0");
  h->GetXaxis()->SetBinLabel(6,"xtkClsMatchD1");
  h->GetXaxis()->SetBinLabel(7,"xtkClsMatchD0_&&_D1");
  h->GetXaxis()->SetBinLabel(8,"no-match_D0&&_D1");
  h->GetXaxis()->SetBinLabel(9,"xtkStubMatchC0");
  
  new TH1I("effVtdc_num",";TDC;#Events",17,-0.5,16.5);
  new TH1I("effVtdc_den",";TDC;#Events",17,-0.5,16.5);

  new TH1F("deltaXPos_trkfei4", "Difference in matched Track impact and Hit X Position", 40000, -20.0, 20.0);
  new TH1F("deltaYPos_trkfei4", "Difference in matched Track Impact and Hit Y Position", 40000, -20.0, 20.0);

}

void Histogrammer::bookTelescopeAnalysisHistograms() {
  fout_->cd();
  fout_->mkdir("TelescopeAnalysis");
  fout_->cd("TelescopeAnalysis");
  new TH1D("nTrack","#Tracks Telescope;#tracks;#events",30,-0.5,29.5);
  new TH1I("nhitsFei4","#Hits in FeI4;#hitsFEI4;#Events",30,-0.5,29.5);

  new TH1F("TkXPos", "Track Impact XPos", 96, -12.0, 12.0);
  new TH1F("TkYPos", "Trach Impact YPos", 240, -12.0, 12.0);
  new TH1F("HtColumn", "Hit Column", 85, -0.5, 84.5);
  new TH1F("HtRow", "Hit Row", 340, -0.5, 339.5);
  new TH1F("HtXPos", "Hit XPos", 96, -12.0, 12.0);
  new TH1F("HtYPos", "Hit YPos", 240, -12.0, 12.0);
  new TH1F("deltaXPos", "Difference in Track impact and Hit X Position", 40000, -20.0, 20.0);
  new TH1F("deltaYPos", "Difference in Track Impact and Hit Y Position", 40000, -20.0, 20.0);
  new TH2F("tkXPosVsHtXPos", "tkXPosVsHtXPos;Xpos of FeI4-Hit(mm);Xpos of Track Impact(mm)", 96, -12.0, 12.0, 96, -12.0, 12.0);
  new TH2F("tkYPosVsHtYPos", "tkYPosVsHtYPos;Ypos of FeI4-Hit(mm);Ypos of Track Impact(mm)", 240, -12.0, 12.0, 240, -12.0, 12.0);
  TH2F* h = dynamic_cast<TH2F*>(Utility::getHist2D("tkXPosVsHtXPos"));
  h->SetOption("colz");
  h = dynamic_cast<TH2F*>(Utility::getHist2D("tkYPosVsHtYPos"));
  h->SetOption("colz");
  new TH1F("deltaXPos_trkfei4", "Difference in Track impact and Hit X Position after alignment", 40000, -20.0, 20.0);
  new TH1F("deltaYPos_trkfei4", "Difference in Track Impact and Hit Y Position after alignment", 40000, -20.0, 20.0);
  new TH1F("deltaXPos_trkfei4M", "Difference in matched Track impact and Hit X Position", 40000, -20.0, 20.0);
  new TH1F("deltaYPos_trkfei4M", "Difference in matched Track Impact and Hit Y Position", 40000, -20.0, 20.0);

}

void Histogrammer::bookTrackFitHistograms(float zMin, float zStep, int zNsteps){
  fout_->cd();
  fout_->mkdir("TrackFit");
  fout_->cd("TrackFit");

  new TH1I("d0_1tk1Hit_diffX","X_{TkAtDUT}-X_{DUT}, d0",10000,-10,10);
  new TH1I("d1_1tk1Hit_diffX","X_{TkAtDUT}-X_{DUT}, d1",10000,-10,10);

  for (int iz=0; iz<zNsteps; iz++){
    new TH1I(Form("d0_1tk1Hit_diffX_iz%i", iz),"X_{TkAtDUT}-X_{DUT}, d0",10000,-10,10);
    new TH1I(Form("d1_1tk1Hit_diffX_iz%i", iz),"X_{TkAtDUT}-X_{DUT}, d0",10000,-10,10);
  }

  float zMax = zMin + ((float)zNsteps) * zStep;
  float shift = zStep/2.;

  new TH1F("d0_offsetVsZ", "x_{DUT} offset vs injected z_{DUT}, d0", zNsteps, zMin-shift, zMax-shift);
  new TH1F("d1_offsetVsZ", "x_{DUT} offset vs injected z_{DUT},d1", zNsteps, zMin-shift, zMax-shift);

  new TH1F("d0_chi2VsZ","chi2 vs injected z_{DUT}, d0", zNsteps, zMin-shift, zMax-shift);
  new TH1F("d1_chi2VsZ","chi2 vs injected z_{DUT}, d1", zNsteps, zMin-shift, zMax-shift);

  new TH1I("d0_1tk1Hit_diffX_aligned","X_{TkAtDUT}-X_{DUT}, d0",10000,-10,10);
  new TH1I("d1_1tk1Hit_diffX_aligned","X_{TkAtDUT}-X_{DUT}, d1",10000,-10,10);
}

TH1* Histogrammer::GetHistoByName(const char* dir, const char* hname){
  fout_->cd(dir);
  return Utility::getHist1D(hname);
}

TH1* Histogrammer::GetHistoByName(const std::string& dir, const std::string& hname) {
  fout_->cd(dir.c_str());
  return Utility::getHist1D(hname);

}

void Histogrammer::FillAlignmentOffsetVsZ(const char* det, const char* histo, int iz, float z, float x, float x_err){

  fout_->cd("TrackFit");
  char histname[50];
  strcpy( histname, det );
  strcat( histname, histo );
  TH1* h = (TH1*) gDirectory->Get(histname);
  h->Fill(z, x);
  h->SetBinError(iz+1, x_err);
}

void Histogrammer::closeFile() { 
  fout_->cd();
  fout_->Write();
  fout_->Close();
  isFileopen_=false;  
}
void Histogrammer::fillClusterHistograms( const char* det, std::vector<tbeam::cluster>& cvec, 
                                          const char* col) {
  fout_->cd(det);
  TString c(col);
  Utility::fillHist1D( "ncluster" + c, cvec.size() );
  for( unsigned int i =0; i<cvec.size(); i++ ) {
    Utility::fillHist1D("clusterWidth" + c,cvec[i].size);
    Utility::fillHist1D("clusterPos" + c,cvec[i].x);
    Utility::fillHistProfile("clusterWidthVsPosProf" + c,cvec[i].x,cvec[i].size);
    Utility::fillHist2D("clusterWidthVsPos2D" + c,cvec[i].x,cvec[i].size);
  } 
}

Histogrammer::~Histogrammer() {
  if(isFileopen_)  {
    std::cout << "You forgot to close the output file!!!Closing it now" << std::endl;  
    closeFile();
   }
  delete fout_; 
}
