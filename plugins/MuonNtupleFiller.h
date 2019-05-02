#ifndef RecoMuon_MuonIdentification_MuonNtupleFiller_H
#define RecoMuon_MuonIdentification_MuonNtupleFiller_H

/** \class MuonNtupleFiller
 *  Dump muon information into a simple ntuple
 *
 *  $Date: 2011/04/06 09:56:30 $
 *  $Revision: 1.4 $
 *  \author P. Traczyk   
 */

// Base Class Headers
#include "FWCore/Framework/interface/EDAnalyzer.h"
//#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"

#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/TrackExtraFwd.h"
#include "DataFormats/TrackReco/interface/TrackToTrackMap.h"
#include "DataFormats/MuonReco/interface/MuonTimeExtra.h"
#include "DataFormats/MuonReco/interface/MuonTimeExtraMap.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "DataFormats/RPCRecHit/interface/RPCRecHitCollection.h"
#include "DataFormats/DTRecHit/interface/DTRecSegment4DCollection.h"
#include "DataFormats/DTRecHit/interface/DTRecHitCollection.h"
#include "DataFormats/CSCRecHit/interface/CSCSegmentCollection.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonShower.h"
#include "DataFormats/L1Trigger/interface/Muon.h"
#include "DataFormats/L1Trigger/interface/L1MuonParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"

#include <TROOT.h>
#include <TSystem.h>
#include "Riostream.h"
#include "TObject.h"
#include <TTree.h>

namespace edm {
  class ParameterSet;
  //  class Event;
  class EventSetup;
  class InputTag;
}

class TFile;
class TTree;

using namespace std;
using namespace edm;
using namespace reco;

class MuonNtupleFiller : public edm::EDAnalyzer {
//class MuonNtupleFiller : public edm::stream::EDProducer<> {
public: 

  explicit MuonNtupleFiller(const edm::ParameterSet&);
  ~MuonNtupleFiller();
  
private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  double iMass(reco::TrackRef imuon, reco::TrackRef iimuon);
  vector<int> countRPChits(reco::TrackRef muon, const edm::Event& iEvent);
  vector<int> countDTsegs(const edm::Event& event, reco::MuonRef muon);
  vector<int> countCSCsegs(const edm::Event& event, reco::MuonRef muon);

  // ----------member data ---------------------------

  edm::ConsumesCollector *iC;
  edm::InputTag TKtrackTags_; 
  edm::InputTag MuonTags_; 
  edm::InputTag TimeTags_; 
  edm::InputTag SIMtrackTags_; 

  string out, open;
  bool debug_;
  bool doSim;
  double theAngleCut;
  double thePtCut;

  edm::EDGetTokenT<reco::BeamSpot> beamSpotToken_;
  edm::EDGetTokenT<reco::TrackCollection> trackToken_;
  edm::EDGetTokenT<reco::MuonCollection> muonToken_;
  edm::EDGetTokenT<l1t::MuonBxCollection> muCollToken_;
  edm::EDGetTokenT<edm::ValueMap<reco::MuonShower>> muons_muonShowerInformation_token_;
  edm::EDGetTokenT<reco::VertexCollection> vertexToken_;
  edm::EDGetTokenT<reco::MuonTimeExtraMap> timeMapCmbToken_;
  edm::EDGetTokenT<reco::MuonTimeExtraMap> timeMapDTToken_;
  edm::EDGetTokenT<reco::MuonTimeExtraMap> timeMapCSCToken_;
  edm::EDGetTokenT<GenParticleCollection> genParticleToken_;
  edm::EDGetTokenT<TrackingParticleCollection> trackingParticleToken_;
  edm::EDGetTokenT<vector<l1extra::L1MuonParticle>> l1extraToken_;
  edm::EDGetTokenT<RPCRecHitCollection> rpcRecHitToken_;
  edm::EDGetTokenT<CSCSegmentCollection> cscSegmentToken_;
  edm::EDGetTokenT<DTRecSegment4DCollection> dtSegmentToken_;

  Handle<reco::MuonCollection> MuCollection;
  Handle<reco::MuonCollection> MuCollectionT;
  Handle<reco::TrackCollection> TKTrackCollection;
  Handle<reco::TrackCollection> STATrackCollection;
  Handle<reco::TrackCollection> GLBTrackCollection;
  Handle<reco::TrackCollection> PMRTrackCollection;
  Handle<reco::TrackCollection> GMRTrackCollection;
  Handle<reco::TrackCollection> FMSTrackCollection;
  Handle<reco::TrackCollection> SLOTrackCollection;
  Handle<edm::SimTrackContainer> SIMTrackCollection;

  Handle<reco::MuonTimeExtraMap> timeMap1;
  Handle<reco::MuonTimeExtraMap> timeMap2;
  Handle<reco::MuonTimeExtraMap> timeMap3;
  
  // ROOT Pointers
  TFile* hFile;
  TTree* t;

  // ***** Tree structure *******
// generator info (if available)
  bool hasSim;
  int genCharge;
  float genPt, genPhi, genEta;
  int genBX;

  bool hasL1;
  int l1Qual[10];
  float l1Pt[10], l1Phi[10], l1Eta[10];
  int l1BX[10];

// event info
  unsigned int event_run;
  unsigned int event_lumi;
  unsigned int event_event;
  int n_vtx;
  double weight;
  bool isCosmic;
  bool isCollision;

// muon ID
  bool isSTA;
  bool isGLB;
  bool isLoose;
  bool isTight;

// muon kinematics and track fit parameters
  int charge;
  float pt, phi, eta, glbpt;
  float dPt;
  float dz;
  float dxy;
  float tkiso;
  
  int nhits[4];
  int nrpchits[4];
  int nsegs[4];
  int nmatches[4];
  
// muon timing
  int muNdof;
  float muTime;
  float muTimeErr;
  int dtNdof;
  float dtTime;
  int cscNdof;
  float cscTime;
  int rpcNdof;
  float rpcTime;
  float rpcTimeErr;

};
#endif
