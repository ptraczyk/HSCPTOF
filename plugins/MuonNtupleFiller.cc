// -*- C++ -*-
//
// Package:    MuonNtupleFiller
// Class:      MuonNtupleFiller
// 
/**\class MuonNtupleFiller MuonNtupleFiller.cc 

 Description: Fill muon timing information histograms 

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Piotr Traczyk
//         Created:  Wed Sep 27 14:54:28 EDT 2006
// $Id: MuonNtupleFiller.cc,v 1.5 2011/04/06 11:44:29 ptraczyk Exp $
//
//

#include "MuonNtupleFiller.h"

// system include files
#include <memory>
#include <string>
#include <iostream>
#include <fstream>
#include <iostream>
#include <iomanip>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/Common/interface/Ref.h"

#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/DTGeometry/interface/DTGeometry.h"
#include "Geometry/DTGeometry/interface/DTLayer.h"
#include "Geometry/DTGeometry/interface/DTSuperLayer.h"
#include "DataFormats/DTRecHit/interface/DTSLRecSegment2D.h"
#include "RecoLocalMuon/DTSegment/src/DTSegmentUpdator.h"
#include "RecoLocalMuon/DTSegment/src/DTSegmentCleaner.h"
#include "RecoLocalMuon/DTSegment/src/DTHitPairForFit.h"

#include <Geometry/CSCGeometry/interface/CSCLayer.h>
#include <DataFormats/MuonDetId/interface/CSCDetId.h>
#include <DataFormats/CSCRecHit/interface/CSCRecHit2D.h>
#include <DataFormats/CSCRecHit/interface/CSCRangeMapAccessor.h>

#include "DataFormats/RPCRecHit/interface/RPCRecHit.h"

#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "DataFormats/MuonReco/interface/MuonTimeExtra.h"
#include "DataFormats/MuonReco/interface/MuonTimeExtraMap.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"

#include "Geometry/Records/interface/GlobalTrackingGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/GlobalTrackingGeometry.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "SimDataFormats/CrossingFrame/interface/CrossingFrame.h"
#include "SimDataFormats/CrossingFrame/interface/MixCollection.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertexContainer.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "DataFormats/L1Trigger/interface/L1MuonParticle.h"

#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TFrame.h>

using namespace muon;

//
// constructors and destructor
//
MuonNtupleFiller::MuonNtupleFiller(const edm::ParameterSet& iConfig) 
  :
  TKtrackTags_(iConfig.getUntrackedParameter<edm::InputTag>("TKtracks")),
  MuonTags_(iConfig.getUntrackedParameter<edm::InputTag>("Muons")),
  TimeTags_(iConfig.getUntrackedParameter<edm::InputTag>("Timing")),
  out(iConfig.getParameter<string>("out")),
  open(iConfig.getParameter<string>("open")),
  debug_(iConfig.getParameter<bool>("debug")),
  doSim(iConfig.getParameter<bool>("mctruthMatching")),
  theAngleCut(iConfig.getParameter<double>("angleCut")),
  thePtCut(iConfig.getParameter<double>("PtCut"))
{
  edm::ConsumesCollector collector(consumesCollector());

  beamSpotToken_ = consumes<reco::BeamSpot>(edm::InputTag("offlineBeamSpot"));
  vertexToken_ = consumes<reco::VertexCollection>(edm::InputTag("offlinePrimaryVertices"));
  trackToken_ = consumes<reco::TrackCollection>(TKtrackTags_);

  muonToken_ = consumes<reco::MuonCollection>(MuonTags_);
  muons_muonShowerInformation_token_  = consumes<edm::ValueMap<reco::MuonShower>>(edm::InputTag("muons", "muonShowerInformation", "RECO"));
  muCollToken_ = consumes<l1t::MuonBxCollection>(edm::InputTag("gmtStage2Digis","Muon"));

  timeMapCmbToken_ = consumes<reco::MuonTimeExtraMap>(edm::InputTag(TimeTags_.label(),"combined"));
  timeMapDTToken_ = consumes<reco::MuonTimeExtraMap>(edm::InputTag(TimeTags_.label(),"dt"));
  timeMapCSCToken_ = consumes<reco::MuonTimeExtraMap>(edm::InputTag(TimeTags_.label(),"csc"));

  genParticleToken_ = consumes<GenParticleCollection>(edm::InputTag("genParticles"));
  trackingParticleToken_ = consumes<TrackingParticleCollection>(edm::InputTag("mix","MergedTrackTruth"));

  rpcRecHitToken_ = consumes<RPCRecHitCollection>(edm::InputTag("rpcRecHits")) ;
  cscSegmentToken_ = consumes<CSCSegmentCollection>(edm::InputTag("cscSegments"));
  dtSegmentToken_ = consumes<DTRecSegment4DCollection>(edm::InputTag("dt4DSegments"));
}


MuonNtupleFiller::~MuonNtupleFiller()
{
  if (hFile!=0) {
    hFile->Close();
    delete hFile;
  }
}


//
// member functions
//


bool isNewHighPtMuon(const reco::Muon& muon, const reco::Vertex& vtx){
  if(!muon.isGlobalMuon()) return false;
  
  bool muValHits = ( muon.globalTrack()->hitPattern().numberOfValidMuonHits()>0 ||
                         muon.tunePMuonBestTrack()->hitPattern().numberOfValidMuonHits()>0 );
  bool muMatchedSt = muon.numberOfMatchedStations()>1;
  if(!muMatchedSt) {
    if( muon.numberOfMatchedStations()==1 && muon.isTrackerMuon() && (
        muon.expectedNnumberOfMatchedStations()<2 ||      // for cracks (ExpectedMSC)
        !(muon.stationMask()==1 || muon.stationMask()==16) ||    // for segment reco fail (ZprimeMSC)
        muon.numberOfMatchedRPCLayers()>2) )              // for 1st stations (ZprimeMSC)
      muMatchedSt = true;
  }

  bool ptErr = muon.tunePMuonBestTrack()->ptError()/muon.tunePMuonBestTrack()->pt() < 0.3;

  bool muID = muValHits && muMatchedSt;
  bool hits = muon.innerTrack()->hitPattern().trackerLayersWithMeasurement() > 5 && muon.innerTrack()->hitPattern().numberOfValidPixelHits() > 0; 
//  bool momQuality = muon.tunePMuonBestTrack()->ptError()/muon.tunePMuonBestTrack()->pt() < 0.3;
  bool ip = fabs(muon.innerTrack()->dxy(vtx.position())) < 0.2 && fabs(muon.innerTrack()->dz(vtx.position())) < 0.5;

//  return muID && hits && momQuality && ip;
  return muID && hits && ip && ptErr;
}


// ------------ method called to for each event  ------------
void
MuonNtupleFiller::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //using namespace edm;
  using reco::TrackCollection;
  using reco::MuonCollection;

  bool tpart=false;

  event_run = iEvent.id().run();
  event_lumi = iEvent.id().luminosityBlock();
  event_event = iEvent.id().event();

  if (debug_)
    cout << endl << " Event: " << iEvent.id() << "  Orbit: " << iEvent.orbitNumber() << "  BX: " << iEvent.bunchCrossing() << endl;

  weight = 1.;
  if( !iEvent.isRealData() ) {
    //---- Generator weights
    edm::Handle<GenEventInfoProduct> gen_ev_info;
    iEvent.getByLabel(edm::InputTag("generator"), gen_ev_info);
    if (gen_ev_info.isValid()) weight = gen_ev_info->weight();
  }

  // count the vertices in the event and store the parameters of the PV
  reco::Vertex::Point posVtx;
  reco::Vertex::Error errVtx;
  edm::Handle<reco::VertexCollection> recVtxs;
  iEvent.getByToken(vertexToken_,recVtxs);
  unsigned int theIndexOfThePrimaryVertex = 999;
  n_vtx=0;
  for (unsigned int ind=0; ind<recVtxs->size(); ++ind) 
    if ( (*recVtxs)[ind].isValid()) {
      if (theIndexOfThePrimaryVertex == 999) theIndexOfThePrimaryVertex = ind;
      n_vtx++;
    }

  if (theIndexOfThePrimaryVertex<100) {
    posVtx = ((*recVtxs)[theIndexOfThePrimaryVertex]).position();
    errVtx = ((*recVtxs)[theIndexOfThePrimaryVertex]).error();
  }

  reco::BeamSpot beamSpot;
  edm::Handle<reco::BeamSpot> beamSpotHandle;
  iEvent.getByToken(beamSpotToken_, beamSpotHandle);
  if (beamSpotHandle.isValid()) beamSpot = *beamSpotHandle;
    else {
      cout << "No beam spot available from EventSetup." << endl;
      return;
    }
  math::XYZPoint beamspot(beamSpot.x0(),beamSpot.y0(), beamSpot.z0());

  const reco::Vertex pvertex(posVtx,errVtx);
  posVtx = beamSpot.position();
  errVtx(0,0) = beamSpot.BeamWidthX();
  errVtx(1,1) = beamSpot.BeamWidthY();
  errVtx(2,2) = beamSpot.sigmaZ();

  edm::Handle<TrackingParticleCollection>  TruthTrackContainer ;
  iEvent.getByToken(trackingParticleToken_, TruthTrackContainer );
  if (!TruthTrackContainer.isValid()) {
    if (debug_) cout << " No TrackingParticle data in the Event" << endl;
  } else tpart=true;

  const TrackingParticleCollection *tPC=0;
  if (tpart) 
    tPC = TruthTrackContainer.product();
  
  edm::Handle<reco::TrackCollection> trackc;
  iEvent.getByToken( trackToken_, trackc);
  const reco::TrackCollection trackC = *(trackc.product());
//  if (trackC.size()>2) isCollision=1;
//    else isCollision=0;

  // Generated particle collection
  Handle<GenParticleCollection> genParticles;
  if (doSim)   
    iEvent.getByToken(genParticleToken_, genParticles);

  iEvent.getByToken(muonToken_,MuCollection);
  const reco::MuonCollection muonC = *(MuCollection.product());
  if (debug_) cout << " Muon collection size: " << muonC.size() << endl;
  if (!muonC.size()) return;
  MuonCollection::const_iterator imuon,iimuon;

  edm::Handle<edm::ValueMap<reco::MuonShower> > muonShowerInformationValueMapH_;
  iEvent.getByToken(muons_muonShowerInformation_token_, muonShowerInformationValueMapH_);

  double maxpt=0;

  float angle=0;
  isCosmic=0;
  for(imuon = muonC.begin(); imuon != muonC.end(); ++imuon) {

//    if (muon::isHighPtMuon(*imuon, pvertex )) 
//      if (imuon->tunePMuonBestTrack()->pt()>maxpt) 
//        maxpt=imuon->tunePMuonBestTrack()->pt();
    if (imuon->pt()>maxpt && muon::isLooseMuon(*imuon)) maxpt=imuon->pt();
      
    // check for back-to-back dimuons
    if (muonC.size()>1) {
      if ((imuon->isGlobalMuon() || imuon->isTrackerMuon()) && (imuon->track().isNonnull())) {
        for(iimuon = imuon+1; iimuon != muonC.end(); ++iimuon) 
          if ((iimuon->isGlobalMuon() || iimuon->isTrackerMuon()) && (iimuon->track().isNonnull())) {
            double cross = imuon->track()->momentum().Dot(iimuon->track()->momentum());
            angle = acos(-cross/iimuon->track()->p()/imuon->track()->p());
            if (angle<theAngleCut) isCosmic=1;
          }
      }
    }
  }

  // only store events with a good quality high pT muon
  if (maxpt<thePtCut) {
    if (debug_) cout << "max pT below threshold at : " << maxpt << " GeV. Aborting" << endl;
    return;
  }

  // Analyze L1 information
  edm::Handle<l1t::MuonBxCollection> muColl;
  iEvent.getByToken(muCollToken_, muColl);

  iEvent.getByToken(timeMapCmbToken_,timeMap1);
//  const reco::MuonTimeExtraMap & timeMapCmb = *timeMap1;
  iEvent.getByToken(timeMapDTToken_,timeMap2);
  const reco::MuonTimeExtraMap & timeMapDT = *timeMap2;
  iEvent.getByToken(timeMapCSCToken_,timeMap3);
  const reco::MuonTimeExtraMap & timeMapCSC = *timeMap3;

  int imucount=0;

  // ---------------------------------------------------------------------------------------------
  // ----------------------- main loop over muons ------------------------------------------------
  // ---------------------------------------------------------------------------------------------

  for(imuon = muonC.begin(); imuon != muonC.end(); ++imuon){
    
    reco::MuonRef muonR(MuCollection,imucount);
    imucount++;    

    if (debug_) 
      cout << endl << "   Found muon. Pt: " << imuon->pt() << "   eta: " << imuon->eta() << endl;

    reco::TrackRef glbTrack = imuon->combinedMuon();
    reco::TrackRef trkTrack = imuon->track();
    reco::TrackRef staTrack = imuon->standAloneMuon();
    reco::MuonShower muonShowerInformation = (*muonShowerInformationValueMapH_)[muonR];
    
    // MuonTimeExtra timec = timeMapCmb[muonR];
    MuonTimeExtra timedt = timeMapDT[muonR];
    MuonTimeExtra timecsc = timeMapCSC[muonR];
    reco::MuonTime timerpc = imuon->rpcTime();
    reco::MuonTime timemuon = imuon->time();

    hasSim = 0;
    isSTA = staTrack.isNonnull();
    isGLB = glbTrack.isNonnull();
    isLoose = muon::isLooseMuon(*imuon);
    isTight = isNewHighPtMuon(*imuon, pvertex );
    
    // fill muon kinematics
    pt = imuon->tunePMuonBestTrack()->pt();
    glbpt=pt;
    if (isGLB) glbpt=glbTrack->pt();
    dPt = imuon->tunePMuonBestTrack()->ptError();
    eta = imuon->tunePMuonBestTrack()->eta();
    phi = imuon->tunePMuonBestTrack()->phi();
    charge = imuon->tunePMuonBestTrack()->charge();
    dxy = imuon->tunePMuonBestTrack()->dxy(pvertex.position());
    dz = imuon->tunePMuonBestTrack()->dz(pvertex.position());
    tkiso=imuon->isolationR03().sumPt/pt;

    // remove some junk, the files are getting too big
    if (!isSTA) continue;
    if (pt < 5) continue;

//    vector<int> rpchits={0,0,0,0};
    vector<int> segments_all={0,0,0,0};
    if (isSTA) {
      vector<int> segments_csc={0,0,0,0};
      segments_all=countDTsegs(iEvent,muonR);
      segments_csc=countCSCsegs(iEvent,muonR);
      for (int i=0;i<4;i++) 
        segments_all[i]+=segments_csc[i];
    }
    
//    double detaphi=999;
    int l1idx=0;
    for (int i=0;i<10;i++) l1Pt[i]=0;
    genPt=0;
    // get L1 information
    for (int ibx=muColl->getFirstBX(); ibx<=muColl->getLastBX(); ibx++)
      for (auto it = muColl->begin(ibx); it != muColl->end(ibx); it++){
        l1t::MuonRef l1muon(muColl, distance(muColl->begin(muColl->getFirstBX()),it) );
        double deta=fabs(l1muon->eta()-eta);
        double dphi=fabs(reco::deltaPhi(l1muon->phi(),phi));
//        double dr=sqrt(deta*deta+dphi*dphi);
        
        // insert a protection against cross-matching L1 candidates for close-by muons
        // If the L1 candidate dR w.r.t. a different Loose muon is smaller then don't look at the candidate here
        // (it's blocked by blowing up the dr so that it's not accepted by the dr<0.15 cut)
        for(iimuon = muonC.begin(); iimuon != muonC.end(); ++iimuon){
          if (iimuon!=imuon && iimuon->pt()>thePtCut && muon::isLooseMuon(*iimuon)) {
            double deta2=fabs(l1muon->eta()-iimuon->tunePMuonBestTrack()->eta());
            double dphi2=fabs(reco::deltaPhi(l1muon->phi(),iimuon->tunePMuonBestTrack()->phi()));
//            double dr2=sqrt(deta2*deta2+dphi2*dphi2);
            if (deta2<deta && dphi2<dphi) {
              deta=1.;
              dphi=1.;
              break;
            }
          }
        }

        // any L1 match falling inside the cone is saved
        // NEW: tight matching in phi and loose in eta
        if (dphi<0.1 && deta<0.4) {
          hasL1=1;
          l1Pt[l1idx]=l1muon->pt();
          l1Eta[l1idx]=l1muon->eta();
          l1Phi[l1idx]=l1muon->phi();
          l1Qual[l1idx]=l1muon->hwQual();
          l1BX[l1idx]=ibx;
          if (l1idx==9) cout << " Too many L1 matches..." << endl;
            else l1idx++;
        }
    }
    
    if (debug_) cout << " found " << l1idx << " L1 matches." << endl;

    muNdof = timemuon.nDof;
    muTime = timemuon.timeAtIpInOut;
    muTimeErr = timemuon.timeAtIpInOutErr;

    rpcNdof = timerpc.nDof;
    rpcTime = timerpc.timeAtIpInOut;
    rpcTimeErr = timerpc.timeAtIpInOutErr;

    dtNdof = timedt.nDof();
    dtTime = timedt.timeAtIpInOut();
    cscNdof = timecsc.nDof();
    cscTime = timecsc.timeAtIpInOut();

    // read muon shower information
    for (int i=0; i<4; i++) { // Loop on stations
      nhits[i]  = (muonShowerInformation.nStationHits).at(i);        // number of all the muon RecHits per chamber crossed by a track (1D hits)
      nsegs[i] = segments_all.at(i);
    }

    bool matched=false;

    if (doSim) 
      for (const auto &iTrack : *genParticles) 
        if (fabs(iTrack.pdgId())==13 && iTrack.p4().Pt()>2) {
          if (debug_) cout << iTrack.p4().Pt() << endl;
          if (trkTrack.isNonnull())
            if ((fabs(iTrack.p4().eta()-imuon->track()->momentum().eta())<0.05) &&
                (fabs(reco::deltaPhi(iTrack.p4().phi(),imuon->track()->momentum().phi()))<0.05)) 
              matched=true;
          if (staTrack.isNonnull())
            if ((fabs(iTrack.p4().eta()-staTrack->momentum().eta())<0.05) &&
                (fabs(reco::deltaPhi(iTrack.p4().phi(),staTrack->momentum().phi()))<0.05))
              matched=true;

          if (matched==true) {
            hasSim=1;
            genPt=iTrack.p4().Pt();
            genEta=iTrack.p4().Eta();
            genPhi=iTrack.p4().Phi();
            genCharge=iTrack.pdgId()/13;
            break;
          }
      }

    if (tpart) {
      for (const auto &iTrack : *tPC)
        if (fabs(iTrack.pdgId())==13 && iTrack.p4().Pt()>2) {
          if (debug_) cout << iTrack.p4().Pt() << endl;
          if (trkTrack.isNonnull())
            if ((fabs(iTrack.p4().eta()-imuon->track()->momentum().eta())<0.05) &&
                (fabs(reco::deltaPhi(iTrack.p4().phi(),imuon->track()->momentum().phi()))<0.05)) 
              matched=true;
          if (staTrack.isNonnull())
            if ((fabs(iTrack.p4().eta()-staTrack->momentum().eta())<0.05) &&
                (fabs(reco::deltaPhi(iTrack.p4().phi(),staTrack->momentum().phi()))<0.05))
              matched=true;

          if (matched==true) {
            hasSim=1;
            genPt=iTrack.p4().Pt();
            genEta=iTrack.p4().Eta();
            genPhi=iTrack.p4().Phi();
            genBX=iTrack.eventId().bunchCrossing();
            genCharge=iTrack.pdgId()/13;
            break;
          }
        }
    }

    t->Fill();
  }

}


// ------------ method called once each job just before starting event loop  ------------
void 
MuonNtupleFiller::beginJob()
{
   hFile = new TFile( out.c_str(), open.c_str() );
   hFile->cd();

   t = new TTree("MuTree", "MuTree");
   t->Branch("hasSim", &hasSim, "hasSim/O");
//   t->Branch("genCharge", &genCharge, "genCharge/I");
//   t->Branch("genPt", &genPt, "genPt/F");
//   t->Branch("genPhi", &genPhi, "genPhi/F");
//   t->Branch("genEta", &genEta, "genEta/F");
//   t->Branch("genBX", &genBX, "genBX/I");

   t->Branch("hasL1", &hasL1, "hasL1/O");
   t->Branch("l1Qual", &l1Qual, "l1Qual[10]/I");
   t->Branch("l1Pt", &l1Pt, "l1Pt[10]/F");
   t->Branch("l1Phi", &l1Phi, "l1Phi[10]/F");
   t->Branch("l1Eta", &l1Eta, "l1Eta[10]/F");
   t->Branch("l1BX", &l1BX, "l1BX[10]/I");

   t->Branch("event_run", &event_run, "event_run/i");
   t->Branch("event_lumi", &event_lumi, "event_lumi/i");
   t->Branch("event_event", &event_event, "event_event/i");
   t->Branch("nVtx", &n_vtx, "n_vtx/i");
//   t->Branch("weight", &weight, "weight/F");
   t->Branch("isCosmic", &isCosmic, "isCosmic/O");
//   t->Branch("isCollision", &isCollision, "isCollision/O");

   t->Branch("isSTA", &isSTA, "isSTA/O");
   t->Branch("isGLB", &isGLB, "isGLB/O");
   t->Branch("isLoose", &isLoose, "isLoose/O");
   t->Branch("isTight", &isTight, "isTight/O");

   t->Branch("charge", &charge, "charge/I");
   t->Branch("pt", &pt, "pt/F");
   t->Branch("glbpt", &glbpt, "glbpt/F");
   t->Branch("phi", &phi, "phi/F");
   t->Branch("eta", &eta, "eta/F");
   t->Branch("dPt", &dPt, "dPt/F");
//   t->Branch("dz", &dz, "dz/F");
//   t->Branch("dxy", &dxy, "dxy/F");
   t->Branch("tkiso", &tkiso, "tkiso/F");

   t->Branch("nhits", &nhits, "nhits[4]/I");
//   t->Branch("nrpchits", &nrpchits, "nrpchits[4]/I");
   t->Branch("nsegs", &nsegs, "nsegs[4]/I");
//   t->Branch("nmatches", &nmatches, "nmatches[4]/I");

//   t->Branch("muNdof", &muNdof, "muNdof/I");
//   t->Branch("muTime", &muTime, "muTime/F");
//   t->Branch("muTimeErr", &muTimeErr, "muTimeErr/F");
   t->Branch("dtNdof", &dtNdof, "dtNdof/I");
   t->Branch("dtTime", &dtTime, "dtTime/F");
//   t->Branch("cscNdof", &cscNdof, "cscNdof/I");
//   t->Branch("cscTime", &cscTime, "cscTime/F");
   t->Branch("rpcNdof", &rpcNdof, "rpcNdof/I");
   t->Branch("rpcTime", &rpcTime, "rpcTime/F");
   t->Branch("rpcTimeErr", &rpcTimeErr, "rpcTimeErr/F");
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MuonNtupleFiller::endJob() {

  hFile->cd();
  t->Write();
  hFile->Write();
  delete t;  
}

double MuonNtupleFiller::iMass(reco::TrackRef imuon, reco::TrackRef iimuon) {
  double energy1 = sqrt(imuon->p() * imuon->p() + 0.011163691);
  double energy2 = sqrt(iimuon->p() * iimuon->p() + 0.011163691);
  math::XYZTLorentzVector ip4(imuon->px(),imuon->py(),imuon->pz(),energy1);
  math::XYZTLorentzVector iip4(iimuon->px(),iimuon->py(),iimuon->pz(),energy2);
  math::XYZTLorentzVector psum = ip4+iip4;
  double mmumu2 = psum.Dot(psum);
  return sqrt(mmumu2);
}

vector<int> MuonNtupleFiller::countRPChits(reco::TrackRef muon, const edm::Event& iEvent) {
  double RPCCut = 30.;

  int layercount[8]={0,0,0,0,0,0,0,0};
  vector<int> stations;
  edm::Handle<RPCRecHitCollection> rpcRecHits;
  iEvent.getByToken(rpcRecHitToken_, rpcRecHits);

  for(RPCRecHitCollection::const_iterator hitRPC = rpcRecHits->begin(); hitRPC != rpcRecHits->end(); hitRPC++) {
    if ( !hitRPC->isValid()) continue;
    
    // only cound in-time RPC hits
    if ( hitRPC->BunchX()!=0) continue;
    RPCDetId myChamber((*hitRPC).geographicalId().rawId());
    LocalPoint posLocalRPC = hitRPC->localPosition();

    for(trackingRecHit_iterator hitC = muon->recHitsBegin(); hitC != muon->recHitsEnd(); ++hitC) {
      if (!(*hitC)->isValid()) continue; 
      if ( (*hitC)->geographicalId().det() != DetId::Muon ) continue; 
      if ( (*hitC)->geographicalId().subdetId() != MuonSubdetId::RPC ) continue;

      // Check that we're in the same RPC chamber
      DetId id = (*hitC)->geographicalId();
      RPCDetId rpcDetIdHit(id.rawId());
      if (rpcDetIdHit!=myChamber) continue;

      // Compare local positions of the muon hit and the hit being considered
      LocalPoint posLocalMuon = (*hitC)->localPosition();
      if((fabs(posLocalMuon.x()-posLocalRPC.x())<RPCCut)) 
        layercount[(rpcDetIdHit.station()-1)*2+rpcDetIdHit.layer()-1]++;
    }
  }
  
  stations.push_back(max(layercount[0],layercount[1]));
  stations.push_back(max(layercount[2],layercount[3]));
  stations.push_back(layercount[4]);
  stations.push_back(layercount[6]);
  
  return stations;
}



vector<int> MuonNtupleFiller::countDTsegs(const edm::Event& iEvent, reco::MuonRef muon) {
  double DTCut = 25.;

  std::vector<int> stations={0,0,0,0};

  edm::Handle<DTRecSegment4DCollection> dtSegments;
  iEvent.getByToken(dtSegmentToken_, dtSegments);

  for (const auto &ch : muon->matches()) {
    if( ch.detector() != MuonSubdetId::DT )  continue;
    DTChamberId DTid( ch.id.rawId() );

    int nsegs_temp = 0;
    std::vector<float> nsegs_x_temp, nsegs_y_temp;

    for (auto seg = dtSegments->begin(); seg!=dtSegments->end(); ++seg) {
      DTChamberId myChamber((*seg).geographicalId().rawId());
      if (!(DTid==myChamber))  continue;
      LocalPoint posLocalSeg = seg->localPosition();

      if( ( posLocalSeg.x()!=0 && ch.x!=0 ) && (fabs(posLocalSeg.x()-ch.x)<DTCut) ) {
        bool found = false;
        for( auto prev_x : nsegs_x_temp) {
          if( fabs(prev_x-posLocalSeg.x()) < 0.1 ) {
            found = true;
            break;
          }
        }
        if( !found )  nsegs_x_temp.push_back(posLocalSeg.x());
      }

      if( ( posLocalSeg.y()!=0 && ch.y!=0 ) && (fabs(posLocalSeg.y()-ch.y)<DTCut) ) {
        bool found = false;
        for( auto prev_y : nsegs_y_temp) {
          if( fabs(prev_y-posLocalSeg.y()) < 0.1 ) {
            found = true;
            break;
          }
        }
        if( !found )  nsegs_y_temp.push_back(posLocalSeg.y());
      }

    }

    nsegs_temp = (int)std::max(nsegs_x_temp.size(), nsegs_y_temp.size());

    //--- subtract best matched segment from given muon
    bool isBestMatched = false;
    for(std::vector<reco::MuonSegmentMatch>::const_iterator matseg = ch.segmentMatches.begin(); matseg != ch.segmentMatches.end(); matseg++) {
      if( matseg->isMask(reco::MuonSegmentMatch::BestInChamberByDR) ) {
        isBestMatched = true;
        break;
      }
    }

    if(isBestMatched) nsegs_temp = nsegs_temp-1;
    if(nsegs_temp>0)  stations[ch.station()-1] += nsegs_temp;

  }

  return stations;
}


vector<int> MuonNtupleFiller::countCSCsegs(const edm::Event& iEvent, reco::MuonRef muon) {
  double CSCCut = 25.;

  std::vector<int> stations={0,0,0,0};

  edm::Handle<CSCSegmentCollection> cscSegments;
  iEvent.getByToken(cscSegmentToken_, cscSegments);  

  for (const auto &ch : muon->matches()) {
    if( ch.detector() != MuonSubdetId::CSC )  continue;
    CSCDetId CSCid( ch.id.rawId() );

    int nsegs_temp = 0;
    std::vector<float> nsegs_phi_temp;
    std::vector<float> nsegs_nhit_temp;

    for (auto seg = cscSegments->begin(); seg!=cscSegments->end(); ++seg) {
      CSCDetId myChamber((*seg).geographicalId().rawId());
      if (!(CSCid==myChamber))  continue;
      LocalPoint posLocalSeg = seg->localPosition();

      if( (posLocalSeg.x()!=0 && posLocalSeg.y()!=0) && (sqrt( (posLocalSeg.x()-ch.x)*(posLocalSeg.x()-ch.x) + (posLocalSeg.y()-ch.y)*(posLocalSeg.y()-ch.y) )<CSCCut) )  { 
        bool found = false;
        for( auto prev_phi : nsegs_phi_temp) {
          if( fabs(prev_phi-posLocalSeg.phi()) < 0.0002 ) {
            for( auto prev_nhit2 : nsegs_nhit_temp) {
              if( fabs(prev_nhit2-seg->nRecHits()) < 2 ) {
                found = true;
                break;
              }
            }
          }
        }
        if( !found ) {
          nsegs_phi_temp.push_back(posLocalSeg.phi());
          nsegs_nhit_temp.push_back(seg->nRecHits());
          nsegs_temp++;
        }
      }
    }

    //--- subtract best matched segment from given muon
    bool isBestMatched = false;
    for(std::vector<reco::MuonSegmentMatch>::const_iterator matseg = ch.segmentMatches.begin(); matseg != ch.segmentMatches.end(); matseg++) {
      if( matseg->isMask(reco::MuonSegmentMatch::BestInChamberByDR) ) {
        isBestMatched = true;
        break;
      }
    }

    if(isBestMatched) nsegs_temp = nsegs_temp-1;
    if(nsegs_temp>0)  stations[ch.station()-1] += nsegs_temp;
  }

  return stations;
}



//define this as a plug-in
DEFINE_FWK_MODULE(MuonNtupleFiller);
