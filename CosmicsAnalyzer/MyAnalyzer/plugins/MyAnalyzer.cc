// -*- C++ -*-
//
// Package:    CosmicsAnalyzer/MyAnalyzer
// Class:      MyAnalyzer
//
/**\class MyAnalyzer MyAnalyzer.cc CosmicsAnalyzer/MyAnalyzer/plugins/MyAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Danyi Zhang
//         Created:  Mon, 24 Jul 2023 22:53:01 GMT
//
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/MuonDetId/interface/MuonSubdetId.h"

// use TFileService
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "Geometry/Records/interface/MuonGeometryRecord.h"

#include "Geometry/DTGeometry/interface/DTGeometry.h"

#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"
//
// class declaration
//

const int kMuonNMax = 10000;
const int kTrackNMax = 10000;
const int kTrackRecHitsNMax = 100;
const float kMuonMass = 0.10566; // GeV

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.

using reco::TrackCollection;

class MyAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit MyAnalyzer(const edm::ParameterSet&);
  ~MyAnalyzer() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;

  // ----------member data ---------------------------
  // edm::EDGetTokenT<TrackCollection> tracksToken_;  //used to select what tracks to read from configuration file
  edm::EDGetTokenT<std::vector<reco::Muon>> muonToken_;
  // edm::EDGetTokenT<std::vector<reco::Track>> tracksToken_;
  // edm::ESGetToken<DTGeometry, MuonGeometryRecord> muonDTGeomToken_;
  const DTGeometry* muonDTGeom;

  int    muon_n_;
  float  muon_pt_[kMuonNMax];
  float  muon_p_[kMuonNMax];
  float  muon_eta_[kMuonNMax];
  float  muon_phi_[kMuonNMax];
  float  muon_track_vx_[kMuonNMax];
  float  muon_track_vy_[kMuonNMax];
  float  muon_track_min_dv_;
  float  muon_track_minDv_vx1_;  // positive phi
  float  muon_track_minDv_vy1_;  // positive phi
  float  muon_track_minDv_vz1_;  // positive phi
  float  muon_track_minDv_vx2_;  // negative phi
  float  muon_track_minDv_vy2_;  // negative phi
  float  muon_track_minDv_vz2_;  // negative phi
  int    diMuon_mass_;

  /*
  int track_n_;
  float  track_vx_[kTrackNMax];
  float  track_vy_[kTrackNMax];
  float  track_vz_[kTrackNMax];
  float  track_phi_[kTrackNMax];
  float  track_py_[kTrackNMax];
  float  track_dtSeg_glbX_[kTrackNMax][kTrackRecHitsNMax];
  float  track_dtSeg_glbY_[kTrackNMax][kTrackRecHitsNMax];
  float  track_dtSeg_glbZ_[kTrackNMax][kTrackRecHitsNMax];
  */

  // TFile* outputFile_;  // outfile 
  TTree* outputTree_;  // outtree
  // TTree* outputTree2_;  // outtree2
  
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
MyAnalyzer::MyAnalyzer(const edm::ParameterSet& iConfig)
    : muonToken_(consumes<std::vector<reco::Muon>>(iConfig.getParameter<edm::InputTag>("muonCollection")))
    // tracksToken_(consumes<std::vector<reco::Track>>(iConfig.getParameter<edm::InputTag>("cosmicMuonCollection"))),
    // muonDTGeomToken_(esConsumes()) 
    {

}

MyAnalyzer::~MyAnalyzer() {
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}

//
// member functions
//

// ------------ method called for each event  ------------
void MyAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;
  using namespace std;

  // Handle: a smart pointer and is used to hold and access data from the event
  edm::Handle<std::vector<reco::Muon>> muons;
  iEvent.getByToken(muonToken_, muons);  // retrieve data products from the event using tokens
  // muonDTGeom = &iSetup.getData(muonDTGeomToken_);

  muon_track_min_dv_ = 999.;
  muon_track_minDv_vx1_ = 999.;
  muon_track_minDv_vy1_ = 999.;
  muon_track_minDv_vz1_ = 999.;
  muon_track_minDv_vx2_ = 999.;
  muon_track_minDv_vy2_ = 999.;
  muon_track_minDv_vz2_ = 999.;
  diMuon_mass_ = 999.;
  TLorentzVector v1;
  TLorentzVector v2;

  muon_n_ = 0;
  for (const auto& muon : *muons) {
    muon_pt_[muon_n_] = muon.pt();
    muon_p_[muon_n_] = muon.p();
    muon_pt_[muon_n_] = muon.pt();
    muon_eta_[muon_n_] = muon.eta();
    muon_phi_[muon_n_] = muon.phi();
    muon_n_ ++;

    if (muon.outerTrack().isNonnull()) {
      reco::TrackRef recoTrack = muon.outerTrack();
      muon_track_vx_[muon_n_] = recoTrack->vx();
      muon_track_vy_[muon_n_] = recoTrack->vy();

      v1.SetPtEtaPhiM(muon.pt(), muon.eta(), muon.phi(), kMuonMass);
      // at least two tracks
      if (muons->size() >= 2) {
        for (const auto& muon2 : *muons) {
          // choose two tracks in opposite direction
          if (muon.phi() * muon2.phi() < 0 && muon2.outerTrack().isNonnull()) {
            reco::TrackRef recoTrack2 = muon2.outerTrack();
            double dx = recoTrack->vx() - recoTrack2->vx();
            double dy = recoTrack->vy() - recoTrack2->vy();
            double dz = recoTrack->vz() - recoTrack2->vz();
            double dxy = sqrt(pow(dx, 2) + pow(dy, 2) + pow(dz, 2));
            // cout << "muon vx = " << recoTrack->vx() << ", vy = " << recoTrack->vy();
            // cout << "muon2 vx = " << recoTrack2->vx() << ", vy = " << recoTrack2->vy();
            // cout << "dxy is "<< dxy << endl; 

            if (dxy < muon_track_min_dv_) {
              muon_track_min_dv_ = dxy;
              muon_track_minDv_vx1_ = (muon.phi() > 0) ? muon.vx() : muon2.vx();
              muon_track_minDv_vy1_ = (muon.phi() > 0) ? muon.vy() : muon2.vy();
              muon_track_minDv_vz1_ = (muon.phi() > 0) ? muon.vz() : muon2.vz();
              muon_track_minDv_vx2_ = (muon.phi() < 0) ? muon.vx() : muon2.vx();
              muon_track_minDv_vy2_ = (muon.phi() < 0) ? muon.vy() : muon2.vy();
              muon_track_minDv_vz2_ = (muon.phi() < 0) ? muon.vz() : muon2.vz();

              v2.SetPtEtaPhiM(muon2.pt(), muon2.eta(), muon2.phi(), kMuonMass);
              diMuon_mass_ = (v1 + v2).M();
            }
          }
        }
      }
    }
  }

  // Fill the TTree with the data for each muon
  outputTree_->Fill();

  /*
  edm::Handle<std::vector<reco::Track>> tracks;
  iEvent.getByToken(tracksToken_, tracks);
  track_n_ = 0;
  for (const auto& track : *tracks) {
    track_vx_[track_n_] = track.vx();
    track_vy_[track_n_] = track.vy();
    track_vz_[track_n_] = track.vz();
    track_phi_[track_n_] = track.phi();
    track_py_[track_n_] = track.py();
    
    int track_recHits_n_ = 0;
    for (trackingRecHit_iterator hit = track.recHitsBegin(); hit != track.recHitsEnd(); ++hit) {
      const TrackingRecHit* recHit = *hit;
      if (recHit->isValid()) {
        DetId hitId = recHit->geographicalId();
        // pick recHits from the DT segment
        if (hitId.det() == DetId::Muon && hitId.subdetId() == MuonSubdetId::DT){
          const GeomDet* dtDet = muonDTGeom->idToDet(hitId);
          GlobalPoint globalPoint = dtDet->toGlobal(recHit->localPosition());
          track_dtSeg_glbX_[track_n_][track_recHits_n_] = globalPoint.x();
          track_dtSeg_glbY_[track_n_][track_recHits_n_] = globalPoint.y();
          track_dtSeg_glbZ_[track_n_][track_recHits_n_] = globalPoint.z();
        }
        track_recHits_n_ ++;
      }
    }

  track_n_ ++;
    
  }

  outputTree2_->Fill();
  */
}

// ------------ method called once each job just before starting event loop  ------------
void MyAnalyzer::beginJob() {

  // Get the TFileService
  edm::Service<TFileService> fs;
  // Create a TTree using the TFileService
  outputTree_ = fs->make<TTree>("MuonTree", "Muon Distributions");
  // outputTree2_ = fs->make<TTree>("CosmicMuonTree", "Cosmic Muon Distributions");
  // Add branches to the TTree

  outputTree_  -> Branch("muon_n",   &muon_n_);
  outputTree_  -> Branch("muon_pt",  muon_pt_, "muon_pt[muon_n]/F");
  outputTree_  -> Branch("muon_p",   muon_p_,  "muon_p[muon_n]/F");
  outputTree_  -> Branch("muon_eta", muon_eta_, "muon_eta[muon_n]/F");
  outputTree_  -> Branch("muon_phi", muon_phi_, "muon_phi[muon_n]/F");
  outputTree_  -> Branch("muon_track_vx", muon_track_vx_, "muon_track_vx[muon_n]/F");
  outputTree_  -> Branch("muon_track_vy", muon_track_vy_, "muon_track_vy[muon_n]/F");
  outputTree_  -> Branch("muon_track_min_dv",   &muon_track_min_dv_);
  outputTree_  -> Branch("muon_track_minDv_vx1",   &muon_track_minDv_vx1_);
  outputTree_  -> Branch("muon_track_minDv_vy1",   &muon_track_minDv_vy1_);
  outputTree_  -> Branch("muon_track_minDv_vz1",   &muon_track_minDv_vz1_);
  outputTree_  -> Branch("muon_track_minDv_vx2",   &muon_track_minDv_vx2_);
  outputTree_  -> Branch("muon_track_minDv_vy2",   &muon_track_minDv_vy2_);
  outputTree_  -> Branch("muon_track_minDv_vz2",   &muon_track_minDv_vz2_);
  outputTree_  -> Branch("diMuon_mass",   &diMuon_mass_);

  /*
  outputTree2_ -> Branch( "track_n",          &track_n_);
  outputTree2_ -> Branch( "track_vx",         track_vx_,         "track_vx[track_n]/F");
  outputTree2_ -> Branch( "track_vy",         track_vy_,         "track_vy[track_n]/F");
  outputTree2_ -> Branch( "track_vz",         track_vz_,         "track_vz[track_n]/F");
  outputTree2_ -> Branch( "track_phi",        track_phi_,        "track_phi[track_n]/F");
  outputTree2_ -> Branch( "track_py",         track_py_,         "track_py[track_n]/F");
  outputTree2_ -> Branch( "track_dtSeg_glbX", track_dtSeg_glbX_, "track_dtSeg_glbX[track_n][100]/F");
  outputTree2_ -> Branch( "track_dtSeg_glbY", track_dtSeg_glbY_, "track_dtSeg_glbY[track_n][100]/F");
  outputTree2_ -> Branch( "track_dtSeg_glbZ", track_dtSeg_glbZ_, "track_dtSeg_glbZ[track_n][100]/F");
  */
}

// ------------ method called once each job just after ending the event loop  ------------
void MyAnalyzer::endJob() {
  // please remove this method if not needed
  // outputTree_->Fill();
  // outputTree2_->Fill();
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void MyAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addWithDefaultLabel(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MyAnalyzer);
