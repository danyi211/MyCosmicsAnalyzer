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

// use TFileService
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "Geometry/Records/interface/MuonGeometryRecord.h"

#include "Geometry/DTGeometry/interface/DTGeometry.h"

#include "TFile.h"
#include "TTree.h"
//
// class declaration
//

const int kTrackNMax = 10000;
const int kTrackRecHitsNMax = 100;

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
  edm::EDGetTokenT<std::vector<reco::Track>> tracksToken_;
  edm::ESGetToken<DTGeometry, MuonGeometryRecord> muonDTGeomToken_;
  const DTGeometry* muonDTGeom;

  std::vector<double> muonPhis_;  // vector of muon phi's to store in outfile
  std::vector<double> muonPts_;
  std::vector<double> muonEtas_;

  int track_n_;
  std::vector<double> track_vx_;
  std::vector<double> track_vy_;
  std::vector<double> track_vz_;
  std::vector<double> track_phi_;
  std::vector<double> track_py_;
  float track_dtSeg_glbX_[kTrackNMax][kTrackRecHitsNMax];
  float track_dtSeg_glbY_[kTrackNMax][kTrackRecHitsNMax];
  float track_dtSeg_glbZ_[kTrackNMax][kTrackRecHitsNMax];

  // TFile* outputFile_;  // outfile 
  TTree* outputTree_;  // outtree
  TTree* outputTree2_;  // outtree2
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
    : muonToken_(consumes<std::vector<reco::Muon>>(iConfig.getParameter<edm::InputTag>("muonCollection"))),
    tracksToken_(consumes<std::vector<reco::Track>>(iConfig.getParameter<edm::InputTag>("cosmicMuonCollection"))),
    muonDTGeomToken_(esConsumes()) {

  //now do what ever initialization is needed

  // Get the TFileService
  edm::Service<TFileService> fs;
  // Create a TTree using the TFileService
  outputTree_ = fs->make<TTree>("MuonTree", "Muon Distributions");
  outputTree2_ = fs->make<TTree>("CosmicMuonTree", "Cosmic Muon Distributions");
  // Add branches to the TTree
  outputTree_->Branch("muonPhi", &muonPhis_);
  outputTree_->Branch("muonPt", &muonPts_);
  outputTree_->Branch("muonEta", &muonEtas_);

  outputTree2_->Branch("track_vx", &track_vx_);
  outputTree2_->Branch("track_vy", &track_vy_);
  outputTree2_->Branch("track_vz", &track_vz_);
  outputTree2_->Branch("track_phi", &track_phi_);
  outputTree2_->Branch("track_py", &track_py_);
  outputTree2_ -> Branch ( "track_n",    &track_n_);
  outputTree2_ -> Branch ( "track_dtSeg_glbX", track_dtSeg_glbX_, "track_dtSeg_glbX[track_n][100]/F");
  outputTree2_ -> Branch ( "track_dtSeg_glbY", track_dtSeg_glbY_, "track_dtSeg_glbY[track_n][100]/F");
  outputTree2_ -> Branch ( "track_dtSeg_glbZ", track_dtSeg_glbZ_, "track_dtSeg_glbZ[track_n][100]/F");
}

MyAnalyzer::~MyAnalyzer() {
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  // 
  // please remove this method altogether if it would be left empty

  // delete outputTree_;

}

//
// member functions
//

// ------------ method called for each event  ------------
void MyAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  // using namespace edm;
  
  using namespace std;

  // Handle: a smart pointer and is used to hold and access data from the event
  edm::Handle<std::vector<reco::Muon>> muons;
  iEvent.getByToken(muonToken_, muons);  // retrieve data products from the event using tokens
  muonDTGeom = &iSetup.getData(muonDTGeomToken_);

  for (const auto& muon : *muons) {
    double phi = muon.phi();
    double pt = muon.pt();
    double eta = muon.eta();
    muonPhis_.push_back(phi);
    muonPts_.push_back(pt);
    muonEtas_.push_back(eta);
    // Fill the TTree with the data for each muon
    outputTree_->Fill();
  }

  edm::Handle<std::vector<reco::Track>> tracks;
  iEvent.getByToken(tracksToken_, tracks);
  track_n_ = 0;
  for (const auto& track : *tracks) {
    track_vx_.push_back(track.vx());
    track_vy_.push_back(track.vy());
    track_vz_.push_back(track.vz());
    track_phi_.push_back(track.phi());
    track_py_.push_back(track.py());
    
    int track_recHits_n_ = 0;
    for (trackingRecHit_iterator hit = track.recHitsBegin(); hit != track.recHitsEnd(); ++hit) {
      const TrackingRecHit* recHit = *hit;
      if (recHit->isValid()) {
        // GlobalPoint globalPos = recHit->globalPosition();
        const GeomDet* dtDet = muonDTGeom->idToDet(recHit->geographicalId());
        GlobalPoint globalPoint = dtDet->toGlobal(recHit->localPosition());
        track_dtSeg_glbX_[track_n_][track_recHits_n_] = globalPoint.x();
        track_dtSeg_glbY_[track_n_][track_recHits_n_] = globalPoint.y();
        track_dtSeg_glbZ_[track_n_][track_recHits_n_] = globalPoint.z();
        track_recHits_n_ ++;
      }
    }

    track_n_ ++;
    outputTree2_->Fill();
  }
  

  // for (const auto& track : iEvent.get(tracksToken_)) {
  //   // do something with track parameters, e.g, plot the charge.
  //   // int charge = track.charge();
  // }
}

// ------------ method called once each job just before starting event loop  ------------
void MyAnalyzer::beginJob() {
  // please remove this method if not needed

}

// ------------ method called once each job just after ending the event loop  ------------
void MyAnalyzer::endJob() {
  // please remove this method if not needed
  outputTree_->Fill();
  outputTree2_->Fill();
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
