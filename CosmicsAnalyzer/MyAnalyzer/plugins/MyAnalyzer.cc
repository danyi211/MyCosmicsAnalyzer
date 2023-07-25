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

// use TFileService
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TFile.h"
#include "TTree.h"
//
// class declaration
//

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
  std::vector<double> muonPhis_;  // vector of muon phi's to store in outfile
  std::vector<double> muonPts_;
  std::vector<double> muonEtas_;
  // TFile* outputFile_;  // outfile 
  TTree* outputTree_;  // outtree
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
    : muonToken_(consumes<std::vector<reco::Muon>>(iConfig.getParameter<edm::InputTag>("muonCollection"))) {

  //now do what ever initialization is needed

  // Get the TFileService
  edm::Service<TFileService> fs;
  // Create a TTree using the TFileService
  outputTree_ = fs->make<TTree>("MuonTree", "Muon Distributions");
  // Add branches to the TTree
  outputTree_->Branch("muonPhi", &muonPhis_);
  outputTree_->Branch("muonPt", &muonPts_);
  outputTree_->Branch("muonEta", &muonEtas_);
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
  
  // Handle: a smart pointer and is used to hold and access data from the event
  edm::Handle<std::vector<reco::Muon>> muons;
  iEvent.getByToken(muonToken_, muons);  // retrieve data products from the event using tokens

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
