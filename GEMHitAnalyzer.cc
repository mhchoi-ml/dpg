#ifndef GEMHitAnalyzer_H
#define GEMHitAnalyzer_H
// cd /cms/ldap_home/iawatson/scratch/GEM/CMSSW_10_1_5/src/ && eval `scramv1 runtime -sh` && eval `scramv1 runtime -sh` && scram b -j 10
// cd ../../.. && source /cvmfs/cms.cern.ch/cmsset_default.sh && eval `scramv1 runtime -sh` && eval `scramv1 runtime -sh` && scram b -j 10
// system include files
#include <memory>
#include <cmath>
#include <iostream>
#include <map>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "RecoMuon/TrackingTools/interface/MuonServiceProxy.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "TrackingTools/PatternTools/interface/Trajectory.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
// GEM
#include "DataFormats/GEMDigi/interface/GEMDigiCollection.h"
#include "DataFormats/GEMRecHit/interface/GEMRecHitCollection.h"
#include "DataFormats/GEMRecHit/interface/GEMSegmentCollection.h"
#include "DataFormats/MuonDetId/interface/GEMDetId.h"
#include "Geometry/GEMGeometry/interface/GEMGeometry.h"
#include "Geometry/GEMGeometry/interface/GEMEtaPartition.h"
#include "Geometry/GEMGeometry/interface/GEMEtaPartitionSpecs.h"
#include "Geometry/CommonTopologies/interface/GEMStripTopology.h"
// Muon
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeedCollection.h"
#include "DataFormats/TrackingRecHit/interface/KfComponentsHolder.h"

#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Run.h"

#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TString.h"
#include "TGraphAsymmErrors.h"
#include "TLorentzVector.h"
#include "TTree.h"

using namespace std;

typedef std::tuple<int, int, int, int> Key4;

class GEMHitAnalyzer : public edm::one::EDAnalyzer<edm::one::WatchRuns> {  
public:
  explicit GEMHitAnalyzer(const edm::ParameterSet&);
  ~GEMHitAnalyzer();

private:
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void beginJob() override;
  virtual void endJob() override;

  virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  virtual void endRun(edm::Run const&, edm::EventSetup const&) override;

  // ----------member data ---------------------------
  edm::Service<TFileService> fs;
  edm::EDGetTokenT<GEMDigiCollection> gemDigis_;
  edm::EDGetTokenT<GEMRecHitCollection> gemRecHits_;
  edm::ESGetToken<GEMGeometry, MuonGeometryRecord> hGEMGeom_; 
  edm::ESGetToken<GEMGeometry, MuonGeometryRecord> hGEMGeomBeginRun_; 

  std::map<Key4, TH2D*> digi_occ_;
  std::map<Key4, TH2D*> rechit_occ_;

  TTree *t_hit;
  int b_region, b_station, b_chamber, b_layer, b_etaPartition, b_firstStrip, b_nStrips;
};

GEMHitAnalyzer::GEMHitAnalyzer(const edm::ParameterSet& iConfig)
  : hGEMGeom_(esConsumes()),
    hGEMGeomBeginRun_(esConsumes<edm::Transition::BeginRun>())
{ 
  gemRecHits_ = consumes<GEMRecHitCollection>(iConfig.getParameter<edm::InputTag>("gemRecHit"));
  gemDigis_ = consumes<GEMDigiCollection>(iConfig.getParameter<edm::InputTag>("gemDigi"));
//  hGEMGeomBegin_ = esConsumes<GEMGeometry, MuonGeometryRecord>(); 
//  hGEMGeom_ = esConsumes<GEMGeometry, MuonGeometryRecord>(); 
}

#endif


GEMHitAnalyzer::~GEMHitAnalyzer(){}

void
GEMHitAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{


  t_hit = fs->make<TTree>("Hit", "Hit");
  t_hit->Branch("region", &b_region, "region/I");
  t_hit->Branch("station", &b_station, "station/I");
  t_hit->Branch("chamber", &b_chamber, "chamber/I");
  t_hit->Branch("layer", &b_layer, "layer/I");
  t_hit->Branch("etaPartition", &b_etaPartition, "etaPartition/I");
  t_hit->Branch("firstStrip", &b_firstStrip, "firstStrip/I");
  t_hit->Branch("nStrips", &b_nStrips, "nStrips/I");

  /* GEM Geometry */
  edm::ESHandle<GEMGeometry> hGEMGeom;
  hGEMGeom = iSetup.getHandle(hGEMGeom_);
//  iSetup.getByToken(hGEMGeom_, hGEMGeom);

//  edm::ESHandle<GEMGeometry> hGEMGeom;
//  iSetup.get<MuonGeometryRecord>().get(hGEMGeom);
  const GEMGeometry* GEMGeometry_ = &*hGEMGeom;

  edm::Handle<GEMRecHitCollection> gemRecHits;
  iEvent.getByToken(gemRecHits_,gemRecHits);

  edm::Handle<GEMDigiCollection> gemDigis;
  iEvent.getByToken(gemDigis_,gemDigis);
  
  for (const GEMEtaPartition* etaPart : GEMGeometry_->etaPartitions()) {
    GEMDetId etaPartId = etaPart->id();

    int region = etaPartId.region();
    int station = etaPartId.station();
    int chamber = etaPartId.chamber();
    int layer = etaPartId.layer();
    int ieta = etaPartId.ieta();

    b_region = region;
    b_station = station;
    b_chamber = chamber;
    b_layer = layer;
    b_etaPartition = ieta;

    Key4 key4(region, station, chamber, ieta);

    auto range = gemRecHits->get(etaPartId);
    for (auto rechit = range.first; rechit != range.second; ++rechit) {
      auto firstStrip = rechit->firstClusterStrip();
      auto clsSize = rechit->clusterSize();

      b_firstStrip = rechit->firstClusterStrip();
      b_nStrips = rechit->clusterSize();

      for (int i =0; i < clsSize; i++) {
        rechit_occ_[key4]->Fill(firstStrip + i, ieta);
      }
    }

    auto digiRange = gemDigis->get(etaPartId); 
    for (auto digi = digiRange.first; digi != digiRange.second; ++digi) {
      auto strip = digi->strip();
      
      digi_occ_[key4]->Fill(strip, ieta);
    }
    t_hit->Fill();
  }
}

void GEMHitAnalyzer::beginJob(){}
void GEMHitAnalyzer::endJob(){}

void GEMHitAnalyzer::beginRun(const edm::Run& run, const edm::EventSetup& iSetup) { 
  /* GEM Geometry */
  edm::ESHandle<GEMGeometry> hGEMGeom;
  hGEMGeom = iSetup.getHandle(hGEMGeomBeginRun_);

//  iSetup.get<MuonGeometryRecord>().get(hGEMGeom);
  const GEMGeometry* GEMGeometry_ = &*hGEMGeom;

  for (const GEMChamber* chamber : GEMGeometry_->chambers()) {
    GEMDetId cId = chamber->id();
    int re = cId.region();
    int st = cId.station();
    int ch = cId.chamber();
    int la = cId.layer();
    if (st != 1) continue;
    Key4 key4(re, st, ch, la);
    int nEta = 8;
    digi_occ_[key4] = fs->make<TH2D>(Form("digi_occ_GE%d_%d_%d_%d", re, st, ch, la),
                                     Form("Occupancy from Track GE%d %d %d %d", re, st, ch, la),
                                     384, -0.5, 383.5,
                                     nEta, 0.5, nEta+0.5);
    rechit_occ_[key4] = fs->make<TH2D>(Form("rechit_occ_GE%d_%d_%d_%d", re, st, ch, la),
                                       Form("Occupancy from Matched RecHit : GE%d %d %d %d", re, st, ch, la),
                                       384, -0.5, 383.5,
                                       nEta, 0.5, nEta+0.5);
  }
}
void GEMHitAnalyzer::endRun(edm::Run const&, edm::EventSetup const&){}
                   
//define this as a plug-in
DEFINE_FWK_MODULE(GEMHitAnalyzer);
