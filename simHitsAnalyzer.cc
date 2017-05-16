
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h" 
#include "FWCore/Framework/interface/EventSetup.h" 
#include "FWCore/Framework/interface/ESHandle.h" 
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h" 
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

//#include "DataFormats/SiPixelDetId/interface/PixelEndCapName.h" 
//#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"
#include "DataFormats/DetId/interface/DetId.h" 
#include "SimDataFormats/TrackingHit/interface/PSimHit.h" 
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h" 
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h" 
#include "Geometry/CommonDetUnit/interface/GeomDetUnit.h" 
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h" 
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"
#include "DataFormats/SiPixelDetId/interface/PXBDetId.h"
#include "DataFormats/SiPixelDetId/interface/PXFDetId.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"


#include <TH1.h>


class simHitsAnalyzer : public edm::EDAnalyzer{
	public: 
		explicit simHitsAnalyzer(const edm::ParameterSet&);
		~simHitsAnalyzer();

	private:
		virtual void beginJob();
		virtual void analyze(const edm::Event&, const edm::EventSetup&);
		virtual void endJob();
		int getLayerId(const PSimHit & simHit, bool isBarrel=true);


		edm::EDGetTokenT<TrackingParticleCollection> tpHitsSrc_;
		edm::EDGetTokenT<edm::PSimHitContainer> hitContainer_;
		const TrackerGeometry * geo_;
		//	edm::ESHandle<ParticleDataTable> pdt;

		enum { 
			BPix1 = 0, BPix2= 1, BPix3=2, BPix4=3,
			FPix1_neg=4, FPix2_neg=5,
			FPix1_pos=6, FPix2_pos=7,
			nLayers = 8
		};



		TH1I *h_hits_1;
		TH1I *h_hits_2;
		std::string SimHitLabel;
};

simHitsAnalyzer::simHitsAnalyzer(const edm::ParameterSet& iConfig):
	tpHitsSrc_(consumes<TrackingParticleCollection>(iConfig.getParameter<edm::InputTag>("tpHitsSrc"))),
	hitContainer_(consumes<edm::PSimHitContainer>(iConfig.getParameter<edm::InputTag>("hitSrc"))){
		edm::Service<TFileService> fs;
		h_hits_1 = fs->make<TH1I>("h_hits_1","Sim Hits",10, -0.5, 9.5);
		h_hits_2 = fs->make<TH1I>("h_hits_2","Sim Hits",10, -0.5, 9.5);
	}

simHitsAnalyzer::~simHitsAnalyzer()
{
}

void simHitsAnalyzer::beginJob(){
}


void simHitsAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){

	using namespace edm;

	edm::ESHandle<TrackerGeometry> tGeo;
	//	iSetup.get<TrackerDigiGeometryRecord>().get(tGeo);
	//
	geo_ = tGeo.product();
	//	//iSetup.getData(pdt);
	//
	edm::Handle<TrackingParticleCollection> tpCollection;
	iEvent.getByToken(tpHitsSrc_, tpCollection);

	std::vector<PSimHit> simHits;
	edm::Handle<PSimHitContainer> PixelBarrelHitsLowTof;
	edm::Handle<PSimHitContainer> PixelBarrelHitsHighTof;
	iEvent.getByToken(hitContainer_,PixelBarrelHitsLowTof);
	iEvent.getByToken(hitContainer_,PixelBarrelHitsHighTof);
	simHits.insert(simHits.end(), PixelBarrelHitsLowTof->begin(), PixelBarrelHitsLowTof->end()); 
	//	simHits.insert(simHits.end(), PixelBarrelHitsHighTof->begin(), PixelBarrelHitsHighTof->end()); 


	struct trkHit{
		unsigned int trackid;
		int hitPid;
		int nhit;
	};
	std::vector<trkHit> track_hits ; 
	long int idIndx;
	for(auto iHit= simHits.begin(); iHit!=simHits.end(); ++iHit){
		idIndx= -1;
		for(int i=0; i<int(track_hits.size()); ++i){ 
			if((*iHit).trackId() == (track_hits.at(i)).trackid \
					&& (*iHit).particleType() == (track_hits.at(i)).hitPid) {
				idIndx = i;
				break;
			}
		}
		if(idIndx ==-1) {
			trkHit newone;
			newone.trackid =(*iHit).trackId();
			newone.hitPid =(*iHit).particleType();
			newone.nhit =1;
			track_hits.push_back(newone);
		}
		else {
			(track_hits.at(idIndx )).nhit ++;
		}
	}
	for(auto i =track_hits.begin(); i!= track_hits.end(); ++i){
		h_hits_1->Fill((*i).nhit);
		std::cout<<(*i).nhit<<std::endl;
		std::cout<<(*i).hitPid<<std::endl;
	}
	simHits.clear();
	track_hits.clear();

	//bool haveFirstLayer = false;
	/*
	   std::vector<std::pair<unsigned int , unsigned int>> track_hits ; 
	   long int idIndx;
	   for(auto iHit= simHits.begin(); iHit!=simHits.end(); ++iHit){
	   idIndx= -1;
	   for(int i=0; i<int(track_hits.size()); ++i){ 
	   if((*iHit).trackId()== track_hits.at(i).first) {
	   idIndx = i;
	   break;
	   }
	   }
	   if(idIndx ==-1) {
	   track_hits.push_back(std::make_pair ( (*iHit).trackId(),1));
	   }
	   else {
	   (track_hits.at(idIndx )).second ++;
	   }
	   }
	   for(auto i =track_hits.begin(); i!= track_hits.end(); ++i){
	   h_hits_1->Fill((*i).second);
	   std::cout<<(*i).second<<std::endl;
	   }
	   simHits.clear();
	   track_hits.clear();
	   */

}

void simHitsAnalyzer::endJob(){
}

int simHitsAnalyzer::getLayerId(const PSimHit & simHit, bool isBarrel){
	unsigned int id = simHit.detUnitId();

	if(isBarrel){
		if(geo_->idToDetUnit(id)->subDetector() ==
				GeomDetEnumerators::PixelBarrel)
		{
			PXBDetId pid(id);
			return pid.layer() - 1; // 0, 1, 2
		}
	}
	else {
		if(geo_->idToDetUnit(id)->subDetector() ==
				GeomDetEnumerators::PixelEndcap)
		{
			PXFDetId pid(id);
			return BPix3 + ((pid.side()-1) << 1) + pid.disk(); // 3 -
		}
	}
	return -1;
}
DEFINE_FWK_MODULE(simHitsAnalyzer);
