
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
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"
#include "DataFormats/SiPixelDetId/interface/PXBDetId.h"
#include "DataFormats/SiPixelDetId/interface/PXFDetId.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
#include "DataFormats/SiStripDetId/interface/StripSubdetector.h"


#include <TH1.h>
#include <TH3F.h>
#include <TNtuple.h>


class simHitsAnalyzer : public edm::EDAnalyzer{
	public: 
		explicit simHitsAnalyzer(const edm::ParameterSet&);
		~simHitsAnalyzer();

	private:
		virtual void beginJob();
		virtual void analyze(const edm::Event&, const edm::EventSetup&);
		virtual void endJob();
		//	int getLayerId(const PSimHit & simHit, bool isBarrel=true);
		//int layer(DetId & id, bool isBarrel=true);
		bool allowDifferentProcessTypeForDifferentDetectors_;


		edm::EDGetTokenT<TrackingParticleCollection> tpHitsSrc_;
		edm::EDGetTokenT<edm::PSimHitContainer> hitHighTofContainer_;
		edm::EDGetTokenT<edm::PSimHitContainer> hitLowTofContainer_;
		edm::EDGetTokenT<edm::PSimHitContainer> hitEndHighTofContainer_;
		edm::EDGetTokenT<edm::PSimHitContainer> hitEndLowTofContainer_;
		edm::EDGetTokenT<edm::PSimHitContainer> hitTECHighTofContainer_;
		edm::EDGetTokenT<edm::PSimHitContainer> hitTECLowTofContainer_;
		edm::EDGetTokenT<edm::PSimHitContainer> hitTIBHighTofContainer_;
		edm::EDGetTokenT<edm::PSimHitContainer> hitTIBLowTofContainer_;
		edm::EDGetTokenT<edm::PSimHitContainer> hitTIDHighTofContainer_;
		edm::EDGetTokenT<edm::PSimHitContainer> hitTIDLowTofContainer_;
		edm::EDGetTokenT<edm::PSimHitContainer> hitTOBHighTofContainer_;
		edm::EDGetTokenT<edm::PSimHitContainer> hitTOBLowTofContainer_;
		//		edm::EDGetTokenT<edm::SimTrackContainer> simTrackContainer_;
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
		TH1I *h_hits_3;
		TH1F *hpt;
		TH3F *hdist;
		TNtuple* tnp;
		std::string SimHitLabel;
};

simHitsAnalyzer::simHitsAnalyzer(const edm::ParameterSet& iConfig):
	allowDifferentProcessTypeForDifferentDetectors_(iConfig.getParameter<bool>("allowDifferentSimHitProcesses") ),
	tpHitsSrc_(consumes<TrackingParticleCollection>(iConfig.getParameter<edm::InputTag>("tpHitsSrc"))),
	hitHighTofContainer_(consumes<edm::PSimHitContainer>(iConfig.getParameter<edm::InputTag>("hitHighTofSrc"))),
	hitLowTofContainer_(consumes<edm::PSimHitContainer>(iConfig.getParameter<edm::InputTag>("hitLowTofSrc"))),
	hitEndHighTofContainer_(consumes<edm::PSimHitContainer>(iConfig.getParameter<edm::InputTag>("hitEndHighTofSrc"))),
	hitEndLowTofContainer_ (consumes<edm::PSimHitContainer>(iConfig.getParameter<edm::InputTag>("hitEndLowTofSrc"))),
	hitTECHighTofContainer_(consumes<edm::PSimHitContainer>(iConfig.getParameter<edm::InputTag>("hitTECHighTofSrc"))),
	hitTECLowTofContainer_ (consumes<edm::PSimHitContainer>(iConfig.getParameter<edm::InputTag>("hitTECLowTofSrc"))),
	hitTIBHighTofContainer_(consumes<edm::PSimHitContainer>(iConfig.getParameter<edm::InputTag>("hitTIBHighTofSrc"))),
	hitTIBLowTofContainer_ (consumes<edm::PSimHitContainer>(iConfig.getParameter<edm::InputTag>("hitTIBLowTofSrc"))),
	hitTIDHighTofContainer_(consumes<edm::PSimHitContainer>(iConfig.getParameter<edm::InputTag>("hitTIDHighTofSrc"))),
	hitTIDLowTofContainer_ (consumes<edm::PSimHitContainer>(iConfig.getParameter<edm::InputTag>("hitTIDLowTofSrc"))),
	hitTOBHighTofContainer_(consumes<edm::PSimHitContainer>(iConfig.getParameter<edm::InputTag>("hitTOBHighTofSrc"))),
	hitTOBLowTofContainer_ (consumes<edm::PSimHitContainer>(iConfig.getParameter<edm::InputTag>("hitTOBLowTofSrc")))
	//simTrackContainer_(consumes<edm::SimTrackContainer>(iConfig.getParameter<edm::InputTag>("simTrackSrc")))
{
	edm::Service<TFileService> fs;
	h_hits_1 = fs->make<TH1I>("h_hits_1","Sim Hits",10, -0.5, 9.5);
	h_hits_2 = fs->make<TH1I>("h_hits_2","Sim Hits",10, -0.5, 9.5);
	h_hits_3 = fs->make<TH1I>("h_hits_3","Sim Hits",10, -0.5, 9.5);
	hpt = fs->make<TH1F>("hpt","tp pT",200, 0,50);
	hdist = fs->make<TH3F>("hdist","Hits",100, -500. , 500., 100,-500,500, 100, -500,500 );
	tnp = fs->make<TNtuple>("tnp","hits_tnp","x:y:z" );
}

simHitsAnalyzer::~simHitsAnalyzer()
{
}

void simHitsAnalyzer::beginJob(){
}


void simHitsAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){

	using namespace edm;

	edm::ESHandle<TrackerGeometry> tGeo;
	//iSetup.get<TrackerDigiGeometryRecord>().get(tGeo);
	//
	geo_ = tGeo.product();
	//iSetup.getData(pdt);
	//
	edm::Handle<TrackingParticleCollection> tpCollection;
	edm::Handle<PSimHitContainer> PixelBarrelHitsHighTof;
	edm::Handle<PSimHitContainer> PixelBarrelHitsLowTof;
	edm::Handle<PSimHitContainer> PixelEndcapHitsHighTof;
	edm::Handle<PSimHitContainer> PixelEndcapHitsLowTof;
	//edm::Handle<PSimHitContainer> SiTECHitsHighTof;
	//edm::Handle<PSimHitContainer> SiTECHitsLowTof;
	//edm::Handle<PSimHitContainer> SiTIBHitsHighTof;
	//edm::Handle<PSimHitContainer> SiTIBHitsLowTof;
	//edm::Handle<PSimHitContainer> SiTIDHitsHighTof;
	//edm::Handle<PSimHitContainer> SiTIDHitsLowTof;
	//edm::Handle<PSimHitContainer> SiTOBHitsHighTof;
	//edm::Handle<PSimHitContainer> SiTOBHitsLowTof;
	//	edm::Handle<SimTrackContainer> simTrackCollection;

	iEvent.getByToken(tpHitsSrc_, tpCollection);
	iEvent.getByToken(hitHighTofContainer_,PixelBarrelHitsHighTof);
	iEvent.getByToken(hitLowTofContainer_, PixelBarrelHitsLowTof);
	iEvent.getByToken(hitEndHighTofContainer_,PixelEndcapHitsHighTof);
	iEvent.getByToken(hitEndLowTofContainer_, PixelEndcapHitsLowTof);
	//iEvent.getByToken(hitTECHighTofContainer_,SiTECHitsHighTof);
	//iEvent.getByToken(hitTECLowTofContainer_, SiTECHitsLowTof);
	//iEvent.getByToken(hitTIBHighTofContainer_,SiTIBHitsHighTof);
	//iEvent.getByToken(hitTIBLowTofContainer_, SiTIBHitsLowTof);
	//iEvent.getByToken(hitTIDHighTofContainer_,SiTIDHitsHighTof);
	//iEvent.getByToken(hitTIDLowTofContainer_, SiTIDHitsLowTof);
	//iEvent.getByToken(hitTOBHighTofContainer_,SiTOBHitsHighTof);
	//iEvent.getByToken(hitTOBLowTofContainer_, SiTOBHitsLowTof);
	//	iEvent.getByToken(simTrackContainer_, simTrackCollection);


	std::vector<PSimHit> simHits;
	simHits.insert(simHits.end(), PixelBarrelHitsHighTof->begin(), PixelBarrelHitsHighTof->end()); 
	simHits.insert(simHits.end(), PixelBarrelHitsLowTof->begin(),  PixelBarrelHitsLowTof->end()); 
	simHits.insert(simHits.end(), PixelEndcapHitsHighTof->begin(), PixelEndcapHitsHighTof->end()); 
	simHits.insert(simHits.end(), PixelEndcapHitsLowTof->begin(),  PixelEndcapHitsLowTof->end()); 
	//simHits.insert(simHits.end(), SiTIBHitsHighTof->begin(), SiTIBHitsHighTof->end()); 
	//simHits.insert(simHits.end(), SiTIBHitsLowTof->begin(),  SiTIBHitsLowTof->end()); 
	//simHits.insert(simHits.end(), SiTOBHitsHighTof->begin(), SiTOBHitsHighTof->end()); 
	//simHits.insert(simHits.end(), SiTOBHitsLowTof->begin(),  SiTOBHitsLowTof->end()); 
	//simHits.insert(simHits.end(), SiTIDHitsHighTof->begin(), SiTIDHitsHighTof->end()); 
	//simHits.insert(simHits.end(), SiTIDHitsLowTof->begin(),  SiTIDHitsLowTof->end()); 
	//simHits.insert(simHits.end(), SiTECHitsHighTof->begin(), SiTECHitsHighTof->end()); 
	//simHits.insert(simHits.end(), SiTECHitsLowTof->begin(),  SiTECHitsLowTof->end()); 

	//	std::vector<trkHit> track_hits ; 
	std::multimap<unsigned int, size_t> trackId_hitsId; 
	//build the PSimHits index
	for(size_t i =0;i<simHits.size();++i){
		trackId_hitsId.insert( std::make_pair( (simHits.at(i)).trackId(),i) );
	}

	//loop over TrackingParticle to pick up the SimTrack in each class
	for(TrackingParticleCollection::size_type i=0; i<tpCollection->size(); i++){
		TrackingParticleRef tpr(tpCollection, i);
		TrackingParticle* tp=const_cast<TrackingParticle*>(tpr.get());

		//TrackingParticle selections:
		if(tp->status() < 0 || tp->charge()==0 || tp->pt()<1 || TMath::Abs(tp->eta())>2.4) continue;

		const std::vector<SimTrack>& simTrackVec = tp->g4Tracks();
		int pdgId = tp->pdgId();

		for(size_t j=0; j<simTrackVec.size(); ++j){
			const SimTrack & simTrack = simTrackVec.at(j);
			if( simTrack.type()!= pdgId)continue;
			//int pdgId = simTrack.type();
			//	size_t matchedHits=0; // As far as I can tell, this is the number of tracker layers with hits on them, i.e. taking account of overlaps.
			size_t numberOfHits=0; // The total number of hits not taking account of overlaps.
			size_t numberOfTrackerHits=0;
			bool init=true;
			int processType=0;
			//int particleType=0;
			double tof =0;
			//int oldLayer = 0;
			//unsigned int newLayer = 0;
			DetId oldDetector;
			DetId newDetector;

			/* find out the hit with the longest tof
			*/
			for( auto iHitIndex=trackId_hitsId.lower_bound( simTrack.trackId() ), end=trackId_hitsId.upper_bound( simTrack.trackId() );
					iHitIndex!=end;
					++iHitIndex ){
/*
*/
				const auto & pSimHit=simHits[ iHitIndex->second ];
				if(pSimHit.particleType() != pdgId) continue;
				if( init ) {
					tof = pSimHit.timeOfFlight();
					processType =pSimHit.processType();
					newDetector=DetId( pSimHit.detUnitId() );
					init = false;
				}
				else if ( tof< pSimHit.timeOfFlight()){
					tof = pSimHit.timeOfFlight();
					processType =pSimHit.processType();
					newDetector=DetId( pSimHit.detUnitId() );
				}
			}
			for( auto iHitIndex=trackId_hitsId.lower_bound( simTrack.trackId() ), end=trackId_hitsId.upper_bound( simTrack.trackId() );
					iHitIndex!=end;
					++iHitIndex ){
				const auto & pSimHit=simHits[ iHitIndex->second ];
				if(pSimHit.particleType() != pdgId) continue;
				//	if( init ) {
				//		processType =pSimHit.processType();
				//		newDetector=DetId( pSimHit.detUnitId() );
				//	}

				oldDetector=newDetector;
				newDetector=DetId( pSimHit.detUnitId() );

				//	if( allowDifferentProcessTypeForDifferentDetectors_ && newDetector.det()!=oldDetector.det() ) processType=pSimHit.processType();

				if( processType==pSimHit.processType() )
				{
					++numberOfHits;
					//		oldLayer=newLayer;
					//		newLayer=0;
					if( newDetector.det() == DetId::Tracker ) {
						++numberOfTrackerHits;

						//	newLayer=tTopo->layer( newDetector );
						//	newLayer=layer( newDetector );

						//	// Count hits using layers for glued detectors
						//	if( (oldLayer!=newLayer || (oldLayer==newLayer && oldDetector.subdetId()!=newDetector.subdetId())) ) ++matchedHits;
					}
				}
				tnp->Fill((pSimHit.entryPoint()).x(), (pSimHit.entryPoint()).y(),(pSimHit.entryPoint()).z());
			}//end loop over PSimHit for this sim Track.
			if(j==0) h_hits_1->Fill(numberOfHits); //1*
		}
		//	h_hits_1->Fill(tp->numberOfHits());
		h_hits_2->Fill(tp->numberOfTrackerHits());
		hpt ->Fill(tp->pt());
	}
std::cout<<tpCollection->size()<<std::endl;
}

void simHitsAnalyzer::endJob(){
}
/*
int simHitsAnalyzer::layer(DetId & id, bool isBarrel){
	if(isBarrel) {
		return int((id.rawId()>>pbVals_.layerStartBit_) & pbVals_.layerMask_);
	}
	else {
		return int((id.rawId()>>pfVals_.diskStartBit_) & pfVals_.diskMask_);
	}
}
*/
DEFINE_FWK_MODULE(simHitsAnalyzer);
