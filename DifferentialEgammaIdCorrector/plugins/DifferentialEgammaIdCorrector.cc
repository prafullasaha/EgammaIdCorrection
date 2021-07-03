// -*- C++ -*-
//
// Package:    EgammaIdCorrection/DifferentialEgammaIdCorrector
// Class:      DifferentialEgammaIdCorrector
// 
/**\class DifferentialEgammaIdCorrector DifferentialEgammaIdCorrector.cc EgammaIdCorrection/DifferentialEgammaIdCorrector/plugins/DifferentialEgammaIdCorrector.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Prafulla Saha
//         Created:  Tue, 21 Jul 2020 11:04:04 GMT
//
//


// system include files
#include <memory>
#include <iostream>
#include "CLHEP/Random/RandomEngine.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "FWCore/Utilities/interface/RandomNumberGenerator.h"

#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "RecoEgamma/EgammaTools/interface/EffectiveAreas.h"
#include "CommonTools/Utils/interface/StringObjectFunction.h"

#include "EgammaIdCorrection/DifferentialEgammaIdCorrector/interface/GlobalVariablesComputer.h"
#include "EgammaIdCorrection/DifferentialEgammaIdCorrector/interface/MVAComputer.h"
#include "EgammaIdCorrection/DataFormats/interface/PhotonCorr.h"
#include "EgammaIdCorrection/DataFormats/interface/ElectronCorr.h"
#include "EgammaIdCorrection/DifferentialEgammaIdCorrector/interface/EGShowerShapeCorrector.h"
//#include "EgammaIdCorrection/DifferentialEgammaIdCorrector/interface/PhotonShowerShapeCorrector.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "TFormula.h"
#include "TTree.h"

using namespace std;
using namespace edm;
using namespace reco;

class DifferentialEgammaIdCorrector : public edm::EDProducer {
   public:
      explicit DifferentialEgammaIdCorrector(const edm::ParameterSet&);
      ~DifferentialEgammaIdCorrector();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
   private:
      	virtual void produce(edm::Event&, const edm::EventSetup&) override;
	edm::EDGetTokenT<edm::View<pat::Electron>> 	electronToken;
        edm::EDGetTokenT<edm::View<pat::Photon>> 	photonToken;
	edm::EDGetTokenT<pat::ElectronCollection> 	electronsColl_Token;
       	edm::EDGetTokenT<pat::PhotonCollection> 	photonsColl_Token;
        edm::EDGetTokenT<double> rhoToken_;
	edm::EDGetTokenT<View<reco::Vertex> > vertexToken_;
	edm::EDGetTokenT<View<pat::PackedCandidate> > pfcandidateToken_;
        edm::EDGetTokenT<EcalRecHitCollection> ecalHitEBToken_;
        edm::EDGetTokenT<EcalRecHitCollection> ecalHitEEToken_;
        edm::EDGetTokenT<EcalRecHitCollection> ecalHitESToken_;
	EffectiveAreas effectiveAreas_;
//	EGShowerShapeCorrector<reco::Photon, pat::Photon> *photon_correction;	
	EGShowerShapeCorrector<reco::Photon, 		PhotonCorr, 	true,	true> 	*photon_correction;
	EGShowerShapeCorrector<reco::GsfElectron, 	ElectronCorr, 	true, 	false>	*electron_correction;
//	CorrectionEGSSforPhoton photon_correction;
//	PhotonShowerShapeCorrector *photon_correction;
	ConsumesCollector cc_;
        GlobalVariablesComputer globalVariablesComputer_;
	static vector<string> showerShapes_;
	bool correctShowerShapes_;
	bool useVtx0ForNeutralIso_;
};

DifferentialEgammaIdCorrector::DifferentialEgammaIdCorrector(const edm::ParameterSet& iConfig):
	rhoToken_		( consumes<double>( iConfig.getParameter<edm::InputTag>( "rhoFixedGridCollection" ) ) ),
	vertexToken_		( consumes<View<reco::Vertex> >( iConfig.getParameter<InputTag> ( "vertexes" ) ) ),
	pfcandidateToken_	( consumes<View<pat::PackedCandidate> >( iConfig.getParameter<InputTag> ( "PFCandidatesTag" ) ) ),
	ecalHitEBToken_		( consumes<EcalRecHitCollection>( iConfig.getParameter<edm::InputTag>( "reducedBarrelRecHitCollection" ) ) ),
        ecalHitEEToken_		( consumes<EcalRecHitCollection>( iConfig.getParameter<edm::InputTag>( "reducedEndcapRecHitCollection" ) ) ),
        ecalHitESToken_		( consumes<EcalRecHitCollection>( iConfig.getParameter<edm::InputTag>( "reducedPreshowerRecHitCollection" ) ) ),
	effectiveAreas_		((iConfig.getParameter<edm::FileInPath>("effAreasConfigFile")).fullPath()),
	cc_( consumesCollector() ),
        globalVariablesComputer_(iConfig.getParameter<edm::ParameterSet>("globalVariables"), cc_)
{
	electronToken          = consumes<edm::View<pat::Electron>> (iConfig.getParameter<edm::InputTag>("Electrons"));
        photonToken            = consumes<edm::View<pat::Photon>> (iConfig.getParameter<edm::InputTag>("Photons"));

	electronsColl_Token          = consumes<pat::ElectronCollection> (iConfig.getParameter<edm::InputTag>("Electrons"));
	photonsColl_Token            = consumes<pat::PhotonCollection> (std::string("slimmedPhotons"));
	produces<std::vector<PhotonCorr> >( "NewPhotonCorrections" ).setBranchAlias( "NewPhotonCorrections" );
	produces<std::vector<ElectronCorr> >( "NewElectronCorrections" ).setBranchAlias( "NewElectronCorrections" );
//	photon_correction = new PhotonShowerShapeCorrector(iConfig, globalVariablesComputer_);
        photon_correction 	=new EGShowerShapeCorrector<reco::Photon, PhotonCorr, true, true>(iConfig, globalVariablesComputer_);
        electron_correction 	=new EGShowerShapeCorrector<reco::GsfElectron, ElectronCorr, true, false>( iConfig, globalVariablesComputer_);
	useVtx0ForNeutralIso_ 	= iConfig.getParameter<bool>( "useVtx0ForNeutralIso" );
}

DifferentialEgammaIdCorrector::~DifferentialEgammaIdCorrector()
{
}

void
DifferentialEgammaIdCorrector::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

	edm::Handle<edm::View<pat::Electron>> pat_electrons;
	iEvent.getByToken(electronToken, pat_electrons);

	edm::Handle<edm::View<pat::Photon>> pat_photons;
	iEvent.getByToken(photonToken, pat_photons);

	edm::Handle<pat::ElectronCollection> electrons;
    	iEvent.getByToken(electronsColl_Token, electrons);

    	edm::Handle<pat::PhotonCollection> photons;
    	iEvent.getByToken(photonsColl_Token, photons);

        edm::Handle<double> rhoHandle;
        iEvent.getByToken( rhoToken_, rhoHandle );
        const double rhoFixedGrd = *( rhoHandle.product() );

	Handle<View<pat::PackedCandidate> > pfcandidates;
        iEvent.getByToken( pfcandidateToken_, pfcandidates );

        Handle<View<reco::Vertex> > vertices;
        iEvent.getByToken( vertexToken_, vertices );
	const reco::Vertex *neutVtx = ( useVtx0ForNeutralIso_ ? &vertices->at( 0 ) : 0 );	

        edm::Service<edm::RandomNumberGenerator> rng;
	CLHEP::HepRandomEngine& engine = rng->getEngine(iEvent.streamID());
	noZS::EcalClusterLazyTools noZsLazyTool( iEvent, iSetup, ecalHitEBToken_, ecalHitEEToken_, ecalHitESToken_ );

	std::unique_ptr<vector<PhotonCorr> > 	phoCorr( new vector<PhotonCorr> );
	std::unique_ptr<vector<ElectronCorr> > 	eleCorr( new vector<ElectronCorr> );
//=========================================PHOTON========================================================================================
        for (const auto & orig_pho : *photons) {
        pat::Photon pho = orig_pho;
	const reco::CaloClusterPtr  seed_clu  = pho.superCluster()->seed();
        const reco::SuperClusterRef super_clu = pho.superCluster();
        std::vector<float> viCov = noZsLazyTool.localCovariances( *seed_clu );
        float sieip = viCov[1];
        float s4 = noZsLazyTool.e2x2( *seed_clu ) / noZsLazyTool.e5x5( *seed_clu );
std::cout<<"R9 from pat photon before correction	= "<<pho.full5x5_r9()<<std::endl;	
	PhotonCorr pho_corr = PhotonCorr( pho );
	pho_corr.setS4(s4);
        pho_corr.setSieip(sieip);
//cout<<"rho----------->"<<rhoFixedGrd<<endl;
//            if( !iEvent.isRealData() ) {                
	PhotonCorr pho_corr_ = photon_correction->corrector(pho_corr, engine, rhoFixedGrd, noZsLazyTool, pfcandidates->ptrs(), *neutVtx);
std::cout<<"R9 from pat photon after correction	 	= "<<pho_corr_.full5x5_r9()<<std::endl;
std::cout<<"pfPhoIso from pat photon BEFORE correction	= "<<pho_corr_.pfPhoIso03()<<std::endl;
std::cout<<"pfPhoIso from pat photon AFTER correction  	= "<<pho_corr_.pfPhoIso03Corr()<<std::endl;
//	    }
	phoCorr->push_back( pho_corr_ );	
        }//end of Photon loop
//=====================================================================================================================================
//==========================================ELECTRON===================================================================================	
	for (const auto & orig_ele : *electrons) {
	const pat::Electron ele = orig_ele;	
        const reco::CaloClusterPtr  seed_clu  = ele.superCluster()->seed();
        const reco::SuperClusterRef super_clu = ele.superCluster();
        std::vector<float> viCov = noZsLazyTool.localCovariances( *seed_clu );
        float sieip = viCov[1];
        float s4 = noZsLazyTool.e2x2( *seed_clu ) / noZsLazyTool.e5x5( *seed_clu );
std::cout<<"R9 from pat photon before correction        = "<<ele.full5x5_r9()<<std::endl;
        ElectronCorr ele_corr = ElectronCorr( ele );
        ele_corr.setS4(s4);
        ele_corr.setSieip(sieip);		
//        ElectronCorr ele_corr_ = electron_correction->SScorrector(ele_corr, engine, rhoFixedGrd, noZsLazyTool, pfcandidates->ptrs(), *neutVtx);  //Fix which variables to correct etc. before running it
//std::cout<<"R9 from pat photon after correction         = "<<ele_corr_.full5x5_r9()<<std::endl;
	eleCorr->push_back( ele_corr );//Fix Me
	}
//========================================================================================================================================

/*
	for (size_t i = 0; i < electrons->size(); ++i) {
	const pat::Electron & ele = (*electrons)[i];
	std::cout<<"R9-------- = "<<ele.r9()<<std::endl;
	}
*/
	iEvent.put( std::move( phoCorr ), "NewPhotonCorrections" ); 
	iEvent.put( std::move( eleCorr ), "NewElectronCorrections" );
	//iEvent.put( std::move(phoCorr ) );
}
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
DifferentialEgammaIdCorrector::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(DifferentialEgammaIdCorrector);
