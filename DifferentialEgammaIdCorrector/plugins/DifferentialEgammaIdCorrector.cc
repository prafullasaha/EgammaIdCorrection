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

//#include "EgammaIdCorrection/DifferentialEgammaIdCorrector/interface/GlobalVariablesComputer.h"
#include "EgammaIdCorrection/DifferentialEgammaIdCorrector/interface/MVAComputer.h"
#include "EgammaIdCorrection/DataFormats/interface/PhotonCorr.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "TFormula.h"
#include "TTree.h"

using namespace std;
using namespace edm;
using namespace reco;
/*
struct eventInfo {
	int nphoton;
};
struct photonInfo {
	float pho_UnCorr_r9;
	float pho_UnCorr_s4;
	float pho_UnCorr_sieie;
	float pho_UnCorr_sieip;
	float pho_UnCorr_etaWidth;
	float pho_UnCorr_phiWidth;
	float pho_Corr_r9;
	float pho_Corr_s4;
	float pho_Corr_sieie;
	float pho_Corr_sieip;
	float pho_Corr_etaWidth;
	float pho_Corr_phiWidth;
	float pho_pt;
	float pho_SCeta;
	float pho_phi;
	float pho_rho;
};
*/

const bool useXGB=false;
class DifferentialEgammaIdCorrector : public edm::EDProducer {
   public:
      explicit DifferentialEgammaIdCorrector(const edm::ParameterSet&);
      ~DifferentialEgammaIdCorrector();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
   private:
//	virtual void beginJob() override;
      	virtual void produce(edm::Event&, const edm::EventSetup&) override;
//	virtual void endJob() override;
	void correctPhoton(PhotonCorr& ph, CLHEP::HepRandomEngine& engine, double rho, noZS::EcalClusterLazyTools& LazyTool);
	edm::EDGetTokenT<edm::View<pat::Electron>> electronToken;
        edm::EDGetTokenT<edm::View<pat::Photon>> photonToken;
	edm::EDGetTokenT<pat::ElectronCollection> electronsColl_Token;
       	edm::EDGetTokenT<pat::PhotonCollection> photonsColl_Token;
        edm::EDGetTokenT<double> rhoToken_;
        edm::EDGetTokenT<EcalRecHitCollection> ecalHitEBToken_;
        edm::EDGetTokenT<EcalRecHitCollection> ecalHitEEToken_;
        edm::EDGetTokenT<EcalRecHitCollection> ecalHitESToken_;
	EffectiveAreas effectiveAreas_;
	ConsumesCollector cc_;
        GlobalVariablesComputer globalVariablesComputer_;
	static vector<string> showerShapes_;
	bool correctShowerShapes_;
	std::map<string, MVAComputer<Photon, StringObjectFunction<Photon, true>, useXGB> > correctionsEB_;
        map<string, TFormula> correctionScalingsEB_;
        std::map<string, MVAComputer<Photon, StringObjectFunction<Photon, true>, useXGB> > correctionsEE_;
        map<string, TFormula> correctionScalingsEE_;

/*	edm::Service<TFileService> fs_;
	void initEventStructure();
	TTree *eventTree;
	TTree *photonTree;
	eventInfo evtInfo;
	photonInfo phoInfo;
*/      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------
};

vector<string> DifferentialEgammaIdCorrector::showerShapes_ = {"r9", "s4", "sieie", "sieip", "etaWidth", "phiWidth"};

DifferentialEgammaIdCorrector::DifferentialEgammaIdCorrector(const edm::ParameterSet& iConfig):
	rhoToken_( consumes<double>( iConfig.getParameter<edm::InputTag>( "rhoFixedGridCollection" ) ) ),
	ecalHitEBToken_( consumes<EcalRecHitCollection>( iConfig.getParameter<edm::InputTag>( "reducedBarrelRecHitCollection" ) ) ),
        ecalHitEEToken_( consumes<EcalRecHitCollection>( iConfig.getParameter<edm::InputTag>( "reducedEndcapRecHitCollection" ) ) ),
        ecalHitESToken_( consumes<EcalRecHitCollection>( iConfig.getParameter<edm::InputTag>( "reducedPreshowerRecHitCollection" ) ) ),
	effectiveAreas_((iConfig.getParameter<edm::FileInPath>("effAreasConfigFile")).fullPath()),
	cc_( consumesCollector() ),
        globalVariablesComputer_(iConfig.getParameter<edm::ParameterSet>("globalVariables"), cc_)
{
	electronToken          = consumes<edm::View<pat::Electron>> (iConfig.getParameter<edm::InputTag>("Electrons"));
        photonToken            = consumes<edm::View<pat::Photon>> (iConfig.getParameter<edm::InputTag>("Photons"));

//	electronsColl_Token          = consumes<edm::View<pat::ElectronCollection>> (iConfig.getParameter<edm::InputTag>("Electrons"));
//   	photonsColl_Token            = consumes<edm::View<pat::PhotonCollection>> (iConfig.getParameter<edm::InputTag>("Photons"));
	electronsColl_Token          = consumes<pat::ElectronCollection> (iConfig.getParameter<edm::InputTag>("Electrons"));
	photonsColl_Token            = consumes<pat::PhotonCollection> (std::string("slimmedPhotons"));
	produces<std::vector<PhotonCorr> >( "NewPhotonCorrections" ).setBranchAlias( "NewPhotonCorrections" );
//	produces<vector<PhotonCorr> >();
//---Load shower shapes corrections
bool correctShowerShapes_ = true;
        if(correctShowerShapes_)
        {
            for(auto& ss_var : showerShapes_)
            {
                //---EB
                auto xgb_config = iConfig.getParameter<edm::ParameterSet>(ss_var+"_corrector_config_EB"); 
                correctionsEB_[ss_var] = MVAComputer<Photon, StringObjectFunction<Photon, true>, useXGB>(xgb_config, &globalVariablesComputer_);
                correctionScalingsEB_[ss_var] = TFormula("", xgb_config.getParameter<string>("regr_output_scaling").c_str());

                //---EE
                xgb_config = iConfig.getParameter<edm::ParameterSet>(ss_var+"_corrector_config_EE");
                correctionsEE_[ss_var] = MVAComputer<Photon, StringObjectFunction<Photon, true>, useXGB>(xgb_config, &globalVariablesComputer_);
                correctionScalingsEE_[ss_var] = TFormula("", xgb_config.getParameter<string>("regr_output_scaling").c_str());

            }

        }

}

    void DifferentialEgammaIdCorrector::correctPhoton(PhotonCorr& pho, CLHEP::HepRandomEngine& engine, double rho, noZS::EcalClusterLazyTools& lazyTool) 
{
        const auto* corrections = std::abs(pho.superCluster()->eta())<1.5 ? &correctionsEB_ : &correctionsEE_;
        const auto* correctionScalings = std::abs(pho.superCluster()->eta())<1.5 ? &correctionScalingsEB_ : &correctionScalingsEE_;

        const reco::CaloClusterPtr  seed_clu = pho.superCluster()->seed();
        const reco::SuperClusterRef super_clu = pho.superCluster();

        std::vector<float> viCov = lazyTool.localCovariances( *seed_clu );
	float sieip = viCov[1];
	float s4 = lazyTool.e2x2( *seed_clu ) / lazyTool.e5x5( *seed_clu );
//std::cout<<"test[2]S4-------- = "<<s4<<std::endl;
if(!(std::isnan(pho.full5x5_r9()) || std::isnan(s4) || std::isnan(pho.full5x5_sigmaIetaIeta()) || std::isnan(sieip) || std::isnan(pho.superCluster()->etaWidth()) || std::isnan(pho.superCluster()->phiWidth())/* || std::isnan(pho.pfPhoIso03()) || std::isnan(pho.pfChgIsoWrtChosenVtx03()) || std::isnan(pho.pfChgIsoWrtWorstVtx03())*/ || std::isinf(pho.full5x5_r9()) || std::isinf(s4) || std::isinf(pho.full5x5_sigmaIetaIeta()) || std::isinf(sieip) || std::isinf(pho.superCluster()->etaWidth()) || std::isinf(pho.superCluster()->phiWidth())/* || std::isinf(pho.pfPhoIso03()) || std::isinf(pho.pfChgIsoWrtChosenVtx03()) || std::isinf(pho.pfChgIsoWrtWorstVtx03())*/))
        {
            if(correctShowerShapes_)
            {
                reco::Photon::ShowerShape correctedShowerShapes = pho.full5x5_showerShapeVariables();
                pho.addUserFloat("uncorr_r9", pho.full5x5_r9());
                pho.addUserFloat("uncorr_s4", s4);
                pho.addUserFloat("uncorr_sieie", pho.full5x5_sigmaIetaIeta());
                pho.addUserFloat("uncorr_sieip", sieip);
                pho.addUserFloat("uncorr_etaWidth", pho.superCluster()->etaWidth());
                pho.addUserFloat("uncorr_phiWidth", pho.superCluster()->phiWidth());
              
                //---Compute corrections
                // R9 (store it inside e3x3) 
                correctedShowerShapes.e3x3 = (pho.full5x5_r9()+correctionScalings->at("r9").Eval(corrections->at("r9")(pho)[0]))*pho.superCluster()->rawEnergy();                            
                //S4
                auto s4_corr = s4+correctionScalings->at("s4").Eval(corrections->at("s4")(pho)[0]);
                // SiEiE
                correctedShowerShapes.sigmaIetaIeta = pho.full5x5_sigmaIetaIeta()+correctionScalings->at("sieie").Eval(corrections->at("sieie")(pho)[0]);
                // SiEiP
                auto sieip_corr = sieip+correctionScalings->at("sieip").Eval(corrections->at("sieip")(pho)[0]);
                // etaWidth
                pho.addUserFloat("etaWidth", (pho.superCluster()->etaWidth()+correctionScalings->at("etaWidth").Eval(corrections->at("etaWidth")(pho)[0])));
                // phiWidth
                pho.addUserFloat("phiWidth", (pho.superCluster()->phiWidth()+correctionScalings->at("phiWidth").Eval(corrections->at("phiWidth")(pho)[0])));
        
                //---set shower shapes
                pho.setS4(s4_corr);
                pho.setSieip(sieip_corr);
                pho.full5x5_setShowerShapeVariables(correctedShowerShapes);
//std::cout<<"Corrected r9 -------> "<<pho.full5x5_r9()<<endl;
/*		phoInfo.pho_Corr_etaWidth = pho.userFloat("etaWidth");
		phoInfo.pho_Corr_phiWidth = pho.userFloat("phiWidth");
*/
                if (pho.full5x5_r9()<0.)
                    {
                        std::cout << "WARNING: R9<0, Original R9: " << pho.userFloat("uncorr_r9") << " corrected R9: " << pho.full5x5_r9() << std::endl;
                    }
            }
	}
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

        edm::Service<edm::RandomNumberGenerator> rng;
	CLHEP::HepRandomEngine& engine = rng->getEngine(iEvent.streamID());
	noZS::EcalClusterLazyTools noZsLazyTool( iEvent, iSetup, ecalHitEBToken_, ecalHitEEToken_, ecalHitESToken_ );

	std::unique_ptr<vector<PhotonCorr> > phoCorr( new vector<PhotonCorr> );

//Fill event tree
//	evtInfo.nphoton = photons->size();
//	eventTree->Fill();
//Filled event tree
        for (const auto & orig_pho : *photons) {
        pat::Photon pho = orig_pho;
//	for (size_t i = 0; i < photons->size(); ++i) {
//        pat::Photon & pho = (*photons)[i];
//	Ptr<pat::Photon> pho = photons->ptrAt( i );

	const reco::CaloClusterPtr  seed_clu = pho.superCluster()->seed();
        const reco::SuperClusterRef super_clu = pho.superCluster();

        std::vector<float> viCov = noZsLazyTool.localCovariances( *seed_clu );
        float sieip = viCov[1];
        float s4 = noZsLazyTool.e2x2( *seed_clu ) / noZsLazyTool.e5x5( *seed_clu );

	

//std::cout<<"R9 from pat photon--------> = "<<pho.full5x5_r9()<<std::endl;
	
	PhotonCorr pho_corr = PhotonCorr( pho );
	pho_corr.setS4(s4);
        pho_corr.setSieip(sieip);

//std::cout<<"R9 from new photon object------> = "<<pho_corr.full5x5_r9()<<std::endl;
//std::cout<<"Uncorrected etaWidth ------> = "<<pho.superCluster()->etaWidth()<<std::endl;
/*
        phoInfo.pho_UnCorr_r9 = pho.full5x5_r9();
	phoInfo.pho_UnCorr_s4 = s4;
        phoInfo.pho_UnCorr_sieip = sieip;
        phoInfo.pho_UnCorr_sieie = pho.full5x5_sigmaIetaIeta();
        phoInfo.pho_UnCorr_etaWidth = pho.superCluster()->etaWidth();
        phoInfo.pho_UnCorr_phiWidth = pho.superCluster()->phiWidth();
	phoInfo.pho_SCeta = pho.superCluster()->eta();
	phoInfo.pho_pt = pho.pt();
	phoInfo.pho_phi = pho.phi();
	phoInfo.pho_rho = rhoFixedGrd;
*/
//cout<<"rho----------->"<<rhoFixedGrd<<endl;
            if( !iEvent.isRealData() ) {
                
                //---compute correction
                correctPhoton(pho_corr, engine, rhoFixedGrd, noZsLazyTool);
	    }
//std::cout<<"SuperCluster eta from new photon object------> = "<<pho_corr.superCluster()->eta()<<std::endl;
//std::cout<<"etaWidth from new photon object------> = "<<pho_corr.userFloat("etaWidth")<<std::endl;
//Fill photon tree-------------------------------------------------------------------------------
/*	if( !iEvent.isRealData() ) {
	phoInfo.pho_Corr_r9 = pho_corr.full5x5_r9();
	phoInfo.pho_Corr_s4 = pho_corr.s4();
	phoInfo.pho_Corr_sieie = pho_corr.full5x5_sigmaIetaIeta();
	phoInfo.pho_Corr_sieip = pho_corr.sieip();
//	phoInfo.pho_Corr_etaWidth = pho_corr.superCluster()->etaWidth();
//	phoInfo.pho_Corr_phiWidth = pho_corr.superCluster()->phiWidth();
	}
	photonTree->Fill();
*/	
//------------------------------------------------------------------------------------------------
	phoCorr->push_back( pho_corr );	
        }//end of Photon loop
	
	for (const auto & orig_ele : *electrons) {
	const pat::Electron ele = orig_ele;
	}
/*
	for (size_t i = 0; i < electrons->size(); ++i) {
	const pat::Electron & ele = (*electrons)[i];
	std::cout<<"R9-------- = "<<ele.r9()<<std::endl;
	}
*/
	iEvent.put( std::move( phoCorr ), "NewPhotonCorrections" ); 
	//iEvent.put( std::move(phoCorr ) );
}
/*
void
DifferentialEgammaIdCorrector::beginJob()
{
//per-event tree----------------------------------------------------------------------------------
	eventTree = fs_->make<TTree>( "eventTree", "per-event tree" );
	eventTree->Branch( "nphoton", &evtInfo.nphoton, "nphoton/I" );
//per-photon tree
	photonTree = fs_->make<TTree>( "photonTree", "per-photon tree" );
	photonTree->Branch( "pho_UnCorr_r9", &phoInfo.pho_UnCorr_r9, "pho_UnCorr_r9/F" );
	photonTree->Branch( "pho_UnCorr_s4", &phoInfo.pho_UnCorr_s4, "pho_UnCorr_s4/F" );
	photonTree->Branch( "pho_UnCorr_sieie", &phoInfo.pho_UnCorr_sieie, "pho_UnCorr_sieie/F" );
	photonTree->Branch( "pho_UnCorr_sieip", &phoInfo.pho_UnCorr_sieip, "pho_UnCorr_sieip/F" );
	photonTree->Branch( "pho_UnCorr_etaWidth", &phoInfo.pho_UnCorr_etaWidth, "pho_UnCorr_etaWidth/F" );
	photonTree->Branch( "pho_UnCorr_phiWidth", &phoInfo.pho_UnCorr_phiWidth, "pho_UnCorr_phiWidth/F" );
	photonTree->Branch( "pho_Corr_r9", &phoInfo.pho_Corr_r9, "pho_Corr_r9/F" );
        photonTree->Branch( "pho_Corr_s4", &phoInfo.pho_Corr_s4, "pho_Corr_s4/F" );
        photonTree->Branch( "pho_Corr_sieie", &phoInfo.pho_Corr_sieie, "pho_Corr_sieie/F" );
        photonTree->Branch( "pho_Corr_sieip", &phoInfo.pho_Corr_sieip, "pho_Corr_sieip/F" );
        photonTree->Branch( "pho_Corr_etaWidth", &phoInfo.pho_Corr_etaWidth, "pho_Corr_etaWidth/F" );
        photonTree->Branch( "pho_Corr_phiWidth", &phoInfo.pho_Corr_phiWidth, "pho_Corr_phiWidth/F" );
	photonTree->Branch( "pho_SCeta", &phoInfo.pho_SCeta, "pho_SCeta/F" );
	photonTree->Branch( "pho_pt", &phoInfo.pho_pt, "pho_pt/F" );
	photonTree->Branch( "pho_phi", &phoInfo.pho_phi, "pho_phi/F" );
	photonTree->Branch( "pho_rho", &phoInfo.pho_rho, "pho_rho/F" );
//End tree ----------------------------------------------------------------------------------------
}

void
DifferentialEgammaIdCorrector::initEventStructure(){
//per-event tree
	evtInfo.nphoton = 0;
//per-photon tree
	phoInfo.pho_UnCorr_r9 = -999;
	phoInfo.pho_UnCorr_s4 = -999;
	phoInfo.pho_UnCorr_sieie = -999;
	phoInfo.pho_UnCorr_sieip = -999;
	phoInfo.pho_UnCorr_etaWidth = -999;
	phoInfo.pho_UnCorr_phiWidth = -999;
	phoInfo.pho_SCeta = -999;
	phoInfo.pho_pt = -999;
	phoInfo.pho_phi = -999;
	phoInfo.pho_rho = -999;
	phoInfo.pho_Corr_r9 = -999;
        phoInfo.pho_Corr_s4 = -999;
        phoInfo.pho_Corr_sieie = -999;
        phoInfo.pho_Corr_sieip = -999;
        phoInfo.pho_Corr_etaWidth = -999;
        phoInfo.pho_Corr_phiWidth = -999;

};
*/
// ------------ method called once each stream before processing any runs, lumis or events  ------------
/*void
DifferentialEgammaIdCorrector::beginStream(edm::StreamID)
{
}
*/
// ------------ method called once each stream after processing all runs, lumis and events  ------------
/*
void
DifferentialEgammaIdCorrector::endStream() {
}
*/
// ------------ method called when starting to processes a run  ------------
/*
void
DifferentialEgammaIdCorrector::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
DifferentialEgammaIdCorrector::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
DifferentialEgammaIdCorrector::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
DifferentialEgammaIdCorrector::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
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
