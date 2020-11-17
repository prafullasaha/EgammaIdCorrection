// -*- C++ -*-
//
// Package:    EgammaIdCorrection/EgammaNtuleWritter
// Class:      EgammaNtuleWritter
//
/**\class EgammaNtuleWritter EgammaNtuleWritter.cc EgammaIdCorrection/EgammaNtuleWritter/plugins/EgammaNtuleWritter.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Prolay Mal
//         Created:  Mon, 14 Sep 2020 14:43:16 GMT
//
//


// system include files
#include <memory>
#include <vector>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
 #include "FWCore/Utilities/interface/InputTag.h"
 #include "DataFormats/TrackReco/interface/Track.h"
 #include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "EgammaIdCorrection/DataFormats/interface/PhotonCorr.h"
#include "TTree.h"


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


class EgammaNtuleWritter : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit EgammaNtuleWritter(const edm::ParameterSet&);
      ~EgammaNtuleWritter();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

  	edm::EDGetTokenT<std::vector<PhotonCorr> > photonCorrToken_;  //used to read from configuration file
  	edm::InputTag photonCorrIT_;
  	edm::Service<TFileService> fs_;
	edm::EDGetTokenT<double> rhoToken_;
	void initEventStructure();
  	TTree *eventTree;
  	TTree *photonTree;
        eventInfo evtInfo;
        photonInfo phoInfo;

};

EgammaNtuleWritter::EgammaNtuleWritter(const edm::ParameterSet& iConfig)
 :
  photonCorrToken_(consumes<std::vector<PhotonCorr> >(iConfig.getParameter<edm::InputTag>("PhoCorrectionContainer"))),
  rhoToken_( consumes<double>( iConfig.getParameter<edm::InputTag>( "rhoFixedGridCollection" ) ) )
{
  photonCorrIT_ = iConfig.getParameter<edm::InputTag>("PhoCorrectionContainer");
  usesResource("TFileService");
  edm::Service<TFileService> m_fileservice;

}


EgammaNtuleWritter::~EgammaNtuleWritter()
{
}

void
EgammaNtuleWritter::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   int n_photon=0;

   edm::Handle<double> rhoHandle;
   iEvent.getByToken( rhoToken_, rhoHandle );
   const double rhoFixedGrd = *( rhoHandle.product() );

   edm::Handle<std::vector<PhotonCorr> > p_vec_photoncorr;
   //iEvent.getByToken(photonCorrToken_, p_vec_photoncorr );
   iEvent.getByLabel(photonCorrIT_, p_vec_photoncorr );
   if (!p_vec_photoncorr.isValid())
     std::cout<<"photonCorr problematic"<<std::endl;
   
   //const auto p_vec_photoncorr = iEvent.get(photonCorrToken_);
   //std::cout<<"p_vec_photoncorr: "<<p_vec_photoncorr.size()<<std::endl;
   //std::vector<PhotonCorr> 
   for(size_t iph=0; iph<(*p_vec_photoncorr).size(); ++iph){
   
   PhotonCorr pho = p_vec_photoncorr->at(iph);

/*        phoInfo.pho_UnCorr_r9 = pho.userFloat("uncorr_r9");
        phoInfo.pho_UnCorr_s4 = pho.userFloat("uncorr_s4");
        phoInfo.pho_UnCorr_sieip = pho.userFloat("uncorr_sieie");
        phoInfo.pho_UnCorr_sieie = pho.userFloat("uncorr_sieip");
        phoInfo.pho_UnCorr_etaWidth = pho.userFloat("uncorr_etaWidth");
        phoInfo.pho_UnCorr_phiWidth = pho.userFloat("uncorr_phiWidth");;
*/        phoInfo.pho_SCeta = pho.superCluster()->eta();
        phoInfo.pho_pt = pho.pt();
        phoInfo.pho_phi = pho.phi();
        phoInfo.pho_rho = rhoFixedGrd;

        if( !iEvent.isRealData() ) {
        phoInfo.pho_Corr_r9 = pho.full5x5_r9();
        phoInfo.pho_Corr_s4 = pho.s4();
        phoInfo.pho_Corr_sieie = pho.full5x5_sigmaIetaIeta();
        phoInfo.pho_Corr_sieip = pho.sieip();
/*        phoInfo.pho_Corr_etaWidth = pho.userFloat("etaWidth");
        phoInfo.pho_Corr_phiWidth = pho.userFloat("phiWidth");
*/        }
        photonTree->Fill();

//std::cout<<"pho etawidth"<<pho.userFloat("etaWidth")<<std::endl;
//std::cout<<"pho uncorr_r9"<<pho.userFloat("uncorr_r9")<<std::endl;
std::cout<<"pho etawidth"<<pho.full5x5_r9()<<std::endl;	
   }
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
}


// ------------ method called once each job just before starting event loop  ------------
void
EgammaNtuleWritter::beginJob()
{
        eventTree = fs_->make<TTree>( "eventTree", "per-event tree" );
        eventTree->Branch( "nphoton", &evtInfo.nphoton, "nphoton/I" );

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

}

void
EgammaNtuleWritter::initEventStructure(){
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
}
// ------------ method called once each job just after ending the event loop  ------------
void
EgammaNtuleWritter::endJob()
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
EgammaNtuleWritter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

}

//define this as a plug-in
DEFINE_FWK_MODULE(EgammaNtuleWritter);
