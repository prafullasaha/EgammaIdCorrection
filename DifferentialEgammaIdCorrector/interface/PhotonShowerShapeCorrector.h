#ifndef PhotonShowerShapeCorrector_h
#define PhotonShowerShapeCorrector_h

#include <iostream>
#include <tuple>
#include <vector>
#include <string>
#include "TFormula.h"
#include "CLHEP/Random/RandomEngine.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "EgammaIdCorrection/DifferentialEgammaIdCorrector/interface/MVAComputer.h"
#include "EgammaIdCorrection/DifferentialEgammaIdCorrector/interface/GlobalVariablesComputer.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"

using namespace std;

//using namespace reco;


class PhotonShowerShapeCorrector
{
public:
	PhotonShowerShapeCorrector(){};
	PhotonShowerShapeCorrector( const edm::ParameterSet &iConfig, GlobalVariablesComputer global_var );
        ~PhotonShowerShapeCorrector();
//        const static bool useXGB = true;
	pat::Photon corrector( pat::Photon &obj, CLHEP::HepRandomEngine& engine, double rho, noZS::EcalClusterLazyTools& LazyTool);
private:
	static vector<string> showerShapes_;
        const static bool useXGB = true;
//	edm::ConsumesCollector cc_;
//        GlobalVariablesComputer globalVariablesComputer_;
	std::map<string, MVAComputer<reco::Photon, StringObjectFunction<reco::Photon, true>, useXGB> > correctionsEB_;
        map<string, TFormula> correctionScalingsEB_;
        std::map<string, MVAComputer<reco::Photon, StringObjectFunction<reco::Photon, true>, useXGB> > correctionsEE_;
        map<string, TFormula> correctionScalingsEE_;
	bool correctShowerShapes_ = true;
};

vector<string> PhotonShowerShapeCorrector::showerShapes_ = {"r9", "s4", "sieie", "sieip", "etaWidth", "phiWidth"};

PhotonShowerShapeCorrector::PhotonShowerShapeCorrector( const edm::ParameterSet &iConfig, GlobalVariablesComputer global_var )
{
        if(correctShowerShapes_)
        {
            for(auto& ss_var : showerShapes_)
            {
		//---EB
                auto xgb_config = iConfig.getParameter<edm::ParameterSet>(ss_var+"_corrector_config_EB");
                correctionsEB_[ss_var] = MVAComputer<reco::Photon, StringObjectFunction<reco::Photon, true>, useXGB>(xgb_config, &global_var);
                correctionScalingsEB_[ss_var] = TFormula("", xgb_config.getParameter<string>("regr_output_scaling").c_str());
		//---EE
                xgb_config = iConfig.getParameter<edm::ParameterSet>(ss_var+"_corrector_config_EE");
                correctionsEE_[ss_var] = MVAComputer<reco::Photon, StringObjectFunction<reco::Photon, true>, useXGB>(xgb_config, &global_var);
                correctionScalingsEE_[ss_var] = TFormula("", xgb_config.getParameter<string>("regr_output_scaling").c_str());

            }
	}
}


PhotonShowerShapeCorrector::~PhotonShowerShapeCorrector() {}

     	pat::Photon PhotonShowerShapeCorrector::corrector( pat::Photon &obj , CLHEP::HepRandomEngine& engine, double rho, noZS::EcalClusterLazyTools& lazyTool)
{
        const auto* corrections = std::abs(obj.superCluster()->eta())<1.5 ? &correctionsEB_ : &correctionsEE_;
        const auto* correctionScalings = std::abs(obj.superCluster()->eta())<1.5 ? &correctionScalingsEB_ : &correctionScalingsEE_;

        const reco::CaloClusterPtr  seed_clu = obj.superCluster()->seed();
        const reco::SuperClusterRef super_clu = obj.superCluster();

        std::vector<float> viCov = lazyTool.localCovariances( *seed_clu );
        float sieip = viCov[1];
        float s4 = lazyTool.e2x2( *seed_clu ) / lazyTool.e5x5( *seed_clu );
if(!(std::isnan(obj.full5x5_r9()) || std::isnan(s4) || std::isnan(obj.full5x5_sigmaIetaIeta()) || std::isnan(sieip) || std::isnan(obj.superCluster()->etaWidth()) || std::isnan(obj.superCluster()->phiWidth())
// || std::isnan(obj.pfPhoIso03()) || std::isnan(obj.pfChgIsoWrtChosenVtx03()) || std::isnan(obj.pfChgIsoWrtWorstVtx03())
 || std::isinf(obj.full5x5_r9()) || std::isinf(s4) || std::isinf(obj.full5x5_sigmaIetaIeta()) || std::isinf(sieip) || std::isinf(obj.superCluster()->etaWidth()) || std::isinf(obj.superCluster()->phiWidth())
// || std::isinf(obj.pfPhoIso03()) || std::isinf(obj.pfChgIsoWrtChosenVtx03()) || std::isinf(obj.pfChgIsoWrtWorstVtx03())
))
        {
            if(correctShowerShapes_)
            {
//                reco::Photon::ShowerShape correctedShowerShapes = obj.full5x5_showerShapeVariables();
		auto correctedShowerShapes = obj.full5x5_showerShapeVariables();
                obj.addUserFloat("uncorr_r9", obj.full5x5_r9());
                obj.addUserFloat("uncorr_s4", s4);
                obj.addUserFloat("uncorr_sieie", obj.full5x5_sigmaIetaIeta());
                obj.addUserFloat("uncorr_sieip", sieip);
                obj.addUserFloat("uncorr_etaWidth", obj.superCluster()->etaWidth());
                obj.addUserFloat("uncorr_phiWidth", obj.superCluster()->phiWidth());

cout<<"[DEBUG]Hi I am here in corrector function"<<endl;

                //---Compute corrections
                // R9 (store it inside e3x3) 
                correctedShowerShapes.e3x3 = (obj.full5x5_r9()+correctionScalings->at("r9").Eval(corrections->at("r9")(obj)[0]))*obj.superCluster()->rawEnergy();
                //S4
                auto s4_corr = s4+correctionScalings->at("s4").Eval(corrections->at("s4")(obj)[0]);
                // SiEiE
                correctedShowerShapes.sigmaIetaIeta = obj.full5x5_sigmaIetaIeta()+correctionScalings->at("sieie").Eval(corrections->at("sieie")(obj)[0]);
                auto sieip_corr = sieip+correctionScalings->at("sieip").Eval(corrections->at("sieip")(obj)[0]);
                obj.addUserFloat("etaWidth", (obj.superCluster()->etaWidth()+correctionScalings->at("etaWidth").Eval(corrections->at("etaWidth")(obj)[0])));
                obj.addUserFloat("phiWidth", (obj.superCluster()->phiWidth()+correctionScalings->at("phiWidth").Eval(corrections->at("phiWidth")(obj)[0])));
//                obj.setS4(s4_corr);
//                obj.setSieip(sieip_corr);
                obj.full5x5_setShowerShapeVariables(correctedShowerShapes);
                if (obj.full5x5_r9()<0.)
                    {
                        std::cout << "WARNING: R9<0, Original R9: " << obj.userFloat("uncorr_r9") << " corrected R9: " << obj.full5x5_r9() << std::endl;
                    }
            }
        }
return obj;
}
typedef PhotonShowerShapeCorrector CorrectionSSforPhoton;     
//typedef EGShowerShapeCorrector<reco::GsfElectron, pat::Electron> CorrectionEGSSforElectron;
#endif
