#ifndef EGShowerShapeCorrector_h
#define EGShowerShapeCorrector_h

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
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"

using namespace std;
using namespace reco;


template<class RECOObj , class PATObj, bool correctShowerShapes, bool correctIsolations>
class EGShowerShapeCorrector
{
public:
	typedef RECOObj reco_obj;
	typedef PATObj pat_obj;
	EGShowerShapeCorrector();
	EGShowerShapeCorrector( const edm::ParameterSet &iConfig, GlobalVariablesComputer& global_var);
        ~EGShowerShapeCorrector();
//        const static bool useXGB = true;
	pat_obj corrector( pat_obj &obj, CLHEP::HepRandomEngine& engine, double rho, noZS::EcalClusterLazyTools& LazyTool, const std::vector<edm::Ptr<pat::PackedCandidate> > &, const reco::Vertex &vtx );
        pat_obj SScorrector( pat_obj &obj, CLHEP::HepRandomEngine& engine, double rho, noZS::EcalClusterLazyTools& LazyTool, const std::vector<edm::Ptr<pat::PackedCandidate> > &, const reco::Vertex &vtx );
private:
	static vector<string> showerShapes_;
        const static bool useXGB = true;
//	edm::ConsumesCollector cc_;
//        GlobalVariablesComputer globalVariablesComputer_;
	std::map<string, MVAComputer<reco_obj, StringObjectFunction<reco_obj, true>, useXGB> > correctionsEB_;
        map<string, TFormula> correctionScalingsEB_;
        std::map<string, MVAComputer<reco_obj, StringObjectFunction<reco_obj, true>, useXGB> > correctionsEE_;
        map<string, TFormula> correctionScalingsEE_;
//	bool correctShowerShapes_;
//	bool correctIsolations_;
};

template<class F, class T, bool correctShowerShapes, bool correctIsolations>
vector<string> EGShowerShapeCorrector< F , T , correctShowerShapes, correctIsolations>::showerShapes_ = {"r9", "s4", "sieie", "sieip", "etaWidth", "phiWidth"};

template<class F, class T, bool correctShowerShapes, bool correctIsolations>
EGShowerShapeCorrector< F , T , correctShowerShapes, correctIsolations>::EGShowerShapeCorrector( const edm::ParameterSet &iConfig, GlobalVariablesComputer& global_var)
{
//	correctShowerShapes_ = correctShowerShapes;
//	correctIsolations_   = correctIsolations;
        if(correctShowerShapes)
        {
            for(auto& ss_var : showerShapes_)
            {
		//---EB
                auto xgb_config = iConfig.getParameter<edm::ParameterSet>(ss_var+"_corrector_config_EB");
                correctionsEB_[ss_var] = MVAComputer<reco_obj, StringObjectFunction<reco_obj, true>, useXGB>(xgb_config, &global_var);
                correctionScalingsEB_[ss_var] = TFormula("", xgb_config.getParameter<string>("regr_output_scaling").c_str());
		//---EE
                xgb_config = iConfig.getParameter<edm::ParameterSet>(ss_var+"_corrector_config_EE");
                correctionsEE_[ss_var] = MVAComputer<reco_obj, StringObjectFunction<reco_obj, true>, useXGB>(xgb_config, &global_var);
                correctionScalingsEE_[ss_var] = TFormula("", xgb_config.getParameter<string>("regr_output_scaling").c_str());

            }
	}
        if(correctIsolations)
        {
	    //---EB
	    // pho iso
            auto iso_config = iConfig.getParameter<edm::ParameterSet>("phoIso_corrector_config_EB"); 
            correctionsEB_["phoIsoClfData"] = MVAComputer<reco_obj, StringObjectFunction<reco_obj, true>, useXGB>(
                iso_config.getParameter<edm::ParameterSet>("clf_data"), 
                &global_var);
            correctionsEB_["phoIsoClfMC"] = MVAComputer<reco_obj, StringObjectFunction<reco_obj, true>, useXGB>(
                iso_config.getParameter<edm::ParameterSet>("clf_mc"), 
                &global_var);
            correctionsEB_["phoIsoPeak2Tail"] = MVAComputer<reco_obj, StringObjectFunction<reco_obj, true>, useXGB>(
                iso_config.getParameter<edm::ParameterSet>("peak2tail"), 
                &global_var);
            correctionsEB_["phoIsoMorphing"] = MVAComputer<reco_obj, StringObjectFunction<reco_obj, true>, useXGB>(
                iso_config.getParameter<edm::ParameterSet>("morphing"), 
                &global_var);
            correctionScalingsEB_["phoIsoMorphing"] = TFormula("", iso_config.getParameter<edm::ParameterSet>("morphing").getParameter<string>("regr_output_scaling").c_str());
            // ch iso
/*            iso_config = iConfig.getParameter<edm::ParameterSet>("chIso_corrector_config_EB"); 
            correctionsEB_["chIsoClfData"] = MVAComputer<reco_obj, StringObjectFunction<reco_obj, true>, useXGB>(
                iso_config.getParameter<edm::ParameterSet>("clf_data"), 
                &global_var);
            correctionsEB_["chIsoClfMC"] = MVAComputer<reco_obj, StringObjectFunction<reco_obj, true>, useXGB>(
                iso_config.getParameter<edm::ParameterSet>("clf_mc"), 
                &global_var);
            correctionsEB_["chIsoPeak2Tail"] = MVAComputer<reco_obj, StringObjectFunction<reco_obj, true>, useXGB>(
                iso_config.getParameter<edm::ParameterSet>("chIso_peak2tail"), 
                &global_var);
            correctionsEB_["chIsoMorphing"] = MVAComputer<reco_obj, StringObjectFunction<reco_obj, true>, useXGB>(
                iso_config.getParameter<edm::ParameterSet>("chIso_morphing"), 
                &global_var);
            correctionScalingsEB_["chIsoMorphing"] = TFormula("", iso_config.getParameter<edm::ParameterSet>("chIso_morphing").getParameter<string>("regr_output_scaling").c_str());            
            correctionsEB_["chIsoWorstPeak2Tail"] = MVAComputer<reco_obj, StringObjectFunction<reco_obj, true>, useXGB>(
                iso_config.getParameter<edm::ParameterSet>("chIsoWorst_peak2tail"), 
                &global_var);
            correctionsEB_["chIsoWorstMorphing"] = MVAComputer<reco_obj, StringObjectFunction<reco_obj, true>, useXGB>(
                iso_config.getParameter<edm::ParameterSet>("chIsoWorst_morphing"), 
                &global_var);
            correctionScalingsEB_["chIsoWorstMorphing"] = TFormula("", iso_config.getParameter<edm::ParameterSet>("chIsoWorst_morphing").getParameter<string>("regr_output_scaling").c_str());            
*/
            //---EE
            // pho iso
            iso_config = iConfig.getParameter<edm::ParameterSet>("phoIso_corrector_config_EE"); 
            correctionsEE_["phoIsoClfData"] = MVAComputer<reco_obj, StringObjectFunction<reco_obj, true>, useXGB>(
                iso_config.getParameter<edm::ParameterSet>("clf_data"), 
                &global_var);
            correctionsEE_["phoIsoClfMC"] = MVAComputer<reco_obj, StringObjectFunction<reco_obj, true>, useXGB>(
                iso_config.getParameter<edm::ParameterSet>("clf_mc"), 
                &global_var);
            correctionsEE_["phoIsoPeak2Tail"] = MVAComputer<reco_obj, StringObjectFunction<reco_obj, true>, useXGB>(
                iso_config.getParameter<edm::ParameterSet>("peak2tail"), 
                &global_var);
            correctionsEE_["phoIsoMorphing"] = MVAComputer<reco_obj, StringObjectFunction<reco_obj, true>, useXGB>(
                iso_config.getParameter<edm::ParameterSet>("morphing"), 
                &global_var);
            correctionScalingsEE_["phoIsoMorphing"] = TFormula("", iso_config.getParameter<edm::ParameterSet>("morphing").getParameter<string>("regr_output_scaling").c_str());
            // ch iso
/*            iso_config = iConfig.getParameter<edm::ParameterSet>("chIso_corrector_config_EE"); 
            correctionsEE_["chIsoClfData"] = MVAComputer<reco_obj, StringObjectFunction<reco_obj, true>, useXGB>(
                iso_config.getParameter<edm::ParameterSet>("clf_data"), 
                &global_var);
            correctionsEE_["chIsoClfMC"] = MVAComputer<reco_obj, StringObjectFunction<reco_obj, true>, useXGB>(
                iso_config.getParameter<edm::ParameterSet>("clf_mc"), 
                &global_var);
            correctionsEE_["chIsoPeak2Tail"] = MVAComputer<reco_obj, StringObjectFunction<reco_obj, true>, useXGB>(
                iso_config.getParameter<edm::ParameterSet>("chIso_peak2tail"), 
                &global_var);
            correctionsEE_["chIsoMorphing"] = MVAComputer<reco_obj, StringObjectFunction<reco_obj, true>, useXGB>(
                iso_config.getParameter<edm::ParameterSet>("chIso_morphing"), 
                &global_var);
            correctionScalingsEE_["chIsoMorphing"] = TFormula("", iso_config.getParameter<edm::ParameterSet>("chIso_morphing").getParameter<string>("regr_output_scaling").c_str());            
            correctionsEE_["chIsoWorstPeak2Tail"] = MVAComputer<reco_obj, StringObjectFunction<reco_obj, true>, useXGB>(
                iso_config.getParameter<edm::ParameterSet>("chIsoWorst_peak2tail"), 
                &global_var);
            correctionsEE_["chIsoWorstMorphing"] = MVAComputer<reco_obj, StringObjectFunction<reco_obj, true>, useXGB>(
                iso_config.getParameter<edm::ParameterSet>("chIsoWorst_morphing"), 
                &global_var);
            correctionScalingsEE_["chIsoWorstMorphing"] = TFormula("", iso_config.getParameter<edm::ParameterSet>("chIsoWorst_morphing").getParameter<string>("regr_output_scaling").c_str());            
*/        }
}

template<class F, class T , bool correctShowerShapes, bool correctIsolations> EGShowerShapeCorrector< F, T , correctShowerShapes, correctIsolations>::~EGShowerShapeCorrector() {}

template<class F, class T, bool correctShowerShapes, bool correctIsolations>
     T EGShowerShapeCorrector< F, T , correctShowerShapes, correctIsolations>::SScorrector( pat_obj &obj , CLHEP::HepRandomEngine& engine, double rho, noZS::EcalClusterLazyTools& lazyTool, const std::vector<edm::Ptr<pat::PackedCandidate> > &pfcands, const reco::Vertex &neut_vtx )
{
        const auto* corrections = std::abs(obj.superCluster()->eta())<1.5 ? &correctionsEB_ : &correctionsEE_;
        const auto* correctionScalings = std::abs(obj.superCluster()->eta())<1.5 ? &correctionScalingsEB_ : &correctionScalingsEE_;

        const reco::CaloClusterPtr  seed_clu = obj.superCluster()->seed();
        const reco::SuperClusterRef super_clu = obj.superCluster();

        std::vector<float> viCov = lazyTool.localCovariances( *seed_clu );
        float sieip = viCov[1];
        float s4 = lazyTool.e2x2( *seed_clu ) / lazyTool.e5x5( *seed_clu );
if(!(std::isnan(obj.full5x5_r9()) || std::isnan(s4) || std::isnan(obj.full5x5_sigmaIetaIeta()) || std::isnan(sieip) || std::isnan(obj.superCluster()->etaWidth()) || std::isnan(obj.superCluster()->phiWidth())
/* || std::isnan(obj.pfPhoIso03()) || std::isnan(obj.pfChgIsoWrtChosenVtx03()) || std::isnan(obj.pfChgIsoWrtWorstVtx03())*/
 || std::isinf(obj.full5x5_r9()) || std::isinf(s4) || std::isinf(obj.full5x5_sigmaIetaIeta()) || std::isinf(sieip) || std::isinf(obj.superCluster()->etaWidth()) || std::isinf(obj.superCluster()->phiWidth())
// || std::isinf(obj.pfPhoIso03()) || std::isinf(obj.pfChgIsoWrtChosenVtx03()) || std::isinf(obj.pfChgIsoWrtWorstVtx03())
 ))
        {
            if(correctShowerShapes)
            {
//                reco::Photon::ShowerShape correctedShowerShapes = obj.full5x5_showerShapeVariables();
                auto correctedShowerShapes = obj.full5x5_showerShape(); //electron based function
                obj.addUserFloat("uncorr_r9", obj.full5x5_r9());
                obj.addUserFloat("uncorr_s4", s4);
                obj.addUserFloat("uncorr_sieie", obj.full5x5_sigmaIetaIeta());
                obj.addUserFloat("uncorr_sieip", sieip);
                obj.addUserFloat("uncorr_etaWidth", obj.superCluster()->etaWidth());
                obj.addUserFloat("uncorr_phiWidth", obj.superCluster()->phiWidth());
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
               obj.setS4(s4_corr);
               obj.setSieip(sieip_corr);
               obj.full5x5_setShowerShape(correctedShowerShapes); //electron based function

               if (obj.full5x5_r9()<0.)
               {
               std::cout << "WARNING: R9<0, Original R9: " << obj.userFloat("uncorr_r9") << " corrected R9: " << obj.full5x5_r9() << std::endl;
               }
           }
        }
return obj;
}
template<class F, class T, bool correctShowerShapes, bool correctIsolations>
     T EGShowerShapeCorrector< F, T , correctShowerShapes, correctIsolations>::corrector( pat_obj &obj , CLHEP::HepRandomEngine& engine, double rho, noZS::EcalClusterLazyTools& lazyTool, const std::vector<edm::Ptr<pat::PackedCandidate> > &pfcands, const reco::Vertex &neut_vtx )
{
        const auto* corrections = std::abs(obj.superCluster()->eta())<1.5 ? &correctionsEB_ : &correctionsEE_;
        const auto* correctionScalings = std::abs(obj.superCluster()->eta())<1.5 ? &correctionScalingsEB_ : &correctionScalingsEE_;

        const reco::CaloClusterPtr  seed_clu = obj.superCluster()->seed();
        const reco::SuperClusterRef super_clu = obj.superCluster();

        std::vector<float> viCov = lazyTool.localCovariances( *seed_clu );
        float sieip = viCov[1];
        float s4 = lazyTool.e2x2( *seed_clu ) / lazyTool.e5x5( *seed_clu );
if(!(std::isnan(obj.full5x5_r9()) || std::isnan(s4) || std::isnan(obj.full5x5_sigmaIetaIeta()) || std::isnan(sieip) || std::isnan(obj.superCluster()->etaWidth()) || std::isnan(obj.superCluster()->phiWidth())
/* || std::isnan(obj.pfPhoIso03()) || std::isnan(obj.pfChgIsoWrtChosenVtx03()) || std::isnan(obj.pfChgIsoWrtWorstVtx03())*/
 || std::isinf(obj.full5x5_r9()) || std::isinf(s4) || std::isinf(obj.full5x5_sigmaIetaIeta()) || std::isinf(sieip) || std::isinf(obj.superCluster()->etaWidth()) || std::isinf(obj.superCluster()->phiWidth())
// || std::isinf(obj.pfPhoIso03()) || std::isinf(obj.pfChgIsoWrtChosenVtx03()) || std::isinf(obj.pfChgIsoWrtWorstVtx03())
))
        {
            if(correctShowerShapes)
            {
//                reco::Photon::ShowerShape correctedShowerShapes = obj.full5x5_showerShapeVariables();
		auto correctedShowerShapes = obj.full5x5_showerShapeVariables();
                obj.addUserFloat("uncorr_r9", obj.full5x5_r9());
                obj.addUserFloat("uncorr_s4", s4);
                obj.addUserFloat("uncorr_sieie", obj.full5x5_sigmaIetaIeta());
                obj.addUserFloat("uncorr_sieip", sieip);
                obj.addUserFloat("uncorr_etaWidth", obj.superCluster()->etaWidth());
                obj.addUserFloat("uncorr_phiWidth", obj.superCluster()->phiWidth());
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
                obj.setS4(s4_corr);
                obj.setSieip(sieip_corr);
                obj.full5x5_setShowerShapeVariables(correctedShowerShapes);

                if (obj.full5x5_r9()<0.)
                    {
                        std::cout << "WARNING: R9<0, Original R9: " << obj.userFloat("uncorr_r9") << " corrected R9: " << obj.full5x5_r9() << std::endl;
                    }
            }
            if(correctIsolations)
            {
		const reco::Vertex *neut_vtx_ = &neut_vtx;
		//pat_obj *obj_ = &obj;
		float pfPhoIso04 = obj.pfCaloIso( obj, pfcands, 0.4, 0.0, 0.070, 0.015, 0.0, 0.0, 0.0, PFCandidate::gamma, neut_vtx_ );
                float pfPhoIso03 = obj.pfCaloIso( obj, pfcands, 0.3, 0.0, 0.070, 0.015, 0.0, 0.0, 0.0, PFCandidate::gamma, neut_vtx_ ); 
		obj.setpfPhoIso04( pfPhoIso04 );
            	obj.setpfPhoIso03( pfPhoIso03 );		
                //--OC-Photon isolation
                obj.addUserFloat("uncorr_pfPhoIso03", obj.pfPhoIso03());

                // peak to tail shift (and viceversa), with conversion from TMVA value [-1,1] to probability [0,1]
                auto p_tail_data = 1./(1.+sqrt(2./(1.+corrections->at("phoIsoClfData")(obj)[0])-1.));
                auto p_tail_mc =  1./(1.+sqrt(2./(1.+corrections->at("phoIsoClfMC")(obj)[0])-1.));
                auto p_peak_data = 1 - p_tail_data;
                auto p_peak_mc = 1 - p_tail_mc;
                auto migration_rnd_value = engine.flat();
            
                double p_move_to_tail = (p_tail_data-p_tail_mc)/p_peak_mc;
                double p_move_to_peak = (p_peak_data-p_peak_mc)/p_tail_mc;
            
                if(obj.pfPhoIso03() == 0 && p_tail_data > p_tail_mc && migration_rnd_value < p_move_to_tail)
                {
                    obj.addUserFloat("peak2tail_rnd", engine.flat()*(0.99-0.01)+0.01);
                    obj.setpfPhoIso03(corrections->at("phoIsoPeak2Tail")(obj)[0]);
                }
                else if(obj.pfPhoIso03() > 0 && p_peak_data > p_peak_mc && migration_rnd_value <= p_move_to_peak)
                    obj.setpfPhoIso03(0.);

                // tail morphing
                if(obj.pfPhoIso03() > 0.)
                {
                    obj.setpfPhoIso03(obj.pfPhoIso03()+correctionScalings->at("phoIsoMorphing").Eval(corrections->at("phoIsoMorphing")(obj)[0]));
                }
                //Make sure pfPhoIso03Corr is same as pfPhoIso03
                obj.setpfPhoIso03Corr(obj.pfPhoIso03());

            }
        }
return obj;
}
//typedef EGShowerShapeCorrector<reco::Photon, pat::Photon> CorrectionEGSSforPhoton;     
//typedef EGShowerShapeCorrector<reco::GsfElectron, pat::Electron> CorrectionEGSSforElectron;
#endif
