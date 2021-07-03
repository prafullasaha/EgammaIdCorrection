#ifndef DataFormats_PhotonCorr_h
#define DataFormats_PhotonCorr_h

#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "FWCore/Utilities/interface/EDMException.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/Common/interface/RefToPtr.h"

#include <map>
#include <string>
#include <set>

    class PhotonCorr : public pat::Photon
    {

    public:
        PhotonCorr();
        PhotonCorr(const pat::Photon & );
	~PhotonCorr();	
	virtual PhotonCorr *clone() const;
	void ZeroVariables();
        void setSipip( float val ) {sipip_ = val;};
        void setSieip( float val ) {sieip_ = val;};
	void setS4( float val ) {S4_ = val;};
        void setpfPhoIso04( float val ) {pfPhoIso04_ = val;};
        void setpfPhoIso03( float val ) {pfPhoIso03_ = val;};
        void setpfPhoIso03Corr( float val ) {pfPhoIso03Cor_ = val;};

	float const sipip() const {return sipip_;};
        float const sieip() const {return sieip_;};
	float const s4() const {return S4_;};
        float const pfPhoIso04() const {return pfPhoIso04_;};
        float const pfPhoIso03() const {return pfPhoIso03_;};
        float const pfPhoIso03Corr() const {return pfPhoIso03Cor_;};

	static bool vetoPackedCand( const pat::Photon &photon, const edm::Ptr<pat::PackedCandidate> &pfcand );	
	float pfCaloIso( const pat::Photon &,
                         const std::vector<edm::Ptr<pat::PackedCandidate> > &,
                         float, float, float, float, float, float, float, reco::PFCandidate::ParticleType, const reco::Vertex* &vtx );	
    private:

	float sipip_;
	float sieip_;
	float S4_;
        float pfPhoIso04_;
        float pfPhoIso03_;
        float pfPhoIso03Cor_;

    };

#endif
