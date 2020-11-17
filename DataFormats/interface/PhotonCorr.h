#ifndef DataFormats_PhotonCorr_h
#define DataFormats_PhotonCorr_h

#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "FWCore/Utilities/interface/EDMException.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

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

	float const sipip() const {return sipip_;};
        float const sieip() const {return sieip_;};
	float const s4() const {return S4_;};
	
    private:

	float sipip_;
	float sieip_;
	float S4_;


    };

#endif
