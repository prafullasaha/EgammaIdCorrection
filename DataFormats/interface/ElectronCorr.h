#ifndef DataFormats_ElectronCorr_h
#define DataFormats_ElectronCorr_h

#include "DataFormats/EgammaCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "FWCore/Utilities/interface/EDMException.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/Common/interface/RefToPtr.h"

#include <map>
#include <string>
#include <set>

    class ElectronCorr : public pat::Electron
    {

    public:
        ElectronCorr();
        ElectronCorr(const pat::Electron & );
	~ElectronCorr();	
	virtual ElectronCorr *clone() const;
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
