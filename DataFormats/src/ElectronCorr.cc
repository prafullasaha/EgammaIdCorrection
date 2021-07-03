#include "EgammaIdCorrection/DataFormats/interface/ElectronCorr.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include <limits>
#include <algorithm>

ElectronCorr::ElectronCorr() : pat::Electron::Electron()
{
	ZeroVariables();
}

void ElectronCorr::ZeroVariables()
{
    sipip_ = 0.;
    sieip_ = 0.;
    S4_ = 0.;
}

// N.B. Other attributes are set using methods in header file
ElectronCorr::ElectronCorr(const pat::Electron &aElectron ) : pat::Electron::Electron( aElectron )
{
    ZeroVariables();
}

ElectronCorr::~ElectronCorr()
{}

ElectronCorr *ElectronCorr::clone() const { return new ElectronCorr( *this ); }

