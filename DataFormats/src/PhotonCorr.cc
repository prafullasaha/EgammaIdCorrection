#include "EgammaIdCorrection/DataFormats/interface/PhotonCorr.h"
#include "FWCore/Utilities/interface/Exception.h"
#include <limits>
#include <algorithm>

PhotonCorr::PhotonCorr() : pat::Photon::Photon()
{
	ZeroVariables();
}

void PhotonCorr::ZeroVariables()
{
    sipip_ = 0.;
    sieip_ = 0.;
    S4_ = 0.;
}

// N.B. Other attributes are set using methods in header file
PhotonCorr::PhotonCorr(const pat::Photon &aPhoton ) : pat::Photon::Photon( aPhoton )
{
    ZeroVariables();
}

PhotonCorr::~PhotonCorr()
{}

PhotonCorr *PhotonCorr::clone() const { return new PhotonCorr( *this ); }
