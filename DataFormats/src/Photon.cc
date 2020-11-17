#include "EgammaIdCorrection/DataFormats/interface/Photon.h"
#include "FWCore/Utilities/interface/Exception.h"
#include <limits>
#include <algorithm>

using namespace corr;

Photon::Photon() : pat::Photon::Photon()
{
	ZeroVariables();
}

void Photon::ZeroVariables()
{
    sipip_ = 0.;
    sieip_ = 0.;
    S4_ = 0.;
}

// N.B. Other attributes are set using methods in header file
Photon::Photon(const pat::Photon &aPhoton ) : pat::Photon::Photon( aPhoton )
{
    ZeroVariables();
}

Photon::~Photon()
{}

Photon *Photon::clone() const { return new Photon( *this ); }
