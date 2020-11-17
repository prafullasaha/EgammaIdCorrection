#include "EgammaIdCorrection/DataFormats/interface/PhotonCorr.h"
#include "EgammaIdCorrection/DataFormats/interface/Photon.h"
#include "DataFormats/Common/interface/Wrapper.h"
#include <vector>
#include <map>

struct dictionary {
	PhotonCorr pCorr;
	std::vector<PhotonCorr> vec_pCorr;
	edm::Wrapper<std::vector<PhotonCorr> > wrp_vec_pCorr;	

};
