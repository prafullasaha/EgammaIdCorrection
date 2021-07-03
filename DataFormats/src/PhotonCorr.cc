#include "EgammaIdCorrection/DataFormats/interface/PhotonCorr.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include <limits>
#include <algorithm>
using namespace reco;

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

bool PhotonCorr::vetoPackedCand( const pat::Photon &photon, const edm::Ptr<pat::PackedCandidate> &pfcand )
{
    edm::RefVector<pat::PackedCandidateCollection> associated =  photon.associatedPackedPFCandidates();

    int nass = 0;
    for( unsigned int ipc = 0; ipc < associated.size(); ipc++ ) {
        edm::Ptr<pat::PackedCandidate> associatedPtr = edm::refToPtr( associated[ipc] );
        if( associatedPtr == pfcand )  { nass++; break; }
    }

    return ( nass > 0 );
}


float PhotonCorr::pfCaloIso( 	const pat::Photon &photon,
                                const std::vector<edm::Ptr<pat::PackedCandidate> > &pfcandidates,
                                float dRMax,
                                float dRVetoBarrel,
                                float dRVetoEndcap,
                                float etaStripBarrel,
                                float etaStripEndcap,
                                float minEnergyBarrel,
                                float minEnergyEndcap,
                                reco::PFCandidate::ParticleType type,
                                const reco::Vertex* &vtx
                              )
{
    // just used to translate particle types to pdgId. Sure: it's very hugly...
    static reco::PFCandidate helper;
    int pdgId = helper.translateTypeToPdgId( type );

    float isovalue = 0;

    float dRVeto = 99;
    float maxetaStrip = 99;
    bool  removeOverlappingCandidates_ = true;

    if( photon.isEB() ) {
        dRVeto        = dRVetoBarrel;
        maxetaStrip  = etaStripBarrel;
    } else if( photon.isEE() ) {
        dRVeto        = dRVetoEndcap;
        maxetaStrip  = etaStripEndcap;
    }

    //// map<float,tuple<edm::Ptr<pat::PackedCandidate>,float,float> > candidates;
    for( size_t ipf = 0; ipf < pfcandidates.size(); ipf++ ) {

        edm::Ptr<pat::PackedCandidate> pfcand = pfcandidates[ipf];

        if( pfcand->pdgId() != pdgId ) { continue; }
        if( photon.isEB() ) if( fabs( pfcand->pt() ) < minEnergyBarrel )     { continue; }
        if( photon.isEE() ) if( fabs( pfcand->energy() ) < minEnergyEndcap ) { continue; }

        if( removeOverlappingCandidates_ && vetoPackedCand( photon, pfcand ) ) { continue; }

        double vx, vy, vz;
        if( vtx ) {
            vx = vtx->x(), vy = vtx->y(), vz = vtx->z();
        } else {
            math::XYZPoint  pfcandvtx = pfcand->vertex();
            vx = pfcandvtx.x(), vy = pfcandvtx.y(), vz = pfcandvtx.z();
        }
        math::XYZVector SCdirectionWrtCandVtx( photon.superCluster()->x() - vx,
                                               photon.superCluster()->y() - vy,
                                               photon.superCluster()->z() - vz
                                             );

        float dEta = fabs( SCdirectionWrtCandVtx.Eta() - pfcand->momentum().Eta() );
        float dR   = deltaR( SCdirectionWrtCandVtx.Eta(), SCdirectionWrtCandVtx.Phi(), pfcand->momentum().Eta(), pfcand->momentum().Phi() );
        //// float dPhi = deltaPhi( SCdirectionWrtCandVtx.Phi(),pfcand->momentum().Phi() );

        if( dEta < maxetaStrip )        { continue; }
        if( dR < dRVeto || dR > dRMax ) { continue; }

        isovalue += pfcand->pt();
        //// candidates[dR] = make_tuple(pfcand,dEta,dPhi);
    }

    return isovalue;
}

