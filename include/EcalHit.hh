#ifndef ECALHIT_HH 
#define ECALHIT_HH 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"


class EcalHit : public G4VHit {
public:

    EcalHit() : fCellID(0), fEdep(0.), fNph_scint(0), fNph_cherenkov(0) { step_x.clear(); step_y.clear(); step_z.clear(); step_E.clear(); }
    EcalHit(G4int cellID) : fCellID(cellID), fEdep(0), fNph_scint(0), fNph_cherenkov(0) { step_x.clear(); step_y.clear(); step_z.clear(); step_E.clear(); }
    virtual ~EcalHit() = default;

    void setCellID(G4int _id) { fCellID = _id; }
    G4int GetCellID() const { return fCellID; }

    void setEdep(G4double _en) { fEdep = _en;}
    void addEdep(G4double _en) { fEdep += _en; }
    G4double GetEdep() const { return fEdep; }

    void setNphScint(G4int _nph) { fNph_scint = _nph;}
    void addNphScint(G4int _nph) { fNph_scint += _nph; }
    G4double GetNphScint() const { return fNph_scint; }    

    void setNphChren(G4int _nph) { fNph_cherenkov = _nph;}
    void addNphChren(G4int _nph) { fNph_cherenkov += _nph; }
    G4double GetNphChren() const { return fNph_cherenkov; }

    void addStep(G4double _posx, G4double _posy, G4double _posz, G4double _en){
      step_x.push_back(_posx);
      step_y.push_back(_posy);
      step_z.push_back(_posz);
      step_E.push_back(_en);
    }
    int GetStepSize() const { return step_E.size(); }
    G4ThreeVector GetStepPos( size_t i ) const {
      if(i>=step_x.size()) return G4ThreeVector(0., 0., 0.);
      G4ThreeVector pos(step_x[i], step_y[i], step_z[i]);
      return pos;
    }
    G4double GetStepE(size_t i) const { return i<step_E.size() ? step_E[i] : -999; }

private:
    G4int fCellID;     
    G4double fEdep;
    G4int fNph_scint;
    G4int fNph_cherenkov;
    std::vector<G4double> step_x;
    std::vector<G4double> step_y;
    std::vector<G4double> step_z;
    std::vector<G4double> step_E;

};

using EcalHitsCollection = G4THitsCollection<EcalHit>;


#endif
