/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.
    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.
    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.
Class
    Foam::radiation::kDistModel
Description
    Class for generation of k-distribution and other related parameters such as
    a-functions and scaling factors
SourceFiles
    kDistModel.C
\*---------------------------------------------------------------------------*/

#ifndef kDistModel_H
#define kDistModel_H

#include "absorptionEmissionModel.H"
#include "fluidThermo.H"

#include "spectralMaths.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{

/*---------------------------------------------------------------------------*\
                  Class kDistModel Declaration
\*---------------------------------------------------------------------------*/

class kDistModel
:
    public absorptionEmissionModel
{
private:

    // Private data
    label nop;
    scalar dlmn[4][4][4];
    scalar almn[4][4][4];
    scalar blmn[3][3][2];

    scalarList *k1,*k2;
    scalar kmin, kmax;

    //Maths Variable 
        label Nq;
        // g is the quadrature point and w is the weight
        scalarList w, g;

    //reference state of the gas
    scalar Tp,XCO2p,XH2Op;
    
    //- Thermo package
    const fluidThermo& thermo_;

    //functions
        //Read Data
        void readData();
        //Set Reference
        void setRefState();
        //kDist Functions
        scalarList aFunction(const scalarList&, const scalarList&, label);
        scalarList fskdistmix(scalar,scalar,scalar,scalar,const scalarList&,label);
        scalar fskDco2(scalar, scalar, scalar);
        scalar fskDh2o(scalar, scalar, scalar, scalar);

public:

    //- Runtime type information
    TypeName("kDistModel");


    // Constructors
    kDistModel(const dictionary& dict, const fvMesh& mesh);


    //- Destructor
    virtual ~kDistModel();


    //- Absorption coefficient for continuous phase
    tmp<volScalarField> aCont(const label bandi = 0) const;

    //- Emission coefficient for continuous phase
    tmp<volScalarField> eCont(const label bandi = 0) const;

    inline bool isGrey() const
    {
        return false;
    }

    //- Correct rays
    void correct
    (
        volScalarField& a,
        PtrList<volScalarField>& aLambda
    ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace radiation
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //