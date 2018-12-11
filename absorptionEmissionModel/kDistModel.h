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

#include "interpolationLookUpTable.H"
#include "absorptionEmissionModel.H"
#include "HashTable.H"
#include "absorptionCoeffs.H"
#include "fluidThermo.H"

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

        //- Absorption model dictionary
        dictionary coeffsDict_;

        //- Hash table with species names
        HashTable<label> speciesNames_;

        //- Indices of species in the look-up table
        FixedList<label, nSpecies_> specieIndex_;

        //- Proportion of the heat released rate emitted
// ??   FixedList<scalar, maxBands_> iEhrrCoeffs_;

        //- Look-up table of species related to ft
        mutable autoPtr<interpolationLookUpTable<scalar>> lookUpTablePtr_;

        //- Thermo package
        const fluidThermo& thermo_;

        //- Pointer list of species being solved involved in the absorption
        UPtrList<volScalarField> Yj_;


public:

    //- Runtime type information
    TypeName("kDistModel");


    // Constructors

        //- Construct from components
        kDistModel(const dictionary& dict, const fvMesh& mesh);


    //- Destructor
    virtual ~kDistModel();


    // Member Functions

        //- Absorption coefficient for continuous phase
        tmp<volScalarField> aCont(const label bandi = 0) const;

        //- Emission coefficient for continuous phase
        tmp<volScalarField> eCont(const label bandi = 0) const;

        //- Emission contribution for continuous phase
        tmp<volScalarField> ECont(const label bandi = 0) const;


        inline bool isGrey() const
        {
            return false;
        }

        //- Correct rays
/* ??   void correct
        (
            volScalarField& a,
            PtrList<volScalarField>& aLambda
        ) const;                                ??  */
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace radiation
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
