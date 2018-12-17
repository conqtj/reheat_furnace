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
\*---------------------------------------------------------------------------*/

#include "WSGGModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace radiation
    {
        defineTypeNameAndDebug(WSGGModel, 0);

        addToRunTimeSelectionTable
        (
            absorptionEmissionModel,
            WSGGModel,
            dictionary
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::WSGGModel::WSGGModel
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    absorptionEmissionModel(dict, mesh),
    thermo_(mesh.lookupObject<fluidThermo>(basicThermo::dictName)),
    Nq(8),
    Tref(1200)
{
    readData();
}

// * * * * * * * * * * * * * Public Functions * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::radiation::WSGGModel::aCont(const label bandi = 0) const
{
    const volScalarField& T = thermo_.T();

    volScalarField H2O
    (
        IOobject
        (
            "H2O",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );
    volScalarField CO2
    (
        IOobject
        (
            "CO2",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    tmp<volScalarField> ta
    (
        new volScalarField
        (
            IOobject
            (
                "a",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedScalar("a", dimless/dimLength, 0)
        )
    );

    scalarField& a = ta.ref().primitiveFieldRef();

    scalarList kg(Nq);
    scalarList aFg(Nq);

    scalar Ts;

    g2 = fskdistmix(XCO2p,XH2Op,Tp,Tp,*k1,2);

    forAll(a, celli)
    {
        scalar XH2O = std::max(1e-06, H2O[celli]);
        scalar XCO2 = std::max(1e-06, CO2[celli]);

        scalar Tg = T[celli];

        kg[0]=1e-12;
        kg[1]=KG[0]*(XCO2+XH2O);
        kg[2]=KG[1]*(XCO2+XH2O);
        kg[3]=KG[2]*(XCO2+XH2O);
        
        for(int i=0;i<3;i++)
        {
            aFg[i+1]=0.0;
            for(int j=0;j<4;j++)
            aFg[i+1]+=bij[i][j]*pow(Tg,j);
        }
        
        //clear gas
        aFg[0]=1.0-aFg[1]-aFg[2]-aFg[3];

        a[celli]
    }

    return ta;
}


Foam::tmp<Foam::volScalarField>
Foam::radiation::WSGGModel::eCont(const label bandi = 0) const
{
    const volScalarField& T = thermo_.T();

    volScalarField H2O
    (
        IOobject
        (
            "H2O",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );
    volScalarField CO2
    (
        IOobject
        (
            "CO2",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    tmp<volScalarField> te
    (
        new volScalarField
        (
            IOobject
            (
                "e",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedScalar("e", dimless/dimLength, 0)
        )
    );

    scalarField& e = te.ref().primitiveFieldRef();

    scalarList kg(Nq);
    scalarList aFg(Nq);

    scalar Ts;

    forAll(e, celli)
    {
        scalar XH2O = std::max(1e-06, H2O[celli]);
        scalar XCO2 = std::max(1e-06, CO2[celli]);

        scalar Tg = T[celli];

        kg[0]=1e-12;
        kg[1]=KG[0]*(XCO2+XH2O);
        kg[2]=KG[1]*(XCO2+XH2O);
        kg[3]=KG[2]*(XCO2+XH2O);
        
        for(int i=0;i<3;i++)
        {
            aFg[i+1]=0.0;
            for(int j=0;j<4;j++)
            aFg[i+1]+=bij[i][j]*pow(Tg,j);
        }
        
        //clear gas
        aFg[0]=1.0-aFg[1]-aFg[2]-aFg[3];

        e[celli]
    }

    return te;
}

// * * * * * * * * * * * * * * Read Data * * * * * * * * * * * * * * * * * //

void Foam::radiation::WSGGModel::readData()
{
    //read coefficients for the WSGG Model 
    string s;

    ifstream dataIn;
    dataIn.open("wsgg.xml");

    if(dataIn)
    {
        getline(dataIn,s);
        getline(dataIn,s);

        for(int i=0;i<3;i++)
            dataIn>>KG[i]>>bij[i][0]>>bij[i][1]>>bij[i][2]>>bij[i][3];


    for(int j=0;j<4;j++)
    {
        for(int i=0;i<3;i++)
            dataIn>>CAijk[i][j][0]>>CAijk[i][j][1]>>CAijk[i][j][2]>>CAijk[i][j][3];
    }

    }

    dataIn.close();

    return;
}

// * * * * * * * * * * * * * Access Functions * * * * * * * * * * * * * * //
Foam::scalarList Foam::radiation::WSGGModel::getw()
{
    return w;
}