/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::heHumidityRhoThermo

Description
    Energy for a mixture based on density

SourceFiles
    heHumidityRhoThermo.C

\*---------------------------------------------------------------------------*/

#ifndef heHumidityRhoThermo_H
#define heHumidityRhoThermo_H

#include "humidityRhoThermo.H"
#include "heThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class heHumidityRhoThermo Declaration
\*---------------------------------------------------------------------------*/

template<class BasicPsiThermo, class MixtureType>
class heHumidityRhoThermo
:
    public heThermo<BasicPsiThermo, MixtureType>
{
    // Private Member Functions

        //- Calculate the thermo variables
        void calculate();

        //- Construct as copy (not implemented)
        heHumidityRhoThermo(const heHumidityRhoThermo<BasicPsiThermo, MixtureType>&);

public:

    //- Runtime type information
    TypeName("heHumidityRhoThermo");


    // Constructors

        //- Construct from mesh and phase name
        heHumidityRhoThermo
        (
            const fvMesh&,
            const word& phaseName
        );


    //- Destructor
    virtual ~heHumidityRhoThermo();


    // Member functions

        //- Update properties
        virtual void correct();

        //- Calculate the saturation pressure of H2O [Pa]
        virtual void pSatH2O();

        //- Calculate partial pressure of water saturation pressure [Pa]
        virtual void partialPressureH2O();

        //- Calculate the absolute water content [kg/m^3]
        virtual void absHumidity();

        //- Calculate the mass of water [g]
        virtual void massH2O();

        //- Conserved transport equation for H2O
        virtual void H2OTransport();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "heHumidityRhoThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
