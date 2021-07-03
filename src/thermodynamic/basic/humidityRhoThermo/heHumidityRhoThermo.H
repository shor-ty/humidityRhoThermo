/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2015-2017 OpenCFD Ltd.
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

#include "rhoThermo.H"
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
        void calculate
        (
            const volScalarField& p,
            volScalarField& T,
            volScalarField& he,
            volScalarField& psi,
            volScalarField& rho,
            volScalarField& mu,
            volScalarField& alpha,
            const bool doOldTimes
        );

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


        //- Construct from mesh,dictionary,phase name with a single temperature
        heHumidityRhoThermo
        (
            const fvMesh&,
            const word& phaseName,
            const word& dictName
        );


    //- Destructor
    virtual ~heHumidityRhoThermo();


    // Member functions

        //- Initialize the thermo:specificHumidity using the thermo:relHum field
        virtual void initialize();

        //- Update properties
        virtual void correct();

        //- Calculate the saturation pressure of H2O [Pa]
        virtual void pSatH2O();

        //- Calculate partial pressure of water saturation pressure [Pa]
        //  based from the relative humidity (only for initialization)
        virtual void partialPressureH2OFromRelHumidity();

        //- Calculate partial pressure of water saturation pressure [Pa]
        virtual void partialPressureH2O();

        //- Calculate the relative humidity [-]
        virtual void relHumidity();

        //- Calculate the water vapor content [kg(water)/m^3(air)]
        virtual void waterVapor();

        //- Calculate the maximum specific humidity (for stabilization)
        virtual void maxSpecificHumidity();

        //- Calculate the water mass [kg(water)]
        virtual void waterMass();

        //- Conserved transport equation for the specific humidity [kg/kg]
        virtual void specificHumidityTransport();

        //- Initialize the specific humidity field by using the relative
        //  humidity field
        virtual void initialSpecificHumidityFromRelHumidity();

        //- Update density field
        virtual void updateRho(volScalarField&);

        //- Limit the specific humidity to the max range
        virtual void limitMax();
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
