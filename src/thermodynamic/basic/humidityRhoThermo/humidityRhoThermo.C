/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "humidityRhoThermo.H"
#include "volFields.H"
#include "fixedHumidityFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(humidityRhoThermo, 0);
    defineRunTimeSelectionTable(humidityRhoThermo, fvMesh);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::humidityRhoThermo::humidityRhoThermo
(
    const fvMesh& mesh,
    const word& phaseName
)
:
    fluidThermo(mesh, phaseName),
    rho_
    (
        IOobject
        (
            phasePropertyName("thermo:rho"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimDensity
    ),

    psi_
    (
        IOobject
        (
            phasePropertyName("thermo:psi"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionSet(0, -2, 2, 0, 0)
    ),

    mu_
    (
        IOobject
        (
            phasePropertyName("thermo:mu"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionSet(1, -1, -1, 0, 0)
    ),

    relHum_
    (
        IOobject
        (
            phasePropertyName("thermo:relHum"),
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimless
    ),

    waterMass_
    (
        IOobject
        (
            phasePropertyName("thermo:waterMass"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimMass
    ),

    waterVapor_
    (
        IOobject
        (
            phasePropertyName("thermo:waterVapor"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimDensity
    ),

    maxWaterVapor_
    (
        IOobject
        (
            phasePropertyName("maxWaterVapor"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimDensity
    ),

    specificHumidity_
    (
        IOobject
        (
            phasePropertyName("thermo:specificHumidity"),
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimless
    ),

    maxSpecificHumidity_
    (
        IOobject
        (
            phasePropertyName("thermo:maxSpecificHumidity"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimless
    ),

    pSatH2O_
    (
        IOobject
        (
            phasePropertyName("thermo:pSatH2O"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimPressure
    ),

    partialPressureH2O_
    (
        IOobject
        (
            phasePropertyName("thermo:partialPressureH2O"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimPressure
    ),

    method_("buck"),

    initWithRelHumidity_(false),

    muEff_
    (
        IOobject
        (
            phasePropertyName("thermo:muEff"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionSet(1,-1,-1,0,0,0,0)
    )
{
    // Read or build the specificHumidity field
    readOrInitSpecificHumidity();

    readMethod();
}


Foam::humidityRhoThermo::humidityRhoThermo
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& phaseName
)
:
    fluidThermo(mesh, dict, phaseName),
    rho_
    (
        IOobject
        (
            phasePropertyName("thermo:rho"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimDensity
    ),

    psi_
    (
        IOobject
        (
            phasePropertyName("thermo:psi"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionSet(0, -2, 2, 0, 0)
    ),

    mu_
    (
        IOobject
        (
            phasePropertyName("thermo:mu"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionSet(1, -1, -1, 0, 0)
    ),

    relHum_
    (
        IOobject
        (
            phasePropertyName("thermo:relHum"),
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimless
    ),

    waterMass_
    (
        IOobject
        (
            phasePropertyName("thermo:waterMass"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimMass
    ),

    waterVapor_
    (
        IOobject
        (
            phasePropertyName("thermo:waterVapor"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimDensity
    ),

    maxWaterVapor_
    (
        IOobject
        (
            phasePropertyName("thermo:maxWaterVapor"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimDensity
    ),

    specificHumidity_
    (
        IOobject
        (
            phasePropertyName("thermo:specificHumidity"),
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimless
    ),

    maxSpecificHumidity_
    (
        IOobject
        (
            phasePropertyName("thermo:maxSpecificHumidity"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimless
    ),

    pSatH2O_
    (
        IOobject
        (
            phasePropertyName("thermo:pSatH2O"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimPressure
    ),

    partialPressureH2O_
    (
        IOobject
        (
            phasePropertyName("thermo:partialPressureH2O"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimPressure
    ),

    method_("buck"),

    initWithRelHumidity_(false),

    muEff_
    (
        IOobject
        (
            phasePropertyName("thermo:muEff"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionSet(1,-1,-1,0,0,0,0)
    )
{
    // Read or build the specificHumidity field
    readOrInitSpecificHumidity();

    // Set the method regarding the calulation of the pSat equation
    readMethod();
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::humidityRhoThermo> Foam::humidityRhoThermo::New
(
    const fvMesh& mesh,
    const word& phaseName
)
{
    return basicThermo::New<humidityRhoThermo>(mesh, phaseName);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::humidityRhoThermo::~humidityRhoThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::humidityRhoThermo::rho() const
{
    return rho_;
}


Foam::tmp<Foam::scalarField>
Foam::humidityRhoThermo::rho(const label patchi) const
{
    return rho_.boundaryField()[patchi];
}


Foam::volScalarField& Foam::humidityRhoThermo::rho()
{
    return rho_;
}


void Foam::humidityRhoThermo::correctRho(const Foam::volScalarField& deltaRho)
{
    rho_ += deltaRho;
}


const Foam::volScalarField& Foam::humidityRhoThermo::psi() const
{
    return psi_;
}


Foam::tmp<Foam::volScalarField> Foam::humidityRhoThermo::mu() const
{
    return mu_;
}


Foam::tmp<Foam::scalarField>
Foam::humidityRhoThermo::mu(const label patchi) const
{
    return mu_.boundaryField()[patchi];
}


void Foam::humidityRhoThermo::readMethod()
{
    //- Change method with respect to fixedHumidity BC
    if (this->rho_.mesh().foundObject<IOList<word>>("methodName"))
    {
        method_ =
            this->rho_.mesh().lookupObject<IOList<word>>("methodName")[0];
    }

    Info<< "Saturation pressure calculation based on "
        << method_ << "\n" << endl;
}


void Foam::humidityRhoThermo::readOrInitSpecificHumidity()
{
    // specificHumidity field is available and was read before
    if (specificHumidity_.typeHeaderOk<volScalarField>())
    {
        Info<< "Initilize humidity by using the thermo:specificHumidity field\n"
            << endl;

        return;
    }

    // Relative humidity field not provided
    if (!relHum_.typeHeaderOk<volScalarField>())
    {
        FatalErrorInFunction
            << "Neither the thermo:specificHumidity or the thermo:relHum "
            << "field was provided in the time-folder"
            << exit(FatalError);
    }
    else
    {
        initWithRelHumidity_ = true;

        Info<< "Initilize humidity by using the thermo:relHum field\n"
            << endl;
    }
}

// ************************************************************************* //
