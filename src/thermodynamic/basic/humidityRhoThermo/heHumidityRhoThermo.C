/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

#include "heHumidityRhoThermo.H"
#include "fvMatricesFwd.H"
#include "fvModels.H"
#include "fvConstraints.H"
#include "fvCFD.H"
#include "bound.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class BasicPsiThermo, class MixtureType>
void Foam::heHumidityRhoThermo<BasicPsiThermo, MixtureType>::calculate()
{
    const scalarField& hCells = this->he();
    const scalarField& pCells = this->p_;

    scalarField& TCells = this->T_.primitiveFieldRef();
    scalarField& CpCells = this->Cp_.primitiveFieldRef();
    scalarField& CvCells = this->Cv_.primitiveFieldRef();
    scalarField& psiCells = this->psi_.primitiveFieldRef();
    scalarField& muCells = this->mu_.primitiveFieldRef();
    scalarField& alphaCells = this->alpha_.primitiveFieldRef();

    forAll(TCells, celli)
    {
        const typename MixtureType::thermoMixtureType& thermoMixture =
            this->cellThermoMixture(celli);

        const typename MixtureType::transportMixtureType& transportMixture =
            this->cellTransportMixture(celli, thermoMixture);

        TCells[celli] = thermoMixture.THE
        (
            hCells[celli],
            pCells[celli],
            TCells[celli]
        );

        CpCells[celli] = thermoMixture.Cp(pCells[celli], TCells[celli]);
        CvCells[celli] = thermoMixture.Cv(pCells[celli], TCells[celli]);
        psiCells[celli] = thermoMixture.psi(pCells[celli], TCells[celli]);
        muCells[celli] = transportMixture.mu(pCells[celli], TCells[celli]);
        alphaCells[celli] =
            transportMixture.kappa(pCells[celli], TCells[celli])
           /thermoMixture.Cp(pCells[celli], TCells[celli]);
    }

    volScalarField::Boundary& pBf =
        this->p_.boundaryFieldRef();

    volScalarField::Boundary& TBf =
        this->T_.boundaryFieldRef();
    volScalarField::Boundary& CpBf =
        this->Cp_.boundaryFieldRef();

    volScalarField::Boundary& CvBf =
        this->Cv_.boundaryFieldRef();

    volScalarField::Boundary& psiBf =
        this->psi_.boundaryFieldRef();

    volScalarField::Boundary& heBf =
        this->he().boundaryFieldRef();

    volScalarField::Boundary& muBf =
        this->mu_.boundaryFieldRef();

    volScalarField::Boundary& alphaBf =
        this->alpha_.boundaryFieldRef();

    forAll(this->T_.boundaryField(), patchi)
    {
        fvPatchScalarField& pp = pBf[patchi];
        fvPatchScalarField& pT = TBf[patchi];
        fvPatchScalarField& pCp = CpBf[patchi];
        fvPatchScalarField& pCv = CvBf[patchi];
        fvPatchScalarField& ppsi = psiBf[patchi];
        fvPatchScalarField& phe = heBf[patchi];
        fvPatchScalarField& pmu = muBf[patchi];
        fvPatchScalarField& palpha = alphaBf[patchi];

        if (pT.fixesValue())
        {
            forAll(pT, facei)
            {
                const typename MixtureType::thermoType& thermoMixture =
                    this->patchFaceThermoMixture(patchi, facei);

                const typename MixtureType::transportMixtureType&
                    transportMixture =
                    this->patchFaceTransportMixture
                    (patchi, facei, thermoMixture);

                phe[facei] = thermoMixture.HE(pp[facei], pT[facei]);
                pCp[facei] = thermoMixture.Cp(pp[facei], pT[facei]);
                pCv[facei] = thermoMixture.Cv(pp[facei], pT[facei]);
                ppsi[facei] = thermoMixture.psi(pp[facei], pT[facei]);
                ppsi[facei] = thermoMixture.psi(pp[facei], pT[facei]);
                pmu[facei] = transportMixture.mu(pp[facei], pT[facei]);
                palpha[facei] =
                    transportMixture.kappa(pp[facei], pT[facei])
                   /thermoMixture.Cp(pp[facei], pT[facei]);
            }
        }
        else
        {
            forAll(pT, facei)
            {
                const typename MixtureType::thermoType& thermoMixture =
                    this->patchFaceThermoMixture(patchi, facei);

                const typename MixtureType::transportMixtureType&
                    transportMixture =
                    this->patchFaceTransportMixture
                    (patchi, facei, thermoMixture);

                pT[facei] = thermoMixture.THE(phe[facei], pp[facei], pT[facei]);
                pCp[facei] = thermoMixture.Cp(pp[facei], pT[facei]);
                pCv[facei] = thermoMixture.Cv(pp[facei], pT[facei]);
                ppsi[facei] = thermoMixture.psi(pp[facei], pT[facei]);
                pmu[facei] = transportMixture.mu(pp[facei], pT[facei]);
                palpha[facei] =
                    transportMixture.kappa(pp[facei], pT[facei])
                   /thermoMixture.Cp(pp[facei], pT[facei]);
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicPsiThermo, class MixtureType>
Foam::heHumidityRhoThermo<BasicPsiThermo, MixtureType>::heHumidityRhoThermo
(
    const fvMesh& mesh,
    const word& phaseName
)
:
    heThermo<BasicPsiThermo, MixtureType>(mesh, phaseName)
{

    if (this->initWithRelHumidity_)
    {
        initialize();
    }


    calculate();

    // First initialisation of the density
    updateRho(this->rho_);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasicPsiThermo, class MixtureType>
Foam::heHumidityRhoThermo<BasicPsiThermo, MixtureType>::~heHumidityRhoThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicPsiThermo, class MixtureType>
void Foam::heHumidityRhoThermo<BasicPsiThermo, MixtureType>::initialize()
{
    pSatH2O();

    partialPressureH2OFromRelHumidity();

    initialSpecificHumidityFromRelHumidity();

    //std::terminate();
}

template<class BasicPsiThermo, class MixtureType>
void Foam::heHumidityRhoThermo<BasicPsiThermo, MixtureType>::correct()
{
    if (debug)
    {
        InfoInFunction << endl;
    }

    calculate();

    //Info<< "   Solve transport equation for specific humidity\n";
    specificHumidityTransport();

    //Info<< "   Calculate the saturation pressure of water\n";
    pSatH2O();

    //Info<< "   Calculate the partial pressure of water\n";
    partialPressureH2O();

    //Info<< "   Calculate the relative humidity\n";
    relHumidity();

    //Info<< "   Calculate the water vapor content\n";
    waterVapor();

    //Info<< "   Calculate the maximum specific humidity\n";
    maxSpecificHumidity();

    //Info<< "   Calculate the water mass\n";
    waterMass();

    //Info<< "   Calculate the density field\n";
    updateRho(this->rho_);

    //- Keep the physical bound of the maximum values
    limitMax();

    if (debug)
    {
        Info<< "    Finished" << endl;
    }
}

template<class BasicPsiThermo, class MixtureType>
void Foam::heHumidityRhoThermo<BasicPsiThermo, MixtureType>::pSatH2O()
{
    const volScalarField theta =
        (this->T_
      - dimensionedScalar("Kelvin", dimensionSet(0,0,0,1,0,0,0), scalar(273.15)))
      * dimensionedScalar("perK", dimensionSet(0,0,0,-1,0,0,0), scalar(1));

    //- Magnus formulation
    //  Valid between -50 to 100 degC and 1013.25 hPa
    if (this->method_ == "magnus")
    {
        const dimensionedScalar pre1
        (
            "pre1",
            dimensionSet(1,-1,-2,0,0,0,0),
            scalar(611.2)
        );

        const dimensionedScalar pre2
        (
            "pre2",
            dimensionSet(0,0,0,0,0,0,0),
            scalar(17.62)
        );

        const dimensionedScalar value1
        (
            "value1",
            dimensionSet(0,0,0,0,0,0,0),
            scalar(243.12)
        );

        this->pSatH2O_ = pre1*exp((pre2*(theta))/(value1+theta));
    }
    //  Buck formula [1996]
    //  Valid between 0 to 100 degC and 1013.25 hPa
    //  Very accurate between 0 degC and 50 degC
    else if (this->method_ == "buck")
    {
        const dimensionedScalar pre1
        (
            "pre1",
            dimensionSet(1,-1,-2,0,0,0,0),
            scalar(611.21)
        );

        const dimensionedScalar value1
        (
            "value1",
            dimensionSet(0,0,0,0,0,0,0),
            scalar(18.678)
        );

        const dimensionedScalar value2
        (
            "value2",
            dimensionSet(0,0,0,0,0,0,0),
            scalar(234.5)
        );

        const dimensionedScalar value3
        (
            "value3",
            dimensionSet(0,0,0,0,0,0,0),
            scalar(257.14)
        );

        this->pSatH2O_ =
            pre1*exp(((value1-(theta/value2))*theta/(value3+theta)));
    }
    else
    {
        FatalErrorInFunction
            << "The specified method to calculate the saturation pressure is "
            << "not supported: " << this->method_ << ". Supported methods are "
            << "'buck' and 'magnus'."
            << exit(FatalError);
    }
}


template<class BasicPsiThermo, class MixtureType>
void Foam::heHumidityRhoThermo<BasicPsiThermo, MixtureType>::
partialPressureH2O()
{
    const dimensionedScalar RSpecificH2O
    (
        "gasConstantH2O",
        dimensionSet(0,2,-2,-1,0,0,0),
        scalar(461.51)
    );

    const dimensionedScalar RSpecificDryAir
    (
        "gasConstantDryAir",
        dimensionSet(0,2,-2,-1,0,0,0),
        scalar(287.058)
    );

    const volScalarField& p = this->p_;
    const volScalarField& sH = this->specificHumidity_;

    volScalarField& pPH2O = this->partialPressureH2O_;

    // Equation 20
    pPH2O =
        p*pow(1 - RSpecificDryAir/RSpecificH2O * (1-pow(sH+VSMALL, -1)), -1);
}


template<class BasicPsiThermo, class MixtureType>
void Foam::heHumidityRhoThermo<BasicPsiThermo, MixtureType>::
partialPressureH2OFromRelHumidity()
{
    const volScalarField& relHum = this->relHum_;
    const volScalarField& pSatH2O = this->pSatH2O_;

    volScalarField& partialPressureH2O = this->partialPressureH2O_;

    partialPressureH2O = relHum * pSatH2O;
}


template<class BasicPsiThermo, class MixtureType>
void Foam::heHumidityRhoThermo<BasicPsiThermo, MixtureType>::relHumidity()
{
    volScalarField& relHum = this->relHum_;

    relHum = this->partialPressureH2O_ / this->pSatH2O_;

    forAll(relHum, cellI)
    {
        if (relHum[cellI] > 1.)
        {
            WarningInFunction
                << "Humidity exeeds 100 percent, condensation occur in cell "
                << cellI << " (not implemented)" << endl;
            break;
        }
    }
}


template<class BasicPsiThermo, class MixtureType>
void Foam::heHumidityRhoThermo<BasicPsiThermo, MixtureType>::waterVapor()
{
    const volScalarField& pPH2O = this->partialPressureH2O_;
    const volScalarField& pSatH2O = this->pSatH2O_;
    const volScalarField& T = this->T_;

    const dimensionedScalar RSpecificH2O
    (
        "gasConstantH2O",
        dimensionSet(0,2,-2,-1,0,0,0),
        scalar(461.51)
    );

    // Water vapor in air
    volScalarField& waterVapor = this->waterVapor_;
    waterVapor = pPH2O / (RSpecificH2O * T);

    // Max water vapor possible
    volScalarField& maxWaterVapor = this->maxWaterVapor_;
    maxWaterVapor = pSatH2O / (RSpecificH2O * T);
}


template<class BasicPsiThermo, class MixtureType>
void Foam::heHumidityRhoThermo<BasicPsiThermo, MixtureType>::
maxSpecificHumidity()
{
    const volScalarField& pSatH2O = this->pSatH2O_;
    const volScalarField& p = this->p_;
    const volScalarField& maxWaterVapor = this->maxWaterVapor_;
    const volScalarField& T = this->T_;

    volScalarField& maxSpecificHumidity = this->maxSpecificHumidity_;

    const dimensionedScalar RSpecificDryAir
    (
        "gasConstantDryAir",
        dimensionSet(0,2,-2,-1,0,0,0),
        scalar(287.058)
    );

    // Calculate the maximum possible specific humidity value equation (23)
    maxSpecificHumidity =
        maxWaterVapor / ((p-pSatH2O)/(RSpecificDryAir*T) + maxWaterVapor);
}


template<class BasicPsiThermo, class MixtureType>
void Foam::heHumidityRhoThermo<BasicPsiThermo, MixtureType>::waterMass()
{
    const volScalarField& waterVapor = this->waterVapor_;
    const scalarField& V = waterVapor.mesh().V();

    scalarField& waterMass = this->waterMass_;

    waterMass = waterVapor*V;

    Info<< "   Total water = " << gSum(waterMass) << " kg" << endl;
}


template<class BasicPsiThermo, class MixtureType>
void Foam::heHumidityRhoThermo<BasicPsiThermo, MixtureType>::
specificHumidityTransport()
{
    volScalarField& specHum = this->specificHumidity_;
    volScalarField& muEff = this->muEff_;

    const volScalarField& mu = this->mu_;

    const volScalarField& rho = this->rho_;

    const surfaceScalarField& phi = this->db().objectRegistry
        ::lookupObject<surfaceScalarField>("phi");

    const IOdictionary& turbProp =
        this->db().objectRegistry
        ::lookupObject<IOdictionary>("momentumTransport");

    const word turbulenceMode = turbProp.lookup("simulationType");

    if (turbulenceMode == "RAS")
    {
        const volScalarField& nut =
            this->db().objectRegistry::lookupObject<volScalarField>("nut");

        muEff = rho*nut + mu;
    }
    else
    {
        muEff = mu;
    }

    // fvOptions has been replaced by fvConstraints and fvMoldels
    const Foam::fvModels& fvModels(Foam::fvModels::New(phi.mesh()));
    const Foam::fvConstraints& fvConstraints
    (
        Foam::fvConstraints::New(phi.mesh())
    );

    fvScalarMatrix specHumEqn
    (
        fvm::ddt(rho, specHum)
      + fvm::div(phi, specHum)
     ==
        fvm::laplacian(muEff, specHum)
      + fvModels.source(rho, specHum)
    );

    specHumEqn.relax();
    fvConstraints.constrain(specHumEqn);
    specHumEqn.solve();
    fvConstraints.constrain(specHum);

    //- To keep physical range
    //  Defined between 0 and max water vaper content based on saturation pressure
    //  Maximum is not limited yet but afterward in the limit() function
    bound
    (
        specHum,
        dimensionedScalar("tmp", dimensionSet(0,0,0,0,0,0,0), scalar(0))
    );
}


template<class BasicPsiThermo, class MixtureType>
void Foam::heHumidityRhoThermo<BasicPsiThermo, MixtureType>::
initialSpecificHumidityFromRelHumidity()
{
    const volScalarField& pPH2O = this->partialPressureH2O_;
    const volScalarField& p = this->p_;
    volScalarField& specHum = this->specificHumidity_;

    // Specific gas constant dry air [J/kg/K]
    const dimensionedScalar RdryAir
    (
        "RdryAir",
        dimensionSet(0,2,-2,-1,0,0,0),
        scalar(287.058)
    );

    // Specific gas constant water vapor [J/kg/K]
    const dimensionedScalar RwaterVapor
    (
        "RwaterVapor",
        dimensionSet(0,2,-2,-1,0,0,0),
        scalar(461.51)
    );

    // Initialize the specific humidity field
    scalarField& specHumCells = specHum.primitiveFieldRef();

    forAll(specHumCells, celli)
    {
        specHumCells[celli] =
            pow
            (
                1
              - (1- p[celli]/(pPH2O[celli]+VSMALL))
              * RwaterVapor.value()/RdryAir.value(),
               -1
            );
    }


    // Assign the same boundary conditions
    const volScalarField::Boundary& relHumBf = this->relHum_.boundaryField();
    volScalarField::Boundary& specHumBf = specHum.boundaryFieldRef();

    // Not the best way
    forAll(specHumBf, patchi)
    {
        if
        (
            relHumBf[patchi].fixesValue()
         && (relHumBf[patchi].type() != "fixedHumidity")
        )
        {
            FatalErrorInFunction
                << "The boundary type '" << relHumBf[patchi].type() << "' "
                << "for the thermo:relHum field is not allowed"
                << "Use the 'fixedHumidity' boundary condition instead"
                << exit(FatalError);
        }
    }

    // Update the boundary data
    forAll(specHumBf, patchi)
    {
        fvPatchScalarField& pspecHum = specHumBf[patchi];

        forAll(pspecHum, facei)
        {
            pspecHum[facei] =
                pow
                (
                    1
                  - (1- p[facei]/(pPH2O[facei]+VSMALL))
                  * RwaterVapor.value()/RdryAir.value(),
                   -1
                );
        }
    }

    specHum.write();
}

template<class BasicPsiThermo, class MixtureType>
void Foam::heHumidityRhoThermo<BasicPsiThermo, MixtureType>::updateRho
(
    volScalarField& rho
)
{
    const volScalarField& T = this->T_;
    const volScalarField& pPH2O = this->partialPressureH2O_;
    const volScalarField& p = this->p_;

    // Specific gas constant dry air [J/kg/K]
    const dimensionedScalar RdryAir
    (
        "RdryAir",
        dimensionSet(0,2,-2,-1,0,0,0),
        scalar(287.058)
    );

    // Specific gas constant water vapor [J/kg/K]
    const dimensionedScalar RwaterVapor
    (
        "RwaterVapor",
        dimensionSet(0,2,-2,-1,0,0,0),
        scalar(461.51)
    );

    scalarField& rhoCells = rho.primitiveFieldRef();

    forAll(rhoCells, celli)
    {
        rhoCells[celli] =
            1/T[celli]
          * (
                (p[celli] - pPH2O[celli])
              / RdryAir.value() + pPH2O[celli]/RwaterVapor.value()
            );
    }

    volScalarField::Boundary& rhoBf = rho.boundaryFieldRef();

    forAll(rhoBf, patchi)
    {
        fvPatchScalarField& prho = rhoBf[patchi];

        forAll(prho, facei)
        {
            prho[facei] =
                1/T[facei]
              * (
                    (p[facei] - pPH2O[facei])
                  / RdryAir.value() + pPH2O[facei]/RwaterVapor.value()
                );
        }
    }
}


template<class BasicPsiThermo, class MixtureType>
void Foam::heHumidityRhoThermo<BasicPsiThermo, MixtureType>::limitMax()
{
    const volScalarField& maxSpecHum = this->maxSpecificHumidity_;

    volScalarField& specHum = this->specificHumidity_;

    label max = 0;

    forAll(specHum, cellI)
    {
        if (specHum[cellI] > maxSpecHum[cellI])
        {
            specHum[cellI] = maxSpecHum[cellI];
            max++;
        }
    }

    if (max > 0)
    {
        Info<< "    Correcting " << max << " cells which were higher than max\n";
    }
}


// ************************************************************************* //
