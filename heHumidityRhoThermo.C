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

\*---------------------------------------------------------------------------*/

#include "heHumidityRhoThermo.H"
#include "fvMatricesFwd.H"
#include "fvCFD.H"
#include "bound.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class BasicPsiThermo, class MixtureType>
void Foam::heHumidityRhoThermo<BasicPsiThermo, MixtureType>::calculate()
{
    const scalarField& hCells = this->he();
    const scalarField& pCells = this->p_;

    scalarField& TCells = this->T_.primitiveFieldRef();
    scalarField& psiCells = this->psi_.primitiveFieldRef();
    scalarField& rhoCells = this->rho_.primitiveFieldRef();
    scalarField& muCells = this->mu_.primitiveFieldRef();
    scalarField& alphaCells = this->alpha_.primitiveFieldRef();

    forAll(TCells, celli)
    {
        const typename MixtureType::thermoType& mixture_ =
            this->cellMixture(celli);

        TCells[celli] = mixture_.THE
        (
            hCells[celli],
            pCells[celli],
            TCells[celli]
        );

        psiCells[celli] = mixture_.psi(pCells[celli], TCells[celli]);
        rhoCells[celli] = mixture_.rho(pCells[celli], TCells[celli]);

        muCells[celli] = mixture_.mu(pCells[celli], TCells[celli]);
        alphaCells[celli] = mixture_.alphah(pCells[celli], TCells[celli]);
    }

    volScalarField::Boundary& pBf =
        this->p_.boundaryFieldRef();

    volScalarField::Boundary& TBf =
        this->T_.boundaryFieldRef();

    volScalarField::Boundary& psiBf =
        this->psi_.boundaryFieldRef();

    volScalarField::Boundary& rhoBf =
        this->rho_.boundaryFieldRef();

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
        fvPatchScalarField& ppsi = psiBf[patchi];
        fvPatchScalarField& prho = rhoBf[patchi];
        fvPatchScalarField& phe = heBf[patchi];
        fvPatchScalarField& pmu = muBf[patchi];
        fvPatchScalarField& palpha = alphaBf[patchi];

        if (pT.fixesValue())
        {
            forAll(pT, facei)
            {
                const typename MixtureType::thermoType& mixture_ =
                    this->patchFaceMixture(patchi, facei);

                phe[facei] = mixture_.HE(pp[facei], pT[facei]);

                ppsi[facei] = mixture_.psi(pp[facei], pT[facei]);
                prho[facei] = mixture_.rho(pp[facei], pT[facei]);
                pmu[facei] = mixture_.mu(pp[facei], pT[facei]);
                palpha[facei] = mixture_.alphah(pp[facei], pT[facei]);
            }
        }
        else
        {
            forAll(pT, facei)
            {
                const typename MixtureType::thermoType& mixture_ =
                    this->patchFaceMixture(patchi, facei);

                pT[facei] = mixture_.THE(phe[facei], pp[facei], pT[facei]);

                ppsi[facei] = mixture_.psi(pp[facei], pT[facei]);
                prho[facei] = mixture_.rho(pp[facei], pT[facei]);
                pmu[facei] = mixture_.mu(pp[facei], pT[facei]);
                palpha[facei] = mixture_.alphah(pp[facei], pT[facei]);
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
    calculate();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasicPsiThermo, class MixtureType>
Foam::heHumidityRhoThermo<BasicPsiThermo, MixtureType>::~heHumidityRhoThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicPsiThermo, class MixtureType>
void Foam::heHumidityRhoThermo<BasicPsiThermo, MixtureType>::correct()
{
    if (debug)
    {
        InfoInFunction << endl;
    }

    calculate();

    Info<< "   Solve transport equation for specific humidity\n";
    specificHumidityTransport();

    Info<< "   Calculate the saturation pressure of water\n";
    pSatH2O();

    Info<< "   Calculate the partial pressure of water\n";
    partialPressureH2O();

    Info<< "   Calculate the relative humidity\n";
    relHumidity();

    Info<< "   Calculate the water content\n";
    waterContent();

    Info<< "   Accounting for density change based on humidity\n";
    densityChange();

    //- Keep physical bounds
    limit();

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
        dimensionedScalar pre1
        (
            "pre1",
            dimensionSet(1,-1,-2,0,0,0,0),
            scalar(611.2)
        );

        dimensionedScalar pre2
        (
            "pre2",
            dimensionSet(0,0,0,0,0,0,0),
            scalar(17.62)
        );

        dimensionedScalar value1
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
        dimensionedScalar pre1
        (
            "pre1",
            dimensionSet(1,-1,-2,0,0,0,0),
            scalar(611.21)
        );

        dimensionedScalar value1
        (
            "value1",
            dimensionSet(0,0,0,0,0,0,0),
            scalar(18.678)
        );

        dimensionedScalar value2
        (
            "value2",
            dimensionSet(0,0,0,0,0,0,0),
            scalar(234.5)
        );

        dimensionedScalar value3
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
void Foam::heHumidityRhoThermo<BasicPsiThermo, MixtureType>::partialPressureH2O()
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
    const volScalarField& T = this->T_;
    const volScalarField& sH = this->specificHumidity_;

    volScalarField& pPH2O = this->partialPressureH2O_;

    pPH2O =
        (-p/RSpecificDryAir/T)
      * pow
        (
            (-1./(RSpecificDryAir*T))
          + ((1-pow(sH+vSmall, -1))* (1/(RSpecificH2O*T))),
           -1
        );
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
void Foam::heHumidityRhoThermo<BasicPsiThermo, MixtureType>::waterContent()
{
    const volScalarField& pPH2O = this->partialPressureH2O_;
    const volScalarField& pSatH2O = this->pSatH2O_;
    const volScalarField& p = this->p_;
    const volScalarField& T= this->T_;

    const dimensionedScalar RSpecificH2O
    (
        "gasConstantH2O",
        dimensionSet(0,2,-2,-1,0,0,0),
        scalar(461.51)
    );

    volScalarField& waterContent = this->waterContent_;
    waterContent = pPH2O / (RSpecificH2O * T);

    volScalarField& maxWaterContent = this->maxWaterContent_;
    maxWaterContent = pSatH2O / (RSpecificH2O * T);

    const dimensionedScalar RSpecificDryAir
    (
        "gasConstantDryAir",
        dimensionSet(0,2,-2,-1,0,0,0),
        scalar(287.058)
    );

    volScalarField& maxSpecificHumidity = this->maxSpecificHumidity_;
    maxSpecificHumidity =
        maxWaterContent / ((p-pSatH2O)/(RSpecificDryAir*T) + maxWaterContent);
}


template<class BasicPsiThermo, class MixtureType>
void Foam::heHumidityRhoThermo<BasicPsiThermo, MixtureType>::
specificHumidityTransport()
{
    volScalarField& specHum = this->specificHumidity_;
    volScalarField& muEff = this->muEff_;

    const volScalarField& mu = this->mu_;

    //- Old density
    const volScalarField& rho =
        this->db().objectRegistry::lookupObject<volScalarField>("rho");

    const surfaceScalarField& phi =
        this->db().objectRegistry::lookupObject<surfaceScalarField>("phi");

    const IOdictionary& turbProp =
        this->db().objectRegistry
        ::lookupObject<IOdictionary>("turbulenceProperties");

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

    fvScalarMatrix specHumEqn
    (
        fvm::ddt(rho, specHum)
      + fvm::div(phi, specHum)
     ==
        fvm::laplacian(muEff, specHum)
    );

    specHumEqn.relax();
    specHumEqn.solve();

    //- To keep physical range
    //  Defined between 0 and max water content based on saturation pressure
    //  Maximum is not limited yet
    bound
    (
        specHum,
        dimensionedScalar("tmp", dimensionSet(0,0,0,0,0,0,0), scalar(0))
    );
}


template<class BasicPsiThermo, class MixtureType>
void Foam::heHumidityRhoThermo<BasicPsiThermo, MixtureType>::densityChange()
{
    this->rho_ += this->waterContent_;
}


template<class BasicPsiThermo, class MixtureType>
void Foam::heHumidityRhoThermo<BasicPsiThermo, MixtureType>::limit()
{
    volScalarField& specHum = this->specificHumidity_;
    const volScalarField& maxSpecHum = this->maxSpecificHumidity_;

    label min = 0;
    label max = 0;

    forAll(specHum, cellI)
    {
        if (specHum[cellI] < 0)
        {
            specHum[cellI] = 0;
            min++;
        }
        else if (specHum[cellI] > maxSpecHum[cellI])
        {
            specHum[cellI] = maxSpecHum[cellI];
            max++;
        }
    }

    if (min > 0)
    {
        Info<< "Corrected " << min << " cells which were lower than 0\n";
    }
    if (max > 0)
    {
        Info<< "Corrected " << max << " cells which were higher than max\n";
    }
}


// ************************************************************************* //
