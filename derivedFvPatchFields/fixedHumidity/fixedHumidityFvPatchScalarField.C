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

#include "fvPatchFieldMapper.H"
#include "volFields.H"
//#include "basicThermo.H"
#include "humidityRhoThermo.C"
#include "addToRunTimeSelectionTable.H"
#include "fixedHumidityFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fixedHumidityFvPatchScalarField::
fixedHumidityFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    mode_("relative"),
    humidity_(0.0)
{}


Foam::fixedHumidityFvPatchScalarField::
fixedHumidityFvPatchScalarField
(
    const fixedHumidityFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    mode_(ptf.mode_),
    humidity_(ptf.humidity_)
{}


Foam::fixedHumidityFvPatchScalarField::
fixedHumidityFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict),
    mode_(dict.lookupOrDefault<word>("mode", "relative")),
    humidity_(readScalar(dict.lookup("humidity")))
{
//    volScalarField& mH2O =
//        this->db().objectRegistry::lookupObjectRef<volScalarField>("massH2O");

    if (mode_ == "absolute")
    {
        //operator==(humidity_);
    }

    else if (mode_ == "relative")
    {

        //operator==(humidity_);
    }

    else
    {
        FatalErrorInFunction
            << "The specified type is not supported '"
            << mode_ << "'. Supported are 'relative' or 'absolute'"
            << exit(FatalError);
    }
}


Foam::fixedHumidityFvPatchScalarField::
fixedHumidityFvPatchScalarField
(
    const fixedHumidityFvPatchScalarField& tppsf
)
:
    fixedValueFvPatchScalarField(tppsf),
    mode_(tppsf.mode_),
    humidity_(tppsf.humidity_)
{}


Foam::fixedHumidityFvPatchScalarField::
fixedHumidityFvPatchScalarField
(
    const fixedHumidityFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(tppsf, iF),
    mode_(tppsf.mode_),
    humidity_(tppsf.humidity_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fixedHumidityFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const basicThermo& thermo = basicThermo::lookupThermo(*this);
    const label patchi = patch().index();

    scalarField massH2O = mH2O(thermo, patchi);

    //const scalarField& pw = thermo.p().boundaryField()[patchi];
    //fvPatchScalarField& Tw =
    //    const_cast<fvPatchScalarField&>(thermo.T().boundaryField()[patchi]);
    //Tw.evaluate();
    //operator==(thermo.he(pw, Tw, patchi));

    fixedValueFvPatchScalarField::updateCoeffs();
}


Foam::scalarField Foam::fixedHumidityFvPatchScalarField::mH2O
(
    const basicThermo& thermo,
    const label patchi
)
{
    //- a) Calc saturation pressure
    const word method = humidityRhoThermo::method_;
    const scalarField& Tpatch = thermo.T().boundaryField()[patchi];
    const scalarField theta = Tpatch - 273.15; //(Tpatch.size(), scalar(0));

    scalarField pSatH2O(Tpatch.size(), scalar(0));

    if (method == "magnus")
    {
        scalar pre1 = 611.2;
        scalar pre2 = 17.62;
        scalar value1 = 243.12;

        forAll(pSatH2O, facei)
        {
            pSatH2O[facei] =
                pre1*exp((pre2*(theta[facei]))/(value1+theta[facei]));
        }
    }
    //  Buck formula [1996]
    //  Valid between 0 to 100 degC and 1013.25 hPa
    //  Very accurate between 0 degC and 50 degC
    else if (method == "buck")
    {
        scalar pre1 = 611.21;
        scalar value1 = 18.678;
        scalar value2 = 234.5;
        scalar value3 = 257.14;

        forAll(pSatH2O, facei)
        {
            scalar TdC = theta[facei];

            pSatH2O[facei] =
                pre1*exp(((value1-TdC)/value2)*TdC/(value3+TdC));
        }
    }

    //- b) Calc partial pressure of water
    scalarField partialPressureH2O(Tpatch.size(), scalar(0));
    partialPressureH2O = humidity_ * pSatH2O;

    //- c) Calc absolute humidity
    scalarField absHumidity(Tpatch.size(), scalar(0));

    scalar RH2O = 461.51;

    absHumidity = partialPressureH2O / RH2O / Tpatch; 

    //- d) Calc water content / mass
    scalarField mH2O(Tpatch.size(), scalar(0));
//    mH2O = absHumidity * V;

    return pSatH2O;

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        fixedHumidityFvPatchScalarField
    );
}

// ************************************************************************* //
