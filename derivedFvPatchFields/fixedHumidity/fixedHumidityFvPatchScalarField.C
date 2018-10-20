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
#include "basicThermo.H"
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
    quantityType_("relative"),
    absolute_(false)
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
    quantityType_(ptf.quantityType_),
    absolute_(ptf.absolute_)
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
    quantityType_(dict.lookupOrDefault<word>("quantity", "relative")),
    absolute_(false)
{
    if (quantityType_ == "absolute")
    {
        absolute_ = true;
    }

    volScalarField& mH2O = this->db().objectRegistry::lookupObject<volScalarField>("massH2O");
}


Foam::fixedHumidityFvPatchScalarField::
fixedHumidityFvPatchScalarField
(
    const fixedHumidityFvPatchScalarField& tppsf
)
:
    fixedValueFvPatchScalarField(tppsf),
    quantityType_(tppsf.quantityType_),
    absolute_(tppsf.absolute_)
{}


Foam::fixedHumidityFvPatchScalarField::
fixedHumidityFvPatchScalarField
(
    const fixedHumidityFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(tppsf, iF),
    quantityType_(tppsf.quantityType_),
    absolute_(tppsf.absolute_)
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

    //const scalarField& pw = thermo.p().boundaryField()[patchi];
    //fvPatchScalarField& Tw =
    //    const_cast<fvPatchScalarField&>(thermo.T().boundaryField()[patchi]);
    //Tw.evaluate();
    //operator==(thermo.he(pw, Tw, patchi));

    fixedValueFvPatchScalarField::updateCoeffs();
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
