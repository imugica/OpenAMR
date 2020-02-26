 /*---------------------------------------------------------------------------*\
   =========                 |
   \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
    \\    /   O peration     | 
     \\  /    A nd           | Copyright (C) 2020 Ibai Mugica
      \\/     M anipulation  |
 -------------------------------------------------------------------------------
 License
     This file is derivative work of OpenFOAM.
 
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

#include "heatExchangerBC.H"
#include "fvPatchFieldMapper.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::heatExchangerPatchScalarField
::heatExchangerPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchField<scalar>(p, iF)
{
    this->refValue() = 0.0;
    this->refGrad() = 0.0;
    this->valueFraction() = 0.0;
}


Foam::heatExchangerPatchScalarField
::heatExchangerPatchScalarField
(
    const heatExchangerPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper)
{}


Foam::heatExchangerPatchScalarField
::heatExchangerPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF)
{
    this->refValue() = 0.0;
    fvPatchScalarField::operator=(this->patchInternalField());
    this->refGrad() = 0.0;
    this->valueFraction() = 0.0;
}


Foam::heatExchangerPatchScalarField
    ::heatExchangerPatchScalarField
(
    const heatExchangerPatchScalarField& ptf
)
:
    mixedFvPatchScalarField(ptf)
{}


Foam::heatExchangerPatchScalarField
    ::heatExchangerPatchScalarField
(
    const heatExchangerPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::heatExchangerPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

        const therCycleControl& caloricCycleControl = db().lookupObject<therCycleControl>("caloricCycleControl");

        const fvsPatchField<scalar>& phip = patch().lookupPatchField<surfaceScalarField, scalar>("phi"+caloricCycleControl.fileN());

	if (patch().name() == caloricCycleControl.coldName)
	{
	        this->refValue() = caloricCycleControl.coldSideT;
	        this->valueFraction() = 1.0;
	}
	else if (patch().name() == caloricCycleControl.hotName)
	{
	        this->refValue() = caloricCycleControl.hotSideT;
	        this->valueFraction() = 1.0;
	}

	if (gSum(phip)>0)
        {
	         this->refGrad() = 0;
	         this->valueFraction() = 0.0;	
        }

        mixedFvPatchScalarField::updateCoeffs();
}


void Foam::heatExchangerPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        heatExchangerPatchScalarField
    );
}

// ************************************************************************* //
