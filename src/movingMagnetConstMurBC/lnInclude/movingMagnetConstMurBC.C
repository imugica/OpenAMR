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

#include "movingMagnetConstMurBC.H"
#include "fvPatchFieldMapper.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "electromagneticConstants.H"
#include "interpolateSplineXY.H"
#include "Time.H"
#include <math.h>  
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::movingMagnetConstMurBC
::movingMagnetConstMurBC
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


Foam::movingMagnetConstMurBC
::movingMagnetConstMurBC
(
    const movingMagnetConstMurBC& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper)
{}


Foam::movingMagnetConstMurBC
::movingMagnetConstMurBC
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


Foam::movingMagnetConstMurBC
    ::movingMagnetConstMurBC
(
    const movingMagnetConstMurBC& ptf
)
:
    mixedFvPatchScalarField(ptf)
{}


Foam::movingMagnetConstMurBC
    ::movingMagnetConstMurBC
(
    const movingMagnetConstMurBC& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::movingMagnetConstMurBC::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    this->valueFraction() =0;

    const fvMesh& mesh = patch().boundaryMesh().mesh();

    const IOdictionary& transDict=mesh.lookupObject<IOdictionary>("magneticProperties");
    const dimensionedScalar& mur=transDict.lookup("mur");

    const Time& time = db().time();
    const IOdictionary& bcDict=mesh.lookupObject<IOdictionary>("magnet");
    const word magnetType=bcDict.lookup("magnetType");

    const vector offSet=bcDict.lookup("offSet");
    const vector Mrem=bcDict.lookup("Mrem");
    const scalarField& pos=bcDict.lookup("Position");
    const scalarField& tpos=bcDict.lookup("time");

    const double relTime=((time.value()/tpos[tpos.size()-1])-floor(time.value()/tpos[tpos.size()-1]))*tpos[tpos.size()-1];
    const double relPos= interpolateSplineXY(relTime,tpos,pos);	

    const vector Lp=bcDict.lookup("Lp");
    const vector Lt=bcDict.lookup("Lt");


    if (magnetType=="Linear-square")
    {
	  	forAll(patch(),cellI)
	  	{
                      condition1_=false;
                      condition2_=false;
                      condition3_=false;
                      condition4_=false;

                      condition1_=(patch().Cf()[cellI] & (Lp/mag(Lp)))< ((offSet & (Lp/mag(Lp))) + relPos + mag(Lp));
                      condition2_=(patch().Cf()[cellI] & (Lp/mag(Lp)))> ((offSet & (Lp/mag(Lp))) + relPos -mag(Lp));
                      condition3_=(patch().Cf()[cellI] & (Lt/mag(Lt)))< ((offSet & (Lt/mag(Lt)))  + mag(Lt));
                      condition4_=(patch().Cf()[cellI] & (Lt/mag(Lt)))> ((offSet & (Lt/mag(Lt)))  - mag(Lt));

                      scalar& localpsiGrad=this->refGrad()[cellI];

			if (condition1_ && condition2_ && condition3_ && condition4_)
			{
				localpsiGrad=(Mrem & patch().Sf()[cellI])/(patch().magSf()[cellI]*mur.value());
			}
			else
			{
				localpsiGrad=0; //ZeroGradient
			}
	  	}
     }

     else if (magnetType=="Linear-square-cyclic")
     {

		double Nx=0;
		double nx=0;

		div_t divresult;

	  	forAll(patch(),cellI)
	  	{
                      condition1_=false;
                      condition2_=false;
                      condition3_=false;
                      condition4_=false;

			Nx=fabs((patch().Cf()[cellI] & (Lp/mag(Lp)))-(offSet & (Lp/mag(Lp))) -relPos )/(2*mag(Lp))-0.5;
			nx=fabs((patch().Cf()[cellI] & (Lp/mag(Lp)))-(offSet & (Lp/mag(Lp)))- relPos)/mag(Lp);

			divresult = div (floor(Nx),2);

			if ((divresult.rem > 0) || (nx<1))
			{
				condition1_=true;
				condition2_=true;
			}
			else
			{
				condition1_=false;
				condition2_=false;
			}
			
			condition3_=(patch().Cf()[cellI] & (Lt/mag(Lt)))< ((offSet & (Lt/mag(Lt)))  + mag(Lt));
			condition4_=(patch().Cf()[cellI] & (Lt/mag(Lt)))> ((offSet & (Lt/mag(Lt)))  - mag(Lt));

			scalar& localpsiGrad=this->refGrad()[cellI];

			if (condition1_ && condition2_ && condition3_ && condition4_)
			{

				localpsiGrad=(Mrem & patch().Sf()[cellI])/(patch().magSf()[cellI]*mur.value());

			}
			else
			{
				localpsiGrad=0; //ZeroGradient
			}
	  	}
	}
	else if (magnetType=="Linear-circle")
	{
	
		double Raux=mag(Lp);

		forAll(patch(),cellI)
		{
                        condition1_=false;

                        double CfProj_p = (patch().Cf()[cellI] & (Lp/mag(Lp))) - (offSet & (Lp/mag(Lp))) -relPos;
                        double CfProj_t = (patch().Cf()[cellI] & (Lt/mag(Lt))) - (offSet & (Lt/mag(Lt)));

			vector CfProj = CfProj_p*(Lp/mag(Lp)) + CfProj_t*(Lt/mag(Lt));

			condition1_=( mag(CfProj) < Raux );
		
			scalar& localpsiGrad=this->refGrad()[cellI];

			if (condition1_)
			{
				localpsiGrad=(Mrem & patch().Sf()[cellI])/(patch().magSf()[cellI]*mur.value());
			}
			else
			{
				localpsiGrad=0; //ZeroGradient
			}
		}

	}

        mixedFvPatchScalarField::updateCoeffs();
}


void Foam::movingMagnetConstMurBC::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
  
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        movingMagnetConstMurBC
    );
}

// ************************************************************************* //
