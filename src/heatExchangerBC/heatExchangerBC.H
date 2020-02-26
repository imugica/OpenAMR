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
 
 Class
     Foam::heatExchangerPatchScalarField
 
 Description
     Boundary condition modeling a perfect heat exchanger. It's operation
     is modelled by a single temperature value, which is calculated as the
     perfect thermal mixture of the fluid that crosses the boundary.
 
 SourceFiles
     heatExchangerBC.H
     heatExchangerBC.C
 
 \*---------------------------------------------------------------------------*/

#ifndef heatExchangerPatchScalarField_H
#define heatExchangerPatchScalarField_H

#include "mixedFvPatchFields.H"
#include "therCycleControl.H"
#include <string>
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
          Class heatExchangerPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class heatExchangerPatchScalarField
:
    public mixedFvPatchScalarField
{

public:

    //- Runtime scalar information
    TypeName("heatExchangerBC");

    // Constructors

        //- Construct from patch and internal field
        heatExchangerPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        heatExchangerPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  heatExchangerPatchScalarField onto a new patch
        heatExchangerPatchScalarField
        (
            const heatExchangerPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        heatExchangerPatchScalarField
        (
            const heatExchangerPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<scalar>> clone() const
        {
            return tmp<fvPatchField<scalar>>
            (
                new heatExchangerPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        heatExchangerPatchScalarField
        (
            const heatExchangerPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<scalar>> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<scalar>>
            (
                new heatExchangerPatchScalarField(*this, iF)
            );
        }


    // Member functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
