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
     Foam::movingMagnetBC
 
 Description
     Boundary condition representing the poles of a moving magnet. The 
     magnetic field of the poles is represented by a Neumann condition 
     for the scalar potential. Magnetic permeability is a field.
 
 SourceFiles
     movingMagnetBC.H
     movingMagnetBC.C
 
 \*---------------------------------------------------------------------------*/

#ifndef movingMagnetBC_H
#define movingMagnetBC_H

#include "mixedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
          Class movingMagnetBC Declaration
\*---------------------------------------------------------------------------*/

class movingMagnetBC
:
    public mixedFvPatchScalarField
{

private:

    bool condition1_=false;
    bool condition2_=false;
    bool condition3_=false;
    bool condition4_=false;

public:

    //- Runtime scalar information
    TypeName("movingMagnetBC");


    // Constructors

        //- Construct from patch and internal field
        movingMagnetBC
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        movingMagnetBC
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  movingMagnetBC onto a new patch
        movingMagnetBC
        (
            const movingMagnetBC&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        movingMagnetBC
        (
            const movingMagnetBC&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<scalar>> clone() const
        {
            return tmp<fvPatchField<scalar>>
            (
                new movingMagnetBC(*this)
            );
        }

        //- Construct as copy setting internal field reference
        movingMagnetBC
        (
            const movingMagnetBC&,
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
                new movingMagnetBC(*this, iF)
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
