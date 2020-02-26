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
     Foam::magneticProp
 
 Description
     Class that contains fields and interpolation tables that dependend on 
     the internal magnetic field.

 SourceFiles
     magneticProp.H
     magneticProp.C
 
 \*---------------------------------------------------------------------------*/
 
 #ifndef magneticProp_H
 #define magneticProp_H
 

 #include "fvMesh.H"
 #include "volFields.H"
 #include "cycleIDvars.H"
 #include "interpolation2DTable.H"

 // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
 namespace Foam
 {
 
 /*---------------------------------------------------------------------------*\
                          Class magneticControl Declaration
 \*---------------------------------------------------------------------------*/
 
 class magneticProp:
     public cycleIDvars 
 {
 private:

         const bool symH_;
         const interpolation2DTable<scalar> CHtable_;
         const interpolation2DTable<scalar> dMdTtable_; //Specific quantity: divided by density [J/(T*kg*K)]
         const interpolation2DTable<scalar> Mtable_;  ////Specific quantity: divided by density [J/(T*kg)]
	 const word firstTime_;
	 const fvMesh& solidMesh_;
         const double rho_;

 public:
 
// Static Data Members
 
     //- Run-time type information
     TypeName("magneticProp");

     // Constructors
     magneticProp(const fvMesh& , const IOdictionary& ,  const word&, const double, const double);
     magneticProp( const fvMesh&, const word&, const double, const double);
     //- Destructor
    virtual ~magneticProp();
 
    //Data Memebers
    volScalarField Hmag;
    volScalarField CH;
    volScalarField CT;
    volScalarField Wmag;

     // Member Functions

    void calcHmag(const volVectorField& );

    void updateCH(const volScalarField&);
    void updateCT(const volScalarField&);
    void updateWmag(const volScalarField&);

 };
 
 
 // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
 } // End namespace Foam
 
 // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
 #endif
 
 // ************************************************************************* //
