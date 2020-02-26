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
     Foam::cycleIDvars
 
 Description
     This class keeps track of the variables that are relative to the
     timing of the thermodynamic cycle. 
 
 SourceFiles
     cycleIDvars.H
     cycleIDvars.C
 
 \*---------------------------------------------------------------------------*/
 
 #ifndef cycleIDvars_H
 #define cycleIDvars_H

#include "fvCFD.H"
#include "scalarField.H"
#include <string>
#include <math.h>   
 
 // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
 namespace Foam
 {
 
 /*---------------------------------------------------------------------------*\
                         Class cycleIDvars Declaration
 \*---------------------------------------------------------------------------*/
 
 class cycleIDvars
 {
 
 private:
          //- Time line
          const scalarField timeLine_;
          double dt_;
          //- Cycle period
          double tau_;
          const word firstTime_;
          //- Relative time (to the cycle period)
          double trel_;
          int fileN_;
	  double tcompareLow_;
	  double tcompareUp_;
          int Nsteps_;

 public:
 
// Static Data Members
     //- Run-time type information
     TypeName("cycleIDvars");
 
     // Constructors
     cycleIDvars(const scalarField&, const word, const double, const bool);
     cycleIDvars();
 
     //- Destructor
     virtual ~cycleIDvars();

     // Member Functions
     void update(double, const bool);

     std::string fileN() const;

     void setDeltaT(scalar);

     int iphiSumRelStep();

     double dt();
     double tau();
     double trel();
     word firstTime() ;

 };

 
 // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
 } // End namespace Foam
 
 // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
 #endif
 
 // ************************************************************************* //
