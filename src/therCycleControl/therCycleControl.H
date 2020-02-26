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
     Foam::therCycleControl
 
 Description
     Class that contains with useful methods to control and write data concerning
     AMR thermodynamic cyles.
 
 SourceFiles
     therCycleControl.H
     therCycleControl.C
 
 \*---------------------------------------------------------------------------*/
 
 #ifndef therCycleControl_H
 #define therCycleControl_H

 #include "fvCFD.H"
 #include <iostream>
 #include <iomanip>
 #include <fstream>
 #include <math.h>   
 #include "fvMesh.H"
 #include "volFields.H"
 #include "surfaceFields.H"
 #include "dimensionedScalar.H"
 #include "cycleIDvars.H"
 #include  "magneticProp.H"
 #include "regIOobject.H"
 
 // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
 namespace Foam
 {
 
 /*---------------------------------------------------------------------------*\
                          Class therCycleControl Declaration
 \*---------------------------------------------------------------------------*/
 
 class therCycleControl:
     public cycleIDvars,
     public regIOobject
 {
 private:
        // Get the starting time to know where tho read phi
         const dictionary gradDict_;
         const dictionary hexDict_;
         std::ofstream integVal_;
         const dimensionedScalar rhof_;
         const dimensionedScalar cf_;
         const word& fluidName_;
         double tolerance_ = GREAT;
         double cycleTol_ = 1E-3;

 public: 
// Static Data Members
 
     //- Run-time type information
     TypeName("therCycleControl");

     // Constructors
     therCycleControl(const IOobject&, const IOdictionary&, const IOdictionary&, const scalarField&, const word&, const double, const word&);
 
     //- Destructor
     virtual ~therCycleControl();
 
    //Data Memebers
      //- Cycle Vars
     bool tolReached = false;
     bool makeGrad;

     double coldSideT = readScalar(gradDict_.lookup("initialTcold"));
     double hotSideT = readScalar(gradDict_.lookup("initialThot"));
     double coldSideQ;
     double hotSideQ;
     double sumWmag;
   
     double cycleAvgColdSideT = 0;
     double cycleAvgHotSideT = 0;
     double cycleColdSideQ = 0;
     double cycleHotSideQ = 0;
     double cycleSumWmag = 0;

     double OldCycleAvgColdSideT = 0;
     double OldCycleAvgHotSideT  = 0;
     double OldCycleColdSideQ = 0;
     double OldCycleHotSideQ  = 0;
     double OldCycleSumWmag = 0;

     word coldName;
     word hotName;
     word units;

     // Member Functions
     void writeCycleVars(const double);

     void initLinearTGrad(volScalarField&, volScalarField&);

     void updateHEx(const volScalarField&, const surfaceScalarField&);

     bool checkCycle();
 
     void setCycleTol(double);

     word solidPhaseName();

     void calcWmag(const magneticProp&);
  
     double tolerance();

     bool writeData(Ostream&) const;

 };
 
 
 // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
 } // End namespace Foam
 
 // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
 #endif
 
 // ************************************************************************* //
