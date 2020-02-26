forAll(fluidRegions, i)
{
      caloricCycleControl[i].update( runTime.value(), false );
}

forAll(solidRegions, i)
{
     magneticVars[i].update( runTime.value(), HDictTimeLine[i].lookup("symmetric"));
}

