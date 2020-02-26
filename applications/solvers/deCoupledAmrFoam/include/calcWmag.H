forAll(fluidRegions, i)
{
           forAll(solidNames, j)
           {
                if (solidNames[j]==caloricCycleControl[i].solidPhaseName())
                {
                      magneticVars[j].updateWmag(Tsolid[j]);
                      caloricCycleControl[i].calcWmag(magneticVars[j]);
                }
           }
}
