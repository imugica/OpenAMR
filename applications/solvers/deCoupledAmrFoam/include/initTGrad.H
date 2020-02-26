forAll(fluidRegions, i)
{
     if (caloricCycleControl[i].makeGrad)
     {
           forAll(solidNames, j)
           {
                if (solidNames[j]==caloricCycleControl[i].solidPhaseName())
                {
                      caloricCycleControl[i].initLinearTGrad(Tfluid[i], Tsolid[j]);
                }
           }
     }

}


