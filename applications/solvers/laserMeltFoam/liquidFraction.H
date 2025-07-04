{
  //correction for gT = 1
  //update deltaH with new temperature

  gT = gT + relax * (cM/Hf) * (T-invLiquidFraction(gT, Tl, Ts));
  
  forAll(mesh.cells(), cellI)
  {
      if( gT[cellI] > 1.0)
      {
          gT[cellI] = 1.0;
      }
      else if( gT[cellI] < 0.0)
      {
          gT[cellI] = 0.0;
      }
  }
  
  //necessary to trunctate gT on boundary field
  gT.correctBoundaryConditions();
  
  /*Info << "After gT update:\t" << "  Min(gT) = " << min(gT).value()
  << "  Max(gT) = " << max(gT).value() << endl;*/

  Info << "After gT update:\t" << "  Min(gT) = " << gMin(gT)
  << "  Max(gT) = " << gMax(gT) << endl;
}

