Info<<"Updating Laser Parameters \n"<<endl;

flb = 0.0;

forAll(cellColumns, colI)
{
    DynamicList<label>& curCol = cellColumns[colI];
    scalar laserFraction = 1.0;
    scalar oldLaserFraction = 0.0;
    
    forAll(curCol,colJ)
    {
        label cellID = curCol[colJ];
        if(pos(laserFraction))
        {
            oldLaserFraction = laserFraction;
            laserFraction -= alphaM[cellID];
            
            if(pos(laserFraction))
            {
                flb[cellID] = alphaM[cellID];
            }
            else
            {
                flb[cellID] = oldLaserFraction;
            }
        }
    }
}


flb.correctBoundaryConditions();
