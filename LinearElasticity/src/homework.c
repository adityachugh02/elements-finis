#include "fem.h"




void geoMeshGenerate() {

    femGeo* theGeometry = geoGetGeometry();

    double w = theGeometry->LxPlate;
    double h = theGeometry->LyPlate;
    
    int ierr;
    int idPlate = gmshModelOccAddRectangle(0, 0, 0, 1.5, 1.2,-1, 0, &ierr);
    int idHole = gmshModelOccAddDisk(0.48, 0.67, 0, 0.1, 0.1, -1, NULL, 0, NULL, 0, &ierr);

    int idSubstractionHole = gmshModelOccAddDisk(0.20, 0.20, 0, 0.20, 0.20, -1, NULL, 0, NULL, 0, &ierr);
    int idSubstractionPlate = gmshModelOccAddRectangle(-1, -1, 0, 2.5, 1.20,-1, 0, &ierr);

    int idSecondSubstractionHole = gmshModelOccAddDisk(0.55, 0.25, 0, 0.17, 0.17, -1, NULL, 0, NULL, 0, &ierr);

    int idSubstractionWedgeCurvedPlate = gmshModelOccAddRectangle(-1, 0.42, 0, 2, 0.5,-1, 0, &ierr);
    int idSecondSubstractionWedgeCurvedPlate = gmshModelOccAddRectangle(-2, -2.5, 0, 2.55, 3,-1, 0, &ierr);
    int idSubstractionWedgeCurvedHole = gmshModelOccAddDisk(0, -0.70, 0, 1.3, 1.3, -1, NULL, 0, NULL, 0, &ierr);

    int idThirdSubstractionHole = gmshModelOccAddDisk(0.55, 0.25, 0, 0.67, 0.67, -1, NULL, 0, NULL, 0, &ierr);
    int idSecondSubstractionPlate = gmshModelOccAddRectangle(0, 0.60, 0, 0.9, 0.8,-1, 0, &ierr);

    int idFourthSubstractionHole = gmshModelOccAddDisk(0.97, 0.92, 0, 0.16, 0.16, -1, NULL, 0, NULL, 0, &ierr);

    int idFifthSubstractionHole = gmshModelOccAddDisk(0.5, 0.67, 0, 0.93, 0.93, -1, NULL, 0, NULL, 0, &ierr);
    int idThirdSubstractionPlate = gmshModelOccAddRectangle(1, 0, 0, 1, 1.5,-1, 0, &ierr);

    int idSixthSubstractionHole = gmshModelOccAddDisk(1.245, 0.98, 0, 0.13, 0.13, -1, NULL, 0, NULL, 0, &ierr);
    int idFourthSubstractionPlate = gmshModelOccAddRectangle(0.90, 0.98, 0, 0.8, 0.5,-1, 0, &ierr);

    int idHookPlateFix = gmshModelOccAddRectangle(0.395, 0.08, 0, 0.2, 0.2,-1, 0, &ierr);

    int plate[] = {2, idPlate};
    int hole[]  = {2, idHole};
    int substractionHole[]  = {2, idSubstractionHole};
    int substractionPlate[] = {2, idSubstractionPlate};
    int secondSubstractionHole[]  = {2, idSecondSubstractionHole};

    int substractionWedgeCurvedPlate[]  = {2, idSubstractionWedgeCurvedPlate};
    int secondSubstractionWedgeCurvedPlate[]  = {2, idSecondSubstractionWedgeCurvedPlate};
    int substractionWedgeCurvedHole[]  = {2, idSubstractionWedgeCurvedHole};

    int thirdSubstractionHole[]  = {2, idThirdSubstractionHole};
    int secondSubstractionPlate[]  = {2, idSecondSubstractionPlate};
    int fourthSubstractionHole[]  = {2, idFourthSubstractionHole};
    int fifthSubstractionHole[]  = {2, idFifthSubstractionHole};
    int thirdSubstractionPlate[]  = {2, idThirdSubstractionPlate};
    int sixthSubstractionHole[]  = {2, idSixthSubstractionHole};
    int fourthSubstractionPlate[]  = {2, idFourthSubstractionPlate};

    int hookPlateFix[]  = {2, idHookPlateFix};

    gmshModelOccCut(substractionPlate, 2, substractionHole, 2, NULL, NULL, NULL, NULL, NULL, -1, 1, 1, &ierr);
    gmshModelOccCut(plate, 2, substractionPlate, 2, NULL, NULL, NULL, NULL, NULL, -1, 1, 1, &ierr);
    gmshModelOccCut(plate, 2, hole, 2, NULL, NULL, NULL, NULL, NULL, -1, 1, 1, &ierr);
    gmshModelOccCut(plate, 2, secondSubstractionHole, 2, NULL, NULL, NULL, NULL, NULL, -1, 1, 1, &ierr);

    gmshModelOccCut(substractionWedgeCurvedHole, 2, substractionWedgeCurvedPlate, 2, NULL, NULL, NULL, NULL, NULL, -1, 1, 1, &ierr);
    gmshModelOccCut(substractionWedgeCurvedHole, 2, secondSubstractionWedgeCurvedPlate, 2, NULL, NULL, NULL, NULL, NULL, -1, 1, 1, &ierr);
    gmshModelOccCut(plate, 2, substractionWedgeCurvedHole, 2, NULL, NULL, NULL, NULL, NULL, -1, 1, 1, &ierr);

    gmshModelOccCut(secondSubstractionPlate, 2, thirdSubstractionHole, 2, NULL, NULL, NULL, NULL, NULL, -1, 1, 1, &ierr);
    gmshModelOccCut(plate, 2, secondSubstractionPlate, 2, NULL, NULL, NULL, NULL, NULL, -1, 1, 1, &ierr);
    gmshModelOccCut(plate, 2, fourthSubstractionHole, 2, NULL, NULL, NULL, NULL, NULL, -1, 1, 1, &ierr);
    gmshModelOccCut(fourthSubstractionPlate, 2, sixthSubstractionHole, 2, NULL, NULL, NULL, NULL, NULL, -1, 1, 1, &ierr);
    gmshModelOccCut(plate, 2, fourthSubstractionPlate, 2, NULL, NULL, NULL, NULL, NULL, -1, 1, 1, &ierr);
    gmshModelOccCut(thirdSubstractionPlate, 2, fifthSubstractionHole, 2, NULL, NULL, NULL, NULL, NULL, -1, 1, 1, &ierr);
    gmshModelOccCut(plate, 2, thirdSubstractionPlate, 2, NULL, NULL, NULL, NULL, NULL, -1, 1, 1, &ierr);

    gmshModelOccCut(plate, 2, hookPlateFix, 2, NULL, NULL, NULL, NULL, NULL, -1, 1, 1, &ierr);
    gmshModelOccSynchronize(&ierr); 
/*
    if (theGeometry->elementType == FEM_QUAD) {
        gmshOptionSetNumber("Mesh.SaveAll",1,&ierr);
        gmshOptionSetNumber("Mesh.RecombineAll",1,&ierr);
        gmshOptionSetNumber("Mesh.Algorithm",11,&ierr);  
        gmshOptionSetNumber("Mesh.SmoothRatio", 21.5, &ierr);  
        gmshOptionSetNumber("Mesh.RecombinationAlgorithm",1.0,&ierr); 
        gmshModelGeoMeshSetRecombine(2,1,45,&ierr);  
        gmshModelMeshGenerate(2,&ierr);  }
*/  
    //if (theGeometry->elementType == FEM_TRIANGLE) {
        gmshOptionSetNumber("Mesh.SaveAll",1,&ierr);
        gmshModelMeshGenerate(2,&ierr);  //}
 
    return;
}


double *femElasticitySolve(femProblem *theProblem)
{

    femFullSystem  *theSystem = theProblem->system;
    femIntegration *theRule = theProblem->rule;
    femDiscrete    *theSpace = theProblem->space;
    femGeo         *theGeometry = theProblem->geometry;
    femNodes       *theNodes = theGeometry->theNodes;
    femMesh        *theMesh = theGeometry->theElements;
    
    
    double x[4],y[4],phi[4],dphidxsi[4],dphideta[4],dphidx[4],dphidy[4];
    int iElem,iInteg,iEdge,i,j,d,map[4],mapX[4],mapY[4];
    
    int nLocal = theMesh->nLocalNode;

    double a   = theProblem->A;
    double b   = theProblem->B;
    double c   = theProblem->C;      
    double rho = theProblem->rho;
    double g   = theProblem->g;
    double **A = theSystem->A;
    double *B  = theSystem->B;
    
    
  //
  //  A faire :-)
  //                
    for (iElem = 0; iElem < theMesh->nElem; iElem++) {
        for (j=0; j < nLocal; j++) {
            map[j] = theMesh->elem[iElem*nLocal+j];
            mapX[j] = 2*map[j];
            mapY[j] = 2*map[j] + 1;
            x[j] = theNodes->X[map[j]];
            y[j] = theNodes->Y[map[j]];}
            for (iInteg=0; iInteg < theRule->n; iInteg++) {
                double xsi = theRule->xsi[iInteg];
                double eta = theRule->eta[iInteg];
                double weight = theRule->weight[iInteg];
                femDiscretePhi2(theSpace,xsi,eta,phi);
                femDiscreteDphi2(theSpace,xsi,eta,dphidxsi,dphideta);
                double dxdxsi = 0.0;
                double dxdeta = 0.0;
                double dydxsi = 0.0;
                double dydeta = 0.0;
                for (i = 0; i < theSpace->n; i++) {
                    dxdxsi += x[i]*dphidxsi[i];
                    dxdeta += x[i]*dphideta[i];
                    dydxsi += y[i]*dphidxsi[i];
                    dydeta += y[i]*dphideta[i]; 
                }
                double jac = fabs(dxdxsi * dydeta - dxdeta * dydxsi);
                for (i = 0; i < theSpace->n; i++) {
                    dphidx[i] = (dphidxsi[i] * dydeta - dphideta[i] * dydxsi) / jac;
                    dphidy[i] = (dphideta[i] * dxdxsi - dphidxsi[i] * dxdeta) / jac; 
                }
                for (i = 0; i < theSpace->n; i++) {
                    for(j = 0; j < theSpace->n; j++) {
                        A[mapX[i]][mapX[j]] += (dphidx[i] * a * dphidx[j] +
                        dphidy[i] * c * dphidy[j]) * jac * weight;
                        A[mapX[i]][mapY[j]] += (dphidx[i] * b * dphidy[j] +
                        dphidy[i] * c * dphidx[j]) * jac * weight;
                        A[mapY[i]][mapX[j]] += (dphidy[i] * b * dphidx[j] +
                        dphidx[i] * c * dphidy[j]) * jac * weight;
                        A[mapY[i]][mapY[j]] += (dphidy[i] * a * dphidy[j] +
                        dphidx[i] * c * dphidx[j]) * jac * weight; 
                    }
                }
                for (i = 0; i < theSpace->n; i++) {
                    B[mapY[i]] -= phi[i] * g * rho * jac * weight; 
                }
            }
        }
        int *theConstrainedNodes = theProblem->constrainedNodes;
        for (int i=0; i < theSystem->size; i++) {
            if (theConstrainedNodes[i] != -1) {
                double value = theProblem->conditions[theConstrainedNodes[i]]->value;
                femFullSystemConstrain(theSystem,i,value); 
            }
        }
    return femFullSystemEliminate(theSystem);
}
