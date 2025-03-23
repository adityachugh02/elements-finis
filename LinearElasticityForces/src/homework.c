#include "fem.h"

double hermiteInterpolation(double d, double d_star, double h0, double h_star) {
    if (d >= d_star) return h_star;
    double t = d / d_star;
    return h0 + (h_star - h0) * (3 * t * t - 2 * t * t * t);
}


double geoSize(double x, double y){

    double h = 0.1;

    // hook inside
    double x0 = 0.55;
    double y0 = 0.25;
    double r0 = 0.17;
    double h0 = 0.02;
    double d0 = 0.5;
  
    // hole
    double x1 = 0.48;
    double y1 = 0.67;
    double r1 = 0.1;
    double h1 = 0.02;
    double d1 = 0.5;

    // hook outside
    /*
    double x2 = 0.20;
    double y2 = 0.20;
    double r2 = 0.20;
    double h2 = 0.01;
    double d2 = 0.2;
    */

//
//     A modifier !
//     
// Your contribution starts here ....
//
    // Calcul des distances aux cercles
    double d_hook_inside = fmax(0, sqrt((x - x0) * (x - x0) + (y - y0) * (y - y0)) - r0);
    double d_hole  = fmax(0, sqrt((x - x1) * (x - x1) + (y - y1) * (y - y1)) - r1);
    //double d_hook_outside = fmax(0, sqrt((x - x2) * (x - x2) + (y - y2) * (y - y2)) - r2);

    // Application de l'interpolation d'Hermite
    double h_hook_inside = hermiteInterpolation(d_hook_inside, d0, h0, h);
    double h_hole  = hermiteInterpolation(d_hole, d1, h1, h);
    //double h_hook_outside = hermiteInterpolation(d_hook_outside, d2, h2, h);

    // Retourner la taille minimale requise pour assurer un bon raffinement
    h = fmin(h_hook_inside, h_hole);
    //h = fmin(h, h_hook_outside);

     
    return h;
    
//   
// Your contribution ends here :-)
//

}

void geoMeshGenerate() {

    
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
        //gmshOptionSetNumber("Mesh.SaveAll",1,&ierr);
        //gmshModelMeshGenerate(2,&ierr);  //}

    geoSetSizeCallback(geoSize);
                                  
    gmshModelOccSynchronize(&ierr);       
    gmshOptionSetNumber("Mesh.SaveAll", 1, &ierr);
    gmshModelMeshGenerate(2, &ierr); 
 
    return;
}

void femElasticityAssembleElements(femProblem *theProblem){
    femFullSystem  *theSystem = theProblem->system;
    femIntegration *theRule = theProblem->rule;
    femDiscrete    *theSpace = theProblem->space;
    femGeo         *theGeometry = theProblem->geometry;
    femNodes       *theNodes = theGeometry->theNodes;
    femMesh        *theMesh = theGeometry->theElements;
    femMesh        *theEdges = theGeometry->theEdges;
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
    
    
    for (iElem = 0; iElem < theMesh->nElem; iElem++) {
        for (j=0; j < nLocal; j++) {
            map[j]  = theMesh->elem[iElem*nLocal+j];
            mapX[j] = 2*map[j];
            mapY[j] = 2*map[j] + 1;
            x[j]    = theNodes->X[map[j]];
            y[j]    = theNodes->Y[map[j]];} 
        
        for (iInteg=0; iInteg < theRule->n; iInteg++) {    
            double xsi    = theRule->xsi[iInteg];
            double eta    = theRule->eta[iInteg];
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
                dydeta += y[i]*dphideta[i]; }
            double jac = fabs(dxdxsi * dydeta - dxdeta * dydxsi);
            
            for (i = 0; i < theSpace->n; i++) {    
                dphidx[i] = (dphidxsi[i] * dydeta - dphideta[i] * dydxsi) / jac;       
                dphidy[i] = (dphideta[i] * dxdxsi - dphidxsi[i] * dxdeta) / jac; }            
            for (i = 0; i < theSpace->n; i++) { 
                for(j = 0; j < theSpace->n; j++) {
                    A[mapX[i]][mapX[j]] += (dphidx[i] * a * dphidx[j] + 
                                            dphidy[i] * c * dphidy[j]) * jac * weight;                                                                                            
                    A[mapX[i]][mapY[j]] += (dphidx[i] * b * dphidy[j] + 
                                            dphidy[i] * c * dphidx[j]) * jac * weight;                                                                                           
                    A[mapY[i]][mapX[j]] += (dphidy[i] * b * dphidx[j] + 
                                            dphidx[i] * c * dphidy[j]) * jac * weight;                                                                                            
                    A[mapY[i]][mapY[j]] += (dphidy[i] * a * dphidy[j] + 
                                            dphidx[i] * c * dphidx[j]) * jac * weight; }}
             for (i = 0; i < theSpace->n; i++) {
                B[mapY[i]] -= phi[i] * g * rho * jac * weight; }}} 
}


void femElasticityAssembleNeumann(femProblem *theProblem)
{
    femFullSystem   *theSystem  = theProblem->system;
    femIntegration  *theRule    = theProblem->ruleEdge;  
    femDiscrete     *theSpace   = theProblem->spaceEdge;
    femGeo          *theGeo     = theProblem->geometry;
    femMesh         *theEdges   = theGeo->theEdges;
    femNodes        *theNodes   = theGeo->theNodes;
    double          *B          = theSystem->B;

    for (int iBnd = 0; iBnd < theProblem->nBoundaryConditions; iBnd++) {

        femBoundaryCondition *theCondition = theProblem->conditions[iBnd];
        femBoundaryType type = theCondition->type;
        double value = theCondition->value;

        if (type == NEUMANN_X || type == NEUMANN_Y) {

            femDomain *theDomain = theCondition->domain;

            for (int iEdgeIndex = 0; iEdgeIndex < theDomain->nElem; iEdgeIndex++) {

                int iEdge = theDomain->elem[iEdgeIndex];

                int node0 = theEdges->elem[2*iEdge + 0];
                int node1 = theEdges->elem[2*iEdge + 1];

                double x0 = theNodes->X[node0];
                double y0 = theNodes->Y[node0];
                double x1 = theNodes->X[node1];
                double y1 = theNodes->Y[node1];

                double dx = x1 - x0;
                double dy = y1 - y0;
                double length = sqrt(dx*dx + dy*dy);
                double jac = length / 2.0;

                int mapU[2];
                if (type == NEUMANN_X) {
                    mapU[0] = 2 * node0;     
                    mapU[1] = 2 * node1; 
                }
                else {
                    mapU[0] = 2 * node0 + 1; 
                    mapU[1] = 2 * node1 + 1; 
                }

                for (int iInteg = 0; iInteg < theRule->n; iInteg++) {
                    double xsi    = theRule->xsi[iInteg];
                    double weight = theRule->weight[iInteg];

                    double phi0 = 0.5*(1.0 - xsi);
                    double phi1 = 0.5*(1.0 + xsi);

                    B[ mapU[0] ] += phi0 * value * jac * weight;
                    B[ mapU[1] ] += phi1 * value * jac * weight;
                }
            }
        }
    }
}

double *femElasticitySolve(femProblem *theProblem)
{
    femFullSystem  *theSystem   = theProblem->system;
    femGeo         *theGeometry = theProblem->geometry;
    femMesh        *theMesh     = theGeometry->theElements;
    femIntegration *theRule     = theProblem->rule;
    femDiscrete    *theSpace    = theProblem->space;
    femNodes       *theNodes    = theGeometry->theNodes;
    double **A = theSystem->A;
    double *B  = theSystem->B;
    
    double x[4], y[4], phi[4], dphidxsi[4], dphideta[4], dphidx[4], dphidy[4];
    int iElem, iInteg, i, j;
    int map[4], mapX[4], mapY[4];
    int nLocal = theMesh->nLocalNode;
    double a   = theProblem->A;
    double b   = theProblem->B;
    double c   = theProblem->C;      
    double rho = theProblem->rho;
    double g   = theProblem->g;

    for (iElem = 0; iElem < theMesh->nElem; iElem++) {
        for (j = 0; j < nLocal; j++) {
            map[j]  = theMesh->elem[iElem * nLocal + j];
            mapX[j] = 2 * map[j];
            mapY[j] = 2 * map[j] + 1;
            x[j]    = theNodes->X[map[j]];
            y[j]    = theNodes->Y[map[j]];
        }
        for (iInteg = 0; iInteg < theRule->n; iInteg++) {
            double xsi    = theRule->xsi[iInteg];
            double eta    = theRule->eta[iInteg];
            double weight = theRule->weight[iInteg];
            femDiscretePhi2(theSpace, xsi, eta, phi);
            femDiscreteDphi2(theSpace, xsi, eta, dphidxsi, dphideta);
            
            double dxdxsi = 0.0, dxdeta = 0.0, dydxsi = 0.0, dydeta = 0.0;
            for (i = 0; i < theSpace->n; i++) {
                dxdxsi += x[i] * dphidxsi[i];
                dxdeta += x[i] * dphideta[i];
                dydxsi += y[i] * dphidxsi[i];
                dydeta += y[i] * dphideta[i];
            }
            double jac = fabs(dxdxsi * dydeta - dxdeta * dydxsi);
            
            for (i = 0; i < theSpace->n; i++) {
                dphidx[i] = (dphidxsi[i] * dydeta - dphideta[i] * dydxsi) / jac;
                dphidy[i] = (dphideta[i] * dxdxsi - dphidxsi[i] * dxdeta) / jac;
            }            
            for (i = 0; i < theSpace->n; i++) {
                for (j = 0; j < theSpace->n; j++) {
                    A[mapX[i]][mapX[j]] += (dphidx[i] * a * dphidx[j] + dphidy[i] * c * dphidy[j]) * jac * weight;
                    A[mapX[i]][mapY[j]] += (dphidx[i] * b * dphidy[j] + dphidy[i] * c * dphidx[j]) * jac * weight;
                    A[mapY[i]][mapX[j]] += (dphidy[i] * b * dphidx[j] + dphidx[i] * c * dphidy[j]) * jac * weight;
                    A[mapY[i]][mapY[j]] += (dphidy[i] * a * dphidy[j] + dphidx[i] * c * dphidx[j]) * jac * weight;
                }
            }
            for (i = 0; i < theSpace->n; i++) {
                B[mapY[i]] -= phi[i] * g * rho * jac * weight;
            }
        }
    }

    femElasticityAssembleNeumann(theProblem);

    int *theConstrainedNodes = theProblem->constrainedNodes;
    for (i = 0; i < theSystem->size; i++) {
        if (theConstrainedNodes[i] != -1) {
            double value = theProblem->conditions[theConstrainedNodes[i]]->value;
            femFullSystemConstrain(theSystem, i, value);
        }
    }
    
    double *tempSol = femFullSystemEliminate(theSystem);
    int size = theSystem->size;

    theProblem->soluce = malloc(sizeof(double) * size);
    for (i = 0; i < size; i++) {
        theProblem->soluce[i] = tempSol[i];
    }
    
    return theProblem->soluce;
}


double *femElasticityForces(femProblem *theProblem)
{
    femFullSystem  *S     = theProblem->system;
    femIntegration *rule  = theProblem->rule;
    femDiscrete    *space = theProblem->space;
    femGeo         *geo   = theProblem->geometry;
    femNodes       *nodes = geo->theNodes;
    femMesh        *mesh  = geo->theElements;
    double *u = theProblem->soluce;
    int size  = S->size;
    
    double *R = theProblem->residuals;

    memset(R, 0, size * sizeof(double));
    
    double a = theProblem->A, b = theProblem->B, c = theProblem->C;
    double rho = theProblem->rho, g = theProblem->g;
    
    int nLocal = mesh->nLocalNode, iElem, iInteg, i, j;
    double x[4], y[4], phi[4], dphidxsi[4], dphideta[4], dphidx[4], dphidy[4];
    int map[4], mapX[4], mapY[4];

    for (iElem = 0; iElem < mesh->nElem; iElem++) {
        for (j = 0; j < nLocal; j++) {
            map[j]  = mesh->elem[iElem * nLocal + j];
            mapX[j] = 2 * map[j];
            mapY[j] = 2 * map[j] + 1;
            x[j]    = nodes->X[map[j]];
            y[j]    = nodes->Y[map[j]];
        }
        for (iInteg = 0; iInteg < rule->n; iInteg++) {
            double xsi = rule->xsi[iInteg], eta = rule->eta[iInteg], w = rule->weight[iInteg];
            femDiscretePhi2(space, xsi, eta, phi);
            femDiscreteDphi2(space, xsi, eta, dphidxsi, dphideta);
            
            double dxdxsi = 0.0, dxdeta = 0.0, dydxsi = 0.0, dydeta = 0.0;
            for (i = 0; i < space->n; i++) {
                dxdxsi += x[i] * dphidxsi[i];
                dxdeta += x[i] * dphideta[i];
                dydxsi += y[i] * dphidxsi[i];
                dydeta += y[i] * dphideta[i];
            }
            double jac = fabs(dxdxsi * dydeta - dxdeta * dydxsi);
            
            for (i = 0; i < space->n; i++) {
                dphidx[i] = (dphidxsi[i] * dydeta - dphideta[i] * dydxsi) / jac;
                dphidy[i] = (dphideta[i] * dxdxsi - dphidxsi[i] * dxdeta) / jac;
            }
            for (i = 0; i < space->n; i++) {
                double sumX = 0.0, sumY = 0.0;
                for (j = 0; j < space->n; j++) {
                    double kxx = (dphidx[i] * a * dphidx[j] + dphidy[i] * c * dphidy[j]) * jac * w;
                    double kxy = (dphidx[i] * b * dphidy[j] + dphidy[i] * c * dphidx[j]) * jac * w;
                    double kyx = (dphidy[i] * b * dphidx[j] + dphidx[i] * c * dphidy[j]) * jac * w;
                    double kyy = (dphidy[i] * a * dphidy[j] + dphidx[i] * c * dphidx[j]) * jac * w;
                    sumX += kxx * u[mapX[j]] + kxy * u[mapY[j]];
                    sumY += kyx * u[mapX[j]] + kyy * u[mapY[j]];
                }
                R[mapX[i]] += sumX;
                R[mapY[i]] += sumY + phi[i] * rho * g * jac * w;
            }
        }
    }
    
    femIntegration *ruleEdge = theProblem->ruleEdge;
    femMesh        *edges    = geo->theEdges;
    int nEdgeLocal = 2;
    double xe[2], ye[2];
    int mapEdge[2], mapU[2];
    
    for (int iBnd = 0; iBnd < theProblem->nBoundaryConditions; iBnd++) {
        femBoundaryCondition *cond = theProblem->conditions[iBnd];
        if (cond->type == NEUMANN_X || cond->type == NEUMANN_Y) {
            femDomain *domain = cond->domain;
            double value = cond->value;
            for (int idx = 0; idx < domain->nElem; idx++) {
                int iEdge = domain->elem[idx];
                mapEdge[0] = edges->elem[2 * iEdge + 0];
                mapEdge[1] = edges->elem[2 * iEdge + 1];
                xe[0] = nodes->X[ mapEdge[0] ];
                xe[1] = nodes->X[ mapEdge[1] ];
                ye[0] = nodes->Y[ mapEdge[0] ];
                ye[1] = nodes->Y[ mapEdge[1] ];
                
                double dx = xe[1] - xe[0], dy = ye[1] - ye[0];
                double length = sqrt(dx * dx + dy * dy);
                double jacEdge = length / 2.0;
                
                if (cond->type == NEUMANN_X) {
                    mapU[0] = 2 * mapEdge[0];
                    mapU[1] = 2 * mapEdge[1];
                } else { 
                    mapU[0] = 2 * mapEdge[0] + 1;
                    mapU[1] = 2 * mapEdge[1] + 1;
                }
                
                for (iInteg = 0; iInteg < ruleEdge->n; iInteg++) {
                    double xsi = ruleEdge->xsi[iInteg], w = ruleEdge->weight[iInteg];
                    double phi0 = 0.5 * (1.0 - xsi);
                    double phi1 = 0.5 * (1.0 + xsi);
                    
                    R[mapU[0]] -= phi0 * value * jacEdge * w;
                    R[mapU[1]] -= phi1 * value * jacEdge * w;
                }
            }
        }
    }
    
    return R;
}