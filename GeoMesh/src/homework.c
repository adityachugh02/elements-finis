#include "fem.h"

double hermiteInterpolation(double d, double d_star, double h0, double h_star) {
    if (d >= d_star) return h_star;
    double t = d / d_star;
    return h0 + (h_star - h0) * (3 * t * t - 2 * t * t * t);
}


double geoSize(double x, double y){

    femGeo* theGeometry = geoGetGeometry();
    
    double h = 0.1;

    // hook inside
    double x0 = 0.55;
    double y0 = 0.25;
    double r0 = 0.17;
    double h0 = 0.01;
    double d0 = 0.5;
  
    // hole
    double x1 = 0.48;
    double y1 = 0.67;
    double r1 = 0.1;
    double h1 = 0.01;
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


#define ___ 0

void geoMeshGenerate() {

    femGeo* theGeometry = geoGetGeometry();

    double w = theGeometry->LxPlate;
    double h = theGeometry->LyPlate;
     
    double x0 = theGeometry->xNotch;
    double y0 = theGeometry->yNotch;
    double r0 = theGeometry->rNotch;
    
    
    double x1 = theGeometry->xHole;
    double y1 = theGeometry->yHole;
    double r1 = theGeometry->rHole;

    printf("x0: %f \n", x0);
    printf("y0: %f \n", y0);
    printf("w: %f \n", w);
    printf("h: %f \n", h);

    printf("x1: %f \n", x1);
    printf("y1: %f \n", y1);
 
//
//  -1- Construction de la géomotrie avec OpenCascade
//      On crée le rectangle
//      On crée les deux cercles
//      On soustrait les cercles du rectangle :-)
//
 
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
//
//  -2- Définition de la fonction callback pour la taille de référence
//      Synchronisation de OpenCascade avec gmsh
//      Génération du maillage (avec l'option Mesh.SaveAll :-)
                  
   
    geoSetSizeCallback(geoSize);
                                  
    gmshModelOccSynchronize(&ierr);       
    gmshOptionSetNumber("Mesh.SaveAll", 1, &ierr);
    gmshModelMeshGenerate(2, &ierr);  
       
//
//  Generation de quads :-)
//
 /**    gmshOptionSetNumber("Mesh.SaveAll", 1, &ierr);
    gmshOptionSetNumber("Mesh.RecombineAll", 1, &ierr);
    gmshOptionSetNumber("Mesh.Algorithm", 8, &ierr);  
    gmshOptionSetNumber("Mesh.RecombinationAlgorithm", 1.0, &ierr); 
    gmshModelGeoMeshSetRecombine(2,1,45,&ierr);  
    gmshModelMeshGenerate(2, &ierr);  **/
   
 
//
//  Plot of Fltk
//
   //gmshFltkInitialize(&ierr);
   //gmshFltkRun(&ierr);  
//
    
}