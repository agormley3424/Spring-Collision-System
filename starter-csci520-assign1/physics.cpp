/*

  USC/Viterbi/Computer Science
  "Jello Cube" Assignment 1 starter code

*/

#include "jello.h"
#include "physics.h"

struct intPoint
{
    int x;
    int y;
    int z;
};

intPoint operator+(intPoint point1, intPoint point2)
{
    point1.x += point2.x;
    point1.y += point2.y;
    point1.z += point2.z;

    return point1;
}

point operator/(point myPoint, double divisor)
{
    myPoint.x /= divisor;
    myPoint.y /= divisor;
    myPoint.z /= divisor;

    return myPoint;
}

point operator+(point point1, point point2)
{
    point1.x += point2.x;
    point1.y += point2.y;
    point1.z += point2.z;

    return point1;
}

point operator*(double floatVal, point& myPoint)
{
    point returnVal;

    returnVal.x = floatVal * myPoint.x;
    returnVal.y = floatVal * myPoint.y;
    returnVal.z = floatVal * myPoint.z;

    return returnVal;
}

point operator*(point& myPoint, double floatVal)
{
    point returnVal;

    returnVal.x = floatVal * myPoint.x;
    returnVal.y = floatVal * myPoint.y;
    returnVal.z = floatVal * myPoint.z;

    return returnVal;
}

// Given myPoint's location, interpolate forcefield values
point TrilinearInterp(point& myPoint, point* forceField, int gridResolution, int gridResolutionInvert)
{
    // Lower left corner of the grid cell the point is in
    int lowerLeftX = floor(myPoint.x * gridResolutionInvert);
    int lowerLeftY = floor(myPoint.y * gridResolutionInvert);
    int lowerLeftZ = floor(myPoint.z * gridResolutionInvert);

    intPoint lowerLeftPoint{ lowerLeftX, lowerLeftY, lowerLeftZ };

    double alpha = (myPoint.x - lowerLeftX) * gridResolutionInvert;
    double beta = (myPoint.y - lowerLeftY) * gridResolutionInvert;
    double gamma = (myPoint.z - lowerLeftZ) * gridResolutionInvert;

    point finalForce = { 0.0, 0.0, 0.0 };

    for (int x = 0; x < 2; ++x)
    {
        for (int y = 0; y < 2; ++y)
        {
            for (int z = 0; z < 2; ++z)
            {
                intPoint thisPoint = lowerLeftPoint + intPoint{x, y, z};

                point& pointForce = forceField[thisPoint.z * gridResolution * gridResolution +
                                               thisPoint.y * gridResolution +
                                               thisPoint.x];

                double alphaVal = (x == 0) ? (1.0 - alpha) : alpha;
                double betaVal = (y == 0) ? (1.0 - beta) : beta;
                double gammaVal = (z == 0) ? (1.0 - gamma) : gamma;

                finalForce = finalForce + (alphaVal * betaVal * gammaVal * pointForce);
            }
        }
    }

    return finalForce;
}

/* Computes acceleration to every control point of the jello cube, 
   which is in state given by 'jello'.
   Returns result in array 'a'. */
void computeAcceleration(struct world * jello, struct point a[8][8][8])
{
    for (int x = 0; x < 8; ++x)
    {
        for (int y = 0; y < 8; ++y)
        {
            for (int z = 0; z < 8; ++z)
            {
                point totalForce = TrilinearInterp(jello->p[x][y][z], jello->forceField,
                    jello->resolution, 1 / jello->resolution);

                a[x][y][z] = totalForce / jello->mass;
            }
        }
    }
}

/* performs one step of Euler Integration */
/* as a result, updates the jello structure */
void Euler(struct world * jello)
{
  int i,j,k;
  point a[8][8][8];

  computeAcceleration(jello, a);
  
  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
        jello->p[i][j][k].x += jello->dt * jello->v[i][j][k].x;
        jello->p[i][j][k].y += jello->dt * jello->v[i][j][k].y;
        jello->p[i][j][k].z += jello->dt * jello->v[i][j][k].z;
        jello->v[i][j][k].x += jello->dt * a[i][j][k].x;
        jello->v[i][j][k].y += jello->dt * a[i][j][k].y;
        jello->v[i][j][k].z += jello->dt * a[i][j][k].z;

      }
}

/* performs one step of RK4 Integration */
/* as a result, updates the jello structure */
void RK4(struct world * jello)
{
  point F1p[8][8][8], F1v[8][8][8], 
        F2p[8][8][8], F2v[8][8][8],
        F3p[8][8][8], F3v[8][8][8],
        F4p[8][8][8], F4v[8][8][8];

  point a[8][8][8];


  struct world buffer;

  int i,j,k;

  buffer = *jello; // make a copy of jello

  computeAcceleration(jello, a);

  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
         pMULTIPLY(jello->v[i][j][k],jello->dt,F1p[i][j][k]);
         pMULTIPLY(a[i][j][k],jello->dt,F1v[i][j][k]);
         pMULTIPLY(F1p[i][j][k],0.5,buffer.p[i][j][k]);
         pMULTIPLY(F1v[i][j][k],0.5,buffer.v[i][j][k]);
         pSUM(jello->p[i][j][k],buffer.p[i][j][k],buffer.p[i][j][k]);
         pSUM(jello->v[i][j][k],buffer.v[i][j][k],buffer.v[i][j][k]);
      }

  computeAcceleration(&buffer, a);

  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
         // F2p = dt * buffer.v;
         pMULTIPLY(buffer.v[i][j][k],jello->dt,F2p[i][j][k]);
         // F2v = dt * a(buffer.p,buffer.v);     
         pMULTIPLY(a[i][j][k],jello->dt,F2v[i][j][k]);
         pMULTIPLY(F2p[i][j][k],0.5,buffer.p[i][j][k]);
         pMULTIPLY(F2v[i][j][k],0.5,buffer.v[i][j][k]);
         pSUM(jello->p[i][j][k],buffer.p[i][j][k],buffer.p[i][j][k]);
         pSUM(jello->v[i][j][k],buffer.v[i][j][k],buffer.v[i][j][k]);
      }

  computeAcceleration(&buffer, a);

  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
         // F3p = dt * buffer.v;
         pMULTIPLY(buffer.v[i][j][k],jello->dt,F3p[i][j][k]);
         // F3v = dt * a(buffer.p,buffer.v);     
         pMULTIPLY(a[i][j][k],jello->dt,F3v[i][j][k]);
         pMULTIPLY(F3p[i][j][k],1.0,buffer.p[i][j][k]);
         pMULTIPLY(F3v[i][j][k],1.0,buffer.v[i][j][k]);
         pSUM(jello->p[i][j][k],buffer.p[i][j][k],buffer.p[i][j][k]);
         pSUM(jello->v[i][j][k],buffer.v[i][j][k],buffer.v[i][j][k]);
      }
         
  computeAcceleration(&buffer, a);


  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
         // F3p = dt * buffer.v;
         pMULTIPLY(buffer.v[i][j][k],jello->dt,F4p[i][j][k]);
         // F3v = dt * a(buffer.p,buffer.v);     
         pMULTIPLY(a[i][j][k],jello->dt,F4v[i][j][k]);

         pMULTIPLY(F2p[i][j][k],2,buffer.p[i][j][k]);
         pMULTIPLY(F3p[i][j][k],2,buffer.v[i][j][k]);
         pSUM(buffer.p[i][j][k],buffer.v[i][j][k],buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],F1p[i][j][k],buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],F4p[i][j][k],buffer.p[i][j][k]);
         pMULTIPLY(buffer.p[i][j][k],1.0 / 6,buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],jello->p[i][j][k],jello->p[i][j][k]);

         pMULTIPLY(F2v[i][j][k],2,buffer.p[i][j][k]);
         pMULTIPLY(F3v[i][j][k],2,buffer.v[i][j][k]);
         pSUM(buffer.p[i][j][k],buffer.v[i][j][k],buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],F1v[i][j][k],buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],F4v[i][j][k],buffer.p[i][j][k]);
         pMULTIPLY(buffer.p[i][j][k],1.0 / 6,buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],jello->v[i][j][k],jello->v[i][j][k]);
      }

  return;  
}
