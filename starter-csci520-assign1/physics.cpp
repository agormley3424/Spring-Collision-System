/*

  USC/Viterbi/Computer Science
  "Jello Cube" Assignment 1 starter code

*/

#include "jello.h"
#include "physics.h"
#include <stdexcept>

double gridLength = 1.0 / 7.0;

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

point operator-(point point1, point point2)
{
    point1.x -= point2.x;
    point1.y -= point2.y;
    point1.z -= point2.z;

    return point1;
}

point operator+(double myFloat, point point)
{
    point.x += myFloat;
    point.y += myFloat;
    point.z += myFloat;

    return point;
}

point operator+(point point, double myFloat)
{
    point.x += myFloat;
    point.y += myFloat;
    point.z += myFloat;

    return point;
}

point operator*(double floatVal, point myPoint)
{
    point returnVal;

    returnVal.x = floatVal * myPoint.x;
    returnVal.y = floatVal * myPoint.y;
    returnVal.z = floatVal * myPoint.z;

    return returnVal;
}

point operator*(point myPoint, double floatVal)
{
    point returnVal;

    returnVal.x = floatVal * myPoint.x;
    returnVal.y = floatVal * myPoint.y;
    returnVal.z = floatVal * myPoint.z;

    return returnVal;
}

// Given myPoint's location, interpolate forcefield values
// myPoint is the point to be interpolated
// forceField is the forcefield vector (3D vector compressed to 1D)
// gridResolution is the size of the grid on one side
// gridLengthInvert is the inverted length between two points on the grid
point TrilinearInterp(point& myPoint, point* forceField, int gridResolution, double gridLengthInvert)
{
    // Speed is good, shouldn't be rotating as much

    point normalizedPoint = myPoint + 2.0;

    // If point is out of bounds
    if (normalizedPoint.x > 4.0 || normalizedPoint.x < 0.0 ||
        normalizedPoint.y > 4.0 || normalizedPoint.y < 0.0 ||
        normalizedPoint.z > 4.0 || normalizedPoint.z < 0.0)
    {
        // printf("Point has escaped bounding box!\n");

        // throw std::runtime_error("Physics Error: TrilinearInterp: Point is out of bounds\n");

        return point{ 0.0, 0.0, 0.0 };
    }

    // Gridlength = 4 / (n - 1)
    // gridLengthInvert = (n - 1) / 4

    // coordPos is the position of our point relative to the force field grid
    // lowerLeft is the position of the lower-left corner of the force field cell our point is in
    double coordPosX = normalizedPoint.x * gridLengthInvert;
    int lowerLeftX = floor(coordPosX);
    double coordPosY = normalizedPoint.y * gridLengthInvert;
    int lowerLeftY = floor(coordPosY);
    double coordPosZ = normalizedPoint.z * gridLengthInvert;
    int lowerLeftZ = floor(coordPosZ);

    intPoint lowerLeftPoint{ lowerLeftX, lowerLeftY, lowerLeftZ };

    // Relative position of our point within the cell
    // No need to normalize, since (coordPos - lowerLeft) is in [0.0, 1.0]
    double alpha = (coordPosX - lowerLeftX);
    double beta = (coordPosY - lowerLeftY);
    double gamma = (coordPosZ - lowerLeftZ);

    point finalForce = { 0.0, 0.0, 0.0 };

    // Find the force at each point in the cell and interpolate to our point
    for (int x = 0; x < 2; ++x)
    {
        for (int y = 0; y < 2; ++y)
        {
            for (int z = 0; z < 2; ++z)
            {
                // A point in the cell
                intPoint thisPoint = lowerLeftPoint + intPoint{x, y, z};


                point pointForce = forceField[thisPoint.x * gridResolution * gridResolution +
                                               thisPoint.y * gridResolution +
                                               thisPoint.z];

                double alphaVal = (x == 0) ? (1.0 - alpha) : alpha;
                double betaVal = (y == 0) ? (1.0 - beta) : beta;
                double gammaVal = (z == 0) ? (1.0 - gamma) : gamma;

                finalForce = finalForce + ((alphaVal * betaVal * gammaVal) * pointForce);
            }
        }
    }

    return finalForce;
}

double PointMagnitude(point& p)
{
    return sqrt(pow(p.x, 2) + pow(p.y, 2) + pow(p.z, 2));
}

double DotProduct(point p1, point& p2)
{
    return (p1.x * p2.x) + (p1.y * p2.y) + (p1.z * p2.z);
}

point PenaltyForce(point& p, point& velocity, double kCoef, double dCoef)
{
    point totalForce = { 0.0, 0.0, 0.0 };
    point direction;
    double penetrationDist;

    point backForceX = { 0.0, 0.0, 0.0 };
    point backDampX = { 0.0, 0.0, 0.0 };
    point backForceY = { 0.0, 0.0, 0.0 };
    point backDampY = { 0.0, 0.0, 0.0 };
    point backForceZ = { 0.0, 0.0, 0.0 };
    point backDampZ = { 0.0, 0.0, 0.0 };
    
    if (p.x > 2.0)
    {
        direction = { -1.0, 0.0, 0.0 };
        penetrationDist = abs(p.x - 2.0);

        backForceX = kCoef * penetrationDist * direction;
        backDampX = (-dCoef) * DotProduct(velocity, direction) * (direction);
    }
    else if (p.x < -2.0)
    {
        direction = { 1.0, 0.0, 0.0 };
        penetrationDist = abs(p.x + 2.0);

        backForceX = kCoef * penetrationDist * direction;
        backDampX = (-dCoef) * DotProduct(velocity, direction) * (direction);
    }

    if (p.y > 2.0)
    {
        direction = { 0.0, -1.0, 0.0 };
        penetrationDist = abs(p.y - 2.0);

        backForceY = kCoef * penetrationDist * direction;
        backDampY = (-dCoef) * DotProduct(velocity, direction) * (direction);
    }
    else if (p.y < -2.0)
    {
        direction = { 0.0, 1.0, 0.0 };
        penetrationDist = abs(p.y + 2.0);

        backForceY = kCoef * penetrationDist * direction;
        backDampY = (-dCoef) * DotProduct(velocity, direction) * (direction);
    }

    if (p.z > 2.0)
    {
        direction = { 0.0, 0.0, -1.0 };
        penetrationDist = abs(p.z - 2.0);

        backForceZ = kCoef * penetrationDist * direction;
        backDampZ = (-dCoef) * DotProduct(velocity, direction) * (direction);
    }
    else if (p.z < -2.0)
    {
        direction = { 0.0, 0.0, 1.0 };
        penetrationDist = abs(p.z + 2.0);

        backForceZ = kCoef * penetrationDist * direction;
        backDampZ = (-dCoef) * DotProduct(velocity, direction) * (direction);
    }

    return backForceX + backDampX +
           backForceY + backDampY +
           backForceZ + backDampZ;
}

point CalcSpringForce(world* jello, double restLength, point pos1, point vel1, point pos2, point vel2)
{
    // Elastic force
    point neighborToMe = pos1 - pos2;
    double distMagnitude = PointMagnitude(neighborToMe);
    double invertDistMagnitude = 1.0 / distMagnitude;

    point hookeForce = -(jello->kElastic) * (distMagnitude - restLength)
        * (neighborToMe * invertDistMagnitude);


    // Damping force
    point dampForce = -(jello->dElastic) * (DotProduct((vel1 - vel2), neighborToMe) * invertDistMagnitude)
        * (invertDistMagnitude * neighborToMe);

    return hookeForce + dampForce;
}

// i, j, and k are the coordinates of the point in the point array
point SpringForce(point& myPoint, world* jello, int i, int j, int k)
{
    // Size of the cube is 1 x 1 x 1 undeformed
    // Grid length is 1 / 8 = 0.125

    point neighbor;

    point totalForce = point{ 0.0, 0.0, 0.0 };


    // Structural springs

    double restLength = gridLength;
    point myVelocity = jello->v[i][j][k];

    // For each dimension
    for (int d = 0; d < 3; ++d)
    {
        // Check neighbor in both directions
        for (int n = -1; n < 2; n += 2)
        {
            int newI = d == 0 ? i + n : i;
            int newJ = d == 1 ? j + n : j;
            int newK = d == 2 ? k + n : k;

            // Only use this neighbor if it has a valid index
            if (newI >= 0 && newI < 8 &&
                newJ >= 0 && newJ < 8 &&
                newK >= 0 && newK < 8)
            {
                neighbor = jello->p[newI][newJ][newK];

                totalForce = totalForce + CalcSpringForce(jello, restLength, myPoint, jello->v[i][j][k],
                    neighbor, jello->v[newI][newJ][newK]);
            }
        }
    }

    // Shear springs (connected planes)

    restLength = sqrt(2) * gridLength;

    // Dimension orthogonal to this point
    for (int d1 = 0; d1 < 3; ++d1)
    {
        // Dimension orthogonal to that dimension
        for (int d2 = 0; d2 < 3; ++d2)
        {
            if (d2 == d1) continue;

            // Which side in the first dimension?
            for (int n1 = -1; n1 < 2; n1 += 2)
            {
                // Which side in the second dimension?
                for (int n2 = -1; n2 < 2; n2 += 2)
                {
                    int newI1 = d1 == 0 ? i + n1 : i;
                    int newJ1 = d1 == 1 ? j + n1 : j;
                    int newK1 = d1 == 2 ? k + n1 : k;

                    int newI2 = d2 == 0 ? newI1 + n2 : newI1;
                    int newJ2 = d2 == 1 ? newJ1 + n2 : newJ1;
                    int newK2 = d2 == 2 ? newK1 + n2 : newK1;

                    // Only use this neighbor if it has a valid index
                    if (newI2 >= 0 && newI2 < 8 &&
                        newJ2 >= 0 && newJ2 < 8 &&
                        newK2 >= 0 && newK2 < 8)
                    {
                        neighbor = jello->p[newI2][newJ2][newK2];

                        totalForce = totalForce + CalcSpringForce(jello, restLength, myPoint, jello->v[i][j][k],
                            neighbor, jello->v[newI2][newJ2][newK2]);
                    }
                }
            }
        }
    }

    // Shear springs (main diagonal neighbors)

    restLength = sqrt(3) * gridLength;

    for (int x = -1; x < 2; x += 2)
    {
        for (int y = -1; y < 2; y += 2)
        {
            for (int z = -1; z < 2; z += 2)
            {
                int newI = i + x;
                int newJ = j + y;
                int newK = k + z;

                // Only use this neighbor if it has a valid index
                if (newI >= 0 && newI < 8 &&
                    newJ >= 0 && newJ < 8 &&
                    newK >= 0 && newK < 8)
                {
                    neighbor = jello->p[newI][newJ][newK];

                    totalForce = totalForce + CalcSpringForce(jello, restLength, myPoint, jello->v[i][j][k],
                        neighbor, jello->v[newI][newJ][newK]);
                }
            }
        }
    }

    // Bend Springs

    restLength = gridLength * 2;

    // For each dimension
    for (int d = 0; d < 3; ++d)
    {
        // Check neighbor in both directions
        for (int n = -2; n < 4; n += 4)
        {
            int newI = d == 0 ? i + n : i;
            int newJ = d == 1 ? j + n : j;
            int newK = d == 2 ? k + n : k;

            // Only use this neighbor if it has a valid index
            if (newI >= 0 && newI < 8 &&
                newJ >= 0 && newJ < 8 &&
                newK >= 0 && newK < 8)
            {
                neighbor = jello->p[newI][newJ][newK];

                totalForce = totalForce + CalcSpringForce(jello, restLength, myPoint, jello->v[i][j][k],
                    neighbor, jello->v[newI][newJ][newK]);
            }
        }
    }

    return totalForce;
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
                int resolution = jello->resolution;

                point& myPoint = jello->p[x][y][z];
                
                point forceFieldForce = TrilinearInterp(myPoint, jello->forceField,
                    resolution, (resolution - 1) / 4.0);

                point penaltyForce = PenaltyForce(myPoint, jello->v[x][y][z], jello->kCollision, jello->dCollision);

                point springForce = SpringForce(myPoint, jello, x, y, z);

                a[x][y][z] = (forceFieldForce + penaltyForce + springForce) / jello->mass;
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
