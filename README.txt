<Please submit this file with your solution.>

CSCI 520, Assignment 1

Aaron Gormley

================

In computeAcceleration, acceleration for each point is calculated as a = F / m, where m is the mass of a point, and F is the total force
acting on it.

Total force is the sum of results from three methods: SpringForce, PenaltyForce, and TrilinearInterp.

TrilinearInterp calculates the force field force acting on the point.
  It does this by first normalizing the point relative to the bounding box boundaries
    (center point is then 2 instead of 0, end is 4 instead of 2, etc.)
  The force field index of the lower left point is computed by dividing the normalized point by the
    length of a grid segment between two points (dependent on the resolution of the force field)
  Alpha, beta, and gamma, the length of the normalized point relative to each axis of the cell, are calculated
  For each of the eight points in the cell, their interpolated force value is calculated with the params from above
  The interpolated force values are added together for the overall interpolated force acting on the point

PenaltyForce determines whether a collision has occurred, and forces the object back if so
  To identify whether a collision has occurred, the point's value in each axis is checked.
    If the value is outside the dimensions of the bounding box, a collision is occurring along that dimension
  This is checked for all dimensions independently, to account for corners
  When a collision is occurring along some dimension, a spring is modeled to drag the point back
  The elastic force of this spring, or springs if a corner was penetrated, pulls the point
    back along the direction it left from, while the damping force reduces its velocity
  The total force returned is the sum of all these forces

SpringForce models structural, shear, and bend springs connecting points to their neighbors
  Structural springs find every neighbor along all viable ends of each dimension relative to the original point
    and connect them
  Shear springs along the same plane find all diagonal neighbors by choosing two different dimensions to move along,
  before connecting them
  Shear springs along the cube diagonals of the up to four connected cells are found by making a step forward or back
  in each dimension relative to the original point
  Bend springs work the same way as structural, except they move two steps along a dimension instead of one
  Invalid neighbors outside of the cube aren't accounted for
  The total force returned is the sum of elastic and damping forces from all these springs

Point positions are calculated by RK4 or Euler, depending on what the world specifies.
  Update calculations are called in showCube, just before the cube render is updated

The only other notable thing I did was add my own vector class for integers (intPoint),
and write a bunch of overload operations for intPoints and points, to make the code
easier to read and to write, vs. calling a bunch of operations every time I wanted to add points
together or multiply them by a float or something
