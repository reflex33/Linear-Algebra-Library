# c# Linear Algebra Library
Library for linear algebra.  Includes sub-libraries for matrices and geometry currently.

## TODO
1. Add a class for a line.
2. Add a class for homogeneous transformation matrix?  Benefits: strong typing.  Detriments: lots of extra code.
3. Add a class for vector/vector3d?  Benefits: strong typing.  Detriments: lots of extra code.

## 5/7/2015
Converted function comments to xml code comments

## 4/30/2015
Added equals_threshold to matrix
Added is_normalized_3d_vector to matrix
Modified is_3d_transformation_matrix to check normalized and orthogonal vectors in addition to other checks
Fixed a bug where normalizing a vector of all zeros was returning a value

## 4/24/2015
Added "distance to" point and plane

## 4/23/2015
Added a bracket operator to the matrix class

## 4/22/2015.2
Exceptionified the matrix class

## 4/22/2015
Started linear algebra rewrite.  Library will now be for concepts of linear algebra.  There will be various sub-libraries, matrices will be one of those.  Added a new sub-library for geometry.  This new library contains planes and 3d points to start.

## 4/17/2015
Added a "Changed" event.  This fires whenever any values of the matrix is changed.  Also added a test program.

## 4/7/2015
Nothing to say here... First release
