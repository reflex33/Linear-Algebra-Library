# c# Linear Algebra Library
Library for linear algebra.  Includes sub-libraries for matrices and geometry currently.

## TODO
1. Add a class for a line.
2. Add a class for homogeneous transformation matrix?  Benefits: strong typing.  Detriments: lots of extra code.
3. Add a class for vector/vector3d?  Benefits: strong typing.  Detriments: lots of extra code.

## Changelog
### 2015.06.16
Added a collection of points class, including centroid and bounding box calculations

### 2015.06.12
Bug fix for point_normal_pair_3d constructors

### 2015.06.02
Added another constructor to point_and_normal_pair_3d class
Added point_and_normal_pair_3d class

### 2015.05.25
Added a rotation submatrix property

### 2015.05.20
Added a copy constructor to the point3d and plane classes

### 2015.05.07
Converted function comments to xml code comments

### 2015.04.30
Added equals_threshold to matrix
Added is_normalized_3d_vector to matrix
Modified is_3d_transformation_matrix to check normalized and orthogonal vectors in addition to other checks
Fixed a bug where normalizing a vector of all zeros was returning a value

### 2015.04.24
Added "distance to" point and plane

### 2015.04.23
Added a bracket operator to the matrix class

### 2015.04.22
Exceptionified the matrix class
Started linear algebra rewrite.  Library will now be for concepts of linear algebra.  There will be various sub-libraries, matrices will be one of those.  Added a new sub-library for geometry.  This new library contains planes and 3d points to start.

### 2015.04.17
Added a "Changed" event.  This fires whenever any values of the matrix is changed.  Also added a test program.

### 2015.04.07
Nothing to say here... First release
