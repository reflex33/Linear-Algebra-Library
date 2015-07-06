using matrix_library;

namespace geometry_library
{
    /// <summary>
    /// A 3D Point
    /// </summary>
    public class point3d
    {
        /// <summary>
        /// Gets/Sets the 'x' of the point
        /// </summary>
        public double x { get; set; }
        /// <summary>
        /// Gets/Sets the 'y' of the point
        /// </summary>
        public double y { get; set; }
        /// <summary>
        /// Gets/Sets the 'z' of the point
        /// </summary>
        public double z { get; set; }
        /// <summary>
        /// Sets the values of the point using 3x1 matrix form
        /// </summary>
        public matrix matrix_representation
        {
            set
            {
                if (value.rows != 3 || value.cols != 1)
                    throw new System.ArgumentException("3D Point must be a 3x1 matrix/vector!");

                x = value[0, 0];
                y = value[1, 0];
                z = value[2, 0];
            }
        }
        /// <summary>
        /// Gets the values of the point using 3x1 matrix form
        /// </summary>
        public matrix matrix_representation_copy
        {
            get
            {
                matrix m = new matrix(3, 1);
                m[0, 0] = x;
                m[1, 0] = y;
                m[2, 0] = z;

                return m;
            }
        }


        /// <summary>
        /// Determines if the point is front of a plane (i.e. is it in the direction of the plane normal)
        /// </summary>
        /// <param name="p">The plane</param>
        /// <returns>True if it is in front, False if it is on or behind the plane</returns>
        public bool is_in_front_of(plane p)
        {
            if (distance_to(p) > 0)
                return true;
            else
                return false;
        }
        /// <summary>
        /// Determines if the point is near a plane using a threshold
        /// </summary>
        /// <param name="p">The plane</param>
        /// <param name="threshold">The threshold</param>
        /// <returns>True if the point is within the distance indicated with the threshold, False otherwise</returns>
        public bool is_near(plane p, double threshold)
        {
            double result = distance_to(p);
            if (result >= -threshold && result <= threshold)
                return true;
            else
                return false;
        }
        /// <summary>
        /// Determines if the point is behind a plane (i.e. is it in the opposite direction of the plane normal)
        /// </summary>
        /// <param name="p">The plane</param>
        /// <returns>True if it is behind, False if it is on or in front of the plane</returns>
        public bool is_behind(plane p)
        {
            if (distance_to(p) < 0)
                return true;
            else
                return false;
        }
        /// <summary>
        /// Determines the distance of the point to a plane
        /// </summary>
        /// <param name="p">The plane</param>
        /// <returns>The distance</returns>
        public double distance_to(plane p)
        {
            return (p.a * x + p.b * y + p.c * z + p.d) / System.Math.Sqrt(p.a * p.a + p.b * p.b + p.c * p.c);
        }


        /// <summary>
        /// Creates a point where x=y=z=0
        /// </summary>
        public point3d()
        {
            x = 0;
            y = 0;
            z = 0;
        }
        /// <summary>
        /// Creates a point that is a deep copy of another point
        /// </summary>
        /// <param name="point_to_copy">The point to copy</param>
        public point3d(point3d point_to_copy)
        {
            x = point_to_copy.x;
            y = point_to_copy.y;
            z = point_to_copy.z;
        }
    }


    /// <summary>
    /// A collection of points
    /// </summary>
    /// <typeparam name="T">Either point3d or point_and_normal_pair_3d</typeparam>
    public class point_3d_collection<T> : System.Collections.ObjectModel.Collection<T> where T : point3d
    {
        /// <summary>
        /// The centroid of the collection points
        /// </summary>
        public point3d centroid
        {
            get
            {
                if (Count <= 0)
                    throw new System.InvalidOperationException("Can't find the centroid without any points!");

                double min_x = double.MaxValue;
                double min_y = double.MaxValue;
                double min_z = double.MaxValue;
                double max_x = double.MinValue;
                double max_y = double.MinValue;
                double max_z = double.MinValue;

                // Find the object centroid
                foreach (point3d point in this)
                {
                    if (point.x < min_x)
                        min_x = point.x;
                    if (point.y < min_y)
                        min_y = point.y;
                    if (point.z < min_z)
                        min_z = point.z;
                    if (point.x > max_x)
                        max_x = point.x;
                    if (point.y > max_y)
                        max_y = point.y;
                    if (point.z > max_z)
                        max_z = point.z;
                }
                point3d centroid = new point3d();
                centroid.x = (max_x + min_x) / 2.0;
                centroid.y = (max_y + min_y) / 2.0;
                centroid.z = (max_z + min_z) / 2.0;

                return centroid;
            }
        }
        /// <summary>
        /// The bounding box of the collection of points.  Returns a Tuple in the following
        /// format: min_x, max_x, min_y, max_y, min_z, max_z
        /// </summary>
        public System.Tuple<double, double, double, double, double, double> bounding_box
        {
            get
            {
                double min_x = double.MaxValue;
                double min_y = double.MaxValue;
                double min_z = double.MaxValue;
                double max_x = double.MinValue;
                double max_y = double.MinValue;
                double max_z = double.MinValue;

                // Find the object centroid
                foreach (point3d point in this)
                {
                    if (point.x < min_x)
                        min_x = point.x;
                    if (point.y < min_y)
                        min_y = point.y;
                    if (point.z < min_z)
                        min_z = point.z;
                    if (point.x > max_x)
                        max_x = point.x;
                    if (point.y > max_y)
                        max_y = point.y;
                    if (point.z > max_z)
                        max_z = point.z;
                }

                return new System.Tuple<double, double, double, double, double, double>(min_x, max_x, min_y, max_y, min_z, max_z);
            }
        }
    }


    /// <summary>
    /// A 3D point with a matched normal
    /// </summary>
    public class point_and_normal_pair_3d : point3d
    {
        private matrix _normal_data;
        /// <summary>
        /// Gets/Sets the 'x' of the normal
        /// </summary>
        public double normal_x
        {
            get { return _normal_data[0, 0]; }
            set { _normal_data[0, 0] = value; }
        }
        /// <summary>
        /// Gets/Sets the 'y' of the normal
        /// </summary>
        public double normal_y
        {
            get { return _normal_data[1, 0]; }
            set { _normal_data[1, 0] = value; }
        }
        /// <summary>
        /// Gets/Sets the 'z' of the normal
        /// </summary>
        public double normal_z
        {
            get { return _normal_data[2, 0]; }
            set { _normal_data[2, 0] = value; }
        }
        /// <summary>
        /// Sets the values of the normal using 3x1 matrix form
        /// </summary>
        public matrix normal_matrix_representation
        {
            set
            {
                if (value.rows != 3 || value.cols != 1)
                    throw new System.ArgumentException("Normal vector must be a 3x1 matrix/vector!");

                normal_x = value[0, 0];
                normal_y = value[1, 0];
                normal_z = value[2, 0];
            }
        }
        /// <summary>
        /// Gets the values of the normal using 3x1 matrix form
        /// </summary>
        public matrix normal_matrix_representation_copy
        {
            get
            {
                matrix m = new matrix(3, 1);
                m[0, 0] = normal_x;
                m[1, 0] = normal_y;
                m[2, 0] = normal_z;

                return m;
            }
        }


        /// <summary>
        /// Creates a point/normal pair where x=y=z=0 and the normal is 0,0,0
        /// </summary>
        public point_and_normal_pair_3d() : base()
        {
            _normal_data = new matrix(3, 1);
        }
        /// <summary>
        /// Creates a point/normal pair that is a deep copy of another point/normal
        /// </summary>
        /// <param name="point_to_copy">The point to copy</param>
        public point_and_normal_pair_3d(point_and_normal_pair_3d point_to_copy) : this()
        {
            x = point_to_copy.x;
            y = point_to_copy.y;
            z = point_to_copy.z;
            normal_x = point_to_copy.normal_x;
            normal_y = point_to_copy.normal_y;
            normal_z = point_to_copy.normal_z;
        }
        /// <summary>
        /// Creates a point/normal pair that is a combination and deep copy of another point and normal
        /// </summary>
        /// <param name="point_to_copy">The point to copy</param>
        /// <param name="normal_to_copy">The normal to copy</param>
        public point_and_normal_pair_3d(point3d point_to_copy, matrix normal_to_copy) : this()
        {
            if (normal_to_copy.rows != 3 || normal_to_copy.cols != 1)
                throw new System.ArgumentException("Normal vector must be a 3x1 matrix/vector!");

            x = point_to_copy.x;
            y = point_to_copy.y;
            z = point_to_copy.z;
            normal_x = normal_to_copy[0, 0];
            normal_y = normal_to_copy[1, 0];
            normal_z = normal_to_copy[2, 0];
        }
    }
}
