using matrix_library;

namespace geometry_library
{
    /// <summary>
    /// A 3D Point
    /// </summary>
    public class point3d
    {
        private void matrix_changed(object sender, System.EventArgs e)
        {
            matrix m = sender as matrix;
            if (m == null)  // If the sender isn't a matrix, throw an exception.  Note: this should never happen
                throw new System.ArgumentException("Unexpected sender to event matrix_changed!");

            if (m.rows != 3 || m.cols != 1)
                throw new System.ArgumentException("3D Point must be a 3x1 matrix/vector!");

            x = m[0, 0];
            y = m[1, 0];
            z = m[2, 0];
        }


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
        /// Gets/Sets the values of the point using 3x1 matrix form
        /// Note:  You can change the values by giving a complete 3x1 vector, or you may change individual values of the resulting matrix
        /// form.  If you change values in the matrix form, it WILL be reflected in the point.
        /// </summary>
        public matrix matrix_representation
        {
            get
            {
                matrix m = new matrix(3, 1);
                m[0, 0] = x;
                m[1, 0] = y;
                m[2, 0] = z;
                m.Changed += new matrix.MatrixChangedEventHandler(matrix_changed);

                return m;
            }
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
    /// A 3D point with a matched normal
    /// </summary>
    public class point_and_normal_pair_3d : point3d
    {
        private void normal_matrix_changed(object sender, System.EventArgs e)
        {
            matrix m = sender as matrix;
            if (m == null)  // If the sender isn't a matrix, throw an exception.  Note: this should never happen
                throw new System.ArgumentException("Unexpected sender to event normal_matrix_changed!");

            if (m.rows != 3 || m.cols != 1)
                throw new System.ArgumentException("Normal vector must be a 3x1 matrix/vector!");

            normal_x = m[0, 0];
            normal_y = m[1, 0];
            normal_z = m[2, 0];
        } 


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
        /// Gets/Sets the values of the normal using 3x1 matrix form
        /// Note:  You can change the values by giving a complete 3x1 vector, or you may change individual values of the resulting matrix
        /// form.  If you change values in the matrix form, it WILL be reflected in the normal.
        /// </summary>
        public matrix normal_matrix_representation
        {
            get
            {
                matrix m = new matrix(3, 1);
                m[0, 0] = normal_x;
                m[1, 0] = normal_y;
                m[2, 0] = normal_z;
                m.Changed += new matrix.MatrixChangedEventHandler(normal_matrix_changed);

                return m;
            }
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
        /// Creates a point where x=y=z=0 and the normal is 0,0,0
        /// </summary>
        public point_and_normal_pair_3d() : base()
        {
            _normal_data = new matrix(3, 1);
        }
        /// <summary>
        /// Creates a point that is a deep copy of another point
        /// </summary>
        /// <param name="point_to_copy">The point to copy</param>
        public point_and_normal_pair_3d(point_and_normal_pair_3d point_to_copy)
        {
            x = point_to_copy.x;
            y = point_to_copy.y;
            z = point_to_copy.z;
            normal_x = point_to_copy.normal_x;
            normal_y = point_to_copy.normal_y;
            normal_z = point_to_copy.normal_z;
        }
    }
}
