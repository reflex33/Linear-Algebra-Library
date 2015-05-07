using matrix_library;

namespace geometry_library
{
    /// <summary>
    /// A 3D Plane
    /// </summary>
    public class plane
    {
        private void matrix_changed(object sender, System.EventArgs e)
        {
            matrix m = sender as matrix;
            if (m == null)  // If the sender isn't a matrix, throw an exception.  Note: this should never happen
                throw new System.ArgumentException("Unexpected sender to event matrix_changed!");

            if (m.rows != 4 || m.cols != 1)
                throw new System.ArgumentException("Plane equation must be a 4x1 matrix/vector!");

            a = m.get(0, 0);
            b = m.get(1, 0);
            c = m.get(2, 0);
            d = m.get(3, 0);
        }

        
        /// <summary>
        /// Gets/Sets the 'a' of the plane
        /// </summary>
        public double a { get; set; }
        /// <summary>
        /// Gets/Sets the 'b' of the plane
        /// </summary>
        public double b { get; set; }
        /// <summary>
        /// Gets/Sets the 'c' of the plane
        /// </summary>
        public double c { get; set; }
        /// <summary>
        /// Gets/Sets the 'd' of the plane
        /// </summary>
        public double d { get; set; }
        /// <summary>
        /// Gets/Sets the values of the plane using 4x1 matrix form
        /// Note:  You can change the values by giving a complete 4x1 vector, or you may change individual values of the resulting matrix
        /// form.  If you change values in the matrix form, it WILL be reflected in the plane.
        /// </summary>
        public matrix matrix_representation
        {
            get
            {
                matrix m = new matrix(4, 1);
                m.set(0, 0, a);
                m.set(1, 0, b);
                m.set(2, 0, c);
                m.set(3, 0, d);
                m.Changed += new matrix.MatrixChangedEventHandler(matrix_changed);
                
                return m;
            }
            set
            {
                if (value.rows != 4 || value.cols != 1)
                    throw new System.ArgumentException("Plane equation must be a 4x1 matrix/vector!");

                a = value.get(0, 0);
                b = value.get(1, 0);
                c = value.get(2, 0);
                d = value.get(3, 0);
            }
        }


        /// <summary>
        /// Transforms the plane by a homogeneous transformation matrix (m*p)
        /// </summary>
        /// <param name="m">The homogeneous transformation matrix</param>
        /// <returns>New plane</returns>
        public plane transform(matrix m)
        {
            if (!m.is_3d_transformation_matrix)
                throw new System.ArgumentException("Transformation matrix must be 4x4!");

            plane result = new plane();
            result.matrix_representation = m.transpose().inverse() * matrix_representation;

            return result;
        }
        /// <summary>
        /// Determines the distance of the plane to a point
        /// </summary>
        /// <param name="p">The point</param>
        /// <returns>The distance</returns>
        public double distance_to(point3d p)
        {
            return p.distance_to(this);
        }


        /// <summary>
        /// Creates a plane where a=b=c=d=0
        /// </summary>
        public plane()
        {
            a = 0;
            b = 0;
            c = 0;
            d = 0;
        }
    }
}
