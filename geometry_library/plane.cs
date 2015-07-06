using matrix_library;

namespace geometry_library
{
    /// <summary>
    /// A 3D Plane
    /// </summary>
    public class plane
    {
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
        /// Sets the values of the plane using 4x1 matrix form
        /// </summary>
        public matrix matrix_representation
        {
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
        /// Gets the values of the plane using 4x1 matrix form
        /// </summary>
        public matrix matrix_representation_copy
        {
            get
            {
                matrix m = new matrix(4, 1);
                m.set(0, 0, a);
                m.set(1, 0, b);
                m.set(2, 0, c);
                m.set(3, 0, d);

                return m;
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
            result.matrix_representation = m.transpose().inverse() * matrix_representation_copy;

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
        /// <summary>
        /// Creates a plane that is a deep copy of another plane
        /// </summary>
        /// <param name="plane_to_copy">The plane to copy</param>
        public plane(plane plane_to_copy)
        {
            a = plane_to_copy.a;
            b = plane_to_copy.b;
            c = plane_to_copy.c;
            d = plane_to_copy.d;
        }
    }
}
