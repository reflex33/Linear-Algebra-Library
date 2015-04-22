using matrix_library;

namespace geometry_library
{
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

        public double a { get; set; }
        public double b { get; set; }
        public double c { get; set; }
        public double d { get; set; }
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

        public plane transform(matrix m)
        {
            if (!m.is_3d_transformation_matrix)
                throw new System.ArgumentException("Transformation matrix must be 4x4!");

            plane result = new plane();
            result.matrix_representation = m.transpose().inverse() * matrix_representation;

            return result;
        }

        public plane()
        {
            a = 0;
            b = 0;
            c = 0;
            d = 0;
        }
    }
}
