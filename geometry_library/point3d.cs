using matrix_library;

namespace geometry_library
{
    public class point3d
    {
        private void matrix_changed(object sender, System.EventArgs e)
        {
            matrix m = sender as matrix;
            if (m == null)  // If the sender isn't a matrix, throw an exception.  Note: this should never happen
                throw new System.ArgumentException("Unexpected sender to event matrix_changed!");

            if (m.rows != 3 || m.cols != 1)
                throw new System.ArgumentException("3D Point must be a 3x1 matrix/vector!");

            x = m.get(0, 0);
            y = m.get(1, 0);
            z = m.get(2, 0);
        }

        public double x { get; set; }
        public double y { get; set; }
        public double z { get; set; }
        public matrix matrix_representation
        {
            get
            {
                matrix m = new matrix(3, 1);
                m.set(0, 0, x);
                m.set(1, 0, y);
                m.set(2, 0, z);
                m.Changed += new matrix.MatrixChangedEventHandler(matrix_changed);

                return m;
            }
            set
            {
                if (value.rows != 3 || value.cols != 1)
                    throw new System.ArgumentException("3D Point must be a 3x1 matrix/vector!");

                x = value.get(0, 0);
                y = value.get(1, 0);
                z = value.get(2, 0);
            }
        }

        public bool is_in_front_of(plane p)
        {
            if (p.a * x + p.b * y + p.c * z + p.d > 0)
                return true;
            else
                return false;
        }
        public bool is_near(plane p, double threshold)
        {
            double result = p.a * x + p.b * y + p.c * z + p.d;
            if (result >= -threshold && result <= threshold)
                return true;
            else
                return false;
        }
        public bool is_behind(plane p)
        {
            if (p.a * x + p.b * y + p.c * z + p.d < 0)
                return true;
            else
                return false;
        }

        public point3d()
        {
            x = 0;
            y = 0;
            z = 0;
        }
    }
}
