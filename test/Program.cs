using System;
using matrix_library;

namespace test
{
    class Program
    {
        private static void something_changed(object sender, EventArgs e)
        {
            Console.WriteLine("stuff in the matrix changed");
        }

        static void Main(string[] args)
        {
            matrix m = new matrix(5);
            m.Changed += new matrix.MatrixChangedEventHandler(something_changed);

            for (int i = 0; i < m.rows; ++i)
            {
                for (int j = 0; j < m.cols; ++j)
                    Console.Write(m.get(i, j));
                Console.WriteLine();
            }
            Console.WriteLine();

            m.set(0, 2, 8);

            for (int i = 0; i < m.rows; ++i)
            {
                for (int j = 0; j < m.cols; ++j)
                    Console.Write(m.get(i, j));
                Console.WriteLine();
            }
            Console.WriteLine();
        }
    }
}
