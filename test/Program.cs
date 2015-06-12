using System;
using matrix_library;
using geometry_library;

namespace test
{
    class Program
    {
        static void Main(string[] args)
        {
            //matrix trans = new matrix(4);
            //trans[0, 0] = 1; trans[0, 1] = 0; trans[0, 2] = 1; trans[0, 3] = 0;
            //trans[1, 0] = 0; trans[1, 1] = 1; trans[1, 2] = 0; trans[1, 3] = 10;
            //trans[2, 0] = 0; trans[2, 1] = 0; trans[2, 2] = 1; trans[2, 3] = 0;
            //trans[3, 0] = 0; trans[3, 1] = 0; trans[3, 2] = 0; trans[3, 3] = 1;
            matrix trans = new matrix(3, 1);
            trans[0, 0] = 1.05;
            trans[1, 0] = 0.05;
            trans[2, 0] = 0;

            Console.WriteLine(trans.is_3d_vector);
            Console.WriteLine(trans.is_normalized_3d_vector);


            point_and_normal_pair_3d p1 = new point_and_normal_pair_3d();
            p1.x = 5;
            p1.y = 6;
            p1.z = 9;
            p1.normal_x = -1.0;
            p1.normal_y = 0.0;
            p1.normal_z = 0.0;
            point_and_normal_pair_3d p2 = new point_and_normal_pair_3d(p1);
            Console.Write(p2.x + " " + p2.y + " " + p2.z + " " + p2.normal_x + " " + p2.normal_y + " " + p2.normal_z);
        }
    }
}
