using System;
using matrix_library;
using geometry_library;

namespace test
{
    class Program
    {
        static void Main(string[] args)
        {
            matrix trans = new matrix(4);
            trans[0, 0] = 1; trans[0, 1] = 0; trans[0, 2] = 0; trans[0, 3] = 0;
            trans[1, 0] = 0; trans[1, 1] = 1; trans[1, 2] = 0; trans[1, 3] = 10;
            trans[2, 0] = 0; trans[2, 1] = 0; trans[2, 2] = 1; trans[2, 3] = 0;
            trans[3, 0] = 0; trans[3, 1] = 0; trans[3, 2] = 0; trans[3, 3] = 1;

            plane p = new plane();
            p.b = 1;
            p = p.transform(trans);

            Console.WriteLine(p.a);
            Console.WriteLine(p.b);
            Console.WriteLine(p.c);
            Console.WriteLine(p.d);
            Console.WriteLine();

            point3d check = new point3d();
            check.y = 15;
            Console.WriteLine(check.distance_to(p));
        }
    }
}
