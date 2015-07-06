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
            Console.WriteLine();


            point_and_normal_pair_3d collection_test_1 = new point_and_normal_pair_3d();
            collection_test_1.x = -1;
            collection_test_1.y = -1;
            collection_test_1.z = -1;
            point_and_normal_pair_3d collection_test_2 = new point_and_normal_pair_3d();
            collection_test_2.x = 1;
            collection_test_2.y = -1;
            collection_test_2.z = -1;
            point_and_normal_pair_3d collection_test_3 = new point_and_normal_pair_3d();
            collection_test_3.x = -1;
            collection_test_3.y = 1;
            collection_test_3.z = -1;
            point_and_normal_pair_3d collection_test_4 = new point_and_normal_pair_3d();
            collection_test_4.x = 1;
            collection_test_4.y = 1;
            collection_test_4.z = -1;
            point_and_normal_pair_3d collection_test_5 = new point_and_normal_pair_3d();
            collection_test_5.x = -1;
            collection_test_5.y = -1;
            collection_test_5.z = 1;
            point_and_normal_pair_3d collection_test_6 = new point_and_normal_pair_3d();
            collection_test_6.x = 1;
            collection_test_6.y = -1;
            collection_test_6.z = 1;
            point_and_normal_pair_3d collection_test_7 = new point_and_normal_pair_3d();
            collection_test_7.x = -1;
            collection_test_7.y = 1;
            collection_test_7.z = 1;
            point_and_normal_pair_3d collection_test_8 = new point_and_normal_pair_3d();
            collection_test_8.x = 1;
            collection_test_8.y = 1;
            collection_test_8.z = 1;
            point_3d_collection<point_and_normal_pair_3d> collection = new point_3d_collection<point_and_normal_pair_3d>();
            collection.Add(collection_test_1);
            collection.Add(collection_test_2);
            collection.Add(collection_test_3);
            collection.Add(collection_test_4);
            collection.Add(collection_test_5);
            collection.Add(collection_test_6);
            collection.Add(collection_test_7);
            collection.Add(collection_test_8);
            point3d result_c = collection.centroid;
            Tuple<double, double, double, double, double, double> result_b = collection.bounding_box;
            Console.Write(result_c.x + " " + result_c.y + " " + result_c.z);
            Console.WriteLine();
            Console.Write(result_b.Item1 + " " + result_b.Item2 + " " + result_b.Item3 + " " + result_b.Item4 + " " + result_b.Item5 + " " + result_b.Item6);
            Console.WriteLine();


            Console.WriteLine();
            matrix copy_test1 = new matrix(4);
            matrix copy_test2 = new matrix(3, 1);
            copy_test2[0, 0] = 0; copy_test2[1, 0] = 0; copy_test2[2, 0] = 1;
            copy_test1.x_axis_vector = copy_test2;
            for (int x = 0; x < copy_test1.rows; ++x)
            {
                for (int y = 0; y < copy_test1.cols; ++y)
                    Console.Write(copy_test1[x, y] + " ");
                Console.WriteLine();
            }
        }
    }
}
