namespace matrix_library
{
    /// <summary>
    /// Class for reprsenting a matrix
    /// </summary>
    public class matrix
    {
        // Matrix data
        protected double[,] the_matrix = null;

        // Matrix creation
        /// <summary>
        /// Creates an empty matrix
        /// </summary>
        public matrix()
        {
        }
        /// <summary>
        /// Constructor which reads matrix from a file
        /// </summary>
        /// <param name="file_name">Path and filename</param>
        public matrix(string file_name)
        {
            set_from_file(file_name);
        }
        /// <summary>
        /// Constructor which makes a SIZExSIZE identity matrix
        /// </summary>
        /// <param name="identity_size">Number of rows and columns for the identity matrix</param>
        public matrix(int identity_size)
        {
            if (identity_size <= 0)
                throw new System.ArgumentException("Invalid identity size");

            the_matrix = new double[identity_size, identity_size];

            // Assign diagnols to 1
            for (int i = 0; i < identity_size; ++i)
                the_matrix[i, i] = 1;
        }
        /// <summary>
        /// Constructor which makes a ROWSxCOLS matrix with all 0's
        /// </summary>
        /// <param name="num_rows">Number of rows</param>
        /// <param name="num_cols">Number of columns</param>
        public matrix(int num_rows, int num_cols)
        {
            if (num_rows == 0 && num_cols == 0)
                return;  // Make an empty matrix
            else if (num_rows > 0 && num_cols > 0)
                the_matrix = new double[num_rows, num_cols];
            else
                throw new System.ArgumentException("Invalid matrix size");
        }
        /// <summary>
        /// Constructor which deep copies a matrix
        /// </summary>
        /// <param name="matrix_to_copy">The matrix to deep copy</param>
        public matrix(matrix matrix_to_copy)
        {
            if (matrix_to_copy.is_empty)
                return;

            the_matrix = new double[matrix_to_copy.rows, matrix_to_copy.cols];
            for (int i = 0; i < matrix_to_copy.rows; ++i)
                for (int j = 0; j < matrix_to_copy.cols; ++j)
                    the_matrix[i, j] = matrix_to_copy.the_matrix[i, j];
        }
        /// <summary>
        /// Resets the matrix data from file
        /// Note:  this overwrites the current data of 'this' matrix
        /// </summary>
        /// <param name="file_name">Path and filename</param>
        public void set_from_file(string file_name)
        {
            System.IO.StreamReader filter_file = null;
            try
            {
                filter_file = new System.IO.StreamReader(file_name);
                string size_line = filter_file.ReadLine();
                string whole_file = filter_file.ReadToEnd();  // Read the rest of the file
                string[] lines = whole_file.Split('\n');  // Split file by line

                // Get matrix size
                string[] values = size_line.Split(' ');  // Split line by individual values
                int num_rows = System.Convert.ToInt32(values[0]);
                int num_cols = System.Convert.ToInt32(values[1]);
                if (num_rows <= 0 || num_cols <= 0)  // Empty matrix
                    throw new System.FormatException("Invalid number of rows/columns in input file");

                // Allocate space
                the_matrix = new double[num_rows, num_cols];

                // Populate matrix
                for (int i = 0; i < num_rows; ++i)
                {
                    values = lines[i].Split(' ');  // Split line by individual values
                    for (int j = 0; j < num_cols; ++j)
                        the_matrix[i, j] = System.Convert.ToDouble(values[j]);
                }
            }
            catch
            {
                the_matrix = null;
                throw new System.FormatException("Input file improperly formated");
            }
            finally
            {
                try
                {
                    filter_file.Close();
                }
                catch { }  // Nothing can be done

                if (Changed != null) Changed(this, System.EventArgs.Empty);  // Inform the user that the matrix has changed
            }
        }


        // Mathmatical functions
        /// <summary>
        /// Determines if this == rhs
        /// </summary>
        /// <param name="rhs">Matrix to check if 'this' is equal to</param>
        /// <returns>True if the matrices contain the same rows, columns, and values.  False otherwise.</returns>
        public bool Equals(matrix rhs)
        {
            return Equals_threshold(rhs, 0);
        }
        /// <summary>
        /// Determines if this == rhs using thresholding
        /// </summary>
        /// <param name="rhs">Matrix to check if 'this' is equal to</param>
        /// <param name="threshold">Threshold to determine equality</param>
        /// <returns>True if the matrices contain the same rows, columns, and values (within a threshold).  False otherwise.</returns>
        public bool Equals_threshold(matrix rhs, double threshold)
        {
            // Check if dimensions are the same
            if (rows != rhs.rows || cols != rhs.cols)
                return false;

            // Deep check values of the matrix
            for (int i = 0; i < rows; ++i)
                for (int j = 0; j < cols; ++j)
                    if (the_matrix[i, j] < rhs.the_matrix[i, j] - threshold || the_matrix[i, j] > rhs.the_matrix[i, j] + threshold)
                        return false;

            return true;
        }
        /// <summary>
        /// Performs element-wise matrix addition
        /// </summary>
        /// <param name="lhs">This matrix</param>
        /// <param name="rhs">Right-hand-side of the + operator</param>
        /// <returns>New matrix</returns>
        public static matrix operator +(matrix lhs, matrix rhs)
        {
            // Check if either matrix is empty
            if (lhs.is_empty || rhs.is_empty)
                throw new System.ArgumentException("Matrix operation dimension mismatch");

            // Check dimensions for proper addition
            if (lhs.rows != rhs.rows || lhs.cols != rhs.cols)
                throw new System.ArgumentException("Matrix operation dimension mismatch");

            // Add the matrices
            matrix new_matrix = new matrix(lhs.rows, lhs.cols);
            for (int i = 0; i < lhs.rows; ++i)
                for (int j = 0; j < lhs.cols; ++j)
                    new_matrix.the_matrix[i, j] = lhs.the_matrix[i, j] + rhs.the_matrix[i, j];

            return new_matrix;
        }
        /// <summary>
        /// Performs element-wise x + m
        /// </summary>
        /// <param name="x">The number to add to each element</param>
        /// <returns>New matrix</returns>
        public matrix element_add(double x)
        {
            // Add
            matrix new_matrix = new matrix(rows, cols);
            for (int i = 0; i < rows; ++i)
                for (int j = 0; j < cols; ++j)
                    new_matrix.the_matrix[i, j] = the_matrix[i, j] + x;

            return new_matrix;
        }
        /// <summary>
        /// Performs element-wise matrix subtraction
        /// </summary>
        /// <param name="lhs">This matrix</param>
        /// <param name="rhs">Right-hand-side of the - operator</param>
        /// <returns>New matrix</returns>
        public static matrix operator -(matrix lhs, matrix rhs)
        {
            // Check dimensions for proper addition
            if (lhs.rows != rhs.rows || lhs.cols != rhs.cols)
                throw new System.ArgumentException("Matrix operation dimension mismatch");

            // Subtract the matrices
            matrix new_matrix = new matrix(lhs.rows, lhs.cols);
            for (int i = 0; i < lhs.rows; ++i)
                for (int j = 0; j < lhs.cols; ++j)
                    new_matrix.the_matrix[i, j] = lhs.the_matrix[i, j] - rhs.the_matrix[i, j];

            return new_matrix;
        }
        /// <summary>
        /// Performs element-wise m - x
        /// </summary>
        /// <param name="x">Value to subtract from each element</param>
        /// <returns>New matrix</returns>
        public matrix element_subtract_by(double x)
        {
            // Subtract
            matrix new_matrix = new matrix(rows, cols);
            for (int i = 0; i < rows; ++i)
                for (int j = 0; j < cols; ++j)
                    new_matrix.the_matrix[i, j] = the_matrix[i, j] - x;

            return new_matrix;
        }
        /// <summary>
        /// Performs element-wise x - m
        /// </summary>
        /// <param name="x">Value to subtract each element from</param>
        /// <returns>New matrix</returns>
        public matrix element_subtract_from(double x)
        {
            // Subtract
            matrix new_matrix = new matrix(rows, cols);
            for (int i = 0; i < rows; ++i)
                for (int j = 0; j < cols; ++j)
                    new_matrix.the_matrix[i, j] = x - the_matrix[i, j];

            return new_matrix;
        }
        /// <summary>
        /// Performs matrix multiplication
        /// </summary>
        /// <param name="lhs">This matrix</param>
        /// <param name="rhs">Right-hand-side of the * operator</param>
        /// <returns>New matrix</returns>
        public static matrix operator *(matrix lhs, matrix rhs)
        {
            // Check inner dimension for proper multiplication
            if (lhs.cols != rhs.rows)
                throw new System.ArgumentException("Matrix operation dimension mismatch");

            // Multiply the matrices
            double sum;
            matrix new_matrix = new matrix(lhs.rows, rhs.cols);
            for (int i = 0; i < lhs.rows; ++i)
                for (int j = 0; j < rhs.cols; ++j)
                {
                    sum = 0;
                    for (int k = 0; k < lhs.cols; ++k)
                        sum += lhs.the_matrix[i, k] * rhs.the_matrix[k, j];
                    new_matrix.the_matrix[i, j] = sum;
                }

            return new_matrix;
        }
        /// <summary>
        /// Performs m .* rhs
        /// </summary>
        /// <param name="rhs">>Right-hand-side of the .* operator</param>
        /// <returns>New matrix</returns>
        public matrix element_multiply(matrix rhs)
        {
            // Check dimensions for proper addition
            if (rows != rhs.rows || cols != rhs.cols)
                throw new System.ArgumentException("Matrix operation dimension mismatch");

            // Multiply
            matrix new_matrix = new matrix(rows, cols);
            for (int i = 0; i < rows; ++i)
                for (int j = 0; j < cols; ++j)
                    new_matrix.the_matrix[i, j] = the_matrix[i, j] * rhs.the_matrix[i, j];

            return new_matrix;
        }
        /// <summary>
        /// Performs element-wise m * x
        /// </summary>
        /// <param name="x">Value to multiply each element by</param>
        /// <returns>New matrix</returns>
        public matrix element_multiply(double x)
        {
            // Multiply
            matrix new_matrix = new matrix(rows, cols);
            for (int i = 0; i < rows; ++i)
                for (int j = 0; j < cols; ++j)
                    new_matrix.the_matrix[i, j] = the_matrix[i, j] * x;

            return new_matrix;
        }
        /// <summary>
        /// Performs element-wise m / x
        /// </summary>
        /// <param name="x">Value to divide each element by</param>
        /// <returns>New matrix</returns>
        public matrix element_divide_by(double x)
        {
            // Divide
            matrix new_matrix = new matrix(rows, cols);
            for (int i = 0; i < rows; ++i)
                for (int j = 0; j < cols; ++j)
                    new_matrix.the_matrix[i, j] = the_matrix[i, j] / x;

            return new_matrix;
        }
        /// <summary>
        /// Performs element wise x / m
        /// </summary>
        /// <param name="x">Value to have each element divide by</param>
        /// <returns>New matrix</returns>
        public matrix element_divide_denom(double x)
        {
            // Divide
            matrix new_matrix = new matrix(rows, cols);
            for (int i = 0; i < rows; ++i)
                for (int j = 0; j < cols; ++j)
                    new_matrix.the_matrix[i, j] = x / the_matrix[i, j];

            return new_matrix;
        }
        /// <summary>
        /// Performs element wise m^x
        /// </summary>
        /// <param name="x">Value to raise each element to</param>
        /// <returns>New matrix</returns>
        public matrix element_raise_to_power(double x)
        {
            // Raise to power
            matrix new_matrix = new matrix(rows, cols);
            for (int i = 0; i < rows; ++i)
                for (int j = 0; j < cols; ++j)
                    new_matrix.the_matrix[i, j] = System.Math.Pow(the_matrix[i, j], x);

            return new_matrix;
        }
        /// <summary>
        /// Performs element wise x^m
        /// </summary>
        /// <param name="x">Value to raise by each element</param>
        /// <returns>New matrix</returns>
        public matrix element_power_to_raise(double x)
        {
            // Power to raise
            matrix new_matrix = new matrix(rows, cols);
            for (int i = 0; i < rows; ++i)
                for (int j = 0; j < cols; ++j)
                    new_matrix.the_matrix[i, j] = System.Math.Pow(x, the_matrix[i, j]);

            return new_matrix;
        }
        /// <summary>
        /// Performs element wise e^m
        /// </summary>
        /// <returns>New matrix</returns>
        public matrix element_exp()
        {
            // e^m
            matrix new_matrix = new matrix(rows, cols);
            for (int i = 0; i < rows; ++i)
                for (int j = 0; j < cols; ++j)
                    new_matrix.the_matrix[i, j] = System.Math.Exp(the_matrix[i, j]);

            return new_matrix;
        }
        /// <summary>
        /// Performs element wise tanh
        /// </summary>
        /// <returns>New matrix</returns>
        public matrix element_tanh()
        {
            // tanh
            matrix new_matrix = new matrix(rows, cols);
            for (int i = 0; i < rows; ++i)
                for (int j = 0; j < cols; ++j)
                    new_matrix.the_matrix[i, j] = System.Math.Tanh(the_matrix[i, j]);

            return new_matrix;
        }


        // Maxtrix functions
        /// <summary>
        /// Transposes the matrix
        /// </summary>
        /// <returns>New matrix</returns>
        public matrix transpose()
        {
            // Check if matrix is empty
            if (is_empty)
                return new matrix();

            // Transpose
            matrix new_matrix = new matrix(cols, rows);
            for (int i = 0; i < cols; ++i)
                for (int j = 0; j < rows; ++j)
                    new_matrix.the_matrix[i, j] = the_matrix[j, i];

            return new_matrix;
        }
        /// <summary>
        /// Removes a row from the matrix
        /// </summary>
        /// <param name="row_to_remove">Row to remove.  0 indexed.</param>
        /// <returns>New matrix</returns>
        public matrix remove_row(int row_to_remove)
        {
            // Check for valid row to remove
            if (row_to_remove < 0 || row_to_remove >= rows)
                throw new System.ArgumentException("Invalid row to remove");

            // Check if the current matrix is empty or only has one row
            if (is_empty || rows == 1)
                return new matrix();

            // Remove row
            int new_row = 0;
            matrix new_matrix = new matrix(rows - 1, cols);
            for (int i = 0; i < rows; ++i)
            {
                if (i == row_to_remove)
                    continue;
                for (int j = 0; j < cols; ++j)
                    new_matrix.the_matrix[new_row, j] = the_matrix[i, j];
                ++new_row;
            }

            return new_matrix;
        }
        /// <summary>
        /// Removes a column from the matrix
        /// </summary>
        /// <param name="col_to_remove">Column to remove.  0 indexed.</param>
        /// <returns>New matrix</returns>
        public matrix remove_col(int col_to_remove)
        {
            // Check for valid row to remove
            if (col_to_remove < 0 || col_to_remove >= cols)
                throw new System.ArgumentException("Invalid column to remove");

            // Check if the current matrix is empty or only has one row
            if (is_empty || cols == 1)
                return new matrix();

            // Remove col
            int new_col;
            matrix new_matrix = new matrix(rows, cols - 1);
            for (int i = 0; i < rows; ++i)
            {
                new_col = 0;
                for (int j = 0; j < cols; ++j)
                {
                    if (j == col_to_remove)
                        continue;
                    new_matrix.the_matrix[i, new_col] = the_matrix[i, j];
                    ++new_col;
                }
            }

            return new_matrix;
        }
        /// <summary>
        /// Removes a row and a column from the matrix
        /// </summary>
        /// <param name="row_to_remove">Row to remove.  0 indexed.</param>
        /// <param name="col_to_remove">Column to remove.  0 indexed.</param>
        /// <returns>New matrix</returns>
        public matrix remove_row_and_col(int row_to_remove, int col_to_remove)
        {
            matrix new_matrix = remove_row(row_to_remove);
            new_matrix = new_matrix.remove_col(col_to_remove);
            return new_matrix;
        }
        /// <summary>
        /// Creates a sub-matrix
        /// </summary>
        /// <param name="start_row">Starting row to keep.  0 indexed.</param>
        /// <param name="stop_row">Last row to keep.  0 indexed.</param>
        /// <param name="start_col">Starting column to keep.  0 indexed.</param>
        /// <param name="stop_col">Last column to keep.  0 indexed.</param>
        /// <returns>New matrix</returns>
        public matrix sub_matrix(int start_row, int stop_row, int start_col, int stop_col)
        {
            // Check if capable of making submatrix
            if (start_row < 0 || stop_row >= rows || start_row > stop_row || start_col < 0 || stop_col >= cols || start_col > stop_col)
                throw new System.ArgumentException("Invalid values to construct a sub-matrix");

            // Create submatrix
            matrix new_matrix = new matrix(stop_row - start_row + 1, stop_col - start_col + 1);
            int new_row = 0;
            for (int i = 0; i < rows; ++i)
            {
                if (i < start_row || i > stop_row)
                    continue;
                int new_col = 0;
                for (int j = 0; j < cols; ++j)
                {
                    if (j < start_col || j > stop_col)
                        continue;
                    new_matrix.the_matrix[new_row, new_col] = the_matrix[i, j];
                    ++new_col;
                }
                ++new_row;
            }

            return new_matrix;
        }
        /// <summary>
        /// Finds the determinant of the matrix
        /// </summary>
        /// <returns>The determinant</returns>
        public double determinant()
        {
            // Check if the current matrix is empty
            if (is_empty)
                return 0;
            // Check if the current matrix is square
            if (rows != cols)
                return 0;

            // Calculate determinant
            double sum = 0;
            if (cols == 1)
                return the_matrix[0, 0];
            else if (cols == 2)
                return (the_matrix[0, 0] * the_matrix[1, 1] - the_matrix[0, 1] * the_matrix[1, 0]);
            else
            {
                for (int k = 0; k < cols; ++k)
                    sum += System.Math.Pow(-1.0, k) * the_matrix[0, k] * minor(0, k);

                return sum;
            }
        }
        /// <summary>
        /// Calculates the minor of the matrix.  That is the determinant after removing a row and column.
        /// </summary>
        /// <param name="row_to_exclude">Row to exclude.  0 indexed.</param>
        /// <param name="col_to_exclude">Column to exclude.  0 indexed</param>
        /// <returns>The minor</returns>
        public double minor(int row_to_exclude, int col_to_exclude)
        {
            // Check for valid rows and cols to remove
            if (row_to_exclude < 0 || col_to_exclude < 0 || row_to_exclude >= rows || col_to_exclude >= cols)
                throw new System.ArgumentException("Invalid row/column to remove");

            // Check if current matrix is empty
            if (is_empty)
                return 0;
            // Check if current matrix is square
            if (rows != cols)
                return 0;
            // Check if the current matrix has only one row or col
            if (rows == 1 || cols == 1)
                return 0;

            // Calculate minor
            return (remove_row_and_col(row_to_exclude, col_to_exclude)).determinant();
        }
        /// <summary>
        /// Calculates the lower (L) and upper (U) triangle matrices and the permutation matrix (P).  Only works for square matrices
        /// </summary>
        /// <param name="L">Matrix for holding L</param>
        /// <param name="U">Matrix for holding U</param>
        /// <param name="P">Matrix for holding P</param>
        public void LU_decomposition(out matrix L, out matrix U, out matrix P)
        {
            if (!is_square)
                throw new System.ArgumentException("Can't perform LU Decomposition on a non-square matrix");

            L = new matrix(rows);
            U = new matrix(this);
            P = new matrix(rows);

            for (int pivot_row = 0; pivot_row < rows - 1; ++pivot_row)
            {
                // Pivoting:  Swap rows so we use the largest possible pivot element
                int pivot_col = pivot_row;  // equal because the pivot element follows the diagnol
                matrix search_pivot_col = this.sub_matrix(pivot_row, rows - 1, pivot_col, pivot_col);
                double max_value;
                int largest_row, largest_col;
                max_value = search_pivot_col.abs_max(out largest_row, out largest_col);
                largest_row += pivot_row;
                U = U.swap_rows(pivot_row, largest_row);  // Swap the largest pivot row with the current row
                P = P.swap_rows(pivot_row, largest_row);  // Swap the same rows in the P matrix
                for (int i = 0; i < pivot_col; ++i)  // Swap the factors in the same rows in the L matrix
                {
                    double temp_val = L.the_matrix[pivot_row, i];
                    L.the_matrix[pivot_row, i] = L.the_matrix[largest_row, i];
                    L.the_matrix[largest_row, i] = temp_val;
                }

                double pivot_element = U.the_matrix[pivot_row, pivot_col];
                for (int subtract_row = pivot_row + 1; subtract_row < rows; ++subtract_row)
                {
                    double scale_factor = U.the_matrix[subtract_row, pivot_col] / pivot_element;
                    for (int nonzero_cols = pivot_col; nonzero_cols < cols; ++nonzero_cols)
                        U.the_matrix[subtract_row, nonzero_cols] = U.the_matrix[subtract_row, nonzero_cols] - scale_factor * U.the_matrix[pivot_row, nonzero_cols];
                    L.the_matrix[subtract_row, pivot_col] = scale_factor;
                }
            }
        }
        /// <summary>
        /// Calculates the row-echelon form of the matrix
        /// </summary>
        /// <returns>New matrix</returns>
        public matrix row_echelon_form()
        {
            // Nested function for finding the next pivot
            System.Func<int, matrix, System.Tuple<int, int>> find_next_pivot = new System.Func<int, matrix, System.Tuple<int, int>>((current_row, current_matrix) =>
            {
                for (int pivot_col = current_row; pivot_col < cols; ++pivot_col)  // Start on the diagonal
                    for (int pivot_row = current_row; pivot_row < rows; ++pivot_row)  // Start with the current row
                        if (current_matrix.the_matrix[pivot_row, pivot_col] != 0)
                            return new System.Tuple<int, int>(pivot_row, pivot_col);

                // No suitable pivot found
                return new System.Tuple<int, int>(-1, -1);
            });

            // Make a copy of the matrix
            matrix result = new matrix(this);

            for (int current_row = 0; current_row < result.rows && current_row < result.cols; ++current_row)
            {
                // Find the pivot
                System.Tuple<int, int> next_pivot = find_next_pivot(current_row, result);
                if (next_pivot.Item1 == -1 && next_pivot.Item2 == -1)  // No pivot left, we're done
                    break;

                // Swap rows
                result = result.swap_rows(current_row, next_pivot.Item1);

                // Multiply this row to reduce the leading number to 1
                double factor = result.the_matrix[current_row, next_pivot.Item2];  // Save the multiplication factor, otherwise it will get overwritten
                for (int i = next_pivot.Item2; i < result.cols; ++i)
                    result.the_matrix[current_row, i] = result.the_matrix[current_row, i] / factor;

                // Add to all subsequent rows so the column has all 0's in it
                for (int subsequent_row = current_row + 1; subsequent_row < result.rows; ++subsequent_row)
                {
                    factor = result.the_matrix[subsequent_row, next_pivot.Item2];  // Save the addition factor, otherwise it will get overwritten
                    for (int i = next_pivot.Item2; i < result.cols; ++i)
                        result.the_matrix[subsequent_row, i] = result.the_matrix[subsequent_row, i] + -factor * result.the_matrix[current_row, i];
                }
            }

            return result;
        }
        /// <summary>
        /// Calculates the reduced row-echelon form of the matrix
        /// </summary>
        /// <returns>New matrix</returns>
        public matrix reduced_row_echelon_form()
        {
            // Get the matrix in row echelon form
            matrix result = row_echelon_form();

            // Find the last row that has a '1' as its leading number
            for (int current_row = result.rows - 1; current_row >= 0; --current_row)
            {
                int current_col = 0;
                while (current_col < result.cols)
                {
                    if (result.the_matrix[current_row, current_col] == 1)  // Found the pivot
                        break;
                    else
                        ++current_col;
                }

                if (current_col == result.cols)  // This row has all zeros, nothing to do
                    continue;

                // Add to all rows above this one so that the column has all zeros except for the current row (which should be 1)
                for (int previous_row = current_row - 1; previous_row >= 0; --previous_row)
                {
                    double factor = result.the_matrix[previous_row, current_col];  // Save the addition factor, otherwise it will get overwritten
                    for (int i = 0; i < result.cols; ++i)
                        result.the_matrix[previous_row, i] = result.the_matrix[previous_row, i] + -factor * result.the_matrix[current_row, i];
                }
            }

            return result;
        }
        /// <summary>
        /// Finds the maximum value contained within the matrix
        /// </summary>
        /// <returns>The maximum value</returns>
        public double max()
        {
            int row, col;
            return max(out row, out col);
        }
        /// <summary>
        /// Finds the maximum value contained within the matrix
        /// </summary>
        /// <param name="row">Row where the value is</param>
        /// <param name="col">Column where the value is</param>
        /// <returns>The maximum value</returns>
        public double max(out int row, out int col)
        {
            // Check if current matrix is empty
            if (is_empty)
                throw new System.ArgumentException("No maximum value for an empty matrix");

            // Find max
            row = 0;
            col = 0;
            double max = the_matrix[0, 0];
            for (int i = 0; i < rows; ++i)
                for (int j = 0; j < cols; ++j)
                    if (the_matrix[i, j] > max)
                    {
                        row = i;
                        col = j;
                        max = the_matrix[i, j];
                    }

            return max;
        }
        /// <summary>
        /// Finds the maximum absolute value in the matrix
        /// </summary>
        /// <returns>The maximum absolute value</returns>
        public double abs_max()
        {
            int row, col;
            return abs_max(out row, out col);
        }
        /// <summary>
        /// Finds the maximum absolute value in the matrix
        /// </summary>
        /// <param name="row">Row where the value is</param>
        /// <param name="col">Column where the value is</param>
        /// <returns>The maximum absolute value</returns>
        public double abs_max(out int row, out int col)
        {
            // Check if current matrix is empty
            if (is_empty)
                throw new System.ArgumentException("No maximum value for an empty matrix");

            // Find max
            row = 0;
            col = 0;
            double max = the_matrix[0, 0];
            for (int i = 0; i < rows; ++i)
                for (int j = 0; j < cols; ++j)
                    if (System.Math.Abs(the_matrix[i, j]) > System.Math.Abs(max))
                    {
                        row = i;
                        col = j;
                        max = the_matrix[i, j];
                    }

            return max;
        }
        /// <summary>
        /// Finds the minimum value contained within the matrix
        /// </summary>
        /// <returns>The minimum value</returns>
        public double min()
        {
            int row, col;
            return min(out row, out col);
        }
        /// <summary>
        /// Finds the minimum value contained within the matrix
        /// </summary>
        /// <param name="row">Row where the value is</param>
        /// <param name="col">Column where the value is</param>
        /// <returns>The minimum value</returns>
        public double min(out int row, out int col)
        {
            // Check if current matrix is empty
            if (is_empty)
                throw new System.ArgumentException("No mimimum value for an empty matrix");

            // Find min
            row = 0;
            col = 0;
            double min = the_matrix[0, 0];
            for (int i = 0; i < rows; ++i)
                for (int j = 0; j < cols; ++j)
                    if (the_matrix[i, j] < min)
                    {
                        row = i;
                        col = j;
                        min = the_matrix[i, j];
                    }

            return min;
        }
        /// <summary>
        /// Calculates the inverse of the matrix
        /// </summary>
        /// <returns>New matrix</returns>
        public matrix inverse()
        {
            // Check if current matrix is empty
            if (is_empty)
                throw new System.ArgumentException("Can't find an inverse of an empty matrix");

            // Check if current matrix has an inverse
            double det = determinant();
            if (det == 0)
                throw new System.ArithmeticException("This matrix does not have an inverse");

            // Calculate cofactor matrix
            matrix cofactor_matrix = new matrix(rows, cols);
            for (int i = 0; i < rows; ++i)
                for (int j = 0; j < cols; ++j)
                    cofactor_matrix.the_matrix[i, j] = System.Math.Pow(-1.0, i) * System.Math.Pow(-1.0, j) * minor(i, j);
            cofactor_matrix = cofactor_matrix.transpose();

            // Calculate inverse
            matrix new_matrix = new matrix(rows, cols);
            for (int i = 0; i < rows; ++i)
                for (int j = 0; j < cols; ++j)
                    new_matrix.the_matrix[i, j] = (1.0 / det) * cofactor_matrix.the_matrix[i, j];

            return new_matrix;
        }
        /// <summary>
        /// alculates the pseudoinverse of the matrix regarless of whether it has a perfect inverse, uses singular value decomposition
        /// </summary>
        /// <returns>New matrix</returns>
        public matrix pseudoinverse()
        {
            // Calculate the singular value decomposition
            matrix u, w, v;
            svd(out u, out w, out v);

            // Calculate the tolerance for the reciprocal calculation
            double DBL_EPSILON = 2.2204460492503131e-16;  // No DBL_EPSILON in c#
            double tolerance = DBL_EPSILON * System.Math.Max(rows, cols) * w.max();

            // Transpose and replace all non-zero values on the diagnol of 'w' with its reciprocal
            matrix w_plus = w.transpose();
            int smallest_dimension = System.Math.Min(rows, cols);
            for (int i = 0; i < smallest_dimension; ++i)
                if (w_plus.the_matrix[i, i] > tolerance)
                    w_plus.the_matrix[i, i] = 1.0 / w_plus.the_matrix[i, i];

            // Calculate pseudoinverse
            return v * w_plus * u.transpose();
        }
        /// <summary>
        /// Adds a new row to the matrix, where each element of the new row has some value
        /// </summary>
        /// <param name="x">The value each element will have</param>
        /// <returns>New matrix</returns>
        public matrix add_row(double x)
        {
            // Check if the matrix is empty
            if (is_empty)
                throw new System.ArgumentException("Can't add a row to an empty matrix");

            // Copy old matrix data
            matrix new_matrix = new matrix(rows + 1, cols);
            for (int i = 0; i < rows; ++i)
                for (int j = 0; j < cols; ++j)
                    new_matrix.the_matrix[i, j] = the_matrix[i, j];

            // Add extra row
            for (int j = 0; j < new_matrix.cols; ++j)
                new_matrix.the_matrix[new_matrix.rows - 1, j] = x;

            return new_matrix;
        }
        /// <summary>
        /// Adds a new column to the matrix, where each element of the new column has some value
        /// </summary>
        /// <param name="x">The value each element will have</param>
        /// <returns>New matrix</returns>
        public matrix add_col(double x)
        {
            // Check if the matrix is empty
            if (is_empty)
                throw new System.ArgumentException("Can't add a column to an empty matrix");

            // Copy old matrix data
            matrix new_matrix = new matrix(rows, cols + 1);
            for (int i = 0; i < rows; ++i)
                for (int j = 0; j < cols; ++j)
                    new_matrix.the_matrix[i, j] = the_matrix[i, j];

            // Add extra col
            for (int i = 0; i < new_matrix.rows; ++i)
                new_matrix.the_matrix[i, new_matrix.cols - 1] = x;

            return new_matrix;
        }
        /// <summary>
        /// Copies a row in a vector form matrix (nx1)
        /// </summary>
        /// <param name="row">Which row to copy.  0 indexed</param>
        /// <returns>New matrix</returns>
        public matrix get_row(int row)
        {
            // Check for valid row
            if (row < 0 || row >= rows)
                throw new System.ArgumentException("Invalid row to return");
            else
                return sub_matrix(row, row, 0, cols - 1).transpose();
        }
        /// <summary>
        /// Copies a column in a vector form matrix (nx1)
        /// </summary>
        /// <param name="col">Which column to copy</param>
        /// <returns>New matrix</returns>
        public matrix get_col(int col)
        {
            // Check for valid col
            if (col < 0 || col >= cols)
                throw new System.ArgumentException("Invalid column to remove");
            else
            {
                return sub_matrix(0, rows - 1, col, col);
            }
        }
        /// <summary>
        /// Swaps two rows of the matrix
        /// </summary>
        /// <param name="row1">First row.  0 indexed.</param>
        /// <param name="row2">Second row.  0 indexed.</param>
        /// <returns>New matrix</returns>
        public matrix swap_rows(int row1, int row2)
        {
            if (row1 < 0 || row1 >= rows || row2 < 0 || row2 >= rows)
                throw new System.ArgumentException("Invalid rows to swap");

            matrix new_matrix = new matrix(this);

            matrix row_1 = get_row(row1);
            matrix row_2 = get_row(row2);

            // Swap
            for (int i = 0; i < cols; ++i)
            {
                new_matrix.the_matrix[row1, i] = row_2.the_matrix[i, 0];
                new_matrix.the_matrix[row2, i] = row_1.the_matrix[i, 0];
            }

            return new_matrix;
        }
        /// <summary>
        /// Swaps two columns of the matrix
        /// </summary>
        /// <param name="col1">First column.  0 indexed.</param>
        /// <param name="col2">Second column.  0 indexed.</param>
        /// <returns>New matrix</returns>
        public matrix swap_cols(int col1, int col2)
        {
            if (col1 < 0 || col1 >= cols || col2 < 0 || col2 >= cols)
                throw new System.ArgumentException("Invalid columns to sawp");

            matrix new_matrix = new matrix(this);

            matrix col_1 = get_col(col1);
            matrix col_2 = get_col(col2);

            // Swap
            for (int i = 0; i < rows; ++i)
            {
                new_matrix.the_matrix[i, col1] = col_2.the_matrix[i, 0];
                new_matrix.the_matrix[i, col2] = col_1.the_matrix[i, 0];
            }

            return new_matrix;
        }


        // Matrix status functions/properties
        /// <summary>
        /// Returns TRUE if the matrix is empty
        /// </summary>
        public bool is_empty
        {
            get
            {
                if (the_matrix == null)
                    return true;
                else
                    return false;
            }

        }
        /// <summary>
        /// Returns TRUE if the matrix is square
        /// Note: returns false if the matrix is empty
        /// </summary>
        public bool is_square
        {
            get
            {
                if (is_empty)
                    return false;
                else
                    return (rows == cols);
            }
        }
        /// <summary>
        /// Returns TRUE if the matrix is singular (that is if the matrix determinant is 0)
        /// </summary>
        public bool is_singular
        {
            get
            {
                if (determinant() == 0)
                    return true;
                else
                    return false;
            }
        }
        /// <summary>
        ///  Returns TRUE if the matrix is symmetric
        ///  Note: returns false if the matrix is empty
        /// </summary>
        public bool is_symmetric
        {
            get
            {
                if (is_square)
                {
                    for (int i = 0; i < rows; ++i)
                        for (int j = 0; j <= i; ++j)
                            if (the_matrix[i, j] != the_matrix[j, i])
                                return false;

                    return true;
                }
                else
                    return false;
            }
        }
        /// <summary>
        /// Returns TRUE if matrix is 4x4 transformation matrix
        /// Direction vectors must be orthognal and normalized, matrix must be 4x4, and last row must be 0,0,0,1
        /// </summary>
        public bool is_3d_transformation_matrix
        {
            get
            {
                // Check if current matrix is correct size
                if (rows != 4 || cols != 4)
                    return false;

                // Check if last row is correct
                if (the_matrix[3, 0] != 0 || the_matrix[3, 1] != 0 || the_matrix[3, 2] != 0 || the_matrix[3, 3] != 1)
                    return false;

                matrix x_vec = this.sub_matrix(0, 2, 0, 0);
                matrix y_vec = this.sub_matrix(0, 2, 1, 1);
                matrix z_vec = this.sub_matrix(0, 2, 2, 2);

                // Check if vectors are normalized
                if (!x_vec.is_normalized_3d_vector || !y_vec.is_normalized_3d_vector || !z_vec.is_normalized_3d_vector)
                    return false;

                // Check if vectors are orthogonal
                if (!z_vec.Equals_threshold(matrix.cross_product(x_vec, y_vec), 0.01))
                    return false;

                return true;
            }
        }
        /// <summary>
        /// Returns TRUE if matrix contains all 0's
        /// Note: returns false if the matrix is empty
        /// </summary>
        public bool is_null_matrix
        {
            get
            {
                if (is_empty)
                    return false;

                for (int i = 0; i < rows; ++i)
                    for (int j = 0; j < cols; ++j)
                        if (the_matrix[i, j] != 0)  // Contains a non-zero
                            return false;

                return true;  // Contains all 0's
            }
        }
        /// <summary>
        /// Returns TRUE if matrix is a diagonal matrix
        /// Note: returns false if the matrix is empty
        /// </summary>
        public bool is_diagonal_matrix
        {
            get
            {
                if (is_empty)
                    return false;

                for (int i = 0; i < rows; ++i)
                    for (int j = 0; j < cols; ++j)
                    {
                        if (i == j && the_matrix[i, j] == 0)  // Test for diagonal not equal to 0
                            return false;
                        else if (i != j && the_matrix[i, j] != 0)  // Test for non-diagonal equal to 0
                            return false;
                    }

                return true;  // Passed all tests
            }
        }
        /// <summary>
        /// Returns TRUE if matrix is a scaler matrix
        /// Note: returns false if the matrix is empty
        /// </summary>
        public bool is_scaler_matrix
        {
            get
            {
                if (is_diagonal_matrix)
                {
                    int num_of_elements = System.Math.Min(rows, cols);  // Get the number of elements that will be in the vector
                    for (int i = 0; i < num_of_elements; ++i)
                        if (the_matrix[i, i] != the_matrix[0, 0])
                            return false;

                    return true;  // Passed all tests
                }

                return false;  // Not a scaler matrix
            }
        }
        /// <summary>
        /// Returns TRUE if matrix is an identity matrix
        /// Note: returns false if the matrix is empty
        /// </summary>
        public bool is_identiy_matrix
        {
            get
            {
                if (is_square && is_scaler_matrix && the_matrix[0, 0] == 1)
                    return true;
                else
                    return false;
            }
        }
        /// <summary>
        /// Returns TRUE if matrix is nx1
        /// </summary>
        public bool is_vector
        {
            get
            {
                if (!is_empty && cols == 1)
                    return true;
                else
                    return false;
            }
        }
        /// <summary>
        /// Returns TRUE if matrix is 3x1
        /// </summary>
        public bool is_3d_vector
        {
            get
            {
                if (rows == 3 && cols == 1)
                    return true;
                else
                    return false;
            }
        }
        /// <summary>
        /// Returns true if the matrix is 3x1 and normalized
        /// </summary>
        public bool is_normalized_3d_vector
        {
            get
            {
                if (rows != 3 && cols != 1)
                    return false;
                if (!this.Equals_threshold(matrix.normalize(this), 0.01))
                    return false;

                return true;
            }
        }
        /// <summary>
        /// Returns the number of rows in the matrix
        /// </summary>
        public int rows
        {
            get
            {
                if (is_empty)  // Matrix is empty
                    return 0;
                else
                    return the_matrix.GetLength(0);
            }
        }
        /// <summary>
        /// Returns the number of columns in the matrix
        /// </summary>
        public int cols
        {
            get
            {
                if (is_empty)  // Matrix is empty
                    return 0;
                else
                    return the_matrix.GetLength(1);
            }
        }
        /// <summary>
        /// Returns the value of a specfic cell
        /// </summary>
        /// <param name="row">The row of the desired value</param>
        /// <param name="col">The column of the desired value</param>
        /// <returns>The value</returns>
        public double get(int row, int col)
        {
            return this[row, col];
        }
        /// <summary>
        /// Sets a specific cell to a value
        /// </summary>
        /// <param name="row">The row of the value to set</param>
        /// <param name="col">The column of the value to set</param>
        /// <param name="value">The value to store in the matrix</param>
        public void set(int row, int col, double value)
        {
            this[row, col] = value;
        }
        /// <summary>
        /// The [] operator of the matrix
        /// Gets/Sets values in the matrix
        /// </summary>
        /// <param name="row">Row of the value to get/set</param>
        /// <param name="col">Column of the value to get/set</param>
        /// <returns>The value</returns>
        public double this[int row, int col]
        {
            get
            {
                if (row < 0 || row >= rows || col < 0 || col >= cols)
                    throw new System.ArgumentException("Invalid row/column to get");
                else
                    return the_matrix[row, col];
            }
            set
            {
                if (row < 0 || row >= rows || col < 0 || col >= cols)
                    throw new System.ArgumentException("Invalid row/column to set");

                the_matrix[row, col] = value;
                if (Changed != null) Changed(this, System.EventArgs.Empty);  // Inform the user that the matrix has changed
            }
        }
        /// <summary>
        /// Gets/Sets the x axis vector (3x1 matrix) of a transformation matrix
        /// Note:  You can only set this by giving a complete vector.  If you try to set an individual value in the returned vector, that value will NOT be reflected in the matrix.
        /// </summary>
        public matrix x_axis_vector
        {
            get
            {
                // Check if current matrix is a transformation matrix
                if (!is_3d_transformation_matrix)
                    throw new System.ArgumentException("Can't get an axis vector of a non 3d transformation matrix");

                return sub_matrix(0, 2, 0, 0);
            }
            set
            {
                // Check if current matrix is a transformation matrix
                if (!is_3d_transformation_matrix)
                    throw new System.ArgumentException("Can't set an axis vector of a non 3d transformation matrix");

                // Check if input has correct dimensions
                if (!value.is_normalized_3d_vector)
                    throw new System.ArgumentException("Input axis vector is not a 3d vector");

                // Copy data
                the_matrix[0, 0] = value.get(0, 0);
                the_matrix[1, 0] = value.get(1, 0);
                the_matrix[2, 0] = value.get(2, 0);

                if (Changed != null) Changed(this, System.EventArgs.Empty);  // Inform the user that the matrix has changed
            }
        }
        /// <summary>
        /// Gets/Sets the y axis vector (3x1 matrix) of a transformation matrix
        /// Note:  You can only set this by giving a complete vector.  If you try to set an individual value in the returned vector, that value will NOT be reflected in the matrix.
        /// </summary>
        public matrix y_axis_vector
        {
            get
            {
                // Check if current matrix is a transformation matrix
                if (!is_3d_transformation_matrix)
                    throw new System.ArgumentException("Can't get an axis vector of a non 3d transformation matrix");

                return sub_matrix(0, 2, 1, 1);
            }
            set
            {
                // Check if current matrix is a transformation matrix
                if (!is_3d_transformation_matrix)
                    throw new System.ArgumentException("Can't set an axis vector of a non 3d transformation matrix");

                // Check if input has correct dimensions
                if (!value.is_normalized_3d_vector)
                    throw new System.ArgumentException("Input axis vector is not a 3d vector");

                // Copy data
                the_matrix[0, 1] = value.get(0, 0);
                the_matrix[1, 1] = value.get(1, 0);
                the_matrix[2, 1] = value.get(2, 0);

                if (Changed != null) Changed(this, System.EventArgs.Empty);  // Inform the user that the matrix has changed
            }
        }
        /// <summary>
        /// Gets/Sets the z axis vector (3x1 matrix) of a transformation matrix
        /// Note:  You can only set this by giving a complete vector.  If you try to set an individual value in the returned vector, that value will NOT be reflected in the matrix.
        /// </summary>
        public matrix z_axis_vector
        {
            get
            {
                // Check if current matrix is a transformation matrix
                if (!is_3d_transformation_matrix)
                    throw new System.ArgumentException("Can't get an axis vector of a non 3d transformation matrix");

                return sub_matrix(0, 2, 2, 2);
            }
            set
            {
                // Check if current matrix is a transformation matrix
                if (!is_3d_transformation_matrix)
                    throw new System.ArgumentException("Can't set an axis vector of a non 3d transformation matrix");

                // Check if input has correct dimensions
                if (!value.is_normalized_3d_vector)
                    throw new System.ArgumentException("Input axis vector is not a 3d vector");

                // Copy data
                the_matrix[0, 2] = value.get(0, 0);
                the_matrix[1, 2] = value.get(1, 0);
                the_matrix[2, 2] = value.get(2, 0);

                if (Changed != null) Changed(this, System.EventArgs.Empty);  // Inform the user that the matrix has changed
            }
        }
        /// <summary>
        /// Gets/Sets the position vector (3x1 matrix) of a transformation matrix
        /// Note:  You can only set this by giving a complete vector.  If you try to set an individual value in the returned vector, that value will NOT be reflected in the matrix.
        /// </summary>
        public matrix position_vector
        {
            get
            {
                // Check if current matrix is a transformation matrix
                if (!is_3d_transformation_matrix)
                    throw new System.ArgumentException("Can't get position vector of a non 3d transformation matrix");

                return sub_matrix(0, 2, 3, 3);
            }
            set
            {
                // Check if current matrix is a transformation matrix
                if (!is_3d_transformation_matrix)
                    throw new System.ArgumentException("Can't set position vector of a non 3d transformation matrix");

                // Check if input has correct dimensions
                if (!value.is_3d_vector)
                    throw new System.ArgumentException("Input position vector is not a 3d vector");

                // Copy data
                the_matrix[0, 3] = value.get(0, 0);
                the_matrix[1, 3] = value.get(1, 0);
                the_matrix[2, 3] = value.get(2, 0);

                if (Changed != null) Changed(this, System.EventArgs.Empty);  // Inform the user that the matrix has changed
            }
        }
        /// <summary>
        /// Gets vector that represents the main diagonal of the matrix (all elements m[i,j] where i==j)
        /// </summary>
        public matrix main_diagonal
        {
            get
            {
                int num_of_elements = System.Math.Min(rows, cols);  // Get the number of elements that will be in the vector

                // Copy data from the matrix
                matrix new_matrix = new matrix(num_of_elements, 1);
                for (int i = 0; i < num_of_elements; ++i)
                    new_matrix.the_matrix[i, 0] = the_matrix[i, i];

                return new_matrix;
            }
        }
        /// <summary>
        /// Gets the trace of the matrix
        /// </summary>
        public double trace
        {
            get
            {
                double sum = 0;
                if (is_square)
                {
                    // Get the main diagonal of the square matrix
                    matrix main_diag = main_diagonal;

                    // Sum the values
                    for (int i = 0; i < main_diag.rows; ++i)
                        sum += main_diag.the_matrix[i, 0];
                }
                return sum;
            }
        }
        /// <summary>
        /// Gets the rank of the matrix
        /// </summary>
        public int rank
        {
            get
            {
                matrix rref = reduced_row_echelon_form();
                int result = 0;
                for (int i = 0; i < rref.rows; ++i)
                {
                    bool all_zeros = true;
                    for (int j = 0; j < rref.cols; ++j)
                        if (rref.the_matrix[i, j] != 0)  // Found a non-zero
                        {
                            all_zeros = false;
                            break;
                        }

                    if (all_zeros == false)  // This was a non-zero row
                        ++result;
                }

                return result;
            }
        }


        // Events
        /// <summary>
        /// Delegate for when the matrix is changed
        /// </summary>
        /// <param name="sender"></param>
        /// <param name="e"></param>
        public delegate void MatrixChangedEventHandler(object sender, System.EventArgs e);
        /// <summary>
        /// Event triggered when any values in the matrix are changed
        /// </summary>
        public event MatrixChangedEventHandler Changed;


        // SVD functions
        /// <summary>
        /// Singular value decomposition of the matrix
        /// </summary>
        /// <param name="U">Output matrix for U</param>
        /// <param name="W">Output matrix for W</param>
        /// <param name="V">Output matrix for V</param>
        public void svd(out matrix U, out matrix W, out matrix V)
        {
            // Algorithm obtained from:
            // ???

            int m = rows;
            int n = cols;
            int flag, i, its, j, jj, k;
            int nm = 0;
            int l = 2;
            double anorm, c, f, g, h, s, scale, x, y, z;
            double[] rv1;

            double[,] a = new double[m + 1, n + 1];
            for (int row = 1; row <= m; ++row)
                for (int col = 1; col <= n; ++col)
                    a[row, col] = the_matrix[row - 1, col - 1];
            double[] w;
            w = new double[n + 1];
            double[,] v = new double[n + 1, n + 1];

            rv1 = new double[n + 1];
            g = scale = anorm = 0.0;
            for (i = 1; i <= n; i++)
            {
                l = i + 1;
                rv1[i] = scale * g;
                g = s = scale = 0.0;
                if (i <= m)
                {
                    for (k = i; k <= m; k++)
                        scale += System.Math.Abs(a[k, i]);
                    if (scale != 0)
                    {
                        for (k = i; k <= m; k++)
                        {
                            a[k, i] /= scale;
                            s += a[k, i] * a[k, i];
                        }
                        f = a[i, i];
                        g = -sign(System.Math.Sqrt(s), f);
                        h = f * g - s;
                        a[i, i] = f - g;
                        for (j = l; j <= n; j++)
                        {
                            for (s = 0.0, k = i; k <= m; k++)
                                s += a[k, i] * a[k, j];
                            f = s / h;
                            for (k = i; k <= m; k++)
                                a[k, j] += f * a[k, i];
                        }
                        for (k = i; k <= m; k++)
                            a[k, i] *= scale;
                    }
                }
                w[i] = scale * g;
                g = s = scale = 0.0;
                if (i <= m && i != n)
                {
                    for (k = l; k <= n; k++)
                        scale += System.Math.Abs(a[i, k]);
                    if (scale != 0)
                    {
                        for (k = l; k <= n; k++)
                        {
                            a[i, k] /= scale;
                            s += a[i, k] * a[i, k];
                        }
                        f = a[i, l];
                        g = -sign(System.Math.Sqrt(s), f);
                        h = f * g - s;
                        a[i, l] = f - g;
                        for (k = l; k <= n; k++)
                            rv1[k] = a[i, k] / h;
                        for (j = l; j <= m; j++)
                        {
                            for (s = 0.0, k = l; k <= n; k++)
                                s += a[j, k] * a[i, k];
                            for (k = l; k <= n; k++)
                                a[j, k] += s * rv1[k];
                        }
                        for (k = l; k <= n; k++)
                            a[i, k] *= scale;
                    }
                }
                anorm = dmax(anorm, (System.Math.Abs(w[i]) + System.Math.Abs(rv1[i])));
            }
            for (i = n; i >= 1; i--)
            {
                if (i < n)
                {
                    if (g != 0)
                    {
                        for (j = l; j <= n; j++)
                            v[j, i] = (a[i, j] / a[i, l]) / g;
                        for (j = l; j <= n; j++)
                        {
                            for (s = 0.0, k = l; k <= n; k++)
                                s += a[i, k] * v[k, j];
                            for (k = l; k <= n; k++)
                                v[k, j] += s * v[k, i];
                        }
                    }
                    for (j = l; j <= n; j++)
                        v[i, j] = v[j, i] = 0.0;
                }
                v[i, i] = 1.0;
                g = rv1[i];
                l = i;
            }
            for (i = imin(m, n); i >= 1; i--)
            {
                l = i + 1;
                g = w[i];
                for (j = l; j <= n; j++)
                    a[i, j] = 0.0;
                if (g != 0)
                {
                    g = 1.0 / g;
                    for (j = l; j <= n; j++)
                    {
                        for (s = 0.0, k = l; k <= m; k++)
                            s += a[k, i] * a[k, j];
                        f = (s / a[i, i]) * g;
                        for (k = i; k <= m; k++)
                            a[k, j] += f * a[k, i];
                    }
                    for (j = i; j <= m; j++)
                        a[j, i] *= g;
                }
                else
                    for (j = i; j <= m; j++)
                        a[j, i] = 0.0;
                ++a[i, i];
            }
            for (k = n; k >= 1; k--)
            {
                for (its = 1; its <= 30; its++)
                {
                    flag = 1;
                    for (l = k; l >= 1; l--)
                    {
                        nm = l - 1;
                        if ((System.Math.Abs(rv1[l]) + anorm) == anorm)
                        {
                            flag = 0;
                            break;
                        }
                        if ((System.Math.Abs(w[nm]) + anorm) == anorm)
                            break;
                    }
                    if (flag != 0)
                    {
                        c = 0.0;
                        s = 1.0;
                        for (i = l; i <= k; i++)
                        {
                            f = s * rv1[i];
                            rv1[i] = c * rv1[i];
                            if ((System.Math.Abs(f) + anorm) == anorm)
                                break;
                            g = w[i];
                            h = pythag(f, g);
                            w[i] = h;
                            h = 1.0 / h;
                            c = g * h;
                            s = -f * h;
                            for (j = 1; j <= m; j++)
                            {
                                y = a[j, nm];
                                z = a[j, i];
                                a[j, nm] = y * c + z * s;
                                a[j, i] = z * c - y * s;
                            }
                        }
                    }
                    z = w[k];
                    if (l == k)
                    {
                        if (z < 0.0)
                        {
                            w[k] = -z;
                            for (j = 1; j <= n; j++)
                                v[j, k] = -v[j, k];
                        }
                        break;
                    }
                    if (its == 30)
                    {
                        U = new matrix();
                        W = new matrix();
                        V = new matrix();

                        return;
                    }
                    x = w[l];
                    nm = k - 1;
                    y = w[nm];
                    g = rv1[nm];
                    h = rv1[k];
                    f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
                    g = pythag(f, 1.0);
                    f = ((x - z) * (x + z) + h * ((y / (f + sign(g, f))) - h)) / x;
                    c = s = 1.0;
                    for (j = l; j <= nm; j++)
                    {
                        i = j + 1;
                        g = rv1[i];
                        y = w[i];
                        h = s * g;
                        g = c * g;
                        z = pythag(f, h);
                        rv1[j] = z;
                        c = f / z;
                        s = h / z;
                        f = x * c + g * s;
                        g = g * c - x * s;
                        h = y * s;
                        y *= c;
                        for (jj = 1; jj <= n; jj++)
                        {
                            x = v[jj, j];
                            z = v[jj, i];
                            v[jj, j] = x * c + z * s;
                            v[jj, i] = z * c - x * s;
                        }
                        z = pythag(f, h);
                        w[j] = z;
                        if (z != 0)
                        {
                            z = 1.0 / z;
                            c = f * z;
                            s = h * z;
                        }
                        f = c * g + s * y;
                        x = c * y - s * g;
                        for (jj = 1; jj <= m; jj++)
                        {
                            y = a[jj, j];
                            z = a[jj, i];
                            a[jj, j] = y * c + z * s;
                            a[jj, i] = z * c - y * s;
                        }
                    }
                    rv1[l] = 0.0;
                    rv1[k] = f;
                    w[k] = x;
                }
            }

            U = new matrix(m, n);
            W = new matrix(n, n);
            V = new matrix(n, n);
            for (int i1 = 0; i1 < m; ++i1)
                for (int j1 = 0; j1 < n; ++j1)
                    U.the_matrix[i1, j1] = a[i1 + 1, j1 + 1];
            for (int i1 = 0; i1 < n; ++i1)
                W.the_matrix[i1, i1] = w[i1 + 1];
            for (int i1 = 0; i1 < n; ++i1)
                for (int j1 = 0; j1 < n; ++j1)
                    V.the_matrix[i1, j1] = v[i1 + 1, j1 + 1];

            return;
        }
        private double sign(double a, double b)
        {
            if (b >= 0.0)
                return System.Math.Abs(a);
            else
                return -System.Math.Abs(a);
        }
        private double dmax(double a, double b)
        {
            if (a > b)
                return a;
            else
                return b;
        }
        private int imin(int a, int b)
        {
            if (a < b)
                return a;
            else
                return b;
        }
        private double pythag(double a, double b)
        {
            double absa, absb, temp;
            absa = System.Math.Abs(a);
            absb = System.Math.Abs(b);
            if (absa > absb)
            {
                temp = absb / absa;
                return absa * System.Math.Sqrt(1.0 + temp * temp);
            }
            else
            {
                if (absb == 0.0)
                    return 0.0;
                temp = absa / absb;
                return absb * System.Math.Sqrt(1.0 + temp * temp);
            }
        }


        // Related functions
        /// <summary>
        /// Calculates the cross product of two vectors, only for 3D vectors stored in 'matrix' class form
        /// </summary>
        /// <param name="a">First vector</param>
        /// <param name="b">Second vector</param>
        /// <returns>New matrix</returns>
        public static matrix cross_product(matrix a, matrix b)
        {
            // Check for proper input data sizes
            if (!a.is_3d_vector || !b.is_3d_vector)
                throw new System.ArgumentException("Can't perform cross product of non 3d vectors");

            // Find cross product
            matrix new_vector = new matrix(3, 1);
            new_vector.the_matrix[0, 0] = a.get(1, 0) * b.get(2, 0) - a.get(2, 0) * b.get(1, 0);
            new_vector.the_matrix[1, 0] = a.get(2, 0) * b.get(0, 0) - a.get(0, 0) * b.get(2, 0);
            new_vector.the_matrix[2, 0] = a.get(0, 0) * b.get(1, 0) - a.get(1, 0) * b.get(0, 0);

            return new_vector;
        }
        /// <summary>
        /// Calculates the dot product of two vectors
        /// </summary>
        /// <param name="a">First vector</param>
        /// <param name="b">Second vector</param>
        /// <returns>New matrix</returns>
        public static double dot_product(matrix a, matrix b)
        {
            // Check for proper dimensions
            if (!a.is_vector || !b.is_vector)
                throw new System.ArgumentException("Can't perform dot product on non vectors");
            if (a.rows != b.rows)
                throw new System.ArgumentException("Can't perform dot product on vectors of different size");

            double sum = 0;
            for (int i = 0; i < a.rows; ++i)
                sum += a.get(i, 0) * b.get(i, 0);

            return sum;
        }
        /// <summary>
        /// Calculates the magnitude of a vector stored in 'matrix' class form
        /// </summary>
        /// <param name="vec">The vector</param>
        /// <returns>The magnitude</returns>
        public static double magnitude(matrix vec)
        {
            // Check for vector
            if (!vec.is_vector)
                throw new System.ArgumentException("Can't perform magnitude on non vectors");

            // Find magnitude
            double result = 0;
            for (int i = 0; i < vec.rows; ++i)
                result += vec.get(i, 0) * vec.get(i, 0);

            return System.Math.Sqrt(result);
        }
        /// <summary>
        /// Calculates the normalized version of a 3D vector stored in 'matrix' class form
        /// </summary>
        /// <param name="vec">The vector</param>
        /// <returns>New matrix</returns>
        public static matrix normalize(matrix vec)
        {
            // Check for proper input data sizes
            if (!vec.is_3d_vector)
                throw new System.ArgumentException("Can't normalize non 3d vectors");

            double mag = magnitude(vec);
            if (mag == 0)
                throw new System.ArgumentException("Can't normalize a vector with no magnitude");

            return vec.element_divide_by(mag);
        }
        /// <summary>
        /// Calculates the rounded value of 'value' using the symmetric method
        /// </summary>
        /// <param name="value">The value to round</param>
        /// <returns>Symmetric rounded value</returns>
        public static int symmetric_round(double value)
        {
            if (value > 0.0)
                return (int)System.Math.Floor(value + 0.5);
            else
                return (int)System.Math.Ceiling(value - 0.5);
        }
        /// <summary>
        /// Solves the linear system described by 'a' and 'b', must be in the form of aX=b
        /// </summary>
        /// <param name="a">a</param>
        /// <param name="b">b</param>
        /// <returns>New matrix</returns>
        public static matrix solve_linear_system(matrix a, matrix b)
        {
            // Can we solve this system
            if (a.cols != b.rows || b.cols != 1 || a.inverse().is_empty)
                throw new System.ArgumentException("System of equations is not solvable");

            return a.inverse() * b;
        }
    }
}
