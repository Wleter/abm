use nalgebra::DMatrix;

pub struct Diagonalization {
    eigvalues: Vec<f64>,
    eigvectors: DMatrix<f64>,
    rank: usize,
}

impl Diagonalization {
    pub fn new(matrix: DMatrix<f64>, zero_rank: Option<usize>) -> Self {
        if Self::is_symmetric(&matrix) == false {
            println!("Warning: Matrix is not fully symmetric")
        }
        let dim = matrix.nrows();
        let eigen = matrix.symmetric_eigen();
        let eigvalues = eigen.eigenvalues.as_slice().to_vec();

        Self {
            eigvalues,
            eigvectors: eigen.eigenvectors,
            rank: dim - zero_rank.unwrap_or(0),
        }
    }

    pub fn continuos_change(&mut self, matrix: DMatrix<f64>) {
        let eigen = matrix.symmetric_eigen();
        let eigvalues = eigen.eigenvalues.as_slice().to_vec();
        let eigvectors = eigen.eigenvectors;

        let scalar_products = self.eigvectors.transpose() * &eigvectors;

        let mut row_to_old_row: Vec<usize> = Vec::with_capacity(eigvectors.nrows());
        let mut taken_rows = vec![false; eigvectors.nrows()];
        for column in scalar_products.row_iter() {
            let mut max = 0.0;
            let mut max_index = 0;
            for (j, value) in column.iter().enumerate() {
                if value.abs() > max && taken_rows[j] == false {
                    max = value.abs();
                    max_index = j;
                }
            }
            taken_rows[max_index] = true;
            row_to_old_row.push(max_index);
        }

        let mut new_eigvectors = DMatrix::zeros(eigvectors.nrows(), eigvectors.ncols());
        for (i, row) in row_to_old_row.iter().enumerate() {
            new_eigvectors.set_row(i, &eigvectors.row(*row));
        }
        self.eigvectors = new_eigvectors;

        self.eigvalues = row_to_old_row.iter().map(|i| eigvalues[*i]).collect();
    }

    pub fn diagonalize(&mut self, matrix: DMatrix<f64>) {
        let dim = matrix.nrows();
        let eigen = matrix.symmetric_eigen();
        let mut eigvalues = eigen.eigenvalues.as_slice().to_vec();
        let eigvectors = eigen.eigenvectors;

        if dim != self.rank {
            eigvalues.sort_by(|a, b| a.abs().partial_cmp(&b.abs()).unwrap());
            eigvalues = eigvalues.into_iter().skip(dim - self.rank).collect();
        }
        eigvalues.sort_by(|a, b| a.partial_cmp(b).unwrap());

        self.eigvalues = eigvalues;
        self.eigvectors = eigvectors;
    }

    pub fn eigvalues(&self) -> &Vec<f64> {
        &self.eigvalues
    }

    pub fn eigvectors(&self) -> &DMatrix<f64> {
        &self.eigvectors
    }

    fn is_symmetric(matrix: &DMatrix<f64>) -> bool {
        if matrix.nrows() != matrix.ncols() {
            panic!("Matrix is not square")
        }
        let mut is_sym = true;

        // SAFETY: matrix is square and in bounds
        unsafe {
            for i in 0..matrix.nrows() {
                for j in 0..(i+1) {
                    let lower = matrix.get_unchecked((i, j)).abs();
                    let upper = matrix.get_unchecked((j, i)).abs();
                    if lower > upper * 1.00001 || 1.00001 * lower < upper {
                        print!("({i},{j}) {:.3e}, ", lower);
                        println!("({j},{i}) {:.3e}", upper);
                        is_sym = false;
                    }
                }
            }
        }

        is_sym
    }
}
