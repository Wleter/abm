use nalgebra::DMatrix;
use pyo3::prelude::*;

pub fn assert_almost_eq(a: f64, b: f64, tol: f64) {
    assert!(a.abs() < b.abs() * (1.0 + tol));
    assert!(a.abs() > b.abs() * (1.0 - tol));
}

pub fn from_block_diagonal(blocks: &Vec<DMatrix<f64>>, total_dim: usize) -> DMatrix<f64> {
    let mut matrix = DMatrix::<f64>::zeros(total_dim, total_dim);

    let mut row = 0;
    let mut col = 0;
    for block in blocks {
        set_sub_matrix(&mut matrix, block, row, col);
        row += block.nrows();
        col += block.ncols();
    }

    matrix
}

pub fn set_sub_matrix(
    matrix: &mut DMatrix<f64>,
    sub_matrix: &DMatrix<f64>,
    row: usize,
    col: usize,
) {
    let (sub_rows, sub_cols) = sub_matrix.shape();
    for i in 0..sub_rows {
        for j in 0..sub_cols {
            matrix[(row + i, col + j)] = sub_matrix[(i, j)];
        }
    }
}

pub fn clebsch_gordan(dj1: usize, dm1: isize, dj2: usize, dm2: isize, dj: usize, dm: isize) -> f64 {
    if dm != dm1 + dm2 {
        return 0.0;
    }

    let code = include_str!("clebsch_gordan.py");
    let mut result: f64 = 0.0;

    Python::with_gil(|py| -> PyResult<()> {
        let fun: Py<PyAny> = PyModule::from_code(py, code, "", "")?
            .getattr("clebsch_gordan")?
            .into();

        result = fun.call1(py, (dj1, dm1, dj2, dm2, dj, dm))?.extract(py)?;

        Ok(())
    })
    .unwrap();

    result
}
