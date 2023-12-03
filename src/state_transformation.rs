use nalgebra::DMatrix;

use crate::{matrix_builder::MatrixBuilder, basis_element::BasisElement};

pub struct StateTransformation<'a, F: Fn(&[BasisElement], &[BasisElement]) -> f64> {
    builder_from: &'a MatrixBuilder,
    builder_to: &'a MatrixBuilder,
    // transform_from: Vec<&'static str>,
    // transform_into: Vec<&'static str>,
    transformation: F,
}

impl<'a, F: Fn(&[BasisElement], &[BasisElement]) -> f64> StateTransformation<'a, F> {
    pub fn new(
        builder_from: &'a MatrixBuilder,
        builder_to: &'a MatrixBuilder,
        // transform_from: Vec<&'static str>,
        // transform_into: Vec<&'static str>,
        transformation: F,
    ) -> Self {
        Self {
            builder_from,
            builder_to,
            // transform_from,
            // transform_into,
            transformation,
        }
    }

    pub fn matrix_transformation(&self) -> DMatrix<f64> {
        let dim_from = self.builder_from.dim();
        let dim_to = self.builder_to.dim();

        let mut matrix = DMatrix::zeros(dim_to, dim_from);

        for (i, state_from) in self.builder_from.basis_states().iter().enumerate() {
            for (j, state_to) in self.builder_to.basis_states().iter().enumerate() {
                matrix[(j, i)] = (self.transformation)(state_from, state_to);
            }
        }

        matrix
    }
}
