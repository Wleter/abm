use std::fmt;

use crate::{irreducible_state::IrreducibleState, basis_element::BasisElement};

pub trait State: fmt::Debug {
    fn name(&self) -> &str;

    fn dim(&self) -> usize;

    fn basis(&self) -> Vec<&BasisElement>;

    fn id(&self) -> Option<usize>;

    fn set_id(&mut self, id: usize);

    fn get_irreducible(&self, index: usize) -> Option<&IrreducibleState>;
}
