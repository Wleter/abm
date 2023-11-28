use std::fmt;

use crate::{irreducible_state::IrreducibleState, state_base::StateBase};

pub trait State: fmt::Debug {
    fn name(&self) -> &str;

    fn dim(&self) -> usize;

    fn basis(&self) -> Vec<&StateBase>;

    fn id(&self) -> Option<usize>;

    fn set_id(&mut self, id: usize);

    fn get_irreducible(&self, index: usize) -> Option<&IrreducibleState>;
}
