use crate::{state::State, state_base::StateBase};

#[derive(Clone, Debug)]
pub struct IrreducibleState {
    name: String,
    algebra: String,
    q_number: isize,
    basis: Vec<StateBase>,
    id: Option<usize>,
}

impl IrreducibleState {
    pub fn new(
        name: &str,
        algebra: &str,
        q_number: isize,
        basis: Vec<StateBase>,
    ) -> IrreducibleState {
        IrreducibleState {
            name: name.to_string(),
            algebra: algebra.to_string(),
            q_number,
            basis,
            id: None,
        }
    }

    pub fn algebra(&self) -> &str {
        &self.algebra
    }

    pub fn q_number(&self) -> isize {
        self.q_number
    }

    pub fn set_id_from_parent(&mut self, parent_id: usize) {
        self.basis
            .iter_mut()
            .for_each(|b| b.set_state_id(parent_id));
    }

    pub fn set_irreducible_state_id(&mut self, id: usize) {
        self.id = Some(id);
        self.basis
            .iter_mut()
            .for_each(|b| b.set_irreducible_state_id(id));
    }
}

impl State for IrreducibleState {
    fn name(&self) -> &str {
        &self.name
    }

    fn dim(&self) -> usize {
        self.basis.len()
    }

    fn basis(&self) -> Vec<&StateBase> {
        self.basis.iter().collect()
    }

    fn id(&self) -> Option<usize> {
        self.id
    }

    fn set_id(&mut self, id: usize) {
        self.id = Some(id);
        self.basis.iter_mut().for_each(|b| {
            b.set_state_id(id);
            b.set_irreducible_state_id(id);
        });
    }

    fn get_irreducible(&self, id: usize) -> Option<&IrreducibleState> {
        if id == self.id.unwrap() {
            Some(self)
        } else {
            None
        }
    }
}
