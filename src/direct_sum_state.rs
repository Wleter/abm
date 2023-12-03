use crate::{irreducible_state::IrreducibleState, state::State, basis_element::BasisElement};

#[derive(Debug, Clone)]
pub struct DirectSumState {
    name: String,
    states: Vec<IrreducibleState>,
    dims: Vec<usize>,
    id: Option<usize>,
}

impl DirectSumState {
    pub fn new(name: &str, states: Vec<IrreducibleState>) -> DirectSumState {
        let dims = states
            .iter()
            .map(|state| state.dim())
            .collect::<Vec<usize>>();

        DirectSumState {
            name: name.to_string(),
            states,
            dims,
            id: None,
        }
    }

    pub fn direct_add(&mut self, state: IrreducibleState) {
        self.dims.push(state.dim());
        self.states.push(state);
    }
}

impl State for DirectSumState {
    fn name(&self) -> &str {
        &self.name
    }

    fn dim(&self) -> usize {
        self.dims.iter().sum()
    }

    fn basis(&self) -> Vec<&BasisElement> {
        self.states
            .iter()
            .map(|state| state.basis())
            .flatten()
            .collect()
    }

    fn id(&self) -> Option<usize> {
        self.id
    }

    fn set_id(&mut self, id: usize) {
        self.id = Some(id);

        self.states.iter_mut().enumerate().for_each(|(i, state)| {
            state.set_irreducible_state_id(i);
            state.set_id_from_parent(id);
        });
    }

    fn get_irreducible(&self, id: usize) -> Option<&IrreducibleState> {
        Some(&self.states[id])
    }
}
