use std::ops::Index;

use crate::{irreducible_state::IrreducibleState, state::State, basis_element::BasisElement};

#[derive(Debug)]
pub struct CompositeState {
    name: String,
    states: Vec<Box<dyn State>>,
    dims: Vec<usize>,
}

impl CompositeState {
    pub fn new(name: &str, mut states: Vec<Box<dyn State>>) -> CompositeState {
        let dims = states
            .iter()
            .map(|state| state.dim())
            .collect::<Vec<usize>>();

        states.iter_mut().enumerate().for_each(|(i, state)| {
            state.set_id(i);
        });

        CompositeState { name: name.to_string(), states, dims }
    }

    pub fn add_state(&mut self, mut state: Box<dyn State>) {
        self.dims.push(state.dim());

        state.set_id(self.states.len());
        self.states.push(state);
    }

    pub fn dim(&self) -> usize {
        self.dims.iter().product()
    }

    pub fn state_index(&self, name: &str) -> usize {
        let index = self.states.iter().position(|state| state.name() == name);

        if let Some(index) = index {
            index
        } else {
            panic!("State with name {} not found in {}", name, self.name);
        }
    }

    pub fn state_names(&self) -> Vec<&str> {
        self.states.iter().map(|state| state.name()).collect()
    }

    pub fn state_from_base(&self, base: &BasisElement) -> &IrreducibleState {
        &self.states[base.state_id().unwrap()]
            .get_irreducible(base.irreducible_state_id().unwrap())
            .unwrap()
    }

    pub fn states_from_index(&self, index: usize) -> Vec<BasisElement> {
        let mut indices = Vec::with_capacity(self.states.len());

        let mut index = index;

        for dim in self.dims.iter() {
            indices.push(index % dim);
            index /= dim;
        }

        indices
            .iter()
            .zip(self.states.iter())
            .map(|(index, state)| *state.basis()[*index])
            .collect()
    }

    pub fn name(&self) -> &str {
        &self.name
    }
}

impl Index<usize> for CompositeState {
    type Output = Box<dyn State>;

    fn index(&self, index: usize) -> &Self::Output {
        &self.states[index]
    }
}
