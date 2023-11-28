use std::{fmt::Display, rc::Rc};

use nalgebra::DMatrix;

use crate::{composite_state::CompositeState, operator::Operator, state_base::StateBase};

pub struct MatrixBuilder {
    states: Rc<CompositeState>,
    basis_states: Vec<Vec<StateBase>>,
}

impl MatrixBuilder {
    pub fn new(states: Rc<CompositeState>) -> Self {
        let dim = states.dim();

        let mut quantum_numbers = Vec::with_capacity(dim);
        for i in 0..dim {
            quantum_numbers.push(states.states_from_index(i));
        }

        MatrixBuilder {
            states,
            basis_states: quantum_numbers,
        }
    }

    pub fn states(&self) -> &CompositeState {
        &self.states
    }

    pub fn filter_states(&mut self, func: impl Fn(&[StateBase]) -> bool) {
        self.basis_states = self
            .basis_states
            .iter()
            .filter(|basis| func(&basis))
            .map(|state| state.clone())
            .collect();
    }

    pub fn unfilter_states(&mut self) {
        let dim = self.states.dim();

        self.basis_states = Vec::with_capacity(dim);
        for i in 0..dim {
            self.basis_states.push(self.states.states_from_index(i));
        }
    }

    pub fn basis_states(&self) -> &[Vec<StateBase>] {
        &self.basis_states
    }

    pub fn dim(&self) -> usize {
        self.basis_states.len()
    }

    pub fn from_operator(&self, operator: &Operator) -> DMatrix<f64> {
        assert!(self.states.name() == operator.states().name());
        let dim = self.dim();
        let mut matrix = DMatrix::zeros(dim, dim);

        for i in 0..dim {
            for j in 0..dim {
                matrix[(i, j)] = operator.values(&self.basis_states[i], &self.basis_states[j]);
            }
        }

        matrix
    }
}

impl Display for MatrixBuilder {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut string = String::new();

        for name in self.states.state_names() {
            string.push_str(&format!("|{}>", name));
        }
        string.push_str(":\n");

        for basis in self.basis_states.iter() {
            for state_base in basis.iter() {
                string.push_str(&format!(
                    "|{},{}>",
                    self.states.state_from_base(state_base).q_number(),
                    state_base.q_number()
                ));
            }
            string.push_str("\n");
        }

        write!(f, "{}", string)
    }
}
