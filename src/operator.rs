use std::rc::Rc;

use crate::{braket::{Braket, Brakets}, composite_state::CompositeState, state_base::StateBase};

pub struct Operator<'a> {
    states: Rc<CompositeState>,
    vec_spaces: Vec<Vec<&'a str>>,
    op_funcs: Vec<Box<dyn Fn(Brakets) -> f64 + 'a>>,
}

impl<'a> Operator<'a> {
    pub fn new(states: Rc<CompositeState>) -> Self {
        Self {
            states,
            vec_spaces: Vec::new(),
            op_funcs: Vec::new(),
        }
    }

    pub fn states(&self) -> &CompositeState {
        &self.states
    }

    pub fn add_operator(
        &mut self,
        spaces: Vec<&'a str>,
        op: impl Fn(Brakets) -> f64 + 'a,
    ) -> &mut Self {
        self.vec_spaces.push(spaces);
        self.op_funcs.push(Box::new(op));

        self
    }

    pub fn reset_operators(&mut self) {
        self.vec_spaces.clear();
        self.op_funcs.clear();
    }

    pub fn values(&self, q_right: &[StateBase], q_left: &[StateBase]) -> f64 {
        let brakets = q_right
            .iter()
            .zip(q_left.iter())
            .map(|(right, left)| Braket::new(right, left))
            .collect::<Vec<Braket>>();

        self.op_funcs
            .iter()
            .zip(self.vec_spaces.iter())
            .map(|(op, spaces)| {
                let indexes = spaces
                    .iter()
                    .map(|space| self.states.state_index(space))
                    .collect::<Vec<usize>>();

                let brakets_red = indexes
                    .iter()
                    .map(|index| brakets[*index])
                    .collect::<Vec<Braket>>();

                let brakets_red = Brakets::new(brakets_red, spaces);

                let remaining_brakets_different = brakets
                    .iter()
                    .enumerate()
                    .filter(|(index, _)| !indexes.contains(index))
                    .fold(false, |different, (_, braket)| {
                        different || !braket.same_state()
                    });

                if !remaining_brakets_different {
                    op(brakets_red)
                } else {
                    0.0
                }
            })
            .sum()
    }
}
