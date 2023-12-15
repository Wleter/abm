use crate::{composite_state::CompositeState, state::State};

#[derive(Clone, Copy, PartialEq, Debug)]
pub struct BasisElement {
    value: f64,
    q_number: isize,
    state_id: Option<usize>,
    irr_state_id: Option<usize>,
}

impl BasisElement {
    pub fn new(value: f64, quant_number: isize) -> Self {
        Self {
            value,
            q_number: quant_number,
            state_id: None,
            irr_state_id: None,
        }
    }

    pub fn value(&self) -> f64 {
        self.value
    }

    pub fn q_number(&self) -> isize {
        self.q_number
    }

    pub fn state_id(&self) -> Option<usize> {
        self.state_id
    }

    pub fn irreducible_state_id(&self) -> Option<usize> {
        self.irr_state_id
    }

    pub fn set_state_id(&mut self, state_id: usize) {
        self.state_id = Some(state_id);
    }

    pub fn set_irreducible_state_id(&mut self, irr_state_id: usize) {
        self.irr_state_id = Some(irr_state_id);
    }

    pub fn same_as(
        &self,
        self_states: &CompositeState,
        other: &Self,
        other_states: &CompositeState,
    ) -> bool {
        let self_state = self_states.state_from_base(&self);
        let other_state = other_states.state_from_base(&other);

        self.q_number == other.q_number
            && self_state.q_number() == other_state.q_number()
            && self_state.name() == other_state.name()
    }
}

impl Eq for BasisElement {}
