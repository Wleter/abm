use std::{mem::swap, rc::Rc};

use crate::{braket::Braket, composite_state::CompositeState, basis_element::BasisElement};

pub struct BoundOperators {
    composite_state: Rc<CompositeState>,
    triplets: Vec<f64>,
    singlets: Vec<f64>,
    fc_factors: Vec<f64>,
}

impl BoundOperators {
    pub fn new(
        composite_state: Rc<CompositeState>,
        bound_states_name: &str,
        triplets: Vec<f64>,
        singlets: Vec<f64>,
        fc_factors: Vec<f64>,
    ) -> Self {
        let b_state_index = composite_state.state_index(bound_states_name);
        let b_state = composite_state[b_state_index]
            .get_irreducible(b_state_index)
            .unwrap();

        assert!(
            b_state.q_number() as usize == triplets.len(),
            "Triplet bounds length is not equal to the number of triplet states"
        );
        assert!(
            b_state.q_number() as usize == singlets.len(),
            "Singlet bounds length is not equal to the number of singlet states"
        );
        assert!(
            b_state.q_number().pow(2) as usize == fc_factors.len(),
            "Wrong number of Franck-Condon factors"
        );

        Self {
            composite_state,
            triplets,
            singlets,
            fc_factors,
        }
    }

    pub fn id(&self, braket: &Braket) -> f64 {
        assert!(self.is_bound_state(braket.bra));
        assert!(self.is_bound_state(braket.ket));

        if braket.ket.q_number() == braket.bra.q_number() {
            1.0
        } else {
            0.0
        }
    }

    pub fn bound_state(&self, state_type: StateType, braket: &Braket) -> f64 {
        assert!(self.is_bound_state(braket.bra));
        assert!(self.is_bound_state(braket.ket));

        if braket.ket.q_number() == braket.bra.q_number() {
            match state_type {
                StateType::Triplet => self.triplets[braket.ket.q_number() as usize],
                StateType::Singlet => self.singlets[braket.bra.q_number() as usize],
            }
        } else {
            0.0
        }
    }

    pub fn state_type_flip(&self, from_state: StateType, braket: &Braket) -> f64 {
        assert!(self.is_bound_state(braket.bra));
        assert!(self.is_bound_state(braket.ket));

        let mut q_right = braket.ket.q_number();
        let mut q_left = braket.bra.q_number();
        if from_state == StateType::Triplet {
            swap(&mut q_right, &mut q_left)
        }

        let dim = self.triplets.len();
        let fc_index = q_right as usize * dim + q_left as usize;

        self.fc_factors[fc_index]
    }

    fn is_bound_state(&self, base: &BasisElement) -> bool {
        self.composite_state.state_from_base(base).algebra() == "bound"
    }
}

#[derive(PartialEq)]
pub enum StateType {
    Triplet,
    Singlet,
}
