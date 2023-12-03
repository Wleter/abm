use crate::{
    irreducible_state::IrreducibleState, spin_operators::projection_range, basis_element::BasisElement, direct_sum_state::DirectSumState, state::State,
};

/// Creates spin state with basis of projection of the spin, with values of the projection
/// and quantum number as double the projection value.
pub fn create_spin(name: &str, double_spin: usize) -> IrreducibleState {
    let basis = projection_range(double_spin)
        .map(|m| BasisElement::new(m as f64 / 2.0, m))
        .collect();

    IrreducibleState::new(name, "spin", double_spin as isize, basis)
}

pub fn create_spin_sum(name: &str, spin1: &IrreducibleState, spin2: &IrreducibleState) -> DirectSumState {
    assert!(spin1.algebra() == "spin" && spin2.algebra() == "spin", "Spin sum must be spin states");

    if spin1.q_number() == 1 && spin2.q_number() == 1 {
        return DirectSumState::new(name, vec![create_spin("triplet", 2), create_spin("singlet", 0)])
    }
    
    let largest_dspin = (spin1.q_number() + spin2.q_number()) as usize;
    let smallest_dspin = (spin1.q_number() - spin2.q_number()).abs() as usize;

    let mut states = Vec::with_capacity(spin1.dim() * spin2.dim());

    for (i, double_spin) in (smallest_dspin..=largest_dspin).rev().step_by(2).enumerate() {

        let name = if i % 2 == 0 {
            format!("symmetric {}", double_spin)
        } else {
            format!("anti-symmetric {}", double_spin)
        };

        states.push(create_spin(&name, double_spin));
    }

    DirectSumState::new(name, states)
}