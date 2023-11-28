use std::rc::Rc;

use abm::{
    composite_state::CompositeState, diagonalization::Diagonalization,
    matrix_builder::MatrixBuilder, operator::Operator, spin_operators::SpinOperators,
    state_factory::create_spin,
};
use quantum::{saving::save_param_change, utility::linspace};

pub struct HifiSingle;

impl HifiSingle {
    fn prepare() -> CompositeState {
        let s1 = Box::new(create_spin("s1", 1));
        let i1 = Box::new(create_spin("i1", 2));

        CompositeState::new("spins_states", vec![s1, i1])
    }

    pub fn run() {
        println!("Solving hyperfine one atom...");

        let a_hifi_1 = 2.0;
        let gamma_e = 1.0e-2;
        let gamma_i1 = -2.0e-4;

        let spin_states = Rc::new(Self::prepare());
        let matrix_builder = MatrixBuilder::new(spin_states.clone());
        let spin_ops = SpinOperators::new(spin_states.clone());

        let mag_fields = linspace(0.0, 1000.0, 1000);
        let mut values = Vec::new();

        let gammas = [gamma_e, gamma_i1];
        let mut zeeman_prop = Operator::new(spin_states.clone());
        spin_ops.zeeman_operator(&mut zeeman_prop, &["s1", "i1"], &gammas);
        let zeeman_prop_matrix = matrix_builder.from_operator(&zeeman_prop);

        let mut hifi = Operator::new(spin_states.clone());
        hifi.add_operator(vec!["s1", "i1"], |brakets| {
            a_hifi_1 * &spin_ops.dot([&brakets["s1"], &brakets["i1"]])
        });
        let hifi_matrix = matrix_builder.from_operator(&hifi);

        let mut diagonalization = Diagonalization::new(hifi_matrix.clone(), None);
        for mag_field in &mag_fields {
            let hamiltonian = &hifi_matrix + &zeeman_prop_matrix * *mag_field;
            diagonalization.continuos_change(hamiltonian);

            values.push(diagonalization.eigvalues().clone());
        }

        save_param_change(
            "hifi_single",
            mag_fields,
            values,
            vec!["magnetic_field", "energies"],
        )
        .unwrap()
    }
}
