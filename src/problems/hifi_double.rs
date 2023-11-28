use std::rc::Rc;

use abm::{
    composite_state::CompositeState, diagonalization::Diagonalization,
    direct_sum_state::DirectSumState, matrix_builder::MatrixBuilder, operator::Operator,
    spin_operators::SpinOperators, state_factory::create_spin,
};
use quantum::{saving::save_param_change, utility::linspace};

pub struct HifiDouble;

impl HifiDouble {
    fn prepare() -> CompositeState {
        let s1 = Box::new(create_spin("s1", 1));
        let i1 = Box::new(create_spin("i1", 1));
        let s2 = Box::new(create_spin("s2", 1));
        let i2 = Box::new(create_spin("i2", 3));

        CompositeState::new("spin", vec![s1, i1, s2, i2])
    }

    pub fn run() {
        println!("Solving hyperfine two atoms...");

        let a_hifi_1 = 2.0;
        let a_hifi_2 = 2.2;
        let gamma_e = 1.0e-2;
        let gamma_i1 = -2.2e-4;
        let gamma_i2 = -1.8e-4;

        let total_double_m = 0;

        let spin_states = Rc::new(Self::prepare());
        let mut matrix_builder = MatrixBuilder::new(spin_states.clone());
        matrix_builder.filter_states(|q_states| {
            q_states.iter().map(|s| s.q_number()).sum::<isize>() == total_double_m
        });
        let spin_ops = SpinOperators::new(spin_states.clone());

        let gammas = [gamma_e, gamma_i1, gamma_e, gamma_i2];
        let mut zeeman_prop = Operator::new(spin_states.clone());
        spin_ops.zeeman_operator(&mut zeeman_prop, &["s1", "i1", "s2", "i2"], &gammas);
        let zeeman_prop_matrix = matrix_builder.from_operator(&zeeman_prop);

        let mut hifi = Operator::new(spin_states.clone());
        hifi.add_operator(vec!["s1", "i1"], |brakets| {
            a_hifi_1 * spin_ops.dot([&brakets["s1"], &brakets["i1"]])
        })
        .add_operator(vec!["s2", "i2"], |brakets| {
            a_hifi_2 * spin_ops.dot([&brakets["s2"], &brakets["i2"]])
        });
        let hifi_matrix = matrix_builder.from_operator(&hifi);

        let mag_fields = linspace(0.0, 1000.0, 1000);
        let mut values = Vec::new();
        let mut diagonalization = Diagonalization::new(hifi_matrix.clone(), None);
        for mag_field in &mag_fields {
            let hamiltonian = &hifi_matrix + &zeeman_prop_matrix * *mag_field;
            diagonalization.diagonalize(hamiltonian);

            values.push(diagonalization.eigvalues().clone());
        }

        save_param_change(
            "hifi_double",
            mag_fields,
            values,
            vec!["magnetic_field", "energies"],
        )
        .unwrap()
    }

    pub fn run_different() {
        println!("Solving hyperfine two atoms...");

        let a_hifi_1 = 2.0;
        let a_hifi_2 = 2.2;
        let gamma_e = 1.0e-2;
        let gamma_i1 = -2.2e-4;
        let gamma_i2 = -1.8e-4;

        let total_double_m = 0;

        let spin_states = Rc::new({
            let s_tot = Box::new(DirectSumState::new(
                "S_tot",
                vec![create_spin("triplet", 2), create_spin("singlet", 0)],
            ));
            let i1 = Box::new(create_spin("i1", 1));
            let i2 = Box::new(create_spin("i2", 3));

            CompositeState::new("spin", vec![s_tot, i1, i2])
        });

        let mut matrix_builder = MatrixBuilder::new(spin_states.clone());
        matrix_builder.filter_states(|q_states| {
            q_states.iter().map(|s| s.q_number()).sum::<isize>() == total_double_m
        });
        let spin_ops = SpinOperators::new(spin_states.clone());

        let gammas = [gamma_e, gamma_i1, gamma_i2];
        let mut zeeman_prop = Operator::new(spin_states.clone());
        spin_ops.zeeman_operator(&mut zeeman_prop, &["S_tot", "i1", "i2"], &gammas);
        let zeeman_prop_matrix = matrix_builder.from_operator(&zeeman_prop);

        let mut hifi = Operator::new(spin_states.clone());
        hifi.add_operator(vec!["S_tot", "i1", "i2"], |brakets| {
            a_hifi_1 / 2.0 * spin_ops.dot([&brakets["S_tot"], &brakets["i1"]])
                + a_hifi_2 / 2.0 * spin_ops.dot([&brakets["S_tot"], &brakets["i2"]])
        })
        .add_operator(vec!["S_tot", "i1", "i2"], |brakets| {
            a_hifi_1 / 2.0 * spin_ops.dot_minus([&brakets["S_tot"], &brakets["i1"]])
                - a_hifi_2 / 2.0 * spin_ops.dot_minus([&brakets["S_tot"], &brakets["i2"]])
        });
        let hifi_matrix = matrix_builder.from_operator(&hifi);

        let mag_fields = linspace(0.0, 1000.0, 1000);
        let mut values = Vec::new();
        let mut diagonalization = Diagonalization::new(hifi_matrix.clone(), None);
        for mag_field in &mag_fields {
            let hamiltonian = &hifi_matrix + &zeeman_prop_matrix * *mag_field;
            diagonalization.diagonalize(hamiltonian);

            values.push(diagonalization.eigvalues().clone());
        }

        save_param_change(
            "hifi_double",
            mag_fields,
            values,
            vec!["magnetic_field", "energies"],
        )
        .unwrap()
    }
}
