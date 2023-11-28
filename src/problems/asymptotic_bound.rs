use std::rc::Rc;

use abm::{
    composite_state::CompositeState, diagonalization::Diagonalization,
    direct_sum_state::DirectSumState, matrix_builder::MatrixBuilder, operator::Operator,
    spin_operators::SpinOperators, state::State, state_factory::create_spin,
    state_transformation::StateTransformation, utils::clebsch_gordan,
};
use quantum::{saving::save_param_change, utility::linspace};

pub struct AsymptoticBound;

impl AsymptoticBound {
    pub fn run() {
        println!("Solving abm for two atoms...");

        let a_hifi_1 = 1.0;
        let a_hifi_2 = 1.0;
        let gamma_e = 2e-2;
        let gamma_i1 = -1.2e-5;
        let gamma_i2 = -1.2e-5;

        let triplet_state = -5.0;
        let singlet_state = -10.0;
        let gordon = 0.3;

        let total_double_m = 0;

        let i1 = create_spin("i1", 2);
        let i2 = create_spin("i2", 2);
        let s_tot = DirectSumState::new(
            "S_tot",
            vec![create_spin("triplet", 2), create_spin("singlet", 0)],
        );
        let i_tot = DirectSumState::new(
            "I_tot",
            vec![
                create_spin("symmetric 4", 4),
                create_spin("anti-symmetric 2", 2),
                create_spin("symmetric 0", 0),
            ],
        );

        let tot_s_states = Rc::new(CompositeState::new(
            "total_s_spin_states",
            vec![
                Box::new(s_tot.clone()),
                Box::new(i1.clone()),
                Box::new(i2.clone()),
            ],
        ));
        let tot_si_states = Rc::new(CompositeState::new(
            "total_si_spin_states",
            vec![Box::new(s_tot.clone()), Box::new(i_tot.clone())],
        ));

        let mut matrix_builder_si = MatrixBuilder::new(tot_si_states.clone());
        matrix_builder_si.filter_states(|q_states| {
            q_states.iter().map(|s| s.q_number()).sum::<isize>() == total_double_m
        });
        matrix_builder_si.filter_states(|q_states| {
            let state_spin_e = tot_si_states.state_from_base(&q_states[0]);
            let state_spin_n = tot_si_states.state_from_base(&q_states[1]);

            (state_spin_e.name() == "singlet" && !state_spin_n.name().contains("anti-symmetric"))
                || (state_spin_e.name() == "triplet"
                    && state_spin_n.name().contains("anti-symmetric"))
        });

        let mut matrix_builder_s = MatrixBuilder::new(tot_s_states.clone());
        matrix_builder_s.filter_states(|q_states| {
            q_states.iter().map(|s| s.q_number()).sum::<isize>() == total_double_m
        });

        let transformation =
            StateTransformation::new(&matrix_builder_s, &matrix_builder_si, |from, to| {
                let d_spin1 = tot_s_states.state_from_base(&from[1]).q_number() as usize;
                let d_spin2 = tot_s_states.state_from_base(&from[2]).q_number() as usize;
                let tot_dspin = tot_si_states.state_from_base(&to[1]).q_number() as usize;

                if from[0].same_as(&tot_s_states, &to[0], &tot_si_states) {
                    clebsch_gordan(
                        d_spin1,
                        from[1].q_number(),
                        d_spin2,
                        from[2].q_number(),
                        tot_dspin,
                        to[1].q_number(),
                    )
                } else {
                    0.0
                }
            });
        let transformation_matrix = transformation.matrix_transformation();

        let mut bound_states = Operator::new(tot_si_states.clone());
        bound_states.add_operator(vec!["S_tot"], |brakets| {
            if brakets["S_tot"].ket == brakets["S_tot"].bra {
                let state_q_number = tot_si_states.state_from_base(brakets["S_tot"].ket).q_number();
                if state_q_number == 0 {
                    return singlet_state;
                } else {
                    assert!(state_q_number == 2);
                    return triplet_state;
                }
            } else {
                0.0
            }
        });
        let bound_states_matrix = matrix_builder_si.from_operator(&bound_states);

        let spin_ops = SpinOperators::new(tot_s_states.clone());

        let gammas = [gamma_e, gamma_i1, gamma_i2];
        let mut zeeman_prop = Operator::new(tot_s_states.clone());
        spin_ops.zeeman_operator(&mut zeeman_prop, &["S_tot", "i1", "i2"], &gammas);
        let zeeman_prop_matrix = matrix_builder_s.from_operator(&zeeman_prop);
        let zeeman_prop_matrix =
            &transformation_matrix * zeeman_prop_matrix * transformation_matrix.transpose();

        let mut hifi = Operator::new(tot_s_states.clone());
        hifi.add_operator(vec!["S_tot", "i1", "i2"], |brakets| {
            // hifi+
            a_hifi_1 / 2.0 * spin_ops.dot([&brakets["S_tot"], &brakets["i1"]])
                + a_hifi_2 / 2.0 * spin_ops.dot([&brakets["S_tot"], &brakets["i2"]])
        })
        .add_operator(vec!["S_tot", "i1", "i2"], |brakets| {
            // hifi-
            gordon
                * (a_hifi_1 / 2.0 * spin_ops.dot_minus([&brakets["S_tot"], &brakets["i1"]])
                    - a_hifi_2 / 2.0 * spin_ops.dot_minus([&brakets["S_tot"], &brakets["i2"]]))
        });
        let hifi_matrix = matrix_builder_s.from_operator(&hifi);
        let hifi_matrix = &transformation_matrix * hifi_matrix * transformation_matrix.transpose();

        let mag_fields = linspace(0.0, 500.0, 1000);
        let mut values = Vec::new();
        let mut diagonalization = Diagonalization::new(hifi_matrix.clone(), None);

        for mag_field in &mag_fields {
            let hamiltonian = &bound_states_matrix + &hifi_matrix + &zeeman_prop_matrix * *mag_field;

            diagonalization.diagonalize(hamiltonian);
            values.push(diagonalization.eigvalues().clone());
        }

        save_param_change(
            "simple_abm",
            mag_fields,
            values,
            vec!["magnetic_field", "energies"],
        )
        .unwrap()
    }
}
