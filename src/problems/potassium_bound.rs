use std::{collections::VecDeque, rc::Rc};

use abm::{
    bound_operators::{BoundOperators, StateType},
    composite_state::CompositeState,
    diagonalization::Diagonalization,
    direct_sum_state::DirectSumState,
    irreducible_state::IrreducibleState,
    matrix_builder::MatrixBuilder,
    operator::Operator,
    spin_operators::SpinOperators,
    state::State,
    basis_element::BasisElement,
    state_factory::{create_spin, create_spin_sum},
    state_transformation::StateTransformation,
    utils::clebsch_gordan,
};
use quantum::{
    problem_selector::ProblemSelector,
    saving::save_param_change,
    units::{convert_data_units, energy_units::{Energy, MHz, GHz}, Au},
    utility::linspace,
};

use crate::defaults::Defaults;

pub struct PotassiumBound;

impl ProblemSelector for PotassiumBound {
    const NAME: &'static str = "potassium-potassium";

    fn list() -> Vec<&'static str> {
        vec![
            "hyperfine 39-41",
            "potassium bound 39-41",
            "hyperfine 39-39",
            "potassium bound 39-39",
        ]
    }

    fn methods(number: &str, _args: &mut VecDeque<String>) {
        match number {
            "0" => Self::hifi_39_41(),
            "1" => Self::abm_39_41(),
            "2" => Self::hifi_39_39(),
            "3" => Self::abm_39_39(),
            _ => println!("Invalid problem number"),
        }
    }
}

impl PotassiumBound {
    const HIFI_K39_MHZ: f64 = 230.8595;
    const HIFI_K41_MHZ: f64 = 127.007;
    const GAMMA_I_K39: f64 = 0.0; //1.989344e-4;
    const GAMMA_I_K41: f64 = 0.0; //1.091921e-4;

    pub fn hifi_39_41() {
        println!("Solving hyperfine for K39 - K41...");

        let a_hifi_1 = Energy::new(Self::HIFI_K39_MHZ, MHz).to_au();
        let a_hifi_2 = Energy::new(Self::HIFI_K41_MHZ, MHz).to_au();
        let gamma_e = -2.0 * Defaults::BOHR_MAG;
        let gamma_i1 = Energy::new(Self::GAMMA_I_K39, MHz).to_au() * 1e-4;
        let gamma_i2 = Energy::new(Self::GAMMA_I_K41, MHz).to_au() * 1e-4;

        let total_double_m = 4;

        let s1 = Box::new(create_spin("s1", 1));
        let i1 = Box::new(create_spin("i1", 3));
        let s2 = Box::new(create_spin("s2", 1));
        let i2 = Box::new(create_spin("i2", 3));
        let states = Rc::new(CompositeState::new("spins_states", vec![s1, i1, s2, i2]));

        let mut matrix_builder = MatrixBuilder::new(states.clone());
        matrix_builder.filter_states(|q_states| {
            q_states.iter().map(|s| s.q_number()).sum::<isize>() == total_double_m
        });
        let spin_ops = SpinOperators::new(states.clone());

        let gammas = [gamma_e, gamma_i1, gamma_e, gamma_i2];
        let mut zeeman_prop = Operator::new(states.clone());
        spin_ops.zeeman_operator(&mut zeeman_prop, &["s1", "i1", "s2", "i2"], &gammas);
        let zeeman_prop_matrix = matrix_builder.from_operator(&zeeman_prop);

        let mut hifi = Operator::new(states.clone());
        hifi.add_operator(vec!["s1", "i1"], |brakets| {
            a_hifi_1 * spin_ops.dot([&brakets["s1"], &brakets["i1"]])
        })
        .add_operator(vec!["s2", "i2"], |brakets| {
            a_hifi_2 * spin_ops.dot([&brakets["s2"], &brakets["i2"]])
        });
        let hifi_matrix = matrix_builder.from_operator(&hifi);

        let mag_fields = linspace(0.0, 600.0, 1000);
        let mut values = Vec::new();
        let mut diagonalization = Diagonalization::new(hifi_matrix.clone(), None);
        for mag_field in &mag_fields {
            let hamiltonian = &hifi_matrix + &zeeman_prop_matrix * *mag_field;
            diagonalization.diagonalize(hamiltonian);

            let energies = convert_data_units(diagonalization.eigvalues(), |x| {
                Energy::new(x, Au).to(GHz).value
            });
            values.push(energies);
        }

        save_param_change(
            "potassium/hifi_39_41",
            mag_fields,
            values,
            vec!["magnetic_field", "energies"],
        )
        .unwrap()
    }

    fn add_bound_states(states: &mut CompositeState) {
        let bound_space = Box::new(IrreducibleState::new(
            "b",
            "bound",
            2,
            vec![BasisElement::new(0.0, 0), BasisElement::new(1.0, 1)],
        ));

        states.add_state(bound_space);
    }

    pub fn abm_39_41() {
        println!("Solving abm for K39 K41...");

        let a_hifi_1 = Energy::new(Self::HIFI_K39_MHZ, MHz).to_au();
        let a_hifi_2 = Energy::new(Self::HIFI_K41_MHZ, MHz).to_au();
        let gamma_e = -2.0 * Defaults::BOHR_MAG;
        let gamma_i1 = Energy::new(Self::GAMMA_I_K39, MHz).to_au() * 1e-4;
        let gamma_i2 = Energy::new(Self::GAMMA_I_K41, MHz).to_au() * 1e-4;

        let triplet_states = vec![Energy::new(-8.33, MHz), Energy::new(-1282.5, MHz)];
        let singlet_states = vec![Energy::new(-32.1, MHz), Energy::new(-1698.1, MHz)];
        let fc_factors = vec![0.9180, 0.0895, 0.0463, 0.9674];

        let total_double_m = 4;
        let states = Rc::new({
            let s_tot = Box::new(DirectSumState::new(
                "S_tot",
                vec![create_spin("triplet", 2), create_spin("singlet", 0)],
            ));
            let i1 = Box::new(create_spin("i1", 3));
            let i2 = Box::new(create_spin("i2", 3));
            let mut spin_states = CompositeState::new("spins_states", vec![s_tot, i1, i2]);

            Self::add_bound_states(&mut spin_states);
            spin_states
        });

        let mut matrix_builder = MatrixBuilder::new(states.clone());
        matrix_builder.filter_states(|q_states| {
            let bound_index = states.state_index("b");
            q_states.iter().map(|s| s.q_number()).sum::<isize>()
                == total_double_m + q_states[bound_index].q_number()
        });
        let spin_ops = SpinOperators::new(states.clone());
        let bound_ops = BoundOperators::new(
            states.clone(),
            "b",
            triplet_states,
            singlet_states,
            fc_factors,
        );

        let mut bound_states = Operator::new(states.clone());
        bound_states.add_operator(vec!["S_tot", "b"], |brakets| {
            if brakets["S_tot"].ket == brakets["S_tot"].bra {
                let state_q_number = states.state_from_base(brakets["S_tot"].ket).q_number();
                if state_q_number == 0 {
                    return bound_ops.bound_state(StateType::Singlet, &brakets["b"]);
                } else {
                    assert!(state_q_number == 2);
                    return bound_ops.bound_state(StateType::Triplet, &brakets["b"]);
                }
            } else {
                0.0
            }
        });
        let bound_states_matrix = matrix_builder.from_operator(&bound_states);

        let gammas = [gamma_e, gamma_i1, gamma_i2];
        let mut zeeman_prop = Operator::new(states.clone());
        spin_ops.zeeman_operator(&mut zeeman_prop, &["S_tot", "i1", "i2"], &gammas);
        let zeeman_prop_matrix = matrix_builder.from_operator(&zeeman_prop);

        let mut hifi = Operator::new(states.clone());
        hifi.add_operator(vec!["S_tot", "i1", "i2"], |brakets| {
            a_hifi_1 / 2.0 * spin_ops.dot([&brakets["S_tot"], &brakets["i1"]])
                + a_hifi_2 / 2.0 * spin_ops.dot([&brakets["S_tot"], &brakets["i2"]])
        })
        .add_operator(vec!["S_tot", "i1", "i2", "b"], |brakets| {
            let hifi_minus = a_hifi_1 / 2.0 * spin_ops.dot_minus([&brakets["S_tot"], &brakets["i1"]])
                - a_hifi_2 / 2.0 * spin_ops.dot_minus([&brakets["S_tot"], &brakets["i2"]]);

            let state_q_number_ket = states.state_from_base(brakets["S_tot"].ket).q_number();
            let state_q_number_bra = states.state_from_base(brakets["S_tot"].bra).q_number();

            let gordon_factors = if state_q_number_ket == 0 && state_q_number_bra == 2 {
                bound_ops.state_type_flip(StateType::Singlet, &brakets["b"])
            } else if state_q_number_ket == 2 && state_q_number_bra == 0 {
                bound_ops.state_type_flip(StateType::Triplet, &brakets["b"])
            } else {
                0.0
            };

            hifi_minus * gordon_factors
        });
        let hifi_matrix = matrix_builder.from_operator(&hifi);

        let mag_fields = linspace(0.0, 600.0, 1000);
        let mut values = Vec::new();
        let mut solver = Diagonalization::new(hifi_matrix.clone(), None);
        for mag_field in &mag_fields {
            let hamiltonian = &bound_states_matrix + &hifi_matrix + &zeeman_prop_matrix * *mag_field;

            solver.diagonalize(hamiltonian);

            let energies = convert_data_units(solver.eigvalues(), |x| Energy::new(x, Au).to(GHz).value);
            values.push(energies);
        }

        save_param_change(
            "potassium/abm_39_41",
            mag_fields,
            values,
            vec!["magnetic_field", "energies"],
        )
        .unwrap()
    }

    pub fn hifi_39_39() {
        println!("Solving hyperfine for K39 - K39...");

        let a_hifi_1 = Energy::new(Self::HIFI_K39_MHZ, MHz).to_au();
        let a_hifi_2 = Energy::new(Self::HIFI_K39_MHZ, MHz).to_au();
        let gamma_e = -2.0 * Defaults::BOHR_MAG;
        let gamma_i1 = Energy::new(Self::GAMMA_I_K39, MHz).to_au() * 1e-4;
        let gamma_i2 = Energy::new(Self::GAMMA_I_K39, MHz).to_au() * 1e-4;

        let total_double_m = 4;

        let i1 = Box::new(create_spin("i1", 3));
        let i2 = Box::new(create_spin("i2", 3));
        let s_tot = Box::new(DirectSumState::new(
            "S_tot",
            vec![create_spin("triplet", 2), create_spin("singlet", 0)],
        ));

        let states_s_tot = Rc::new(
            CompositeState::new("states_s_tot", vec![s_tot.clone(), i1.clone(), i2.clone()])
        );
        let states_si_tot = Rc::new({
            let i_tot = Box::new(create_spin_sum("I_tot", &i1, &i2));
            CompositeState::new("states_si_tot", vec![s_tot, i_tot])
        });

        let mut matrix_builder_si = MatrixBuilder::new(states_si_tot.clone());
        matrix_builder_si.filter_states(|q_states| {
            q_states.iter().map(|s| s.q_number()).sum::<isize>() == total_double_m
        });
        matrix_builder_si.filter_states(|q_states| {
            let state_spin_e = states_si_tot.state_from_base(&q_states[0]);
            let state_spin_n = states_si_tot.state_from_base(&q_states[1]);

            (state_spin_e.name() == "singlet" && state_spin_n.name().contains("anti-symmetric"))
                || (state_spin_e.name() == "triplet"
                    && !state_spin_n.name().contains("anti-symmetric"))
        });

        let mut matrix_builder_s = MatrixBuilder::new(states_s_tot.clone());
        matrix_builder_s.filter_states(|q_states| {
            q_states.iter().map(|s| s.q_number()).sum::<isize>() == total_double_m
        });

        let transformation =
            StateTransformation::new(&matrix_builder_s, &matrix_builder_si, |from, to| {
                let d_spin1 = states_s_tot.state_from_base(&from[1]).q_number() as usize;
                let d_spin2 = states_s_tot.state_from_base(&from[2]).q_number() as usize;
                let tot_dspin = states_si_tot.state_from_base(&to[1]).q_number() as usize;

                if from[0].same_as(&states_s_tot, &to[0], &states_si_tot) {
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

        let spin_ops_s = SpinOperators::new(states_s_tot.clone());

        let gammas = [gamma_e, gamma_i1, gamma_i2];
        let mut zeeman_prop = Operator::new(states_s_tot.clone());
        spin_ops_s.zeeman_operator(&mut zeeman_prop, &["S_tot", "i1", "i2"], &gammas);
        let zeeman_prop_matrix = &transformation_matrix
            * matrix_builder_s.from_operator(&zeeman_prop)
            * transformation_matrix.transpose();

        let mut hifi = Operator::new(states_s_tot.clone());
        hifi.add_operator(vec!["S_tot", "i1", "i2"], |brakets| {
            a_hifi_1 / 2.0 * spin_ops_s.dot([&brakets["S_tot"], &brakets["i1"]])
                + a_hifi_2 / 2.0 * spin_ops_s.dot([&brakets["S_tot"], &brakets["i2"]])
        })
        .add_operator(vec!["S_tot", "i1", "i2"], |brakets| {
            a_hifi_1 / 2.0 * spin_ops_s.dot_minus([&brakets["S_tot"], &brakets["i1"]])
                - a_hifi_2 / 2.0 * spin_ops_s.dot_minus([&brakets["S_tot"], &brakets["i2"]])
        });
        let hifi_matrix = &transformation_matrix
            * matrix_builder_s.from_operator(&hifi)
            * transformation_matrix.transpose();

        let mag_fields = linspace(0.0, 600.0, 1000);
        let mut values = Vec::new();
        let mut diagonalization = Diagonalization::new(hifi_matrix.clone(), None);
        for mag_field in &mag_fields {
            let hamiltonian = &hifi_matrix + &zeeman_prop_matrix * *mag_field;

            diagonalization.diagonalize(hamiltonian);

            let energies = convert_data_units(diagonalization.eigvalues(), |x| {
                Energy::new(x, Au).to(GHz).value
            });
            values.push(energies);
        }

        save_param_change(
            "potassium/hifi_39_39",
            mag_fields,
            values,
            vec!["magnetic_field", "energies"],
        )
        .unwrap()
    }

    pub fn abm_39_39() {
        println!("Solving abm for K39 K39...");

        let a_hifi_1 = Energy::new(Self::HIFI_K39_MHZ, MHz).to_au();
        let a_hifi_2 = Energy::new(Self::HIFI_K39_MHZ, MHz).to_au();
        let gamma_e = -2.0 * Defaults::BOHR_MAG;
        let gamma_i1 = Energy::new(Self::GAMMA_I_K39, MHz).to_au() * 1e-4;
        let gamma_i2 = Energy::new(Self::GAMMA_I_K39, MHz).to_au() * 1e-4;

        let triplet_states = vec![Energy::new(-8.33, MHz), Energy::new(-1282.5, MHz)];
        let singlet_states = vec![Energy::new(-32.1, MHz), Energy::new(-1698.1, MHz)];
        let fc_factors = vec![0.9180, 0.0895, 0.0463, 0.9674];

        let total_double_m = 4;

        let i1 = Box::new(create_spin("i1", 3));
        let i2 = Box::new(create_spin("i2", 3));
        let s_tot = Box::new(DirectSumState::new(
            "S_tot",
            vec![create_spin("triplet", 2), create_spin("singlet", 0)],
        ));

        let states_s_tot = Rc::new({
            let mut spin_states = CompositeState::new("states_s_tot", vec![s_tot.clone(), i1.clone(), i2.clone()]);

            Self::add_bound_states(&mut spin_states);
            spin_states
        });

        let states_si_tot = Rc::new({
            let i_tot = Box::new(create_spin_sum("I_tot", &i1, &i2));
            let mut spin_states = CompositeState::new("states_si_tot", vec![s_tot, i_tot]);

            Self::add_bound_states(&mut spin_states);
            spin_states
        });

        let mut matrix_builder_si = MatrixBuilder::new(states_si_tot.clone());
        matrix_builder_si.filter_states(|q_states| {
            let bound_index = states_si_tot.state_index("b");
            q_states.iter().map(|s| s.q_number()).sum::<isize>()
                == total_double_m + q_states[bound_index].q_number()
        });
        matrix_builder_si.filter_states(|q_states| {
            let state_spin_e = states_si_tot.state_from_base(&q_states[0]);
            let state_spin_n = states_si_tot.state_from_base(&q_states[1]);

            (state_spin_e.name() == "singlet" && state_spin_n.name().contains("anti-symmetric"))
                || (state_spin_e.name() == "triplet"
                    && !state_spin_n.name().contains("anti-symmetric"))
        });

        let mut matrix_builder_s = MatrixBuilder::new(states_s_tot.clone());
        matrix_builder_s.filter_states(|q_states| {
            let bound_index = states_s_tot.state_index("b");
            q_states.iter().map(|s| s.q_number()).sum::<isize>()
                == total_double_m + q_states[bound_index].q_number()
        });

        let transformation =
            StateTransformation::new(&matrix_builder_s, &matrix_builder_si, |from, to| {
                let d_spin1 = states_s_tot.state_from_base(&from[1]).q_number() as usize;
                let d_spin2 = states_s_tot.state_from_base(&from[2]).q_number() as usize;
                let tot_dspin = states_si_tot.state_from_base(&to[1]).q_number() as usize;
                if &from[3].q_number() != &to[2].q_number() {
                    return 0.0;
                }

                if from[0].same_as(&states_s_tot, &to[0], &states_si_tot) {
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

        let spin_ops_s = SpinOperators::new(states_s_tot.clone());
        let bound_ops_si = BoundOperators::new(
            states_si_tot.clone(),
            "b",
            triplet_states.clone(),
            singlet_states.clone(),
            fc_factors.clone(),
        );
        let bound_ops_s = BoundOperators::new(
            states_s_tot.clone(),
            "b",
            triplet_states,
            singlet_states,
            fc_factors,
        );

        let mut bound_states = Operator::new(states_si_tot.clone());
        bound_states.add_operator(vec!["S_tot", "b"], |brakets| {
            if brakets["S_tot"].ket == brakets["S_tot"].bra {
                let state_q_number = states_si_tot.state_from_base(brakets["S_tot"].ket).q_number();
                if state_q_number == 0 {
                    return bound_ops_si.bound_state(StateType::Singlet, &brakets["b"]);
                } else {
                    assert!(state_q_number == 2);
                    return bound_ops_si.bound_state(StateType::Triplet, &brakets["b"]);
                }
            } else {
                0.0
            }
        });
        let bound_states_matrix = matrix_builder_si.from_operator(&bound_states);

        let gammas = [gamma_e, gamma_i1, gamma_i2];
        let mut zeeman_prop = Operator::new(states_s_tot.clone());
        spin_ops_s.zeeman_operator(&mut zeeman_prop, &["S_tot", "i1", "i2"], &gammas);
        let zeeman_prop_matrix = &transformation_matrix
            * matrix_builder_s.from_operator(&zeeman_prop)
            * transformation_matrix.transpose();

        let mut hifi = Operator::new(states_s_tot.clone());
        hifi.add_operator(vec!["S_tot", "i1", "i2"], |brakets| {
            a_hifi_1 / 2.0 * spin_ops_s.dot([&brakets["S_tot"], &brakets["i1"]])
                + a_hifi_2 / 2.0 * spin_ops_s.dot([&brakets["S_tot"], &brakets["i2"]])
        })
        .add_operator(vec!["S_tot", "i1", "i2", "b"], |brakets| {
            let hifi_minus = a_hifi_1 / 2.0 * spin_ops_s.dot_minus([&brakets["S_tot"], &brakets["i1"]])
                - a_hifi_2 / 2.0 * spin_ops_s.dot_minus([&brakets["S_tot"], &brakets["i2"]]);

            let state_q_number_ket = states_s_tot.state_from_base(brakets["S_tot"].ket).q_number();
            let state_q_number_bra = states_s_tot.state_from_base(brakets["S_tot"].bra).q_number();

            let gordon_factors = if state_q_number_ket == 0 && state_q_number_bra == 2 {
                bound_ops_s.state_type_flip(StateType::Singlet, &brakets["b"])
            } else if state_q_number_ket == 2 && state_q_number_bra == 0 {
                bound_ops_s.state_type_flip(StateType::Triplet, &brakets["b"])
            } else {
                0.0
            };

            hifi_minus * gordon_factors
        });
        let hifi_matrix = &transformation_matrix
            * matrix_builder_s.from_operator(&hifi)
            * transformation_matrix.transpose();

        let mag_fields = linspace(0.0, 600.0, 1000);
        let mut values = Vec::new();
        let mut diagonalization = Diagonalization::new(hifi_matrix.clone(), None);
        for mag_field in &mag_fields {
            let hamiltonian = &bound_states_matrix + &hifi_matrix + &zeeman_prop_matrix * *mag_field;

            diagonalization.diagonalize(hamiltonian);

            let energies = convert_data_units(diagonalization.eigvalues(), |x| {
                Energy::new(x, Au).to(GHz).value
            });
            values.push(energies);
        }

        save_param_change(
            "potassium/abm_39_39",
            mag_fields,
            values,
            vec!["magnetic_field", "energies"],
        )
        .unwrap()
    }
}
