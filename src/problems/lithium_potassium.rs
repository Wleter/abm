use std::{collections::VecDeque, rc::Rc};

use abm::{
    composite_state::CompositeState,
    diagonalization::Diagonalization,
    direct_sum_state::DirectSumState,
    matrix_builder::MatrixBuilder,
    operator::Operator,
    spin_operators::SpinOperators,
    state_factory::create_spin,
};
use quantum::{
    problem_selector::ProblemSelector,
    saving::save_param_change,
    units::{convert_data_units, energy_units::EnergyUnit},
    utility::linspace,
};

use crate::defaults::Defaults;

pub struct LithiumPotassium;

impl ProblemSelector for LithiumPotassium {
    const NAME: &'static str = "potassium-potassium";

    fn list() -> Vec<&'static str> {
        vec![
            "hyperfine",
            "asymptotic bound",
        ]
    }

    fn methods(number: &str, _args: &mut VecDeque<String>) {
        match number {
            "0" => Self::hifi(),
            "1" => Self::abm(),
            _ => println!("Invalid problem number"),
        }
    }
}

impl LithiumPotassium {
    const HIFI_LI6_MHZ: f64 = 228.2 / 1.5;
    const HIFI_K40_MHZ: f64 = -1285.8 / 4.5;

    pub fn hifi() {
        println!("Solving hyperfine for Li6 - K40...");

        // ---------- Li6 ----------
        let a_hifi = EnergyUnit::MHz.to_au(Self::HIFI_LI6_MHZ);
        let gamma_e = -2.0 * Defaults::BOHR_MAG / 2.350517567e9;
        let gamma_i = 0.0;

        let s = Box::new(create_spin("s", 1));
        let i = Box::new(create_spin("i", 2));

        let states = Rc::new(CompositeState::new("spins_states", vec![s, i]));
        let matrix_builder = MatrixBuilder::new(states.clone());
        let spin_ops = SpinOperators::new(states.clone());

        let gammas = [gamma_e, gamma_i];
        let mut zeeman_prop = Operator::new(states.clone());
        spin_ops.zeeman_operator(&mut zeeman_prop, &["s", "i"], &gammas);
        let zeeman_prop_matrix = matrix_builder.from_operator(&zeeman_prop);

        let mut hifi = Operator::new(states.clone());
        hifi.add_operator(vec!["s", "i"], |brakets| {
            a_hifi * spin_ops.dot([&brakets["s"], &brakets["i"]])
        });
        let hifi_matrix = matrix_builder.from_operator(&hifi);

        let mag_fields = linspace(0.0, 1000.0, 1000);
        let mut values = Vec::new();
        let mut solver = Diagonalization::new(hifi_matrix.clone(), None);
        for mag_field in &mag_fields {
            let hamiltonian = &hifi_matrix + &zeeman_prop_matrix * *mag_field;
            solver.diagonalize(hamiltonian);

            let energies = convert_data_units(solver.eigvalues(), |x| EnergyUnit::Au.to_giga_hz(x));
            values.push(energies);
        }

        save_param_change(
            "Li_K/hifi_Li6",
            mag_fields,
            values,
            vec!["magnetic_field", "energies"],
        )
        .unwrap();

        // ---------- K40 ----------
        let a_hifi = EnergyUnit::MHz.to_au(Self::HIFI_K40_MHZ);
        let gamma_e = -2.0 * Defaults::BOHR_MAG / 2.350517567e9;
        let gamma_i = 0.0;

        let s = Box::new(create_spin("s", 1));
        let i = Box::new(create_spin("i", 8));
        let states = Rc::new(CompositeState::new("spins_states", vec![s, i]));

        let matrix_builder = MatrixBuilder::new(states.clone());
        let spin_ops = SpinOperators::new(states.clone());

        let gammas = [gamma_e, gamma_i];
        let mut zeeman_prop = Operator::new(states.clone());
        spin_ops.zeeman_operator(&mut zeeman_prop, &["s", "i"], &gammas);
        let zeeman_prop_matrix = matrix_builder.from_operator(&zeeman_prop);

        let mut hifi = Operator::new(states.clone());
        hifi.add_operator(vec!["s", "i"], |brakets| {
            a_hifi * spin_ops.dot([&brakets["s"], &brakets["i"]])
        });
        let hifi_matrix = matrix_builder.from_operator(&hifi);

        let mag_fields = linspace(0.0, 1000.0, 1000);
        let mut values = Vec::new();
        let mut solver = Diagonalization::new(hifi_matrix.clone(), None);
        for mag_field in &mag_fields {
            let hamiltonian = &hifi_matrix + &zeeman_prop_matrix * *mag_field;
            solver.diagonalize(hamiltonian);

            let energies = convert_data_units(solver.eigvalues(), |x| EnergyUnit::Au.to_giga_hz(x));
            values.push(energies);
        }

        save_param_change(
            "Li_K/hifi_K40",
            mag_fields,
            values,
            vec!["magnetic_field", "energies"],
        )
        .unwrap();
    }

    pub fn abm() {
        println!("Solving abm for K39 K41...");

        let a_hifi_1 = EnergyUnit::MHz.to_au(Self::HIFI_LI6_MHZ);
        let a_hifi_2 = EnergyUnit::MHz.to_au(Self::HIFI_K40_MHZ);
        let gamma_e = -2.0 * Defaults::BOHR_MAG / 2.350517567e9;
        let gamma_i1 = 0.0;
        let gamma_i2 = 0.0;

        let triplet_state = EnergyUnit::MHz.to_au(-427.44);
        let singlet_state = EnergyUnit::MHz.to_au(-720.76);
        let fc_factor = 0.979;

        let total_double_m = -6;
        let states = Rc::new({
            let s_tot = Box::new(DirectSumState::new(
                "S_tot",
                vec![create_spin("triplet", 2), create_spin("singlet", 0)],
            ));
            let i1 = Box::new(create_spin("i1", 2));
            let i2 = Box::new(create_spin("i2", 8));
            CompositeState::new("spins_states", vec![s_tot, i1, i2])
        });

        let mut matrix_builder = MatrixBuilder::new(states.clone());
        matrix_builder.filter_states(|q_states| {
            q_states.iter().map(|s| s.q_number()).sum::<isize>() == total_double_m
        });
        let spin_ops = SpinOperators::new(states.clone());

        let mut bound_states = Operator::new(states.clone());
        bound_states.add_operator(vec!["S_tot"], |brakets| {
            if brakets["S_tot"].ket == brakets["S_tot"].bra {
                let state_q_number = states.state_from_base(brakets["S_tot"].ket).q_number();
                if state_q_number == 0 {
                    singlet_state
                } else {
                    assert!(state_q_number == 2);
                    triplet_state
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
        .add_operator(vec!["S_tot", "i1", "i2"], |brakets| {
            fc_factor * (
                a_hifi_1 / 2.0 * spin_ops.dot_minus([&brakets["S_tot"], &brakets["i1"]])
                    - a_hifi_2 / 2.0 * spin_ops.dot_minus([&brakets["S_tot"], &brakets["i2"]])
            )
        });
        let hifi_matrix = matrix_builder.from_operator(&hifi);

        let mag_fields = linspace(0.0, 400.0, 1000);
        let mut values = Vec::new();
        let mut solver = Diagonalization::new(hifi_matrix.clone(), None);
        for mag_field in &mag_fields {
            let hamiltonian = &bound_states_matrix + &hifi_matrix + &zeeman_prop_matrix * *mag_field;

            solver.diagonalize(hamiltonian);

            let energies = convert_data_units(solver.eigvalues(), |x| EnergyUnit::Au.to_giga_hz(x));
            values.push(energies);
        }

        save_param_change(
            "Li_K/abm",
            mag_fields,
            values,
            vec!["magnetic_field", "energies"],
        )
        .unwrap()
    }
}
