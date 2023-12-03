pub mod bound_operators;
pub mod braket;
pub mod composite_state;
pub mod diagonalization;
pub mod direct_sum_state;
pub mod irreducible_state;
pub mod matrix_builder;
pub mod operator;
pub mod spin_operators;
pub mod state;
pub mod basis_element;
pub mod state_factory;
pub mod state_transformation;
pub mod utils;
pub mod defaults;

extern crate nalgebra;
extern crate quantum;

#[cfg(test)]
mod tests {
    use std::rc::Rc;

    use nalgebra::DMatrix;

    use crate::{
        composite_state::CompositeState,
        direct_sum_state::DirectSumState,
        matrix_builder::MatrixBuilder,
        operator::Operator,
        spin_operators::*,
        state::State,
        state_factory::create_spin,
        state_transformation::StateTransformation,
        utils::{assert_almost_eq, clebsch_gordan},
    };

    fn sample_composite() -> CompositeState {
        let state1 = Box::new(create_spin("spin_1", 1));
        let state2 = Box::new(create_spin("spin_2", 2));
        let state3 = Box::new(create_spin("spin_3", 3));

        CompositeState::new("test_state", vec![state1, state2, state3])
    }

    #[test]
    fn test_composite() {
        let mut composite_state = sample_composite();

        assert_eq!(composite_state.dim(), 24);
        assert_eq!(
            composite_state
                .states_from_index(0)
                .iter()
                .map(|b| b.q_number())
                .collect::<Vec<isize>>(),
            vec![1, 2, 3]
        );
        assert_eq!(
            composite_state
                .states_from_index(1)
                .iter()
                .map(|b| b.q_number())
                .collect::<Vec<isize>>(),
            vec![-1, 2, 3]
        );
        assert_eq!(
            composite_state
                .states_from_index(7)
                .iter()
                .map(|b| b.q_number())
                .collect::<Vec<isize>>(),
            vec![-1, 2, 1]
        );

        let triplet = create_spin("triplet", 2);
        let singlet = create_spin("singlet", 0);
        let sum_state = Box::new(DirectSumState::new("spin_sum", vec![triplet, singlet]));
        composite_state.add_state(sum_state);

        let states = composite_state.states_from_index(2);
        let triplet_base_state = states.last().unwrap();
        assert_eq!(
            composite_state.state_from_base(&triplet_base_state).name(),
            "triplet"
        );

        let states = composite_state.states_from_index(3 * 24 + 2);
        let singlet_base = states.last().unwrap();
        assert_eq!(
            composite_state.state_from_base(&singlet_base).name(),
            "singlet"
        );
    }

    #[test]
    fn test_matrix_builder() {
        let composite_state = Rc::new(sample_composite());
        let mut matrix_builder = MatrixBuilder::new(composite_state);

        assert_eq!(matrix_builder.dim(), 24);

        matrix_builder
            .filter_states(|q_states| q_states.iter().map(|s| s.q_number()).sum::<isize>() == -6);

        assert_eq!(matrix_builder.dim(), 1);

        matrix_builder.unfilter_states();
        matrix_builder
            .filter_states(|q_states| q_states.iter().map(|s| s.q_number()).sum::<isize>() == 0);
        assert_eq!(matrix_builder.dim(), 6);
    }

    #[test]
    fn test_spin_operators() {
        let state_1 = Box::new(create_spin("s1", 1));
        let composite_state = Rc::new(CompositeState::new("test_state", vec![state_1]));
        let matrix_builder = MatrixBuilder::new(composite_state.clone());
        let spin = SpinOperators::new(composite_state.clone());

        let mut proj_z = Operator::new(composite_state.clone());
        proj_z.add_operator(vec!["s1"], |brakets| spin.proj_z(&brakets["s1"]));
        let proj_z_matrix = matrix_builder.from_operator(&proj_z);

        assert_eq!(
            proj_z_matrix,
            DMatrix::from_vec(2, 2, vec![0.5, 0.0, 0.0, -0.5])
        );

        let mut ladder_plus = Operator::new(composite_state.clone());
        ladder_plus.add_operator(vec!["s1"], |brakets| spin.ladder_plus(&brakets["s1"]));
        let ladder_plus_matrix = matrix_builder.from_operator(&ladder_plus);
        assert_eq!(
            ladder_plus_matrix,
            DMatrix::from_vec(2, 2, vec![0.0, 0.0, 1.0, 0.0])
        );

        let mut ladder_minus = Operator::new(composite_state.clone());
        ladder_minus.add_operator(vec!["s1"], |brakets| spin.ladder_minus(&brakets["s1"]));
        let ladder_minus_matrix = matrix_builder.from_operator(&ladder_minus);
        assert_eq!(
            ladder_minus_matrix,
            DMatrix::from_vec(2, 2, vec![0.0, 1.0, 0.0, 0.0])
        );

        let mut summed = Operator::new(composite_state.clone());
        summed.add_operator(vec!["s1"], |brakets| spin.proj_z(&brakets["s1"]));
        summed.add_operator(vec!["s1"], |brakets| spin.ladder_plus(&brakets["s1"]));
        summed.add_operator(vec!["s1"], |brakets| spin.ladder_minus(&brakets["s1"]));
        let summed_matrix = matrix_builder.from_operator(&summed);
        assert_eq!(
            summed_matrix,
            DMatrix::from_vec(2, 2, vec![0.5, 1.0, 1.0, -0.5])
        );
    }

    #[test]
    fn test_dot_minus() {
        let states = Rc::new({
            let s_tot = Box::new(DirectSumState::new(
                "S_tot",
                vec![create_spin("triplet", 2), create_spin("singlet", 0)],
            ));

            CompositeState::new("spin", vec![s_tot])
        });

        let matrix_builder = MatrixBuilder::new(states.clone());
        println!("{}", matrix_builder);
        let spin_ops = SpinOperators::new(states.clone());
        let mut op = Operator::new(states.clone());
        op.add_operator(vec!["S_tot"], |brakets| spin_ops.proj_z_minus(&brakets["S_tot"]));
        let matrix = matrix_builder.from_operator(&op);
        println!("{}", matrix);
        let mut op = Operator::new(states.clone());
        op.add_operator(vec!["S_tot"], |brakets| {
            spin_ops.ladder_minus_minus(&brakets["S_tot"])
        });
        let matrix = matrix_builder.from_operator(&op);
        println!("{}", matrix);
        let mut op = Operator::new(states.clone());
        op.add_operator(vec!["S_tot"], |brakets| {
            spin_ops.ladder_plus_minus(&brakets["S_tot"])
        });
        let matrix = matrix_builder.from_operator(&op);
        println!("{}", matrix);
    }

    #[test]
    fn test_hifi() {
        let states = Rc::new({
            let s_tot = Box::new(DirectSumState::new(
                "S_tot",
                vec![create_spin("triplet", 2), create_spin("singlet", 0)],
            ));
            let i1 = Box::new(create_spin("i1", 1));
            let i2 = Box::new(create_spin("i2", 1));

            CompositeState::new("spin", vec![s_tot, i1, i2])
        });

        let mut matrix_builder = MatrixBuilder::new(states.clone());
        matrix_builder
            .filter_states(|q_states| q_states.iter().map(|s| s.q_number()).sum::<isize>() == 0);

        let spin_ops = SpinOperators::new(states.clone());
        let mut hifi = Operator::new(states.clone());
        hifi.add_operator(vec!["S_tot", "i1", "i2"], |brakets| {
            let value = 0.5 * spin_ops.dot([&brakets["S_tot"], &brakets["i1"]])
                + 0.5 * spin_ops.dot([&brakets["S_tot"], &brakets["i2"]]);

            value
        })
        .add_operator(vec!["S_tot", "i1", "i2"], |brakets| {
            let value = 0.5 * spin_ops.dot_minus([&brakets["S_tot"], &brakets["i1"]])
                - 0.5 * spin_ops.dot_minus([&brakets["S_tot"], &brakets["i2"]]);

            value
        });
        let matrix = matrix_builder.from_operator(&hifi);
        println!("{:.3e}", matrix);

        let states_ssi = Rc::new({
            let s1 = Box::new(create_spin("s1", 1));
            let s2 = Box::new(create_spin("s2", 1));
            let i1 = Box::new(create_spin("i1", 1));
            let i2 = Box::new(create_spin("i2", 1));

            CompositeState::new("spin", vec![s1, s2, i1, i2])
        });
        let mut matrix_builder_ssi = MatrixBuilder::new(states_ssi.clone());
        matrix_builder_ssi
            .filter_states(|q_states| q_states.iter().map(|s| s.q_number()).sum::<isize>() == 0);
        let transf = StateTransformation::new(&matrix_builder_ssi, &matrix_builder, |from, to| {
            let d_spin1 = states_ssi.state_from_base(&from[0]).q_number() as usize;
            let d_spin2 = states_ssi.state_from_base(&from[1]).q_number() as usize;
            let tot_dspin = states.state_from_base(&to[0]).q_number() as usize;

            if from[2].same_as(&states_ssi, &to[1], &states)
                && from[3].same_as(&states_ssi, &to[2], &states)
            {
                clebsch_gordan(
                    d_spin1,
                    from[0].q_number(),
                    d_spin2,
                    from[1].q_number(),
                    tot_dspin,
                    to[0].q_number(),
                )
            } else {
                0.0
            }
        });
        println!("{}", matrix_builder);
        println!("{}", matrix_builder_ssi);
        let transf_matrix = transf.matrix_transformation();
        println!("transformation matrix");
        println!("{:.3e}", transf_matrix);
        println!("transformed matrix");
        println!(
            "{:.3e}",
            transf_matrix.transpose() * matrix * &transf_matrix
        );

        let spin_ops = SpinOperators::new(states_ssi.clone());

        let mut hifi = Operator::new(states_ssi.clone());
        hifi.add_operator(vec!["s1", "i1"], |brakets| {
            spin_ops.dot([&brakets["s1"], &brakets["i1"]])
        })
        .add_operator(vec!["s2", "i2"], |brakets| {
            spin_ops.dot([&brakets["s2"], &brakets["i2"]])
        });
        let hifi_matrix = matrix_builder_ssi.from_operator(&hifi);
        println!("{:.3e}", hifi_matrix);

        let hifi_matrix = &transf_matrix * hifi_matrix * transf_matrix.transpose();
        println!("{:.3e}", hifi_matrix);
    }

    #[test]
    fn test_clebsch_gordan() {
        let result = clebsch_gordan(1, 1, 1, 1, 2, 2);
        assert_almost_eq(result, 1.0, 1e-8);

        let result = clebsch_gordan(1, 1, 1, -1, 2, 0);
        assert_almost_eq(result, 1.0 / 2f64.sqrt(), 1e-8);

        let result = clebsch_gordan(1, -1, 1, 1, 2, 0);
        assert_almost_eq(result, 1.0 / 2f64.sqrt(), 1e-8);

        let result = clebsch_gordan(1, 1, 1, -1, 0, 0);
        assert_almost_eq(result, 1.0 / 2f64.sqrt(), 1e-8);

        let result = clebsch_gordan(1, -1, 1, 1, 0, 0);
        assert_almost_eq(result, -1.0 / 2f64.sqrt(), 1e-8);

        let result = clebsch_gordan(2, -2, 4, 2, 6, 0);
        assert_almost_eq(result, 1.0 / 5f64.sqrt(), 1e-8);
    }
}
