use std::{
    iter::{Rev, StepBy},
    ops::RangeInclusive,
    rc::Rc,
};

use crate::{
    braket::Braket, composite_state::CompositeState, operator::Operator, state_base::StateBase,
    utils::clebsch_gordan,
};

#[derive(Clone)]
pub struct SpinOperators {
    composite_state: Rc<CompositeState>,
}

impl SpinOperators {
    pub fn new(composite_state: Rc<CompositeState>) -> Self {
        Self { composite_state }
    }

    pub fn id(&self, braket: &Braket) -> f64 {
        assert!(self.is_spin_state(braket.bra));
        assert!(self.is_spin_state(braket.ket));
        if !self.same_spin(braket) {
            return 0.0;
        }

        (braket.bra.q_number() == braket.ket.q_number()) as isize as f64
    }

    pub fn proj_z(&self, braket: &Braket) -> f64 {
        assert!(self.is_spin_state(braket.bra));
        assert!(self.is_spin_state(braket.ket));
        if !self.same_spin(braket) {
            return 0.0;
        }

        if braket.ket.q_number() == braket.bra.q_number() {
            braket.ket.value()
        } else {
            0.0
        }
    }

    pub fn ladder_plus(&self, braket: &Braket) -> f64 {
        assert!(self.is_spin_state(braket.bra));
        assert!(self.is_spin_state(braket.ket));
        if !self.same_spin(braket) {
            return 0.0;
        }

        if braket.ket.q_number() + 2 == braket.bra.q_number() {
            let spin = self.composite_state.state_from_base(&braket.ket).q_number() as f64 / 2.0;
            let m = braket.ket.q_number() as f64 / 2.0;
            (spin * (spin + 1.0) - m * (m + 1.0)).sqrt()
        } else {
            0.0
        }
    }

    pub fn ladder_minus(&self, braket: &Braket) -> f64 {
        assert!(self.is_spin_state(braket.bra));
        assert!(self.is_spin_state(braket.ket));
        if !self.same_spin(braket) {
            return 0.0;
        }

        if braket.ket.q_number() - 2 == braket.bra.q_number() {
            let spin = self.composite_state.state_from_base(&braket.ket).q_number() as f64 / 2.0;
            let m = braket.ket.q_number() as f64 / 2.0;
            (spin * (spin + 1.0) - m * (m - 1.0)).sqrt()
        } else {
            0.0
        }
    }

    pub fn spin_sqr(&self, base: &StateBase) -> f64 {
        assert!(self.is_spin_state(base));
        let spin = self.composite_state.state_from_base(&base).q_number() as f64 / 2.0;

        spin * (spin + 1.0)
    }

    pub fn dot(&self, brakets: [&Braket; 2]) -> f64 {
        let val1 = self.proj_z(brakets[0]) * self.proj_z(brakets[1]);
        let val2 = 0.5 * self.ladder_plus(brakets[0]) * self.ladder_minus(brakets[1]);
        let val3 = 0.5 * self.ladder_minus(brakets[0]) * self.ladder_plus(brakets[1]);
        val1 + val2 + val3
    }

    pub fn proj_z_minus(&self, brakets: &Braket) -> f64 {
        let spin_ket = self
            .composite_state
            .state_from_base(&brakets.ket)
            .q_number();
        let spin_bra = self
            .composite_state
            .state_from_base(&brakets.bra)
            .q_number();

        if brakets.ket.q_number() == 0 && brakets.bra.q_number() == 0 {
            if spin_ket != spin_bra {
                1.0
            } else {
                0.0
            }
        } else {
            0.0
        }
    }

    pub fn ladder_plus_minus(&self, brakets: &Braket) -> f64 {
        let spin_ket = self
            .composite_state
            .state_from_base(&brakets.ket)
            .q_number();
        let spin_bra = self
            .composite_state
            .state_from_base(&brakets.bra)
            .q_number();

        if brakets.ket.q_number() + 2 == brakets.bra.q_number() {
            if spin_ket == 0 && spin_bra == 2 {
                -2.0_f64.sqrt()
            } else if spin_ket == 2 && spin_bra == 0 {
                2.0_f64.sqrt()
            } else {
                0.0
            }
        } else {
            0.0
        }
    }

    pub fn ladder_minus_minus(&self, brakets: &Braket) -> f64 {
        let spin_ket = self
            .composite_state
            .state_from_base(&brakets.ket)
            .q_number();
        let spin_bra = self
            .composite_state
            .state_from_base(&brakets.bra)
            .q_number();

        if brakets.ket.q_number() - 2 == brakets.bra.q_number() {
            if spin_ket == 2 && spin_bra == 0 {
                -2.0_f64.sqrt()
            } else if spin_ket == 0 && spin_bra == 2 {
                2.0_f64.sqrt()
            } else {
                0.0
            }
        } else {
            0.0
        }
    }

    pub fn dot_minus(&self, brakets: [&Braket; 2]) -> f64 {
        assert!(self.is_spin_state(brakets[0].bra));
        assert!(self.is_spin_state(brakets[0].ket));
        assert!(self.is_spin_state(brakets[1].bra));
        assert!(self.is_spin_state(brakets[1].ket));

        let val1 = self.proj_z_minus(brakets[0]) * self.proj_z(brakets[1]);
        let val2 = 0.5 * self.ladder_plus_minus(&brakets[0]) * self.ladder_minus(&brakets[1]);
        let val3 = 0.5 * self.ladder_minus_minus(&brakets[0]) * self.ladder_plus(&brakets[1]);

        val1 + val2 + val3
    }

    fn is_spin_state(&self, base: &StateBase) -> bool {
        self.composite_state.state_from_base(&base).algebra() == "spin"
    }

    fn same_spin(&self, braket: &Braket) -> bool {
        let spin_ket = self.composite_state.state_from_base(&braket.ket).q_number();
        let spin_bra = self.composite_state.state_from_base(&braket.bra).q_number();

        spin_ket == spin_bra
    }

    pub fn preserve_projection(&self, brakets: &[&Braket]) -> f64 {
        for braket in brakets {
            assert!(self.is_spin_state(braket.bra));
            assert!(self.is_spin_state(braket.ket));
        }
        let kets_proj = brakets.iter().map(|b| b.ket.q_number()).sum::<isize>();
        let bras_proj = brakets.iter().map(|b| b.bra.q_number()).sum::<isize>();

        if kets_proj == bras_proj {
            1.0
        } else {
            0.0
        }
    }

    pub fn zeeman_operator<'a>(
        &'a self,
        operator: &mut Operator<'a>,
        spins: &'a [&str],
        gamma_factors: &'a [f64],
    ) {
        for (spin, gamma) in spins.iter().zip(gamma_factors.iter()) {
            operator.add_operator(vec![*spin], |brakets| -*gamma * self.proj_z(&brakets.first().unwrap()));
        }
    }

    pub fn symmetrize_spins(&self, brakets: [&Braket; 2]) -> f64 {
        assert!(self.is_spin_state(brakets[0].bra));
        assert!(self.is_spin_state(brakets[0].ket));
        assert!(self.is_spin_state(brakets[1].bra));
        assert!(self.is_spin_state(brakets[1].ket));

        let bases_right = brakets.map(|b| b.ket.q_number());
        let bases_left = brakets.map(|b| b.bra.q_number());

        if bases_right[0] == bases_right[1] && bases_right == bases_left {
            1.0
        } else if bases_right[0] != bases_right[1]
            && (bases_right == bases_left || bases_right == [bases_left[1], bases_left[0]])
        {
            0.5
        } else {
            0.0
        }
    }

    pub fn anti_symmetrize_spins(&self, brakets: [&Braket; 2]) -> f64 {
        assert!(self.is_spin_state(brakets[0].bra));
        assert!(self.is_spin_state(brakets[0].ket));
        assert!(self.is_spin_state(brakets[1].bra));
        assert!(self.is_spin_state(brakets[1].ket));

        let bases_right = brakets.map(|b| b.ket.q_number());
        let bases_left = brakets.map(|b| b.bra.q_number());

        if bases_right[0] != bases_right[1] && bases_right == bases_left {
            0.5
        } else if bases_right[0] != bases_right[1] && bases_right == [bases_left[1], bases_left[0]]
        {
            -0.5
        } else {
            0.0
        }
    }

    pub fn flip_from_singlet(&self, brakets: [&Braket; 2]) -> f64 {
        assert!(self.is_spin_state(brakets[0].bra));
        assert!(self.is_spin_state(brakets[0].ket));
        assert!(self.is_spin_state(brakets[1].bra));
        assert!(self.is_spin_state(brakets[1].ket));

        let bases_right = brakets.map(|b| b.ket.q_number());
        let bases_left = brakets.map(|b| b.bra.q_number());

        let singlet_proj_right = clebsch_gordan(1, bases_right[0], 1, bases_right[1], 0, 0);
        let triplet_proj_left = clebsch_gordan(1, bases_left[0], 1, bases_left[1], 2, -2)
            + clebsch_gordan(1, bases_left[0], 1, bases_left[1], 2, 0)
            + clebsch_gordan(1, bases_left[0], 1, bases_left[1], 2, 2);

        singlet_proj_right * triplet_proj_left
    }

    pub fn flip_from_triplet(&self, brakets: [&Braket; 2]) -> f64 {
        assert!(self.is_spin_state(brakets[0].bra));
        assert!(self.is_spin_state(brakets[0].ket));
        assert!(self.is_spin_state(brakets[1].bra));
        assert!(self.is_spin_state(brakets[1].ket));

        let bases_right = brakets.map(|b| b.ket.q_number());
        let bases_left = brakets.map(|b| b.bra.q_number());

        let triplet_proj_right = clebsch_gordan(1, bases_right[0], 1, bases_right[1], 2, -2)
            + clebsch_gordan(1, bases_right[0], 1, bases_right[1], 2, 0)
            + clebsch_gordan(1, bases_right[0], 1, bases_right[1], 2, 2);

        let singlet_proj_left = clebsch_gordan(1, bases_left[0], 1, bases_left[1], 0, 0);

        triplet_proj_right * singlet_proj_left
    }

    pub fn symmetry_flip(&self, brakets: [&Braket; 2]) -> f64 {
        self.flip_from_singlet(brakets) + self.flip_from_triplet(brakets)
    }

    pub fn symmetry_stay(&self, brakets: [&Braket; 2]) -> f64 {
        let bases_right = brakets.map(|b| b.ket.q_number());
        let bases_left = brakets.map(|b| b.bra.q_number());

        let triplet_proj_right = clebsch_gordan(1, bases_right[0], 1, bases_right[1], 2, -2)
            + clebsch_gordan(1, bases_right[0], 1, bases_right[1], 2, 0)
            + clebsch_gordan(1, bases_right[0], 1, bases_right[1], 2, 2);
        let singlet_proj_right = clebsch_gordan(1, bases_right[0], 1, bases_right[1], 0, 0);

        let triplet_proj_left = clebsch_gordan(1, bases_left[0], 1, bases_left[1], 2, -2)
            + clebsch_gordan(1, bases_left[0], 1, bases_left[1], 2, 0)
            + clebsch_gordan(1, bases_left[0], 1, bases_left[1], 2, 2);
        let singlet_proj_left = clebsch_gordan(1, bases_left[0], 1, bases_left[1], 0, 0);

        triplet_proj_left * triplet_proj_right + singlet_proj_left * singlet_proj_right
    }
}

pub fn projection_range(double_spin: usize) -> StepBy<Rev<RangeInclusive<isize>>> {
    (-(double_spin as isize)..=double_spin as isize)
        .rev()
        .step_by(2)
}
