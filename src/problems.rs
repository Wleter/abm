use std::collections::VecDeque;

use quantum::problem_selector::ProblemSelector;

use self::{
    asymptotic_bound::AsymptoticBound, hifi_double::HifiDouble, hifi_single::HifiSingle,
    potassium_bound::PotassiumBound, lithium_potassium::LithiumPotassium,
};

pub mod asymptotic_bound;
pub mod hifi_double;
pub mod hifi_single;
pub mod potassium_bound;
pub mod lithium_potassium;

pub struct Problems;

impl ProblemSelector for Problems {
    const NAME: &'static str = "ABM implementation";

    fn list() -> Vec<&'static str> {
        vec![
            "hyperfine one atom",
            "hyperfine two atoms",
            "asymptotic bound",
            "potassium bound",
            "lithium potassium"
        ]
    }

    fn methods(number: &str, args: &mut VecDeque<String>) {
        match number {
            "0" => HifiSingle::run(),
            "1" => HifiDouble::run(),
            "2" => AsymptoticBound::run(),
            "3" => PotassiumBound::select(args),
            "4" => LithiumPotassium::select(args),
            _ => println!("Invalid problem number"),
        }
    }
}
