use std::ops::{Index, Deref};

use crate::basis_element::BasisElement;

#[derive(Clone, Copy)]
pub struct Braket<'a> {
    pub bra: &'a BasisElement,
    pub ket: &'a BasisElement,
}

impl<'a> Braket<'a> {
    pub fn new(bra: &'a BasisElement, ket: &'a BasisElement) -> Self {
        Self { bra, ket }
    }

    pub fn same_state(&self) -> bool {
        self.bra == self.ket
    }
}

#[derive(Clone)]
pub struct Brakets<'a, 'b> {
    pub brakets: Vec<Braket<'a>>,
    pub names: &'b [&'a str],
}

impl<'a, 'b> Brakets<'a, 'b> {
    pub fn new(brakets: Vec<Braket<'a>>, names: &'b [&'a str]) -> Self {
        assert!(brakets.len() == names.len());
        Self { 
            brakets, 
            names 
        }
    }
}

impl<'a> Index<&str> for Brakets<'a, '_> {
    type Output = Braket<'a>;

    fn index(&self, index: &str) -> &Self::Output {
        let i = self.names
            .iter()
            .position(|x| *x == index);

        let i = match i {
            Some(i) => i,
            None => panic!("Braket name {} not found", index),
        };

        // SAFETY: index is always in bounds from assert! in Brakets::new
        unsafe {
            self.brakets.get_unchecked(i)
        }
    }
}

impl<'a> Deref for Brakets<'a, '_> {
    type Target = [Braket<'a>];

    fn deref(&self) -> &Self::Target {
        &self.brakets
    }
}