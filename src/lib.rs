// Recommended setting for error_chain.
#![recursion_limit = "1024"]

#[macro_use]
extern crate error_chain;
extern crate itertools;
extern crate linxal;
extern crate ndarray;
extern crate num_complex;
extern crate rayon;
extern crate serde;
#[macro_use]
extern crate serde_derive;
extern crate sxd_document;
extern crate sxd_xpath;

pub mod float;
pub mod vec_util;
pub mod units;

mod model;
pub use model::Model;

pub mod fourier;
pub mod qe;
pub mod w90;

pub mod paths;

pub mod tetra;
pub mod dos;
