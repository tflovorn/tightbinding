extern crate num_complex;
extern crate rulinalg;

mod model;
pub use model::Model;

pub mod fourier;
pub mod w90;

#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {}
}
