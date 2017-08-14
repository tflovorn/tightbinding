extern crate num_complex;
extern crate rulinalg;

mod model;
pub use model::Model;

pub mod fourier;

#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {}
}
