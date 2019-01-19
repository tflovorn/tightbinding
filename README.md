# tightbinding

[![Build Status](https://travis-ci.org/tflovorn/tightbinding.svg?branch=master)](https://travis-ci.org/tflovorn/tightbinding)

Tools for working with tight-binding models (defined by a map from lattice
vectors R to H(R) = <0|H|R>) and models defined directly on k-space as H(k).

Includes an implementation of the linear tetrahedron method with curvature
corrections as defined in Blöchl, Jepsen, and Andersen, PRB 49, 16223 (1994).

## Dependencies

Test data is stored in Git LFS. Before cloning the repository, install Git LFS:

    cd ~
    curl -L -o git-lfs-linux-amd64-2.2.1.tar.gz https://github.com/git-lfs/git-lfs/releases/download/v2.2.1/git-lfs-linux-amd64-2.2.1.tar.gz
    tar -xvzf git-lfs-linux-amd64-2.2.1.tar.gz
    cd git-lfs-2.2.1
    PREFIX=$HOME ./install.sh

This assumes $HOME/bin is on your $PATH. If it is not, add the following to ~/.bashrc:

    export PATH=$HOME/bin:$PATH

Pre-commit hook checks that code installed by package `pre-commit` checks that
code is unchanged after running `cargo fmt`. To get this command:

    rustup component add rustfmt

## Binaries

### `tb_w90_dos`

Extract a tight-binding model from an `hr.dat` file as produced by Wannier90
and an accompanying `data-file.xml` from a self-consistent Quantum Espresso
calculation. Produce the density of states of this model.
