# Dependencies

Test data is stored in Git LFS. Before cloning the repository, install Git LFS:

    cd ~
    curl -L -o git-lfs-linux-amd64-2.2.1.tar.gz https://github.com/git-lfs/git-lfs/releases/download/v2.2.1/git-lfs-linux-amd64-2.2.1.tar.gz
    tar -xvzf git-lfs-linux-amd64-2.2.1.tar.gz
    cd git-lfs-2.2.1
    PREFIX=$HOME ./install.sh

This assumes $HOME/bin is on your $PATH. If it is not, add the following to ~/.bashrc:

    export PATH=$HOME/bin:$PATH

Requires liblapacke. To install on Debian-based distributions:

    sudo apt-get install liblapacke-dev
