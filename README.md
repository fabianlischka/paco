# PACO
======

## Overview

## List of files

Notes
* files containing `main` in the name contain a `main` routine and can be run
* there are three ways to run a problem:
  * `control_main` runs it without domain decomposition
  * `dd_main` runs it with domain decomposition, by simple iteration
  * `dd_gmres` runs it with domain decomposition, with GMRES for the outer loop
* the `_full` version uses a different y_hat.


Length | Name | Comment
-----: | ---- |
       | examples | Directory of examples
  1902 | control.h |
  3556 | control_main.c |
  3411 | control_main_full.c |
  4448 | control_test.c |
 12402 | dd_gmres.c |
 12396 | dd_gmres_full.c |
  6816 | dd_main.c |
  7090 | dd_main_full.c |
  8776 | IO.c |
  1605 | Makefile |
   201 | minitest.c |
 14144 | Problem.c |
   731 | read_parallel.c |
  1480 | Solver.c |
  7930 | test.pbs |

## Useful links

### Theory

* M Gander, F Kwok: "Schwarz Methods for the Time-Parallel Solution of Parabolic Control Problems" [PDF](http://www.math.hkbu.edu.hk/~felix_kwok/docs/Kwok-Felix-Gander-Martin_J.-276.pdf)
* F Kwok: "On the Time-domain Decomposition of Parabolic Optimal Control Problems" [PDF](http://www.math.hkbu.edu.hk/~felix_kwok/docs/kwok_plenary.pdf)
* M Gander, F Kwok, G Wanner: "History of constrained optimization" [PDF]http://www.math.hkbu.edu.hk/~felix_kwok/docs/GanderKwokWanner_OPTPDE.pdf)

### HK Baptist University
* HKBU High Performance Cluster, Sciblade [FAQ](http://site.sci.hkbu.edu.hk/hpccc/sciblade/faq_sciblade.php) and
[installed software](http://site.sci.hkbu.edu.hk/hpccc/sciblade/software.php)

### Software Components

#### FreeFem++

* [FreeFem++](http://www.freefem.org/ff++/)

#### MPI

* [MVAPICH2](http://mvapich.cse.ohio-state.edu/features/#mv2), an
implementation of [MPI-3](http://mpi-forum.org/docs/mpi-3.0/mpi30-report.pdf)
* Prof. Alonso's [lecture notes for CME342](http://adl.stanford.edu/cme342/Lecture_Notes.html) with an introduction
to MPI, including application to Domain Decomposition



## Steps

### Windows

* install [Chocolatey](https://chocolatey.org/install), a Windows package manager
* run (in a PowerShell, as administrator)
  `choco install -y anaconda2 atom git putty rsync vagrant virtualbox`
    * obviously, skip software you have installed already or don't need
    * [Anaconda](https://www.continuum.io/anaconda-overview) is a modern
    Python distribution that contains a lot of scientific libraries,
    such as NumPy, SciPy, Jupyter, etc.
    We choose to download the version with Python 2, as Sciblade is on Python 2.
    * [Atom](https://atom.io/) is a modern, extensible editor
    * [git](https://git-scm.com/) is a modern, distributed version control system
    * [PuTTY](http://www.putty.org/) is a Telnet/SSH client for Windows
    * rsync, required for Vagrant to work
    * [Vagrant](https://www.vagrantup.com/) allows very simple management of virtual machines
    * [VirtualBox]() is the VM "provider" used by Vagrant
      * Note: Vagrant will install VirtualBox by itself, if not installed, but
      it installs a version (5.0.10) that crashes on my Windows 10, so it's
      easiest to just install the latest version with choco.

### virtual CentOS box

Create a virtual machine using CentOS 6, which is (basically) what's running
on Sciblade (as of 2016-10-14).

* create a CentOS box with  `vagrant init centos/6`, then start it with
  `vagrant up`
  * Note: might have to manually install the Guest Additions, see [here](https://www.virtualbox.org/manual/ch04.html)
    * after putting `C:\Program Files\Oracle\VirtualBox\VBoxGuestAdditions.iso`
      into the optical drive in the virtual machine in the VirtualBox GUI,
      log onto the VM using `vagrant ssh` and run
      ```
      sudo yum update
      sudo yum install dkms epel-release kernel-devel
      sudo yum groupinstall 'Development Tools'
      sudo mkdir /mnt/cdrom
      sudo mount /dev/sr0 /mnt/cdrom
      cd /mnt/cdrom
      sh ./VBoxLinuxAdditions.run
      ```
    * if upon `vagrant up`, rsync fails because of a wrong path (`/c/foo/bar` instead
      of `C:/foo/bar`), you'll have to hack a vagrant ruby file,
      `C:\HashiCorp\Vagrant\embedded\gems\gems\vagrant-1.8.6\plugins\synced_folders\rsync\helper.rb`, to replace
      ```
      if Vagrant::Util::Platform.windows?
        # rsync for Windows expects cygwin style paths, always.
        hostpath = Vagrant::Util::Platform.cygwin_path(hostpath)
      end
      ```
      by
      ```
      if Vagrant::Util::Platform.windows?
        # rsync for Windows expects cygwin style paths, always.
        hostpath = "/cygdrive" + Vagrant::Util::Platform.cygwin_path(hostpath)
      end
      ```
      that is, add `/cygdrive` in front of the path. See [issue #3230](https://github.com/mitchellh/vagrant/issues/3230) for more.



* connect to sciblade with `ssh lischka@sciblade.sci.hkbu.edu.hk` or on Windows
after installing PuTTY `putty -ssh lischka@sciblade.sci.hkbu.edu.hk`


### Workflow (old way, for reference)

0. set up stokes.edp as desired.
1. run `FreeFem++-nw stokes.edp`
  * now have files:
      * bay_flux.ps, stokes.msh (can be ignored)
      * mass.txt,  Rih.txt,  stiff.txt, u.txt, v.txt
2. compute matrices A.txt and B.txt, and (hardcoded) C.txt, D.txt
    * pull these into Matlab, generate bay.mat
    * run getA.m, which uses bay.ff, bay.vv, bay.u, bay.v for fvm.m,
      then also bay.mi, bay.mj, bay.ms to compute M, M2,
      then also bay.ki, bay.kj, bay.ks to compute K
      then also bay.bj, bay.bi, ones to compute A,
      then Am (diagonal trans) and Bm
    * these are written out, manually, to A.txt and B.txt
    * C.txt and D.txt are hardcoded
3. run control_main, dd_main, or dd_gmres


### Notes: Installation of FreeFem++

* Installation on Linux [instructions here](http://www.freefem.org/ff++/linux.php)
* Source available [here](http://www.freefem.org/ff++/ftp/freefem++-3.49.tar.gz),
see below wget command
* prerequisites: lapack, arpack, suitesparse, general development tools -
`sudo yum install -y lapack-devel.x86_64 arpack-devel.x86_64 suitesparse-devel.x86_64`
* make like this:
```
wget http://www.freefem.org/ff++/ftp/freefem++-3.49.tar.gz -O ff.tar.gz
tar zxvf ff.tar.gz
cd freefem++-3.49
# autoreconf -i  # will likely fail, needs autoconf version >= 2.69
# instead do the next two lines to download pre-configuration
wget http://www.freefem.org/ff++/ff++/AutoGeneratedFile.tar.gz -O AutoGeneratedFile.tar.gz
tar zxvf AutoGeneratedFile.tar.gz
# we don't need mpi support (though users might want it)
# we definitely don't need graphics support (glut) on the cluster
./configure -without-mpi -without-glut
make
make check   # takes a while
make install
# to test:
cd examples++-tutorial
FreeFem++-nw LaplaceP1.edp
```
* note: for me, distance.cpp in examples++-load failed - just delete it or rename .cpp-skip, and remove all references to distance in examples++-load/Makefile
