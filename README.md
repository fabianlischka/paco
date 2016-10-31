# PACO
======

On [Github](https://github.com/fabianlischka/paco/).
Original code (C, Matlab) written by [Felix Kwok](http://www.math.hkbu.edu.hk/~felix_kwok/), see papers below.

## Overview

This code solves an optimal control problem, where we would like to minimise
the deviation of some state trajectory y(t) from a target trajectory y^hat(t),
and the deviation of the final state y(T) from a target state y^hat_T, by
controlling a function u(t) that influences y(t) via the relationship
```
y_t + Ay = Bu + f
```

In addition to this governing PDE (and initial condition) going forward in time,
the optimal solution satisfies an adjoint PDE (and final condition) going backward
in time, leading to a two-point boundary problem.

The code here employs time domain decomposition to solve this problem in parallel.
(Thus: **Pa** rallel Optimal **Co** ntrol.)

### Motivating Question

Imagine we have a flowing river. Now we inject some color into it, and we are
interested in how the color flows. In particular, we would like it to be in a
certain place at a certain time. When and where do we have to drop the color to
have it as close to the target as possible?

### Motivating Example: Details

We start with an incompressible fluid flow in some 2d geometry, in steady state,
which can be be described by the Stokes equation

laplacian(U) = grad(p)
div(U) = 0

where U =(u(x,y), v(x,y)) is the fluid velocity (a vector field),
and p =p(x,y) is its pressure (a scalar field).

We then solve the advection-diffusion equation

c_t = laplacian(c) - U . grad(c) + R

choosing R to minimise the deviation of c.

Note that FreeFem++ uses the weak formulation of the problem, ie multiply with
a test function (nice and continuous etc., compact support), and integrate by
parts over infinity = the support.

### Quick reminder of vector notations

* gradient of a scalar field s is a vector: (ds/dx, ds/dy, ...)
* divergence of a vector field F=(F1,F2,...) is a scalar: dF1/dx + dF2/dy + ...
* laplacian of a scalar field s is a scalar, namely the divergence of the gradient,
akin to the "curvature". Heat equation (diffusion): ds/dt = c laplacian(s)
* (vector) laplacian of a vector field F is (in cartesian coord's) the vector
of the laplacians of the components of F, laplacian(F) = (laplacian(F1), laplacian(F2), ...)
* product rule: div(sF) = grad(s) . F + s div(F), a scalar

### Formulation

Given f=f(t), with t in (0,T), and matrices A, B,
and y=y(t) (with time derivative y_t = del y/del t), determined by the constraints

(1a) y_t + Ay = Bu + f
(1b) y(0) = y0

and parameters a, b, and matrices C, D and target state y^hat = y^hat(t) and y^hat_T,
this code finds an optimal control u=u(t) with t in (0,T) such that the objective
functional F(u) is minimised, where

(2) F(u) = 1/2 int_0^T |u(t)|^2 dt + a/2 int_0^T |Cy-y^hat|^2 dt + b|Dy(T)-y^hat_T|^2

See F Kwok: "On the Time-domain Decomposition of Parabolic Optimal Control Problems"
in Section "Background" below for details.

### Inputs/Outputs

* The problem as outlined above needs, as inputs:
  * A, B, C, D, a, b, y0, y^hat(t), y^hat_T

The full run proceeds in 3 stages (FreeFem, Python/Matlab, C code).
Stage n takes the inputs/data from stage n-1 to stage n.

1. Stage 1: FreeFem++
  * reads: prefix + 's0_' + **stokes.edp**
  * writes: prefix + 's1_' +
    * u.txt, v.txt, mass.txt, stiff.txt, Rih.txt
    * stokes.msh (mesh, see [Section 5.1.4 "Data Structures and Read/Write Statements for a Mesh"](http://www.freefem.org/ff++/ftp/freefem++doc.pdf#subsection.5.1.4) for format)
    * bay_flux.ps (plot)
1. Stage 2: Python (previously: Matlab)
  * reads: prefix + 's1_' +
    * u.txt, v.txt, mass.txt, stiff.txt, Rih.txt
    * stokes.msh
  * writes: prefix + 's2_' +
    * A.txt, B.txt, C.txt, D.txt
2. Stage 3: C code
  * control_main.c
    * reads: prefix + 's2_' +
      * A.txt, B.txt, C.txt, D.txt
      * params.txt (for beta, Nt, a1, a2, q, pp, qq, max_iter, tol, krylov
      * yhat.txt
    * writes: prefix + 's3_'
      * sol.txt
      * yy.txt (trajectory)
  * control_main_full.c
    * similar to control_main.c, but yhat is computed, not read, and yy.txt is not written
  * dd_main.c
    * reads:
      * like control_main.c
      * existing sol.txt
    * writes:
      * nothing, just outputs the error vis-a-vis the control

## Installation

### Prerequisites

  * Python 2, with numpy
  * FreeFem++ (see installation notes below)

## List of files

Notes
* files containing `main` in the name contain a `main` routine and can be run
* there are three ways to run a problem:
  * `control_main` runs it without domain decomposition
  * `dd_main` runs it with domain decomposition, by simple iteration
  * `dd_gmres` runs it with domain decomposition, with GMRES for the outer loop
* the `_full` version uses a different y_hat.


Length | Name | Comment
-----: | ---- | -------
     . | examples | directory of examples
  1902 | control.h |
  3556 | control_main.c | run this for single-node
  3411 | control_main_full.c | idem, for full variant
  4448 | control_test.c | don't solve, just check inputs
 12402 | dd_gmres.c | run this for multi-node GMRES
 12396 | dd_gmres_full.c | idem, for full variant
  6816 | dd_main.c | run this for multi-node iterative
  7090 | dd_main_full.c | idem, for full variant
  8776 | IO.c | helper routines for reading variables and params
  1605 | Makefile |
   201 | minitest.c |
 14144 | Problem.c | wrapper containing problem variables, params, and solver
   731 | read_parallel.c | util to test file availability on nodes
  1480 | Solver.c | solver for the inner problem
  7930 | test.pbs | batch file for parallel run

## Useful links

### Background

* M Gander, F Kwok: "Schwarz Methods for the Time-Parallel Solution of Parabolic Control Problems" [PDF](http://www.math.hkbu.edu.hk/~felix_kwok/docs/Kwok-Felix-Gander-Martin_J.-276.pdf)
* F Kwok: "On the Time-domain Decomposition of Parabolic Optimal Control Problems" [PDF](http://www.math.hkbu.edu.hk/~felix_kwok/docs/kwok_plenary.pdf)
* M Gander, F Kwok, G Wanner: "History of constrained optimization" [PDF](http://www.math.hkbu.edu.hk/~felix_kwok/docs/GanderKwokWanner_OPTPDE.pdf)

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

If you're developing on a windows box, here is some software you can quickly
install to make life easier, and put a virtual machine with CentOS in place,
so you can test stuff locally before going onto the sciblade cluster.

* install [Chocolatey](https://chocolatey.org/install), a Windows package manager
* install other required software: run (in a PowerShell, as administrator)
  `choco install -y anaconda2 git vagrant virtualbox`
    * obviously, skip software you have installed already or don't need
    * [Anaconda](https://www.continuum.io/anaconda-overview) is a modern
    Python distribution that contains a lot of scientific libraries,
    such as NumPy, SciPy, Jupyter, etc.
    We choose to download the version with Python 2, as Sciblade is on Python 2.
    * [git](https://git-scm.com/) is a modern, distributed version control system
    * [Vagrant](https://www.vagrantup.com/) allows very simple management of virtual machines
    * [VirtualBox]() is the VM "provider" used by Vagrant
      * Note: Vagrant will install VirtualBox by itself, if not installed, but
      it installs a version (5.0.10) that crashes on my Windows 10, so it's
      easiest to just install the latest version with choco.
* install some other useful software as desired: run (in a PowerShell, as administrator)
  `choco install -y atom kdiff3 putty winscp`
    * obviously, skip software you have installed already or don't need
    * [Atom](https://atom.io/) is a free, modern, extensible editor
    * [KDiff3](http://kdiff3.sourceforge.net/) highlights differences between
    (up to 3) text files
    * [PuTTY](http://www.putty.org/) is a free Telnet/SSH client for Windows
    * [WinSCP](https://winscp.net/eng/index.php) is a free SFTP, SCP and FTP client for Windows

### virtual CentOS box

Create a virtual machine using CentOS 6, which is (basically) what's running
on Sciblade (as of 2016-10-14).

* create a CentOS 6 box, see [release details](https://seven.centos.org/2016/10/updated-centos-vagrant-images-available-v1609-01/)
  * go do the directory with Vagrantfile (should be `paco\vm`)
  * add guest additions with `vagrant plugin install vagrant-vbguest`
  * start it with `vagrant up`
  * this will download the box, add the Virtualbox guest additions,
    and provision it, including installation of numpy and FreeFem++

### HKBU SciBlade specific

* connect to sciblade with `ssh lischka@sciblade.sci.hkbu.edu.hk` or on Windows
after installing PuTTY `putty -ssh lischka@sciblade.sci.hkbu.edu.hk`
* Directories:
  * data in `/u1/local/share/felix_grp`
  * FreeFem++ in `/u1/local/ff++`
* ` /u1/local/ff++/bin/FreeFem++-nw stokes.edp`

### Workflow (old way, for reference)

0. set up stokes.edp as desired.
1. run `FreeFem++-nw stokes.edp`
  * now have files:
      * bay_flux.ps, stokes.msh
      * mass.txt,  Rih.txt,  stiff.txt, u.txt, v.txt
2. extract matrices
  * mangle manually, in Matlab:
    * ff, vv from stokes.msh
    * u, v from u.txt, v.txt
    * mi, mj, ms from mass.txt
    * ki, kj, ks from stiff.txt
    * bi, bj, bs from Rih.txt
  * store in, say, bay.mat
3. compute matrices A.txt and B.txt, and (hardcoded) C.txt, D.txt
    * using bay.mat:
    * run getA.m, which uses bay.ff, bay.vv, bay.u, bay.v for fvm.m,
      then also bay.mi, bay.mj, bay.ms to compute M, M2,
      then also bay.ki, bay.kj, bay.ks to compute K
      then also bay.bj, bay.bi, ones to compute A,
      then Am (diagonal trans) and Bm
    * these are written out, manually, with write_csr_matrix.m, to A.txt and B.txt
    * C.txt and D.txt are hardcoded
3. run (on the cluster) control_main, dd_main, or dd_gmres




# OBSOLETE (supserseded notes)

## Vagrant

### Guest Additions for virtualbox

OBSOLETE - Instead, install [vagrant-vbguest](https://github.com/dotless-de/vagrant-vbguest) with `vagrant plugin install vagrant-vbguest`

* Note: might have to manually install the Guest Additions, see [here](https://www.virtualbox.org/manual/ch04.html)
  * after putting `C:\Program Files\Oracle\VirtualBox\VBoxGuestAdditions.iso`
    into the optical drive in the virtual machine in the VirtualBox GUI,
    log onto the VM using `vagrant ssh` and run
    ```
    sudo yum update
    sudo yum install -y dkms epel-release kernel-devel numpy.x86_64
    sudo yum groupinstall -y 'Development Tools'
    sudo mkdir /mnt/cdrom
    sudo mount /dev/sr0 /mnt/cdrom
    cd /mnt/cdrom
    sh ./VBoxLinuxAdditions.run
    ```


### Rsync failure

OBSOLETE - Instead, use "synched-folder" type `SMB` on Windows, or just `virtualbox`, by putting this in your Vagrantfile:
`config.vm.synced_folder ".", "/vagrant", type: "virtualbox"`

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


### Notes: Installation of FreeFem++

OBSOLETE: instead, use in `vm`, the `Vagrantfile`, with its link to `installFreeFem.sh` which does everything below, including patching the makefile
using `removeDistanceExample.patch`

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
