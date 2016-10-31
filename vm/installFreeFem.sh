echo "Installing prerequisites for FreeFrem++"
#    sudo yum install -y dkms epel-release kernel-devel
sudo yum groupinstall -y 'Development Tools'
sudo yum install -y numpy.x86_64 lapack-devel.x86_64 arpack-devel.x86_64 suitesparse-devel.x86_64
echo "Downloading FreeFrem++"
wget -nv http://www.freefem.org/ff++/ftp/freefem++-3.49.tar.gz -O ff.tar.gz
tar zxvf ff.tar.gz
cd freefem++-3.49
# autoreconf -i  # will likely fail, needs autoconf version >= 2.69
# instead do the next two lines to download pre-configuration
wget -nv http://www.freefem.org/ff++/ff++/AutoGeneratedFile.tar.gz -O AutoGeneratedFile.tar.gz
tar zxvf AutoGeneratedFile.tar.gz
# the distance.edp example does not compile for me, with errors
# default: distance.cpp:297: error: reference to 'Tet' is ambiguous
# etc. We apply a patch to just drop it from the Makefile
cd examples++-load/
patch < /vagrant/removeDistanceExample.patch
cd ..
# we don't need mpi support (though users might want it)
# we definitely don't need graphics support (glut) on the cluster
./configure -without-mpi -without-glut
echo "Making FreeFrem++"
make
make check   # takes a while
make install
echo "# to test:"
echo "cd examples++-tutorial"
echo "FreeFem++-nw LaplaceP1.edp"
