GNU package-type make is tried. Namely, you should:

mkdir <builddir>
cd <builddir>

<srcdir>/configure --enable-optimise for NOARCH (single node C++ codes)

or

<srcdir>/configure --build sparc-sun-solaris2.9 --host powerpc-gnu-elf
--enable-target=qcdoc for QCDOC

gmake depend
gmake


environment variable QOS should be set to point where QOS is installed
at compile time.

if you want to build the library or executables from the source tree,
type 'gmake install' from the build tree before doing so. 'gmake
install' copies appropriate Makefile's to the source directory.

Good Luck!!
