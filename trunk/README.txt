Building genomeSIMLA
====================

To build genomeSIMLA, from this directory, type:
./configure
make
make install

You will need to have root permissions to execute "make install" with the 
default configure options.  Additional helpful options that can be passed
to configure are:

--prefix=<directory>
	This defines the directory to install genomeSIMLA into.  genomeSIMLA 
	will be copied into $prefix/bin, and by default, the prefix is /usr/local.
	
--enable-debug
	This option enables all debugging symbols and warnings and turns off
	any optimizations.  This is helpful if you encounter a bug and would like
	to fix it, but running in debug mode will likely result in MUCH slower 
	execution times.
	
Ubuntu Requirements
===================
In addition to the standard build tools, building the software for Ubuntu (tested with
10.04) requires the following packages:
libboost-all-dev
libpng-dev
libfreetype6-dev
	
Font Files
==========

Note that in this new version, we are not installing the font file, but we are
still including it.  This means that you will likely need to change the path
to the FONT directive in the configuration file.  

Also, the default value for the FONT directive has changed to:
/usr/share/fonts/liberation/LiberationMono-Bold.ttf

Test genomeSIMLA
====================

To test with sample inputs, change to sample directory:
cd sample
./testGenomeSIMLA.sh
