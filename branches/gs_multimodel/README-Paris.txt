PARIS Build Instructions

1.0	SOCI and SQLite
1.1	Paris depends on the open source library, SOCI, to connect to the database. For recent versions of Linux and OS X, SQLite is probably already installed. However, SOCI does depend on a specific version of SQLite being installed. I've included instructions for building both libraries. The following directions assume wget and git are both installed. These are pretty standard on modern linux distributions. Because SQLite isn't supported by default, you can't just download the source from the website. Instead, you have to get it from the head via git (I think you can also access it via svn or cvs)
1.2	Building SQLite
1.2.1	Download a recent version of SQlite
1.2.1.1 	Extract and build sources
1.2.1.1.1	wget http://www.sqlite.org/sqlite-amalgamation-3.6.23.1.tar.gz
		Once it's downloaded, extract the tarball to a directory where you can build the library. 
1.2.1.1.2	./configure --prefix=/scratch/soci --enable-shared=no
        	In this example, we are turning off shared libraries and writing someplace other than the
		 default directory. If you have sudo rights, you can leave both of these options off. 
1.2.1.1.3	make install
1.2.2	Download a recent version of SOCI
1.2.2.1		git clone git://soci.git.sourceforge.net/gitroot/soci/soci
1.2.2.1.1	From within the src folder inside the archive, generate the configure script for your platform:
1.2.2.1.2	./autogen.sh
1.2.2.1.3	./configure --enable-backend-sqlite3=yes --with-sqlite3=/scratch/soci --prefix=/scratch/soci --enable-shared=no
                This activates the sqlite3 connectors and points to the path where sqlite can be found. The path is only
		necessary if you don't have sudo rights. By turning off shared, we can avoid having to make the
		library available at a place that is searched. If you are installing directly to the default paths, this 
		isn't necessary.
1.2.2.1.4	make install
1.3	Extract paris sources
1.3.1	Edit the extlibs.make file, set the first line to reflect the path to the directory you installed the soci and 
	SQLite3 libraries to.
1.3.2	make 

The executable will be bin/paris64 or bin/paris depending on if your system is 32 or 64bit linux. 

The software should compile on windows and OS X and any other system for which GCC is available. Windows users
will need to set up MinGW in order to compile the software (www.mingw.org). OS X users need only to ensure that
they have GCC installed (install the developer tools from the installation media).

When compiling for windows, users should use WIN32=1 when building paris (step 1.3.2, i.e. make WIN32=1 )




