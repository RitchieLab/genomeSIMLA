/* Copyright (C) 2002 MySQL AB

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; version 2 of the License.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA */

/*
** example file of UDF (user definable functions) that are dynamicly loaded
** into the standard mysqld core.
**
** The functions name, type and shared library is saved in the new system
** table 'func'.  To be able to create new functions one must have write
** privilege for the database 'mysql'.	If one starts MySQL with
** --skip-grant, then UDF initialization will also be skipped.
**
** Syntax for the new commands are:
** create function <function_name> returns {string|real|integer}
**		  soname <name_of_shared_library>
** drop function <function_name>
**
** Each defined function may have a xxxx_init function and a xxxx_deinit
** function.  The init function should alloc memory for the function
** and tell the main function about the max length of the result
** (for string functions), number of decimals (for double functions) and
** if the result may be a null value.
**
** If a function sets the 'error' argument to 1 the function will not be
** called anymore and mysqld will return NULL for all calls to this copy
** of the function.
**
** All strings arguments to functions are given as string pointer + length
** to allow handling of binary data.
** Remember that all functions must be thread safe. This means that one is not
** allowed to alloc any global or static variables that changes!
** If one needs memory one should alloc this in the init function and free
** this on the __deinit function.
**
** Note that the init and __deinit functions are only called once per
** SQL statement while the value function may be called many times
**
** Function 'metaphon' returns a metaphon string of the string argument.
** This is something like a soundex string, but it's more tuned for English.
**
** Function 'myfunc_double' returns summary of codes of all letters
** of arguments divided by summary length of all its arguments.
**
** Function 'myfunc_int' returns summary length of all its arguments.
**
** Function 'sequence' returns an sequence starting from a certain number.
** On the end is a couple of functions that converts hostnames to ip and
** vice versa.
**
** A dynamicly loadable file should be compiled shared.
** (something like: gcc -shared -o my_func.so myfunc.cc).
** You can easily get all switches right by doing:
** cd sql ; make udf_example.o
** Take the compile line that make writes, remove the '-c' near the end of
** the line and add -shared -o udf_example.so to the end of the compile line.
** The resulting library (udf_example.so) should be copied to some dir
** searched by ld. (/usr/lib ?)
** If you are using gcc, then you should be able to create the udf_example.so
** by simply doing 'make udf_example.so'.
**
** After the library is made one must notify mysqld about the new
** functions with the commands:
**
   CREATE AGGREGATE FUNCTION dpvar RETURNS INTEGER SONAME "dpvar.so"

	create table other.Jonathan_Mitogenes_snps_ldboundaries

	SELECT lower_bound, position, upper_bound FROM
		SELECT ldspline(A.pos2, A.dprime, 1, 0, A.pos1) as upper_bound,
			A.pos1 AS position FROM (
				SELECT * FROM LD.CEU inner join
					(select pos from other.Jonathan_Mitogenes_snps_10KB a inner join LD.index_2_rs b on a.rs_id = b.rs_id)
				as f where LD.CEU.pos1 = f.pos
			) AS A GROUP BY position) AS C
			NATURAL JOIN
		(SELECT B.pos2 AS position, ldspline(B.pos1, B.dprime, 1, 1, B.pos2) as lower_bound FROM (
			SELECT * FROM LD.CEU inner join
			(select pos from other.Jonathan_Mitogenes_snps_10KB a inner join LD.index_2_rs b on a.rs_id = b.rs_id)
			as g where LD.CEU.pos2 = g.pos
	) AS B GROUP BY position) AS D;
**
** After this the functions will work exactly like native MySQL functions.
** Functions should be created only once.
**
** The functions can be deleted by:
**
** DROP FUNCTION dpvar;
**
** The CREATE FUNCTION and DROP FUNCTION update the func@mysql table. All
** Active function will be reloaded on every restart of server
** (if --skip-grant-tables is not given)
**
** If you ge problems with undefined symbols when loading the shared
** library, you should verify that mysqld is compiled with the -rdynamic
** option.
**
** If you can't get AGGREGATES to work, check that you have the column
** 'type' in the mysql.func table.  If not, run 'mysql_fix_privilege_tables'.
**
*/

#include "msql.h"



/* These must be right or mysqld will not find the symbol! */
extern "C" my_bool ldspline_init( UDF_INIT* initid, UDF_ARGS* args, char* message );
extern "C" void ldspline_deinit( UDF_INIT* initid );
extern "C" void ldspline_reset( UDF_INIT* initid, UDF_ARGS* args, char* is_null, char *error );
extern "C" void ldspline_clear( UDF_INIT* initid, char* is_null, char *error );
extern "C" void ldspline_add( UDF_INIT* initid, UDF_ARGS* args, char* is_null, char *error );
extern "C" longlong ldspline( UDF_INIT* initid, UDF_ARGS* args, char* is_null, char *error );



