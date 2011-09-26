#include "ldspline.h"
#include <string.h>

/*
** Average Cost Aggregate Function.
*/
my_bool
ldspline_init( UDF_INIT* initid, UDF_ARGS* args, char* message )
{
 	//Let's make sure that we have all the pieces required to perform the conversion
	if (args->arg_count != 5) {
		strcpy(message, "Invalid arg list: (position, ldvalue, threshold, isDescending, label)");
		return 1;
	}

	//Verify that the types are appropriate!
	if (args->arg_type[0] != INT_RESULT) {
		strcpy(message, "position (arg 0) must be of type integer.");
		return 1;
	}
	if (args->arg_type[1] != REAL_RESULT) {
		strcpy(message, "ldvalue (arg 1 ) must be a float.");
		return 1;
	}
	if (args->arg_type[4] != INT_RESULT) {
		strcpy(message, "ldvalue (arg 4) must be of type integer.");
		return 1;
	}

	//We need to make sure this is a real value, not a string
	args->arg_type[2] = REAL_RESULT;

	//We need to make sure that we interpret the asc/desc variable correctly (0/1)
	args->arg_type[3] = INT_RESULT;

	initid->maybe_null	= 0;		/* The result may be null */
	initid->decimals	= 4;		/* We want 4 decimals in the result */
	initid->max_length	= 20;		/* 6 digits + . + 10 decimals */

	initid->ptr = (char *)new FnData;
  	return 0;
}

void
ldspline_deinit( UDF_INIT* initid ) {
 	FnData* data = (FnData*)initid->ptr;
	delete data;
}


/* This is only for MySQL 4.0 compability */
void
ldspline_reset(UDF_INIT* initid, UDF_ARGS* args, char* is_null, char* message) {
  	ldspline_clear(initid, is_null, message);
  	ldspline_add(initid, args, is_null, message);
}

/* This is needed to get things to work in MySQL 4.1.1 and above */

void
ldspline_clear(UDF_INIT* initid, char* is_null __attribute__((unused)),
              char* message __attribute__((unused)))		
{
 	FnData* data = (FnData*)initid->ptr;
	data->count = 0;	
	data->label = 0;
}


void
ldspline_add(UDF_INIT* initid, UDF_ARGS* args,
            char* is_null __attribute__((unused)),
            char* message __attribute__((unused)))		{
	if (args->args[0] && args->args[1]) 	{
		FnData* data = (FnData*)initid->ptr;
		data->Lock();
	
		//Set the threshold from the first entry
		if (data->count == 0) {
			data->label = *((ulonglong*)args->args[4]);
			if (args->args[2])
				data->threshold = *((double*)args->args[2]);
			if (args->args[3])
				data->sortDescending= *((uint*)args->args[3]) == 1;
		}

		//Set the add the datapoint
		if (data->count < MAX_ROWCOUNT) 
			data->AddData((uint*)args->args[0], (double*)args->args[1]);
		
		data->UnLock();	
	}

}


longlong
ldspline( UDF_INIT* initid, 
		UDF_ARGS* args __attribute__((unused)),
        char* is_null, 
		char* error __attribute__((unused)))	{
	FnData* data = (FnData*)initid->ptr;

 	if (!data->count)   {
		*is_null = 1;
    	return 0;
  	}

	data->Lock();
	
	*is_null = 0;
	double thresh 	= data->threshold;
	int rowCount 	= data->count;
	double posDiff = 0.0;

	data->Sort();

	longlong furthestSnp = data->label;
	longlong nextSnp = furthestSnp;

	bool isFound = true;

	for (int i=0; i<rowCount; i++) {
		furthestSnp = nextSnp;
		nextSnp = data->data[i].position;

		if (data->data[i].ldvalue < thresh) {
			isFound=false;
			break;
		}
		
	}
	if (isFound)
		furthestSnp=nextSnp;
	

	data->UnLock();
	return furthestSnp;
}
