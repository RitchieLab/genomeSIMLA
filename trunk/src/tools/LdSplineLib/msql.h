#ifndef RL_MSQL_H
#define RL_MSQL_H

#define MAX_ROWCOUNT 2048


#include <algo.h>


#ifdef WIN32
	// Insert your headers here
//	#define WIN32_LEAN_AND_MEAN		// Exclude rarely-used stuff from Windows headers

//	#include <windows.h>

	typedef unsigned __int64 ulonglong;	/* Microsofts 64 bit types */
	typedef __int64 		longlong;
	#define MutexType 		CRITICAL_SECTION
	#define InitMutex(M) 	InitializeCriticalSection(&M)
	#define DestroyMutex(M)	DeleteCriticalSection(&M)
	#define LockMutex(M) 	EnterCriticalSection(&M)
	#define UnLockMutex(M)	LeaveCriticalSection(&M)
	#define DllExport		__declspec( dllexport ) 

#ifndef ulonglong
	typedef unsigned long long ulonglong;
#endif //ulonglong

#ifndef longlong
	typedef long long 		longlong;
#endif	//longlong

#else
	#include <pthread.h>


	#define MutexType 		pthread_mutex_t
	#define InitMutex(M) 	pthread_mutex_init(&M, MY_MUTEX_INIT_SLOW)
	#define DestroyMutex(M)	pthread_mutex_destroy(&M)
	#define LockMutex(M) 	pthread_mutex_lock(&M)
	#define UnLockMutex(M)	pthread_mutex_unlock(&M)
	#define DllExport
#endif /*WIN32*/


#ifdef STANDARD
	/* STANDARD is defined, don't use any mysql functions */
	#include <stdlib.h>
	#include <stdio.h>
	#include <string.h>
	#ifdef WIN32
		typedef unsigned __int64 ulonglong;	/* Microsofts 64 bit types */
		typedef __int64 		longlong;
	#else
		typedef unsigned long long ulonglong;
		typedef long long 		longlong;
	#endif /*WIN32*/
#else
	#include <my_global.h>
	#include <my_sys.h>
	#if defined(MYSQL_SERVER)
		#include <m_string.h>		/* To get strmov() */





	#else
		/* when compiled as standalone */
		//#define strmov(a,b) strcpy(a,b)
		//#define bzero(a,b) memset(a,0,b)
		//#define memcpy_fixed(a,b,c) memcpy(a,b,c)
		#include <string.h>
	#endif
#endif
#include <mysql.h>
#include <ctype.h>



struct DataNode {
	double 		ldvalue;			///<The dprime/rsquared for a given SNP pairing
	longlong 	position;			///<The position associated with the other snp

	DataNode() : ldvalue(0.0), position(0) {}
	~DataNode() {}	
};



struct FnData {
	struct ltDataNode {
		bool operator()(const DataNode& l, const DataNode& r) {
			return l.position < r.position;
		}
	};
	
	struct gtDataNode {
		bool operator()(const DataNode& l, const DataNode &r) {
			return l.position > r.position;
		}
	};

	DataNode data[MAX_ROWCOUNT];					
	ulonglong count;				///<The number of records seen
	double threshold;				///<Threshold for restricting which nodes are to be used	
	double stringency;				///<Threshold for window stringency approach
	MutexType dataLock;				///<Let's try using mutexes to prevent order problems
	bool sortDescending;			///<Do we want to sort in descending order?
	ulonglong label;
	
	FnData() : count(0), threshold(0.0), label(0) {
		//Let's set up the lock....this is a mysql setting..not sure what it does
		InitMutex(dataLock);
	};
	~FnData() {
		DestroyMutex(dataLock);
	}
	bool AddData(uint *pos, double *ldvalue) {
		bool success = false;
		if (count < MAX_ROWCOUNT) {
			DataNode &d = data[count++];
			d.ldvalue = *ldvalue;
			d.position = *pos;	
			success = true;
		}
		return success;		
	}

	void Sort() {
		if (sortDescending) {
			sort(data, data+count, gtDataNode());
		} else {
			sort(data, data+count, ltDataNode());
		}
	}

	void Lock() {
		LockMutex(dataLock);
	}
	void UnLock() {
		UnLockMutex(dataLock);
	}

};
#endif //RL_MSQL_H

