//
// C++ Interface: wxgridloci
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef GENOMESIM_GUIWXGRIDLOCI_H
#define GENOMESIM_GUIWXGRIDLOCI_H

#include <wx/grid.h>
#include <vector>
#include "simulation/locus.h"

namespace GenomeSIM {

namespace GUI {


using namespace std;
using namespace Simulation;

/**
&brief interfaces with the grid class to allow the rendering of huge lists of loci. 

	@author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>
*/
class wxGridLoci : public wxGridTableBase {
public:
    wxGridLoci(vector<Locus *> *loci, bool readOnly = false);
    ~wxGridLoci();

    virtual int GetNumberRows();
    virtual int GetNumberCols();
    virtual bool IsEmptyCell( int row, int col );
    virtual wxString GetValue( int row, int col );
    virtual void SetValue( int row, int col, const wxString& value );

    virtual wxString GetColLabelValue( int col );

    virtual wxString GetTypeName( int row, int col );
    virtual bool CanGetValueAs( int row, int col, const wxString& typeName );
    virtual bool CanSetValueAs( int row, int col, const wxString& typeName );

    virtual long GetValueAsLong( int row, int col );
    virtual double GetValueAsDouble( int row, int col );

    virtual void SetValueAsLong( int row, int col, long value );
    virtual void SetValueAsDouble( int row, int col, double value );

	//virtual wxGridCellAttr *GetAttr(int row, int col);

//	void LoadFile(const char *filename);
	void Reset();
	void Refresh();
	void ClearTable(size_t lociCount = 0);

	/**
	 * @brief Returns whether or not the chromosome has been changed
	 */
	bool HasChanged() { return hasChanged; }
	
	/**
	 * @brief Indicates that this is the saved state. Any changes from here on will be recorded
	 */
	void SetStableState() { hasChanged= false; }
protected:
	bool hasChanged;
	vector<Locus *> *loci;
	bool readOnly;
	//wxGridCellAttr *attr;


};


inline
wxGridLoci::wxGridLoci(vector<Locus *> *loci, bool readOnly) : hasChanged(false), loci(loci), readOnly(readOnly)  { 
}

inline
wxGridLoci::~wxGridLoci() { 
	//delete attr;
}

}

}

#endif
