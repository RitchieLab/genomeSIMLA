/////////////////////////////////////////////////////////////////////////////
// Name:        genomesimla.h
// Purpose:     
// Author:      Eric Torstenson
// Modified by: 
// Created:     Thu 29 Nov 2007 02:29:43 PM CST
// RCS-ID:      
// Copyright:   Copyright 2007 Ritchie Lab - See COPYING for License 
// Licence:     
/////////////////////////////////////////////////////////////////////////////

#ifndef _GENOMESIMLA_H_
#define _GENOMESIMLA_H_

#if defined(__GNUG__) && !defined(__APPLE__)
#pragma interface "genomesimla.cpp"
#endif

/*!
 * Includes
 */
#include "appcontroller.h"
////@begin includes
#include "wx/image.h"
#include "mainframe.h"
////@end includes
#include <wx/apptrait.h>
#include <wx/stdpaths.h>

/*!
 * Forward declarations
 */

////@begin forward declarations
////@end forward declarations

/*!
 * Control identifiers
 */

////@begin control identifiers
////@end control identifiers


namespace GenomeSIM {
namespace GUI {
/*!
 * GenomeSIMLAApp class declaration
 */

class GenomeSIMLAApp: public wxApp
{    
    DECLARE_CLASS( GenomeSIMLAApp )
    DECLARE_EVENT_TABLE()

public:
    /// Constructor
    GenomeSIMLAApp();

    /// Initialises the application
    virtual bool OnInit();

    /// Called on exit
    virtual int OnExit();

////@begin GenomeSIMLAApp event handler declarations

////@end GenomeSIMLAApp event handler declarations

////@begin GenomeSIMLAApp member function declarations

////@end GenomeSIMLAApp member function declarations

////@begin GenomeSIMLAApp member variables
////@end GenomeSIMLAApp member variables

protected:
	AppController controller;
};

/*!
 * Application instance declaration 
 */
}
}
////@begin declare app
DECLARE_APP(GenomeSIM::GUI::GenomeSIMLAApp)
////@end declare app

#endif
    // _GENOMESIMLA_H_
