/////////////////////////////////////////////////////////////////////////////
// Name:        genomesimla.cpp
// Purpose:     
// Author:      Eric Torstenson
// Modified by: 
// Created:     Thu 29 Nov 2007 02:29:43 PM CST
// RCS-ID:      
// Copyright:   Copyright 2007 Ritchie Lab - See COPYING for License 
// Licence:     
/////////////////////////////////////////////////////////////////////////////

#if defined(__GNUG__) && !defined(__APPLE__)
#pragma implementation "genomesimla.h"
#endif

// For compilers that support precompilation, includes "wx/wx.h".
#include "wx/wxprec.h"

#ifdef __BORLANDC__
#pragma hdrstop
#endif

#ifndef WX_PRECOMP
#include "wx/wx.h"
#endif

////@begin includes
////@end includes
#include <iostream>
#include <stdio.h>
#include "genomesimla.h"
#include "mainframe.h"
////@begin XPM images
////@end XPM images

using namespace std;
/*!
 * Application instance implementation
 */

////@begin implement app
IMPLEMENT_APP( GenomeSIM::GUI::GenomeSIMLAApp )
////@end implement app

/*!
 * GenomeSIMLAApp type definition
 */
IMPLEMENT_CLASS( GenomeSIM::GUI::GenomeSIMLAApp, wxApp )

/*!
 * GenomeSIMLAApp event table definition
 */
BEGIN_EVENT_TABLE( GenomeSIM::GUI::GenomeSIMLAApp, wxApp )

////@begin GenomeSIMLAApp event table entries
////@end GenomeSIMLAApp event table entries

END_EVENT_TABLE()


namespace GenomeSIM {
namespace GUI {
/*!
 * Constructor for GenomeSIMLAApp
 */

GenomeSIMLAApp::GenomeSIMLAApp()
{
////@begin GenomeSIMLAApp member initialisation
////@end GenomeSIMLAApp member initialisation
}

/*!
 * Initialisation for GenomeSIMLAApp
 */

bool GenomeSIMLAApp::OnInit()
{    
////@begin GenomeSIMLAApp initialisation
	// Remove the comment markers above and below this block
	// to make permanent changes to the code.

#if wxUSE_XPM
	wxImage::AddHandler(new wxXPMHandler);
#endif
#if wxUSE_LIBPNG
	wxImage::AddHandler(new wxPNGHandler);
#endif
#if wxUSE_LIBJPEG
	wxImage::AddHandler(new wxJPEGHandler);
#endif
#if wxUSE_GIF
	wxImage::AddHandler(new wxGIFHandler);
#endif
	string filename = "GenomeSIMLA";
	SetAppName(_("wxGenomeSIMLA"));
cerr<<"Application Name: "<<GetAppName()<<"\n";
	MainFrame* mainWindow = new MainFrame( NULL, ID_GENOMESIMLA_MAINFRAME, _(filename.c_str()));
	//, wxDefaultPosition, wxSize(-1, -1), wxMINIMIZE_BOX | wxMAXIMIZE_BOX | wxRESIZE_BORDER | wxSYSTEM_MENU | wxCAPTION | wxCLIP_CHILDREN );
	mainWindow->Show(true);
////@end GenomeSIMLAApp initialisation
	mainWindow->InitAppController(&controller);
	if (argc > 1) {
		filename = argv[1];

		try {
			controller.SetConfigurationFilename(filename.c_str());
			mainWindow->SetTitle(filename.c_str());
		}
		catch (Exception::FileNotFound& e) {
			stringstream err;
			err<<"The requested file, "<<filename<<", couldn't be opened. Please check that it is spelled correctly, including the correct path. ";

			wxMessageBox(err.str().c_str());
		}
	}

	mainWindow->RefreshSettings();
    return true;
}

/*!
 * Cleanup for GenomeSIMLAApp
 */
int GenomeSIMLAApp::OnExit()
{    
////@begin GenomeSIMLAApp cleanup
	return wxApp::OnExit();
////@end GenomeSIMLAApp cleanup
}
}
}
