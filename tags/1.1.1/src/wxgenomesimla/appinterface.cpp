#include "appinterface.h"
#include "genomesim/timestamp.h"

namespace GenomeSIM {
namespace GUI {
AppInterface::Version::Version() : 
	maj(APPMAJOR),
  min(APPMINOR),
  bug(APPBUGFIX),
  build(BUILD_NUMBER),
  date(BUILD_DATE) { }
}
}
 
