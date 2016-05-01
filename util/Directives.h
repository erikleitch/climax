#ifndef GCP_UTIL_DIRECTIVES_H
#define GCP_UTIL_DIRECTIVES_H

// If we are using corba for any command or data transport, we need to
// include CORBA headers in our code

#define DIR_NEED_CORBA (ANT_CORBA | HAVE_CORR | HAVE_DC | HAVE_DELAY)

// But we can only use CORBA if both the CORBA libs and the carma tree
// exist

#define DIR_USE_CORBA (HAVE_CORBA & HAVE_CARMA & DIR_NEED_CORBA)

// See if we should use CORBA for antenna comms

#define DIR_USE_ANT_CORBA (HAVE_CORBA & HAVE_CARMA & ANT_CORBA)

// See if we should use various CARMA subsystems

#define DIR_USE_CORR (HAVE_CORBA & HAVE_CARMA & HAVE_CORR & HAVE_GCP)

#define DIR_USE_DC (HAVE_CORBA & HAVE_CARMA & HAVE_DC & HAVE_GCP)

#define DIR_USE_DELAY (HAVE_CORBA & HAVE_CARMA & HAVE_DELAY & HAVE_GCP)

#define DIR_USE_STRIP (HAVE_GCP)

#define DIR_USE_WX (HAVE_GCP & HAVE_WX)

#define DIR_HAVE_CARMA (HAVE_CARMA)

#define DIR_HAVE_CORR (HAVE_CARMA)

#define DIR_DEBUG (COMPILE_WITH_DEBUG)

#define DIR_HAVE_RT    (HAVE_RT)
 
#define DIR_HAVE_NUMPY (HAVE_NUMPY)

#endif // End #ifndef GCP_UTIL_DIRECTIVES_H
