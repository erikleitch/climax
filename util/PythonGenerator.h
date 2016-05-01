// $Id: PythonGenerator.h,v 1.1 2007/07/25 00:45:48 eml Exp $

#ifndef GCP_UTIL_PYTHONGENERATOR_H
#define GCP_UTIL_PYTHONGENERATOR_H

/**
 * @file PythonGenerator.h
 * 
 * Tagged: Tue Jul 24 11:13:44 PDT 2007
 * 
 * @version: $Revision: 1.1 $, $Date: 2007/07/25 00:45:48 $
 * 
 * @author Erik Leitch
 */

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

namespace gcp {
  namespace util {

    class PythonGenerator {
    public:

      /**
       * Constructor.
       */
      PythonGenerator(std::string fileName);
      PythonGenerator(std::string fileName, std::string dir);

      void setOutputPrefix(std::string prefix);
      void setOutputCcSuffix(std::string suffix);

      /**
       * Destructor.
       */
      virtual ~PythonGenerator();

      void outputCcFile();

    private:

      std::string outputPrefix_;
      std::string outputCcSuffix_;

      void parseFile();
      std::string caps(std::string inp);

      std::string sourceFileName_;
      std::string sourceFilePrefix_;
      std::string dir_;

      std::vector<std::string> functions_;
      std::vector<std::string> procRetVals_;
      std::vector<std::string> comments_;
      
      bool requiresNumpy_;

      static const bool print_info_on_import=false;

    }; // End class PythonGenerator

  } // End namespace util
} // End namespace gcp



#endif // End #ifndef GCP_UTIL_PYTHONGENERATOR_H
