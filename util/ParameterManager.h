// $Id: $

#ifndef GCP_UTIL_PARAMETERMANAGER_H
#define GCP_UTIL_PARAMETERMANAGER_H

/**
 * @file ParameterManager.h
 * 
 * Tagged: Thu Jul 26 10:03:25 PDT 2012
 * 
 * @version: $Revision: $, $Date: $
 * 
 * @author Erik Leitch
 */
#include "gcp/util/DataType.h"
#include "gcp/util/String.h"

#include <map>

namespace gcp {
  namespace util {

    class ParameterManager {
    public:

      // Struct to encapslate a single parameter

      struct Parameter {

	DataType data_;
	std::string units_;
	std::string comment_;
	bool resizable_;
	ParameterManager* owner_;

	Parameter() {
	  units_ = " ";
	  resizable_ = false;
	  owner_ = 0;
	};

	~Parameter() {
	  deleteResizable();
	};

	void setResizable(bool resizable);
	void deleteResizable();
	void setComment(std::string comment);

	Parameter(const Parameter& param) 
	{
	  *this = param;
	}

	Parameter(Parameter& param) 
	{
	  *this = param;
	}

	void operator=(const Parameter& param)
	{
	  *this = (Parameter&)param;
	}

	void operator=(Parameter& param)
	{
	  data_    = param.data_;
	  units_   = param.units_;
	  comment_ = param.comment_;
	}

      };

      /**
       * Constructor.
       */
      ParameterManager();

      /**
       * Destructor.
       */
      virtual ~ParameterManager();

      //------------------------------------------------------------
      // Methods dealing with parameters
      //------------------------------------------------------------

      // Add a parameter to our map of known parameters

      void addParameter(std::string name, DataType::Type type, std::string comment=" ", bool resizable=false);
      void remParameter(std::string name);

      // Add another parameter manager to our hierarchy of parameters

      void addParameter(ParameterManager& pm);
      void addParameter(std::string name, ParameterManager& pm);
      void updateParameter(std::string name, ParameterManager& pm);

      // Set methods

      virtual void setParameter(std::string name, std::string val, std::string units=" ", bool external=true);
      virtual void incrementParameter(std::string name, std::string val, std::string units=" ", bool external=true);

      // Get methods
      
      Parameter* getParameter(std::string name, bool checkVal=false);
      void listParameters(std::ostringstream& os, bool sort=true);
      void formatParameter(std::ostringstream& os, std::string name, std::string comment, 
			   unsigned maxLen, unsigned lineWidth);
      void appendLineAndWord(std::ostringstream& osLine, unsigned& nInLine, 
			     std::ostringstream& osWord, unsigned& nInWord, 
			     std::ostringstream& filler, std::vector<std::string>& lines,
			     unsigned lineWidth, bool atLineBreak);

      // Copy parameters from another parameter manager.
      //
      // If all == true,   any parameter not excluded will be copied
      // If all == false, only unspecified parameters will be copied
      //
      // We make this virtual so that inheritors that include other
      // ParameterManager objects can, in principle, copy those
      // parameters too when this method is called

      virtual void copyParameters(ParameterManager* pm, std::map<std::string, std::string>& excludedParameters, bool all);
      virtual void setParameters(ParameterManager* pm, std::map<std::string, std::string>& excludedParameters, bool all);

      // Return true if this name matches any parameter (minimum matching allowed)

      bool matches(std::string name);

      // Return true if this name unique matches any parameter (minimum matching allowed)

      bool unique(std::string name);

      // Return the parameter that matches this name

      std::string getMatch(std::string name);

      // Accessor methods for parameter values

      std::string getStringVal(std::string name);
      bool getBoolVal(std::string name);
      unsigned char getUcharVal(std::string name);
      char getCharVal(std::string name);
      unsigned short getUshortVal(std::string name);
      short getShortVal(std::string name);
      unsigned int getUintVal(std::string name);
      int getIntVal(std::string name);
      unsigned long getUlongVal(std::string name);
      long getLongVal(std::string name);
      float getFloatVal(std::string name);
      double getDoubleVal(std::string name);

      std::vector<double> parseRange(String& val);

      // A map of parameters used by this dataset

      std::map<std::string, Parameter*> parameterMap_;
      std::map<Parameter*, std::string> reverseParameterMap_;
      std::vector<Parameter*> parameterVec_;

      std::string name_;

    }; // End class ParameterManager

  } // End namespace util
} // End namespace gcp



#endif // End #ifndef GCP_UTIL_PARAMETERMANAGER_H
