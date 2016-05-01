// $Id: HubbleConstant.h,v 1.1 2012/05/02 23:44:51 eml Exp $

#ifndef GCP_UTIL_HUBBLECONSTANT_H
#define GCP_UTIL_HUBBLECONSTANT_H

/**
 * @file HubbleConstant.h
 * 
 * Tagged: Thu Apr  5 09:34:30 PDT 2012
 * 
 * @version: $Revision: 1.1 $, $Date: 2012/05/02 23:44:51 $
 * 
 * @author Erik Leitch
 */
#include "gcp/util/ConformableQuantity.h"
#include "gcp/util/Length.h"

namespace gcp {
  namespace util {

    class HubbleConstant : public ConformableQuantity {
    public:

      /**
       * Constructor.
       */
      HubbleConstant();

      /**
       * Destructor.
       */
      virtual ~HubbleConstant();

      void setKmPerSecPerMpc(double kmPerSecPerMpc)
      {
	val_ = kmPerSecPerMpc;
      }

      void setH100(double h100)
      {
	setKmPerSecPerMpc(h100 * 100);
      }

      inline double kmPerSecPerMpc() const {
	return val_;
      }

      inline double inverseSeconds() const {
	return val_ / Length::kmPerMpc_;
      }

      inline double h100() const {
	return kmPerSecPerMpc() / 100;
      }

      // Return H as a ratio to the specified kmPerSecPerMpc value

      inline double h(double kmPerSecPerMpcVal) const {
	return kmPerSecPerMpc() / kmPerSecPerMpcVal;
      }

      HubbleConstant operator*(double fac);

      double operator/(const HubbleConstant& hc);
      double operator/(HubbleConstant& hc);

      void operator=(const HubbleConstant& hc);
      void operator=(HubbleConstant& hc);

      friend std::ostream& operator<<(std::ostream& os, const HubbleConstant& hc);
      friend std::ostream& operator<<(std::ostream& os, HubbleConstant& hc);

    }; // End class HubbleConstant

    std::ostream& operator<<(std::ostream& os, const HubbleConstant& hc);
    std::ostream& operator<<(std::ostream& os, HubbleConstant& hc);

  } // End namespace util
} // End namespace gcp



#endif // End #ifndef GCP_UTIL_HUBBLECONSTANT_H
