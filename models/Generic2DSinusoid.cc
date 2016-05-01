#include "gcp/models/Generic2DSinusoid.h"

using namespace std;

using namespace gcp::util;

/**.......................................................................
 * Constructor.
 */
Generic2DSinusoid::Generic2DSinusoid() 
{
  addComponent(wavelength_);
  addComponent(phase_);

  addComponentName(wavelength_, "wave",  "The angular wavelength of the sinusoid");
  addComponentName(phase_,      "phase", "The phase of the sinusoid");

  initializeComponentsToFixed();
}

/**.......................................................................
 * Destructor.
 */
Generic2DSinusoid::~Generic2DSinusoid() {}

/**.......................................................................
 * Evaluate the dimensionless shape of this profile
 */
double Generic2DSinusoid::envelope(unsigned type, double xRad, double yRad)
{
  double waveRad  = wavelength_.radians();
  double phaseRad = phase_.radians();

  if(fabs(waveRad) > 0.0)
    return sin(2*M_PI*(xRad/waveRad) + phaseRad);
  else
    return 1.0;
    
}
