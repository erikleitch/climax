// $Id: $

#ifndef GCP_UTIL_SOLVER_H
#define GCP_UTIL_SOLVER_H

/**
 * @file Solver.h
 * 
 * Tagged: Mon Nov 25 00:20:30 PST 2013
 * 
 * @version: $Revision: $, $Date: $
 * 
 * @author Erik Leitch

 */
#define EVAL_FN(fn) void (*fn)(double x, std::vector<double>& pars, double& y, std::vector<double>& dydpars)

namespace gcp {
  namespace util {

    class Solver {
    public:

      /**
       * Constructor.
       */
      Solver(unsigned nPar);

      /**
       * Destructor.
       */
      virtual ~Solver();

      void setEvalFn(EVAL_FN(*fn));

    private:

      EVAL_FN(*evalFn_);

    }; // End class Solver

  } // End namespace util
} // End namespace gcp



#endif // End #ifndef GCP_UTIL_SOLVER_H
