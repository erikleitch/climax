#include "gcp/fftutil/Generic1DModel.h"
#include "gcp/pgutil/PgUtil.h"

using namespace std;

using namespace gcp::util;

/**.......................................................................
 * Constructor.
 */
Generic1DModel::Generic1DModel() 
{
  dataSetType_ = DataSetType::DATASET_1D;
}

/**.......................................................................
 * Destructor.
 */
Generic1DModel::~Generic1DModel() 
{
  // Delete any thread pool resources that were allocated

  for(unsigned iThread=0; iThread < execData_.size(); iThread++) {
    delete execData_[iThread];
  }
}

/**.......................................................................
 * Evaluate this model at the specified x-coordinate
 */
double Generic1DModel::eval(double x)
{
  ThrowError("No inherited eval1D method has been defined for this model");
}

/**.......................................................................
 * Fill a 1D array with this model
 */
void Generic1DModel::fill1DArray(std::vector<double>& x, std::vector<double>& y)
{
  if(!pool_ || execData_.size() == 0) {
    fill1DArraySingleThread(x, y);
  } else {

    unsigned nThreadExec = initializeExecData(x, y);

    synchronizer_.reset(nThreadExec);

    for(unsigned iThread=0; iThread < nThreadExec; iThread++) {

      ExecData* ed = execData_[iThread];
      
      // Now execute it
      
      synchronizer_.registerPending(iThread);
      pool_->execute(&execFill1DArrayMultiThread, ed);
    }
    
    synchronizer_.wait();
  }
}

/**.......................................................................
 * Fill a 1D array with this model
 */
void Generic1DModel::fill1DArraySingleThread(std::vector<double>& x, std::vector<double>& y)
{
  for(unsigned i=0; i < x.size(); i++)
    y[i] = eval(x[i]);
}

/**.......................................................................
 * Fill a 1D array with this model
 */
void Generic1DModel::fill1DArrayMultiThread(ExecData* ed)
{
  for(ed->i_=ed->iStart_; ed->i_ < ed->iStop_; ed->i_++)
    ed->y_->at(ed->i_) = eval(ed->evalData_, ed->x_->at(ed->i_));
}

EXECUTE_FN(Generic1DModel::execFill1DArrayMultiThread)
{
  ExecData* ed = (ExecData*) args;
  Generic1DModel* model = ed->model_;

  //  COUT("Calling fillImage with start = " << ed->iYStart_ << " stop = " << ed->iYStop_);

  model->fill1DArrayMultiThread(ed);
  model->synchronizer_.registerDone(ed->iSegment_, ed->nSegment_);
}

/**.......................................................................
 * Fill a 1D array with this model, adding
 */
void Generic1DModel::generateFake1DData(std::string fileName, std::vector<double>& x, std::vector<double>& y, double sigma, std::string units)
{
  fill1DArray(x, y);
  
  PgUtil::linePlot(x,y);

  for(unsigned i=0; i < y.size(); i++) {
    y[i] += Sampler::generateGaussianSample(sigma);
  }

  std::ofstream fout;
  fout.open(fileName.c_str(), ios::out);

  if(!fout) {
    ThrowError("Unable to open file: " << fileName);
  }

  for(unsigned i=0; i < y.size(); i++)
    fout << x[i] << " " << y[i] << " " << sigma << std::endl;

  fout.close();
}

unsigned Generic1DModel::initializeExecData(std::vector<double>& x, std::vector<double>& y)
{
  unsigned npt    = x.size();

  unsigned nThreadTotal = pool_->nThread();
  unsigned nThreadExec = nThreadTotal > npt ? npt : nThreadTotal;
  unsigned nPtPerThread = npt / nThreadExec;
  unsigned iStart, iStop;
  
  for(unsigned iThread=0; iThread < nThreadExec; iThread++) {

    iStart = iThread * nPtPerThread;
    iStop  = iStart + nPtPerThread;

    if(iStop > npt)
      iStop = npt;

    ExecData* ed = execData_[iThread];

    ed->x_        = &x;
    ed->y_        = &y;
    ed->iSegment_ = iThread;
    ed->nSegment_ = nThreadExec;
    ed->iStart_   = iStart;
    ed->iStop_    = (iThread == nThreadExec-1 && iStop < npt) ? npt : iStop; // Last thread needs to finish

    initializeEvalData(ed->evalData_);
  }

  return nThreadExec;
}

void Generic1DModel::setThreadPool(ThreadPool* pool)
{
  // Call the base-class method

  Model::setThreadPool(pool);

  // ExecData::ExecData() can throw if inheritor hasn't define an
  // allocateEvalData() method.  This just means that the inheritor
  // doesn't support multi-threaded evaluation.  Quietly catch this
  // and the size of execData_ will be used later to determine if
  // multi-threaded evalulation can be used.

  if(pool) {
    try {
      
      // Now initialize per-thread resources once
      
      COUT("Number of threads = " << pool_->nThread());
      
      for(unsigned iThread=0; iThread < pool_->nThread(); iThread++)
	execData_.push_back(new ExecData(this));
      
    } catch(...) {
      
    }
  }
}

void* Generic1DModel::allocateEvalData()
{
  ThrowError("This inheritor has not defined a method to allocate evaluation data");
}

void Generic1DModel::initializeEvalData(void* evalData)
{
  ThrowError("This inheritor has not defined a method to initialize evaluation data");
}

double Generic1DModel::eval(void* evalData, double x)
{
  ThrowError("No envelope(void* evalData, double x) method has been defined for this inherited class");

}
