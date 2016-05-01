/*
   gcpMatReadMiriad.cc
   Tom Plagge 12/2/09

   Reads a Miriad UV file into Matlab.
 
   12/2/09 -- First committed version

*/

#include <iostream>
#include <math.h>
#include <vector>
#include <string>
#include <algorithm>

#include "mex.h"
#include "matrix.h"
#include "gcp/miriad/miriad.h"

using namespace std;

#define READMIRIAD_DEBUG 0

/* ----------------- 
   Class definitions
   ----------------- */

// We'll consider all Miriad variables
//  to be either strings, integers, or doubles.
//  Really this, like the frame class below,
//  is just a struct.  Maybe later we'll need
//  more functionality, though.
class mirvar {
 public:
  string name;
  vector<string> sval;
  vector<int>    ival;
  vector<double> dval;
};

// A "frame" for our purposes is all the data
//  collected and timestamped with a given mjd.
//  We will keep a list of variables to store
//  on a frame-by-frame basis.  Note that in the
//  Miriad file, these variables can be updated
//  at any time, not necessarily on frame boundaries.
class frame {
 public:
  vector<mirvar> vars;
};

// All the data from a Miriad file is stored in
//  this class. 
class mirfile {
 public:
  mirfile(string filename);
  mxArray *get_mxarray();

 private:
  bool get_miriad_var(vector<mirvar> &v, string name);
  bool fill_miriad_struct();
  bool grok_miriad_file();
  string m_filename;
  vector<mirvar> m_vars;
  vector<frame> m_frames;
  int m_uvh;
  int m_nbase;
  int m_nants;
  int m_nspect;
  int m_nwide;
  int m_nchan;
  int m_nframe;
  vector<int> m_nschan;
  vector<int> m_baselines;
  vector<double> m_mjds;
  int m_nvis;
  vector<vector<vector<vector<double> > > > m_visSpecRe;
  vector<vector<vector<vector<double> > > > m_visSpecIm;
  vector<vector<vector<vector<int> > > > m_flagSpec;
  vector<vector<vector<double> > > m_visWideRe;
  vector<vector<vector<double> > > m_visWideIm;
  vector<vector<vector<int> > > m_flagWide;
  vector<vector<double> > m_u;
  vector<vector<double> > m_v;
  vector<vector<double> > m_w;
};

/* --------------------------- 
   Member function definitions
   --------------------------- */

mirfile::mirfile(string filename) :
  m_filename(filename) {
  fill_miriad_struct();
}

mxArray *mirfile::get_mxarray() {
// Fill the matlab array with the
// data we've read in.

  // Create the struct to return
  int dims[2]={1,1};
  mxArray *retval=mxCreateStructArray(2,dims,0,NULL);

  // Add the wideband visibilities and flags if they exist
  if (m_nwide > 0) {
    mxAddField(retval,"visw");
    mxAddField(retval,"flagw");
    int viswdims[3]={m_visWideRe.size(),
      m_visWideRe[0].size(),m_visWideRe[0][0].size()};
    mxArray *visbrickw = mxCreateNumericArray(3,viswdims,mxDOUBLE_CLASS,mxCOMPLEX);
    mxArray *flagbrickw = mxCreateNumericArray(3,viswdims,mxINT32_CLASS,mxREAL);
    double *realPartW=(double*)mxGetData(visbrickw);
    double *imagPartW=(double*)mxGetImagData(visbrickw);
    int *flagW=(int*)mxGetData(flagbrickw);
    int iiw[3];
    for (iiw[0]=0; iiw[0]<viswdims[0]; iiw[0]++) {
      for (iiw[1]=0; iiw[1]<viswdims[1]; iiw[1]++) {
        for (iiw[2]=0; iiw[2]<viswdims[2]; iiw[2]++) {
          *(realPartW + mxCalcSingleSubscript(visbrickw,3,iiw))=
            m_visWideRe[iiw[0]][iiw[1]][iiw[2]];
          *(imagPartW + mxCalcSingleSubscript(visbrickw,3,iiw))=
            m_visWideIm[iiw[0]][iiw[1]][iiw[2]];
          *(flagW     + mxCalcSingleSubscript(flagbrickw,3,iiw))=
            m_flagWide[iiw[0]][iiw[1]][iiw[2]];
        }
      }
    }
    mxSetField(retval,0,"visw",visbrickw);
    mxSetField(retval,0,"flagw",flagbrickw);
  }

  // Add the spectral visibilities and flags if they exist
  if (m_nspect > 0) {
    mxAddField(retval,"vis");
    mxAddField(retval,"flag");
    int maxschans=0;
    for (int i=0; i<m_visSpecRe[0][0].size(); i++) {
      if (m_visSpecRe[0][0][i].size() > maxschans)
        maxschans=m_visSpecRe[0][0][i].size();
    }
    int visdims[4]={m_visSpecRe.size(),m_visSpecRe[0].size(),m_visSpecRe[0][0].size(),
                    maxschans};
    mxArray *visbrick = mxCreateNumericArray(4,visdims,mxDOUBLE_CLASS,mxCOMPLEX);
    mxArray *flagbrick = mxCreateNumericArray(4,visdims,mxINT32_CLASS,mxREAL);
    double *realPart=(double*)mxGetData(visbrick);
    double *imagPart=(double*)mxGetImagData(visbrick);
    int *flag=(int*)mxGetData(flagbrick);
    int ii[4];
    for (ii[0]=0; ii[0]<visdims[0]; ii[0]++) {
      for (ii[1]=0; ii[1]<visdims[1]; ii[1]++) {
        for (ii[2]=0; ii[2]<visdims[2]; ii[2]++) {
          for (ii[3]=0; ii[3]<visdims[3]; ii[3]++) {
            if (ii[3] < m_visSpecRe[0][0][ii[2]].size()) {
              *(realPart + mxCalcSingleSubscript(visbrick,4,ii))=
                m_visSpecRe[ii[0]][ii[1]][ii[2]][ii[3]];
              *(imagPart + mxCalcSingleSubscript(visbrick,4,ii))=
                m_visSpecIm[ii[0]][ii[1]][ii[2]][ii[3]];
              *(flag     + mxCalcSingleSubscript(flagbrick,4,ii))=
                m_flagSpec[ii[0]][ii[1]][ii[2]][ii[3]];
            } else {
              *(realPart + mxCalcSingleSubscript(visbrick,4,ii))=(double)0.0;
              *(imagPart + mxCalcSingleSubscript(visbrick,4,ii))=(double)0.0;
              *(flag     + mxCalcSingleSubscript(flagbrick,4,ii))=0;
            }
          }
        }
      }
    }
    mxSetField(retval,0,"vis",visbrick);
    mxSetField(retval,0,"flag",flagbrick);

    mxAddField(retval,"u");
    mxAddField(retval,"v");
    mxAddField(retval,"w");
    int uvwdims[2]={m_u.size(),m_u[0].size()};
    mxArray *ubrick = mxCreateNumericArray(2,uvwdims,mxDOUBLE_CLASS,mxREAL);
    mxArray *vbrick = mxCreateNumericArray(2,uvwdims,mxDOUBLE_CLASS,mxREAL);
    mxArray *wbrick = mxCreateNumericArray(2,uvwdims,mxDOUBLE_CLASS,mxREAL);
    double *udata=(double*)mxGetData(ubrick);
    double *vdata=(double*)mxGetData(vbrick);
    double *wdata=(double*)mxGetData(wbrick);
    int iiuvw[2];
    for (iiuvw[0]=0; iiuvw[0]<uvwdims[0]; iiuvw[0]++) {
      for (iiuvw[1]=0; iiuvw[1]<uvwdims[1]; iiuvw[1]++) {
        *(udata + mxCalcSingleSubscript(ubrick,2,iiuvw))=
          m_u[iiuvw[0]][iiuvw[1]];
        *(vdata + mxCalcSingleSubscript(vbrick,2,iiuvw))=
          m_v[iiuvw[0]][iiuvw[1]];
        *(wdata + mxCalcSingleSubscript(wbrick,2,iiuvw))=
          m_w[iiuvw[0]][iiuvw[1]];
      }
    }
    mxSetField(retval,0,"u",ubrick);
    mxSetField(retval,0,"v",vbrick);
    mxSetField(retval,0,"w",wbrick);
  }

  // Add global data
  for (unsigned i=0; i<m_vars.size(); i++) {
    mxAddField(retval,m_vars[i].name.c_str());
    if (m_vars[i].sval.size() != 0) {
      // String-valued variables (no arrays of strings allowed for now)
      dims[1]=1;
      mxSetField(retval,0,m_vars[i].name.c_str(),
        mxCreateString(m_vars[i].sval[0].c_str()));
    } else if (m_vars[i].ival.size() != 0) {
      // Integer-valued variables
      mxArray *thisvar=NULL;
      dims[1]=m_vars[i].ival.size();
      thisvar=mxCreateNumericArray(2,dims,mxINT32_CLASS, mxREAL);
      int *thisptr=(int*)mxGetData(thisvar);
      for (unsigned k=0; k<dims[1]; k++) 
        *(thisptr + k) = (int)m_vars[i].ival[k];
      mxSetField(retval,0,m_vars[i].name.c_str(),thisvar);
    } else if (m_vars[i].dval.size() != 0) {
      // Floating point-valued variables
      mxArray *thisvar=NULL;
      dims[1]=m_vars[i].dval.size();
      thisvar=mxCreateNumericArray(2,dims,mxDOUBLE_CLASS, mxREAL);
      double *thisptr=(double*)mxGetData(thisvar);
      for (unsigned k=0; k<dims[1]; k++)
        *(thisptr + k) = (double)m_vars[i].dval[k];
      mxSetField(retval,0,m_vars[i].name.c_str(),thisvar);
    } else {
      // Empty variables (don't bother with them)
      continue;
    }
  }

  // Create a struct array for frame data (timestamp-by-timestamp variables)
  mxAddField(retval,"frame");
  dims[1]=m_frames.size();
  mxArray *framestruct=mxCreateStructArray(2,dims,0,NULL);

  // Loop over all of the frames in the miriad file
  for (unsigned i=0; i<m_frames.size(); i++) {
    // Loop over all of the variables
    for (unsigned j=0;j<m_frames[i].vars.size(); j++) {
      // If the field name doesn't exist, create it.
      if (!mxGetField(framestruct,0,m_frames[i].vars[j].name.c_str()))
        mxAddField(framestruct,m_frames[i].vars[j].name.c_str());

      if (m_frames[i].vars[j].sval.size() != 0) {
        // String-valued variables (no arrays of strings allowed for now)
        dims[1]=1;
        mxSetField(framestruct,i,m_frames[i].vars[j].name.c_str(),
          mxCreateString(m_frames[i].vars[j].sval[0].c_str()));
      } else if (m_frames[i].vars[j].ival.size() != 0) {
        // Integer-valued variables
        mxArray *thisvar=NULL;
        dims[1]=m_frames[i].vars[j].ival.size();
        thisvar=mxCreateNumericArray(2,dims,mxINT32_CLASS, mxREAL);
        void *thisvptr=mxGetData(thisvar);
        int *thisptr=(int*)thisvptr;
        for (unsigned k=0; k<dims[1]; k++) 
          *(thisptr + k) = (int)m_frames[i].vars[j].ival[k];
        mxSetField(framestruct,i,m_frames[i].vars[j].name.c_str(),thisvar);
      } else if (m_frames[i].vars[j].dval.size() != 0) {
        // Floating point-valued variables
        mxArray *thisvar=NULL;
        dims[1]=m_frames[i].vars[j].dval.size();
        thisvar=mxCreateNumericArray(2,dims,mxDOUBLE_CLASS, mxREAL);
        void *thisvptr=mxGetData(thisvar);
        double *thisptr=(double*)thisvptr;
        for (unsigned k=0; k<dims[1]; k++)
          *(thisptr + k) = (double)m_frames[i].vars[j].dval[k];
        mxSetField(framestruct,i,m_frames[i].vars[j].name.c_str(),thisvar);
      } else {
        // Empty variables (don't bother with them)
        continue;
      }
    }
  }
  // Add the frame data to the returned structure
  mxSetField(retval,0,"frame",framestruct);
  return retval;
}

bool mirfile::grok_miriad_file() {
// Parse the miriad file for sizes, number of channels, and
// crap like that.
  uvopen_c(&(m_uvh), m_filename.c_str(), "old"); 
  uvnext_c(m_uvh);

  // Some day, these values won't be hard-coded in and
  //  you'll be able to request the members you need.  Maybe.
  vector<string> varnames;
  varnames.push_back("epoch"); varnames.push_back("sfreq"); 
  varnames.push_back("freq"); varnames.push_back("wfreq"); 
  varnames.push_back("restfreq"); varnames.push_back("wwidth");
  varnames.push_back("nants"); varnames.push_back("nchan");
  varnames.push_back("npol"); varnames.push_back("nschan");
  varnames.push_back("nspect"); varnames.push_back("nwide");
  varnames.push_back("longitu"); varnames.push_back("latitud");
  varnames.push_back("telescop"); varnames.push_back("version");
  for(vector<string>::iterator i=varnames.begin(); i!=varnames.end(); i++) {
    if (!get_miriad_var(m_vars,*i)) {varnames.erase(i);i--;} 
  }

  // get some relevant parameters
  unsigned inants=find(varnames.begin(),varnames.end(),"nants")-varnames.begin();
  m_nants=m_vars[inants].ival[0];
  int nbase_max=(m_nants-1)*(m_nants)/2;

  unsigned inspect=find(varnames.begin(),varnames.end(),"nspect")-varnames.begin();
  m_nspect=(inspect>0 && inspect<varnames.size())? m_vars[inspect].ival[0] : 0;

  unsigned inwide=(find(varnames.begin(),varnames.end(),"nwide"))-varnames.begin();
  m_nwide=(inwide>0 && inwide<varnames.size())? m_vars[inwide].ival[0] : 0;

  unsigned inschan=(find(varnames.begin(),varnames.end(),"nschan"))-varnames.begin();
  m_nschan=(m_nspect>0)? m_vars[inschan].ival : vector<int>(0);

  unsigned inchan=(find(varnames.begin(),varnames.end(),"nchan"))-varnames.begin();
  m_nchan=(inchan>0 && inchan<varnames.size())? m_vars[inchan].ival[0] : 0;

  // Store the characteristics of the baselines and frames.
  mirvar a1, a2, base,mjds;
  a1.name="ant1";
  a2.name="ant2";
  base.name="baseline";
  mjds.name="mjd";
  float  bufData[2*m_nchan];
  int    bufFlags[m_nchan];
  float  bufSysTemp[m_nants*m_nspect];
  double bufPreamble[5];

  // Parse through the file looking for unique mjd's and baselines.
  uvopen_c(&(m_uvh), m_filename.c_str(), "old"); 
  uvnext_c(m_uvh);
  m_nvis=0;
  while(1) {
    int nread;
    uvset_c(m_uvh, "data",     "channel", 0, 1.0, 1.0, 1.0);
    uvset_c(m_uvh, "preamble", "uvw/time/baseline", 0, 0.0, 0.0, 0.0);
    uvread_c(m_uvh, bufPreamble, bufData, bufFlags, m_nchan, &nread);
    if (!nread) break;
    m_nvis++;
    unsigned thisbase=(unsigned)(bufPreamble[4]);
    unsigned thisa2 = (thisbase & 255)-1;
    unsigned thisa1 = (thisbase-thisa2)/256-1;
    if (find(base.ival.begin(),base.ival.end(),bufPreamble[4])==base.ival.end()) {
      a1.ival.push_back((int)thisa1+1);
      a2.ival.push_back((int)thisa2+1);
      base.ival.push_back((int)bufPreamble[4]);
    }
    if (find(mjds.dval.begin(),mjds.dval.end(),bufPreamble[3])==mjds.dval.end()) {
      mjds.dval.push_back(bufPreamble[3]);
    }
  }
  m_baselines=base.ival;
  m_nbase=base.ival.size();
  m_mjds=mjds.dval;
  m_nframe=mjds.dval.size();
  m_vars.push_back(base);
  m_vars.push_back(a1);
  m_vars.push_back(a2);
  m_vars.push_back(mjds);

  uvflush_c(m_uvh);
  uvclose_c(m_uvh);
  
  // Resize the data arrays
  m_visSpecRe.resize(m_nframe);
  m_visSpecIm.resize(m_nframe);
  m_visWideRe.resize(m_nframe);
  m_visWideIm.resize(m_nframe);
  m_flagSpec.resize(m_nframe);
  m_flagWide.resize(m_nframe);
  m_u.resize(m_nframe);
  m_v.resize(m_nframe);
  m_w.resize(m_nframe);
  m_frames.resize(m_nframe);
  for (int iframe=0;iframe<m_nframe;iframe++) {
    m_visSpecRe[iframe].resize(m_nbase);
    m_visSpecIm[iframe].resize(m_nbase);
    m_visWideRe[iframe].resize(m_nbase);
    m_visWideIm[iframe].resize(m_nbase);
    m_flagSpec[iframe].resize(m_nbase);
    m_flagWide[iframe].resize(m_nbase);
    m_u[iframe].resize(m_nbase);
    m_v[iframe].resize(m_nbase);
    m_w[iframe].resize(m_nbase);
    for (int ibase=0;ibase<m_nbase;ibase++) {
      m_visSpecRe[iframe][ibase].resize(m_nspect);
      m_visSpecIm[iframe][ibase].resize(m_nspect);
      m_visWideRe[iframe][ibase].resize(m_nwide);
      m_visWideIm[iframe][ibase].resize(m_nwide);
      m_flagSpec[iframe][ibase].resize(m_nspect);
      m_flagWide[iframe][ibase].resize(m_nwide);
      for (int ispect=0;ispect<m_nspect;ispect++) {
        m_visSpecRe[iframe][ibase][ispect].resize(m_nschan[ispect]);
        m_visSpecIm[iframe][ibase][ispect].resize(m_nschan[ispect]);
        m_flagSpec[iframe][ibase][ispect].resize(m_nschan[ispect]);
      }
    }
  }
}


bool mirfile::fill_miriad_struct() {
// Read the data from the miriad file into the mirfile class

  // First figure out how big everything is and set it all up.
  grok_miriad_file();
  
  // Now read the data.
  float  bufData[2*m_nchan];
  float  bufWideData[2*m_nwide];
  int    bufFlags[m_nchan];
  int    bufWideFlags[m_nwide];
  float  bufSysTemp[m_nants*m_nspect];
  double bufPreamble[5];
  double curmjd=0.;
  int    curframe=-1;
  int    nread;

  uvopen_c(&m_uvh, m_filename.c_str(), "old"); 
  //uvnext_c(m_uvh);

  // Read the first visibility
  uvset_c(m_uvh, "data",     "channel", 0, 1.0, 1.0, 1.0);
  uvset_c(m_uvh, "preamble", "uvw/time/baseline", 0, 0.0, 0.0, 0.0);
  uvread_c(m_uvh, bufPreamble, bufData, bufFlags, m_nchan, &nread);
  if (!nread) {
    if (READMIRIAD_DEBUG) cout << "Empty file!" << endl;
    return false;
  }
  if (m_nwide) {
    uvset_c(m_uvh, "data",     "wide", 0, 1.0, 1.0, 1.0);
    uvwread_c(m_uvh, bufWideData, bufWideFlags, m_nwide, &nread);
  }

  int nvis=0;
  // Loop until the file has ended
  if (READMIRIAD_DEBUG) cout << "|----------|" << endl << " ";
  while(1) {
    if (READMIRIAD_DEBUG) {if (nvis % m_nvis/10==0) cout << "*";}
    if (bufPreamble[3]!=curmjd) {
      // A new frame has started
      curmjd=bufPreamble[3];
      curframe++;
      vector<string> frame_varnames;
      // Some day, these values won't be hard-coded in and
      //  you'll be able to request the members you need.  Maybe.
      frame_varnames.push_back("airtemp"); frame_varnames.push_back("ambpsys");
      frame_varnames.push_back("antaz"); frame_varnames.push_back("antel");
      frame_varnames.push_back("antpos"); frame_varnames.push_back("axisrms");
      frame_varnames.push_back("cable"); frame_varnames.push_back("chi");
      frame_varnames.push_back("dazim"); frame_varnames.push_back("ddec");
      frame_varnames.push_back("dec"); frame_varnames.push_back("delay");
      frame_varnames.push_back("delev"); frame_varnames.push_back("dra");
      frame_varnames.push_back("epoch"); 
      frame_varnames.push_back("inttime"); frame_varnames.push_back("ischan");
      frame_varnames.push_back("jyperk"); frame_varnames.push_back("jyperka");
      frame_varnames.push_back("lo1"); frame_varnames.push_back("modedesc");
      frame_varnames.push_back("obsdec"); frame_varnames.push_back("obsra");
      frame_varnames.push_back("on"); frame_varnames.push_back("pamatten");
      frame_varnames.push_back("phaselo1"); frame_varnames.push_back("phaselo2");
      frame_varnames.push_back("phasem1"); frame_varnames.push_back("plangle"); 
      frame_varnames.push_back("plmaj"); frame_varnames.push_back("plmin"); 
      frame_varnames.push_back("pltb"); frame_varnames.push_back("pntdec"); 
      frame_varnames.push_back("pntra"); frame_varnames.push_back("pol");
      frame_varnames.push_back("precipmm"); frame_varnames.push_back("pressmb");
      frame_varnames.push_back("psys"); frame_varnames.push_back("psysattn");
      frame_varnames.push_back("purpose"); frame_varnames.push_back("ra");
      frame_varnames.push_back("relhumid"); 
      frame_varnames.push_back("rmspath"); frame_varnames.push_back("sdf");
      frame_varnames.push_back("source");
      frame_varnames.push_back("systemp"); frame_varnames.push_back("tau230");
      frame_varnames.push_back("tcorr"); frame_varnames.push_back("themt"); 
      frame_varnames.push_back("veldop"); frame_varnames.push_back("veltype"); 
      frame_varnames.push_back("vsource");
      frame_varnames.push_back("winddir"); frame_varnames.push_back("windmph");
      frame_varnames.push_back("wsystemp");
      for(vector<string>::iterator i=frame_varnames.begin(); i!=frame_varnames.end(); i++) {
        if (!get_miriad_var(m_frames[curframe].vars,*i)) {frame_varnames.erase(i);i--;} 
      }
      uvgetvr_c(m_uvh, H_REAL, "systemp", (char *)&bufSysTemp,
                m_nants*m_nspect);
      mirvar systempv;
      systempv.name="systemp";
      for (unsigned j=0; j<m_nants*m_nspect; j++) systempv.dval.push_back(bufSysTemp[j]);
      m_frames[curframe].vars.push_back(systempv);
    }

    // baseline number
    unsigned thisbase=(unsigned)(bufPreamble[4]);
    unsigned curbase=find(m_baselines.begin(),m_baselines.end(),thisbase)-
      m_baselines.begin();
    if (curbase >= m_nbase) {
      if (READMIRIAD_DEBUG) cout << "Unexpected baseline " << thisbase << ", skipping" << endl;
      continue;
    }
    // uvw
    m_u[curframe][curbase]=bufPreamble[0];
    m_v[curframe][curbase]=bufPreamble[1];
    m_w[curframe][curbase]=bufPreamble[2];
    // vis
    for(unsigned j=0; j < m_nwide; j++) {
      m_visWideRe[curframe][curbase][j]=bufWideData[2*j];
      m_visWideIm[curframe][curbase][j]=bufWideData[2*j+1];
      m_flagWide[curframe][curbase][j]=bufWideFlags[j];
    }
    int t=0;
    for(unsigned j=0; j < m_nspect; j++) {
      for(unsigned k=0; k<m_nschan[j]; k++) {
        m_visSpecRe[curframe][curbase][j][k]=bufData[2*t+2*k];
        m_visSpecIm[curframe][curbase][j][k]=bufData[2*t+2*k+1];
        m_flagSpec[curframe][curbase][j][k]=bufFlags[t+k];
      }
      t+=m_nschan[j];
    }
    // read next visibility 
    uvset_c(m_uvh, "data",     "channel", 0, 1.0, 1.0, 1.0);
    uvset_c(m_uvh, "preamble", "uvw/time/baseline", 0, 0.0, 0.0, 0.0);
    uvread_c(m_uvh, bufPreamble, bufData, bufFlags, m_nchan, &nread);
    if (!nread) break; // Loop ends if there are no more visibilities.
    if (m_nwide) {
      uvset_c(m_uvh, "data",     "wide", 0, 1.0, 1.0, 1.0);
      uvwread_c(m_uvh, bufWideData, bufWideFlags, m_nwide, &nread);
    }
    nvis++;
  } 
  uvflush_c(m_uvh);
  uvclose_c(m_uvh);
  if (READMIRIAD_DEBUG) cout << endl;
  return true;
}

bool mirfile::get_miriad_var(vector<mirvar> &v, string name) {
  // Get a variable value from the file and add it to
  //  an array of variables for later inclusion in the
  //  matlab structure.
  bool retval=true;
  mirvar var;
  var.name=name;
  char tmpType; int tmpSize, tmpUpdated;
  uvprobvr_c(m_uvh, name.c_str(), &tmpType, &tmpSize, &tmpUpdated);
  if (!tmpSize) {
    if (READMIRIAD_DEBUG) cout << "Could not find variable " << name << endl;
    retval=false;
    return retval;
  }
  float tmpr[tmpSize];
  double tmpd[tmpSize];
  char tmpstr[200]; int tmpstrlen;
  int tmp[tmpSize];
  switch(tmpType) {
    case 'r': //float
      uvgetvr_c(m_uvh,H_REAL,name.c_str(),(char *)tmpr,tmpSize);
      var.dval.resize(tmpSize);
      for (unsigned i=0;i<tmpSize;i++) var.dval[i]=tmpr[i];
      break;
    case 'd': //double
      uvgetvr_c(m_uvh,H_DBLE,name.c_str(),(char *)tmpd,tmpSize);
      var.dval.resize(tmpSize);
      for (unsigned i=0;i<tmpSize;i++) var.dval[i]=tmpd[i];
      break;
    case 'a': //string
      tmpstrlen=200;
      uvgetvr_c(m_uvh,H_BYTE,name.c_str(),tmpstr,tmpstrlen);
      var.sval.push_back(tmpstr);
      break;
    case 'i': //integer
      uvgetvr_c(m_uvh,H_INT,name.c_str(),(char *)tmp,tmpSize);
      var.ival.resize(tmpSize);
      for (unsigned i=0;i<tmpSize;i++) var.ival[i]=tmp[i];
      break;
    default:
      if (READMIRIAD_DEBUG) cout << "Unknown variable type " << tmpType << " for " << name << endl;
      retval=false;
      break;
  }
  v.push_back(var);
  return retval;
}

/* -------------------- 
   Function definitions
   -------------------- */

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[]) {
// This is called from Matlab when you issue the command
//  data=gcpMatReadMiriad('my_filename.mir')

  // Create the matlab structure
  int dims[2]={1,1};

  // Open the miriad file
  string filename = mxArrayToString(prhs[0]);
  mirfile *m=new mirfile(filename);
  plhs[0]=m->get_mxarray();
  nlhs=1;
  delete m;
}

