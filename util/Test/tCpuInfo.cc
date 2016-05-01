#include <iostream>
#include <iomanip>

#include <string>
#include <cmath>

#include "gcp/program/Program.h"

#include "gcp/util/Exception.h"
#include "gcp/util/OsInfo.h"

using namespace std;
using namespace gcp::util;
using namespace gcp::program;

KeyTabEntry Program::keywords[] = {
  { "rad",        "0", "d", "Radians"},
  { END_OF_KEYWORDS}
};

void Program::initializeUsage() {};

void cpuID(unsigned i, unsigned regs[4]) {
#ifdef _WIN32
  __cpuid((int *)regs, (int)i);

#else
  asm volatile
    ("cpuid" : "=a" (regs[0]), "=b" (regs[1]), "=c" (regs[2]), "=d" (regs[3])
     : "a" (i), "c" (0));
  // ECX is set to zero for CPUID function 4
#endif
}


int Program::main()
{
  unsigned regs[4];

  // Get vendor
  char vendor[12];
  cpuID(0, regs);
  ((unsigned *)vendor)[0] = regs[1]; // EBX
  ((unsigned *)vendor)[1] = regs[3]; // EDX
  ((unsigned *)vendor)[2] = regs[2]; // ECX
  string cpuVendor = string(vendor, 12);

  // Get CPU features
  cpuID(1, regs);
  unsigned cpuFeatures = regs[3]; // EDX

  // Logical core count per CPU
  cpuID(1, regs);
  unsigned logical = (regs[1] >> 16) & 0xff; // EBX[23:16]
  cout << " logical cpus: " << logical << endl;
  unsigned cores = logical;

  COUT("Vendor is: " << cpuVendor);

  if (cpuVendor == "GenuineIntel") {
    // Get DCP cache info
    cpuID(4, regs);
    cores = ((regs[0] >> 26) & 0x3f) + 1; // EAX[31:26] + 1

  } else if (cpuVendor == "AuthenticAMD") {
    // Get NC: Number of CPU cores - 1
    cpuID(0x80000008, regs);
    cores = ((unsigned)(regs[2] & 0xff)) + 1; // ECX[7:0] + 1
  }

  cout << "    cpu cores: " << cores << endl;

  // Detect hyper-threads  
  bool hyperThreads = cpuFeatures & (1 << 28) && cores < logical;

  cout << "hyper-threads: " << (hyperThreads ? "true" : "false") << endl;

  return 0;
}

