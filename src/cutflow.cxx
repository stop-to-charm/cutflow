#include <cstdio>
#include "SmartChain.hh"

int main (int narg, const char* argv[]) { 
  SmartChain* chain = new SmartChain("susy"); 
  for (int iii = 1; iii < narg; iii++) { 
    printf("file: %s, %i of %i\n", argv[iii], iii, narg - 1); 
    chain->add(argv[iii]); 
  }

}
