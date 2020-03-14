#include "cphos.h"

void CPhos::input( long double input )
{
  // phosphatase pathway is equal to kinase pathway
  // with a phosphatase operation at the end
  // where the phosphatase operation is approxpimated to be
  // a fractional integration
  CKin::input(input);
  fi.input(expo_sum);
  fi.fIntegrate();
}
