#include "NetworkProjectionDynamic.h"


ProjectionDynamic::ProjectionDynamic() :
    tau_psc_(3.0),
    tau_fac_(0.0),
    tau_rec_(800.0),
    U_(0.5),
    x_(1.0),
    y_(0.0),
    u_(0.0)
  { }