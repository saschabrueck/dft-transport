// This is merely a dummy to facilitate splitting the code structure
// for now
#include "GreenSolver_common.C"

#ifdef MMC
#include "mcalderara/GreenSolver.C"
#else
#include "sbrueck/GreenSolver.C"
#endif
