#ifndef __libcp2k_tools
#define __libcp2k_tools

#include "libcp2k.H"
#include "Types.H"

/**
 *   \brief initializes a transport_parameters type from a cp2k_transport_parameters type
 *   \author Mohammad Hossein Bani-Hashemian
 */
void intit_transport_parameters_from_cp2k(cp2k_transport_parameters&, transport_parameters*);


#endif
