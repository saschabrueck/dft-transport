// A header to abstract away various compilers' C++11 support
#ifndef CP2K_OMENDFT_TRANSPORT_CXX11_CHECK_H_
#define CP2K_OMENDFT_TRANSPORT_CXX11_CHECK_H_

// intel 2011.1.075
#if defined(__ICC) && __ICC == 1110
# define HAS_NO_ENUM_CLASS
# define HAS_NO_2D_VEC_INIT
#endif


#endif
