#ifndef ASSERT_
#define ASSERT_

#ifdef NO_ARG_CHECK
#define Assert(condition, message) ((void)0)
#else /* NO_ARG_CHECK */
#include <iostream>
#include <cstdlib>
#define Assert(condition, message)\
do {\
  if(!(condition)){ \
    std::cerr << "Assertion failed: " << (message) << "\n" \
              << "File: " << __FILE__ << ", Line: " << __LINE__ << std::endl;\
    std::abort(); \
  } \
} while (0)
#endif /* NO_ARG_CHECK */

#ifdef T_DEBUG
#define T_LOG(x) x
#else
#define T_LOG(x)
#endif

#endif /* ASSERT_ */