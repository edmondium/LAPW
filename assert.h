#ifndef ASSERT_
#define ASSERT_

#include <iostream>
#include <source_location>
#include <stdexcept>

// Custom assert functionality using C++23 features
inline void Assert(
    bool condition,
    const std::string& message,
    const std::source_location location = std::source_location::current()) {
    if (!condition) {
        std::cerr << "Assertion failed: " << message << "\n"
                  << "File: " << location.file_name() << "\n"
                  << "Function: " << location.function_name() << "\n"
                  << "Line: " << location.line() << std::endl;
        std::terminate(); // Replace std::abort with std::terminate
    }
}

// Debugging log macro with T_DEBUG support
#ifdef T_DEBUG
#define T_LOG(x) x
#else
#define T_LOG(x)
#endif

#endif // ASSERT_
