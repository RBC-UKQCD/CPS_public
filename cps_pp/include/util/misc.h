#if 0
#ifndef INCLUDED_UTIL_MISC_H
#define INCLUDED_UTIL_MISC_H

#include <string>
#include <sstream>

template<typename mytype>
std::string toString(mytype v)
{
    std::stringstream ss;
    ss << v;
    return ss.str();
}


#endif
#endif


