#ifndef STRING_UTIL_H
#define STRING_UTIL_H

#include <stdio.h>
#include <iostream>
#include <string>
#include <algorithm>

void eraseLastFeed(std::string &s)
{
    if (!s.empty() && s[s.length()-1] == '\n' && s[s.length()-1] == '\r') {
        s.erase(s.length()-1);
    };
};

#endif 