#pragma once
#include "defs.hh"

struct exp_vector
{
    const cvector* data = nullptr;
    exp_vector(const cvector& d) : data(&d) {}
    exp_vector(std::nullptr_t) : data(nullptr) {}

    /* Returns e^(-2πi x/N) -- even if data.size() ≠ N */
    complex operator() (std::size_t i, std::size_t N) const
    {
        //if(!data)
        return complex(std::polar(1.0, -6.2831853071795864769252867665590057683936*i/N));
        //return (*data)[(i % N) * (data->size()/N)];
    }
};

/* Creates e^(-2πi x/N) for x in 0..(N-1) */
exp_vector ExpCircle(std::size_t N, bool force=false);
