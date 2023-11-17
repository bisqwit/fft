#include <cmath>
#include <numbers>
#include <mutex>
#include <unordered_map>

#include "exp.hh"

exp_vector ExpCircle(std::size_t N, bool force)
{
    return nullptr;

    static std::unordered_map<std::size_t, cvector> data;
    static std::mutex lock;

    for(unsigned n=1; n<=10 && (force||n==1); ++n)
        if(auto j = data.find(N*n); j != data.end())
            return j->second;

    double mul = -2.0 * std::numbers::pi_v<double> / N;
    //auto v  = std::views::iota(N) | std::views::transform([=](auto a){ return std::polar(1.0, a*mul); }) | std::views::common;
    //auto vec = std::ranges::to<cvector>(v);
    //cvector vec(v.begin(), v.end());
    // Precalculate e^(-2πi x/N) for x in 0..(N-1)
    cvector vec(N);
    for(std::size_t a=0; a<N; ++a) vec[a] = std::polar(1.0, a*mul);

    std::unique_lock<std::mutex> lk(lock);
    return data.emplace(N, std::move(vec)).first->second;
}
