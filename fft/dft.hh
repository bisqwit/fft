#include "defs.hh"
#include "exp.hh"

cvector DFT(const cvector& input);

void DFT(const complex* input, complex* output,
         exp_vector exps,
         std::size_t N,
         std::size_t instride, std::size_t outstride);
