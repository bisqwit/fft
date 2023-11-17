#include "defs.hh"
#include "exp.hh"

cvector FFT_tukey(const cvector& input);

enum TukeySplit
{
    Tukey_Split_r1,    // @$r_1$@ on pienin tekijä, ei rekursoida koskaan
    Tukey_Split_half,  // Tekijöiden jako tasaisesti, ei rekursoida koskaan
    Tukey_Split_r2,    // @$r_2$@ on pienin tekijä, ei rekursoida koskaan
    Tukey_Split_r1_rec_always,    // Rekursoidaan aina
    Tukey_Split_half_rec_always,
    Tukey_Split_r2_rec_always,
    Tukey_Split_r1_rec_small_r1,  // Rekursoidaan, jos @$r_1$@ koostuu pienistä
    Tukey_Split_half_rec_small_r1,
    Tukey_Split_r2_rec_small_r1,
    Tukey_Split_r1_rec_small_r2,  // Rekursoidaan, jos @$r_2$@ koostuu pienistä
    Tukey_Split_half_rec_small_r2,
    Tukey_Split_r2_rec_small_r2,
    Tukey_Split_r1_rec_smallm1_r1,  // Rek., jos @$r_1$@ on alkuluku ja @$r_1-1$@ koostuu pienistä
    Tukey_Split_half_rec_smallm1_r1,
    Tukey_Split_r2_rec_smallm1_r1,
    Tukey_Split_r1_rec_smallm1_r2,  // Rek., jos @$r_2$@ on alkuluku ja @$r_2-1$@ koostuu pienistä
    Tukey_Split_half_rec_smallm1_r2,
    Tukey_Split_r2_rec_smallm1_r2
};
cvector FFT_tukey_fast(const cvector& input, TukeySplit how);
