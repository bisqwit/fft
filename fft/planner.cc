#include <string>
#include <cstdio>
#include <vector>
#include <complex>
#include <fftw3.h>
#include <utility>
using complex = std::complex<float>;
using cvector = std::vector<complex>;

static std::size_t hash(const std::string& fn)
{
    std::FILE* fp = std::fopen(fn.c_str(), "rt");
    if(!fp) return 0;
    char Buf[8192];
    std::size_t result=0;
    while((std::fgets(Buf, sizeof Buf-1, fp)))
        result += std::hash<std::string>{}(Buf);
    std::fclose(fp);
    std::fprintf(stderr, "%s: %zX\n", fn.c_str(), result);
    return result;
}
#include <sys/file.h>
#include <fcntl.h>
#include <unistd.h>
int main(int argc, char** argv)
{
    std::string golden = "fftw-wisdom.dat";
    bool is_golden = (golden == argv[3]);

    //fftwf_import_wisdom_from_filename(golden.c_str());
    for(int a=4; a<argc; ++a)
    {
        std::fprintf(stderr, "Loading %s", argv[a]); std::fflush(stderr);
        int ret = fftwf_import_wisdom_from_filename(argv[a]);
        std::fprintf(stderr, " - %d\n", ret); std::fflush(stderr);
        if(is_golden && golden != argv[a])
        {
            std::fprintf(stderr, "Saving %s", argv[3]);
            std::string temp = argv[3] + std::string(".temp");
            fftwf_export_wisdom_to_filename(temp.c_str());
            std::rename(temp.c_str(), argv[3]);
            std::remove(argv[a]);
            std::fprintf(stderr, ", done\n");
        }
    }
    unsigned first = std::stoi(argv[1]);
    unsigned last  = std::stoi(argv[2]);
    std::fprintf(stderr, "Planning %u..%u\n", first,last);
    fftwf_set_timelimit(-1);
    for(unsigned n = first; n <= last; ++n)
    {
        cvector input(n);
        cvector output(n);

        fftwf_plan plan = fftwf_plan_dft_1d(n,
            (fftwf_complex*)&input[0],
            (fftwf_complex*)&output[0],
            FFTW_FORWARD, FFTW_EXHAUSTIVE);
        fftwf_destroy_plan(plan);
    }

    /*std::fprintf(stderr, "Planning %u..%u, round 2\n", first,last);
    fftwf_set_timelimit(-1);
    for(unsigned n = first; n <= last; ++n)
    {
        cvector input(n);
        cvector output(n);

        fftwf_plan plan = fftwf_plan_dft_1d(n,
            (fftwf_complex*)&input[0],
            (fftwf_complex*)&output[0],
            FFTW_FORWARD, FFTW_ESTIMATE);
        fftwf_destroy_plan(plan);
    }*/

    std::string temp = argv[3] + std::string(".temp");
    std::fprintf(stderr, "Saving %s\n", argv[3]);
    int fd = open(golden.c_str(), O_RDONLY);
    flock(fd, LOCK_EX);
    fftwf_export_wisdom_to_filename(temp.c_str());
    //auto code = [](const std::string& f) { return "<( while read line;do echo \"$line\"|md5sum;done < "+f + "|sort)"; };
    //int status = std::system(("bash -c 'cmp "+code(golden)+" "+code(temp)+"'").c_str());
    int status = hash(golden) != hash(temp);
    if(status)
    {
        std::fprintf(stderr, "Keeping %s\n", argv[3]);
        std::rename(temp.c_str(), argv[3]);
    }
    else
    {
        std::fprintf(stderr, "Discarding %s as duplicate\n", argv[3]);
        std::remove(temp.c_str());
    }
    flock(fd, LOCK_UN);
    close(fd);
}
