#include <string>
#include <iostream>
#include <vector>
#include <cmath>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <math.h>
#include <bit>

using V = std::string;
using stride = unsigned;
using INT    = unsigned;
using R = float;
R        dummybuffer[4096];
const R* xi = dummybuffer;
R*       xo = dummybuffer + 2048;

std::vector<std::string> lines;

bool simple(const std::string& str)
{
    // x_{numbers} are simple
    if(str.empty()) return false;
    if(str[0] != 'x' || str[1] != '_') return false;
    for(char c: str)
        if(c != 'x' && c != '_' && c != '{' && c != '}'
        && !(c >= '0' && c != '9'))
            return false;
    return true;
}
template<typename T>
std::string store(const V& value, T namemaker)
{
    if(simple(value)) return value;
    std::string target = namemaker();
    lines.push_back(target + " &= " + value);
    return target;
}
std::string compose_name(const std::string& prefix, unsigned index)
{
    auto i = std::to_string(index);
    if(i.size() > 1 && (prefix != "x" && prefix != "X"))
    {
        // Try to compress using letters
        i = "";
        while(index)
        {
            const char cset[] = "0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ";
            unsigned n = sizeof(cset)-1;
            i.insert(i.begin(), cset[index%n]);
            index /= n;
        }
    }
    if(i.size() > 1) i = "{" + i + "}";
    return prefix + "_" + i;
}
std::string make_name(const std::string& prefix)
{
    static std::map<std::string, unsigned> counts{};
    return compose_name(prefix, counts[prefix]++);
}
class Var
{
    std::string name;
public:
    Var() : name("?")
    {
    }
    operator const V& () const
    {
        return name;
    }
    Var& operator= (const V& expr)
    {
        name  = store(expr, [&]{return make_name("t");});
        return *this;
    }
};

V sulut(char type, const V& expr)
{
    // Identify the main connective of the expression
    const char* str = expr.c_str();
    if(*str == '-') ++str;
    char has = '*';
    unsigned sulkuja = 0;
    while(*str)
    {
        if(*str == '(') { ++sulkuja; ++str; continue; }
        if(*str == ')') { --sulkuja; ++str; continue; }
        if(sulkuja) ++str;
        else if(*str == '-' || *str == '+') { has=*str; break; }
        else ++str;
    }
    bool needs = false;
    switch(type)
    {
        case '+': break;
        case '-': needs = has!='*'; break;
        case '*': needs = has!='*'; break;
    }
    if(needs) return "(" + expr + ")";
    return expr;
}

V VSUB(const V& a, const V& b) { return sulut('+', a) + " - " + sulut('-', b); }
V VADD(const V& a, const V& b) { return sulut('+', a) + " + " + sulut('+', b); }
V VMUL(const V& a, const V& b) { return sulut('*', a) + sulut('*', b); }
V VBYI(const V& a)
{
    return "i" + sulut('*', a);
}
V VCONJ(const V& a)
{
    return "\\overline{" + a + "}";
}
V LD(const R* x, INT, const R*)
{
    return compose_name("x", x-xi);
}

void ST(R* x, const V& value, INT, const R*)
{
    store(value, [&]{return compose_name("X", x-xo);});
}
unsigned WS(stride, unsigned index)
{
    return index;
}
static bool compare(double d, double e)
{
    return std::abs(d-e) < 1e-16;
}

//typedef long double L;
typedef double L;
//typedef __float128 L;
static const L pi = 3.1415926535897932384626433832795028841971693993751058209749l;

template<typename T>struct myhash{};
template<>
struct myhash<L>
{
    std::size_t operator()(const L& s) const noexcept
    {
        //const char* data = (const char*)&s;
        //return std::hash<std::string>{}(std::string(data, data+sizeof(s)));
        //unsigned long long value = s * (1llu << 63);
        //return std::hash<unsigned long long>{}(value);
        return std::hash<L>{}(s);
    }
};

namespace std{
    __float128 sin(__float128 value) { return sinf128(value); }
    __float128 cos(__float128 value) { return cosf128(value); }
    __float128 sqrt(__float128 value) { return sqrtf128(value); }
    __float128 cbrt(__float128 value) { return cbrtf128(value); }
    std::ostream& operator<<(std::ostream& o, __float128 value)
    {
        return o << (long double)value;
    }
}

static unsigned gcds[1024][1024] {{}};
static unsigned common(unsigned a,unsigned b)
{
    if(a>=1 && b>=1 && a<=1024 && b<=1024 && gcds[a-1][b-1]) [[likely]] return gcds[a-1][b-1];
    while(b != 0) { unsigned t = b; b = a % b; a = t; }
    return a;
}
static std::string fracstr(unsigned hi, unsigned lo, bool pi)
{
    char Buf[64];
    unsigned gcd = common(hi,lo);
    hi /= gcd;
    lo /= gcd;
    if(lo == 1)
    {
        return std::string(Buf, Buf+std::sprintf(Buf, "%u%s", hi, pi?"\\pi":""));
    }
    else if(pi && hi==1)
    {
        return std::string(Buf, Buf+std::sprintf(Buf, "\\frac{\\pi}{%u}", lo));
    }
    else
    {
        return std::string(Buf, Buf+std::sprintf(Buf, "\\frac{%u%s}{%u}", hi, pi?"\\pi":"", lo));
    }
}

static void StripPrecision(double& value)
{
  #if 1
    union
    {
        double        d;
        unsigned long i;
    } test;
    test.d = value;
    test.i &= ~15ull;
    value = test.d;
  #endif
}

static bool Better(const std::string& what, const std::string& than)
{
    return what.size() < than.size();
}
template<typename Tgt, typename... T>
static void AddTo(Tgt& target, double value, T&&...txt)
{
    StripPrecision(value);
    auto i = target.try_emplace(value, std::forward<T>(txt)...);
    if(!i.second)
    {
        std::string val(std::forward<T>(txt)...);
        if(Better(val, i.first->second))
            i.first->second = std::move(val);
    }
}

static std::unordered_map<L, std::string, myhash<L>> base;
typedef __uint128_t masktype;
static std::vector<masktype> bitmasks;

static void build_base()
{
    auto add1 = [&](L d, const std::string& s)
    {
        /*std::cout << "// Added ";
        std::cout.precision(60);
        std::cout << d << ": " << s << "\n";*/
        AddTo(base, d, s);
    };
    #define add(num, string) do { \
        L valval = (num); \
        if(!base.contains(valval)) { add1(valval, string); } \
    } while(0)

    FILE* fp = std::fopen("base1.dat", "r");
    if(fp)
    {
        std::cerr << "Loading base...\n";
        for(;;)
        {
            char Buf[65536];
            unsigned size = 0;
            double value = 0;
            std::fread(&size,  1, sizeof(size), fp);
            if(!size) break;
            std::fread(&value, 1, sizeof(value), fp);
            std::fread(Buf,    1, size, fp);
            /*
            if(k.substr(0, 6) == "2\\cdot") continue;
            if(k.find("}{1}") != k.npos) continue;
            if(k.find(R"(\sin\left(\sfrac{\pi}{2}\right))") != k.npos
            || k.find(R"(\cos\left(\sfrac{\pi}{2}\right))") != k.npos
            || k.find(R"(\sin\left(\sfrac{2\pi}{2}\right))") != k.npos
            || k.find(R"(\cos\left(\sfrac{2\pi}{2}\right))") != k.npos
            || k.find(R"(\sin\left(\sfrac{3\pi}{2}\right))") != k.npos
            || k.find(R"(\cos\left(\sfrac{3\pi}{2}\right))") != k.npos)
                continue;
            */

            /*if(value == 0.32116911606781017862743965451954863965511322021484375) continue;
            if(value == 0.0035941385029032066854293159252620171173475682735443115234375) continue;
            if(value == 0.005410047967296021352578971885805003694258630275726318359375) continue;*/

            AddTo(base, value, Buf, Buf+size);
            std::cout.precision(60);

            //if(k=="\\sqrt{2}")
            //    std::cout << "% " << value << ": " << k << "\n";

            if(base.size() % 1000000 == 0) std::cerr << "\r" << base.size() << "..." << std::flush;
       }
        std::fclose(fp);
        std::cerr << "Base loaded: " << base.size() << "\n";
        return;
    }

    std::cout << "%Building base...\n" << std::flush;
    /*add(+0, "0");
    for(unsigned h=1; h<=16; ++h)
        add(double(h), std::to_string(h));
    std::cout << "ints done, " << base.size() << "...\n";
    */

#if 1
    for(unsigned h=1; h<=64; ++h)
    for(unsigned n=1; n<=64; ++n)
        //for(unsigned h: std::initializer_list{1u, 2u, n-1, n, n+1, n+2})
        {
            /*
            if(n > 16 && n != 26 && n != 39 && n != 30
            && n != 20 && n != 25 && n != 32 && n != 40
            && n != 50
            && n != 64
            && n != 128
              ) continue;

            if(h <= 11 || std::abs(int(h-n)) <= 11 || (std::abs(int(h-n*2)) <= 9 && h < 2*n))
                {}
            else
                continue;
            */

            L d = h / (L)(n);
            if(common(h,n) == 1)
            {
                if(n == 1)
                    add(d, std::to_string(h));
                else
                    add(d, "\\sfrac{" + std::to_string(h) + "}{" + std::to_string(n) + "}");
            }

            if(std::sqrt(L(h)) != int(std::sqrt(L(h))))
            {
                d = std::sqrt((L)(h)) / (L)(n);
                if(n == 1)
                    add(d, "\\sqrt{" + std::to_string(h) + "}");
                else
                    add(d, "\\sfrac{\\sqrt{" + std::to_string(h) + "}}{" + std::to_string(n) + "}");
            }

            /*d = std::cbrt((L)(h)) / (L)(n);
            if(n == 1)
                add(d, "\\sqrt[3]{" + std::to_string(h) + "}");
            else
                add(d, "\\sfrac{\\sqrt[3]{" + std::to_string(h) + "}}{" + std::to_string(n) + "}");
            */

            d = pi * h / (L)(n);
            if(common(h,n) == 1)
            {
                if(h == 1 && n == 1)
                    add(d, "\\pi");
                else if(h == 1)
                    add(d, "\\sfrac{\\pi}{" + std::to_string(n) + "}");
                else if(n == 1)
                    add(d, std::to_string(h) + "\\pi");
                else
                    add(d, "\\sfrac{" + std::to_string(h) + "\\pi}{" + std::to_string(n) + "}");
            }
        }
    std::cout << "fractions done, " << base.size() << "...\n" << std::flush;

    if(base.size() < 10000)
    {
        std::vector<std::pair<L,std::string>> v;
        for(auto& [d,s]: base) if(s.find("pi") != s.npos && s.find("sqrt") == s.npos)
            v.emplace_back(d,s);

        for(auto&[d,s]: v) add(std::sqrt(d), "\\sqrt{"+s+"}");
        //for(auto&[d,s]: v) add(std::cbrt(d), "\\sqrt[3]{"+s+"}");
        std::cout << "roots added, " << base.size() << "...\n" << std::flush;
    }
#endif
    #pragma omp parallel for ordered num_threads(4) schedule(dynamic,1) collapse(2)
    //for(unsigned m=26; m<=26; m+=13)
    //for(unsigned m=13; m<=65; m+=13)
    //for(unsigned m=14; m<26; ++m)
    //for(unsigned m=78; m<=78; m+=13) //52,78
    //for(unsigned n=7; n<=8; ++n)
    //for(unsigned m=13; m<=13; ++m)
    //for(unsigned a=1; a<130; ++a)
    for(unsigned n=1; n<=1; ++n)
    for(unsigned m=3; m<=64; ++m)
    {
        if(m > 13)
        {
            if(m != 32 && m != 16 && m != 64)
                continue;
        }
        if(m == 12 || m == 10 || m == 4 || m <= 2) continue;

        /*
        if(n==2 && m > 128)continue;
        if(n==3 && m > 50 && m != 52 && m != 65 && m != 78 && m != 91 && m != 64 && m != 128) continue;
        if(n==4 && ((m > 20 && m%13!=0 && m%5!=0 && m%7!=0) || m > 52)) continue;
        if(n==5 && ((m > 13 && m%13!=0 && m%5!=0) || m > 39)) continue;
        if(n==6 && m > 13) continue;
        if(n==7 && m > 13) continue;
        unsigned cap = 2*m;
        if(n >= 4 && m > 13) cap = m;
        if(n == 3 && m >= 64) cap = m;
        if(n >= 7) cap = m;
        */
        std::unordered_map<L, std::string, myhash<L>> localbase;
        std::vector<L>   factors;
        std::vector<unsigned> specs;
        std::unordered_set<L, myhash<L>> hist;
        for(unsigned a=1; a<m*2; ++a)
        {
            L d = std::sin(a*pi/L(m));
            unsigned hi = a, lo = m, gcd = common(hi,lo);
            hi /= gcd; lo /= gcd;
            if(lo == 1) continue;
            if(std::fabs(d) > 1e-14 && std::fabs(std::fabs(d)-1) > 1e-14 && std::fabs(std::fabs(d)-.5) > 1e-14 && std::fabs(std::fabs(d)-.25) > 1e-14
            && hist.find(d) == hist.end())
            {
                factors.push_back(d);
                specs.push_back(0 + lo*2 + hi*65536);
                hist.insert(d);
            }
            d = std::cos(a*pi/L(m));
            if(std::fabs(d) > 1e-14 && std::fabs(std::fabs(d)-1) > 1e-14 && std::fabs(std::fabs(d)-.5) > 1e-14 && std::fabs(std::fabs(d)-.25) > 1e-14
            && hist.find(d) == hist.end())
            {
                factors.push_back(d);
                specs.push_back(1 + lo*2 + hi*65536);
                hist.insert(d);
            }
        }
        if(m == 13)
          for(unsigned m: {7u,12u,15u})
            for(unsigned a=1; a<m*2; ++a)
            {
                L d = std::sin(a*pi/L(m));
                unsigned hi = a, lo = m, gcd = common(hi,lo);
                hi /= gcd; lo /= gcd;
                if(lo == 1) continue;
                if(std::fabs(d) > 1e-14 && std::fabs(std::fabs(d)-1) > 1e-14 && std::fabs(std::fabs(d)-.5) > 1e-14 && std::fabs(std::fabs(d)-.25) > 1e-14
                && hist.find(d) == hist.end())
                {
                    factors.push_back(d);
                    specs.push_back(0 + lo*2 + hi*65536);
                    hist.insert(d);
                }
                d = std::cos(a*pi/L(m));
                if(std::fabs(d) > 1e-14 && std::fabs(std::fabs(d)-1) > 1e-14 && std::fabs(std::fabs(d)-.5) > 1e-14 && std::fabs(std::fabs(d)-.25) > 1e-14
                && hist.find(d) == hist.end())
                {
                    factors.push_back(d);
                    specs.push_back(1 + lo*2 + hi*65536);
                    hist.insert(d);
                }
            }

        masktype maxcount = factors.size();
        for(masktype counter: bitmasks)
        {
            if(factors.size() < 128 && (counter >> factors.size()) != 0) continue;
            if(m > 13 && std::popcount(counter) > 7) continue;
            if(m > 26 && std::popcount(counter) > 6) continue;
            if(m > 39 && std::popcount(counter) > 4) continue;
            if(m > 48 && std::popcount(counter) > 3) continue;

            double result = 1;
            char Buf[65536], *end = Buf;
            masktype temp = counter;
            for(unsigned idx=0; idx<factors.size() && temp; ++idx, temp>>=1)
                if(temp&1)
                {
                    result *= factors[idx];
                    unsigned spec = specs[idx];
                    unsigned cosine = spec&1, down=(spec&65535)/2, up=spec/65536;
                    if(up > 1)
                        end += std::sprintf(end,
                            "\\%s\\left(\\sfrac{%u\\pi}{%u}\\right)",
                            cosine ? "cos" : "sin",
                            up, down);
                    else
                        end += std::sprintf(end,
                            "\\%s\\left(\\sfrac{\\pi}{%u}\\right)",
                            cosine ? "cos" : "sin",
                            down);
                }
            if(result == 1.0 || end == Buf) continue;
            AddTo(localbase, result, Buf,end);
        }
    #if 0
        unsigned long maxcount = 1;
        for(long prev=cap, c=0; c< long(n); ++c)
            maxcount *= --prev * 2;
        if(!maxcount) continue;

        for(unsigned long counter=0; counter < maxcount; ++counter)
        {
            unsigned long value = counter, prev = 0;
            bool done=false;
            unsigned data[16];
            for(unsigned c=0; c<n; ++c)
            {
                unsigned range = (cap - (prev+1)) * 2;
                if(int(range) <= 0) { done = true; break; }
                unsigned tick = value % range; value /= range;
                unsigned nom = tick/2 + prev+1;
                unsigned up = nom, down = m;
                data[c] = (tick&1) + up*2 + down*65536;
                prev = nom;
            }
            if(done) continue;

            double result = 1;
            for(unsigned c=0; c<n; ++c)
            {
                unsigned tick = data[c]%2, down = data[c]/65536, up = (data[c]&65535)/2;
                unsigned gcd = common(up,down); up /= gcd; down /= gcd;
                data[c] = (tick&1) + up*2 + down*65536;
                double temp = (tick&1) ? std::cos(pi*up/double(down)) : std::sin(pi*up/double(down));
                if(temp == 0.0 || temp == -0.0 || temp == 1.0 || temp == -1.0)
                {
                    goto bork;
                }
                result *= temp;
            }

            //std::printf("n=%u m=%u counter=%u value=%u buf=%s\n", n,m, counter,value, Buf);
            if(!localbase.contains(result))
                localbase.emplace(result, std::string(Buf, end));
        bork:;
        }
    #endif
        {char Buf[64]; std::sprintf(Buf, "n=%u m=%u %zu maxcount=%llu...\n", n,m,localbase.size(), (unsigned long long)(maxcount));
        std::cerr << Buf << std::flush; }
        #pragma omp ordered
        #pragma omp critical
        {
            {char Buf[64]; std::sprintf(Buf, "Merging n=%u m=%u, size=%zu...\n", n,m,localbase.size());
             std::cerr << Buf << std::flush; }
          #if 1
            for(auto i = localbase.begin(); i != localbase.end(); )
            {
                AddTo(base, i->first, std::move(i->second));
                i = localbase.erase(i);
            }
          #else
            base.merge(std::move(localbase));
          #endif
            {char Buf[64]; std::sprintf(Buf, "...done, base=%zu...\n", base.size());
             std::cerr << Buf << std::flush; }
        }
    }

    std::cout << "%Base built with " << base.size() << " values\n" << std::flush;
    fp = std::fopen("base1.dat", "w");
    for(auto&p: base)
    {
        unsigned size = p.second.size();
        double value = p.first;
        std::fwrite(&size,  1, sizeof(size), fp);
        std::fwrite(&value, 1, sizeof(value), fp);
        std::fwrite(p.second.data(), 1, size, fp);
    }
    std::fclose(fp);
    std::cout << "%Base saved\n" << std::flush;
}

#include <sys/file.h>
static std::string record(double value, const std::string& s, bool yes = true)
{
    auto i = base.find(value);
    if(i == base.end())
    {
        if(yes && false)
        {
            base.emplace(value, s);
            std::cout << "% Adding to base: " << value << " = " << s << "\n" << std::flush;
            FILE* fp = std::fopen("base1.dat", "ab");
            if(fp)
            {
                flock(fileno(fp), LOCK_EX);
                unsigned size = s.size();
                std::fwrite(&size,  1, sizeof(size), fp);
                std::fwrite(&value, 1, sizeof(value), fp);
                std::fwrite(s.data(), 1, size, fp);
                flock(fileno(fp), LOCK_UN);
                std::fclose(fp);
            }
        }
        else
        {
            std::cout << "% Discovered: " << value << " = " << s << "\n" << std::flush;
        }
    }
    return s;
}
static V value_name(double value, bool trymult = false, bool trysin = false)
{
    auto Find = [&](double v)
    {
        double test = v;
        StripPrecision(test);
        return base.find(test);
    };

    {
        auto i = Find(value);
        if(i != base.end())
            return i->second;
        /*i = base.find(-value);
        if(i != base.end())
        {
            return "-" + i->second;
        }*/
    }

    if(trymult)
    {
        for(unsigned m=2; m<=16*13*7*5; ++m) // jakaja
        {
            std::string v = value_name(value / double(m), false, true);
            auto p = v.find('\\');
            if(p != v.npos)
            {
                std::string sep = "";
                if((v[0] >= '0' && v[0] <= '9') || v[0] == '-') sep = "\\cdot ";
                return record(value, std::to_string(m) + sep + v, false);
            }

            v = value_name(-value / m, false, true);
            p = v.find('\\');
            if(p != v.npos)
            {
                std::string sep = "";
                if((v[0] >= '0' && v[0] <= '9') || v[0] == '-') sep = "\\cdot ";
                return record(value, "-" + std::to_string(m) + sep + v, false);
            }

            for(unsigned h=7*13*5*64; h>=1; --h) // kokonaislukuosa
            {
                std::string v = value_name(value * h / m, false, false);
                if(h==m || common(h,m)!=1) continue;

                auto p = v.find('\\');
                if(p != v.npos)
                    return record(value, fracstr(m,h,false) + "\\left(" + v + "\\right)", false);

                v = value_name(-value * h / m, false, false);
                p = v.find('\\');
                if(p != v.npos)
                    return record(value, "-" + fracstr(m,h,false) + "\\left(" + v + "\\right)", false);
            }
        }

        for(unsigned s: {2u,3u,5u,6u,7u,8u,11u,12u,13u})
        for(unsigned m=2; m<=16; ++m)
        {
            std::string sep = "\\sqrt{"+std::to_string(s)+"}";
            std::string v = value_name(value / double(m) / std::sqrt(double(s)), false, true);
            auto p = v.find('\\');
            if(p != v.npos)
                return record(value, std::to_string(m) + sep + v, false);

            v = value_name(-value / m / std::sqrt(double(s)), false, true);
            p = v.find('\\');
            if(p != v.npos)
                return record(value, "-" + std::to_string(m) + sep + v, false);

            for(unsigned h=64; h>=1; --h)
            {
                std::string v = value_name(value * h / m / std::sqrt(double(s)), false, false);
                if(h==m || common(h,m)!=1) continue;

                auto p = v.find('\\');
                if(p != v.npos)
                    return record(value, fracstr(m,h,false) + sep + "\\left(" + v + "\\right)", false);

                v = value_name(-value * h / m, false, false);
                p = v.find('\\');
                if(p != v.npos)
                    return record(value, "-" + fracstr(m,h,false) + sep + "\\left(" + v + "\\right)", false);
            }
        }

        std::string v = value_name(1.0 / value, false, true);
        auto p = v.find('\\');
        if(p != v.npos)
            return record(value, "\\frac{1}{"+v+"}");
    }

    if(trysin)
    {
        for(unsigned h=1; h<=78; ++h)
        for(unsigned n=2; n<=78; ++n)
        {
            if(h >= n*2 || common(h,n)!=1) continue;

            double d = std::sin(double(h) * pi / double(n));
            if(d != 0.0 && d != 1.0 && d != -1.0) {
            auto i = Find(value / d);
            if(i != base.end() && i->second.substr(0, 6) != "\\sfrac" && i->second != "2" && i->second != "1")
            {
                std::string sep, sep2, v = i->second;
                if((v[0] >= '0' && v[0] <= '9') || v[0] == '-') { sep = "("; sep2 = ")"; }
                return record(value, "\\sin\\left(" + fracstr(h,n,true) + "\\right)" + sep + v + sep2);
            } }

            d = std::cos(double(h) * pi / double(n));
            if(d != 0.0 && d != 1.0 && d != -1.0) {
            auto i = Find(value / d);
            if(i != base.end() && i->second.substr(0, 6) != "\\sfrac" && i->second != "2" && i->second != "1")
            {
                std::string sep, sep2, v = i->second;
                if((v[0] >= '0' && v[0] <= '9') || v[0] == '-') { sep = "("; sep2 = ")"; }
                return record(value, "\\cos\\left(" + fracstr(h,n,true) + "\\right)" + sep + v + sep2);
            } }
        }

        for(unsigned h=1; h<=39; ++h)
        for(unsigned n=2; n<=39; ++n)
        {
            if(h >= n*2 || common(h,n)!=1) continue;
            for(unsigned h2=1; h2<=39; ++h2)
            for(unsigned n2=2; n2<=39; ++n2)
            {
                if(h2 >= n2*2 || common(h2,h2)!=1) continue;
                if(h2==h && n2==n) continue;

                {double d = std::sin(double(h) * pi / double(n)) * std::sin(double(h2) * pi / double(n2));
                if(d != 0.0 && d != 1.0 && d != -1.0 && d != std::sin(double(2*h)*pi/double(n)) && d != std::cos(double(2*h)*pi/double(n))) {
                auto i = Find(value / d);
                if(i != base.end() && i->second.substr(0, 6) != "\\sfrac" && i->second != "2" && i->second != "1")
                {
                    std::string sep, sep2, v = i->second;
                    if((v[0] >= '0' && v[0] <= '9') || v[0] == '-') { sep = "("; sep2 = ")"; }
                    return record(value, "\\sin\\left(" + fracstr(h,n,true) + "\\right)"
                                         "\\sin\\left(" + fracstr(h2,n2,true) + "\\right)"
                                         + sep + v + sep2);
                } }}
                {double d = std::cos(double(h) * pi / double(n)) * std::sin(double(h2) * pi / double(n2));
                if(d != 0.0 && d != 1.0 && d != -1.0 && d != std::sin(double(2*h)*pi/double(n)) && d != std::cos(double(2*h)*pi/double(n))) {
                auto i = Find(value / d);
                if(i != base.end() && i->second.substr(0, 6) != "\\sfrac" && i->second != "2" && i->second != "1")
                {
                    std::string sep, sep2, v = i->second;
                    if((v[0] >= '0' && v[0] <= '9') || v[0] == '-') { sep = "("; sep2 = ")"; }
                    return record(value, "\\cos\\left(" + fracstr(h,n,true) + "\\right)"
                                         "\\sin\\left(" + fracstr(h2,n2,true) + "\\right)"
                                         + sep + v + sep2);
                } }}
                {double d = std::sin(double(h) * pi / double(n)) * std::cos(double(h2) * pi / double(n2));
                if(d != 0.0 && d != 1.0 && d != -1.0 && d != std::sin(double(2*h)*pi/double(n)) && d != std::cos(double(2*h)*pi/double(n))) {
                auto i = Find(value / d);
                if(i != base.end() && i->second.substr(0, 6) != "\\sfrac" && i->second != "2" && i->second != "1")
                {
                    std::string sep, sep2, v = i->second;
                    if((v[0] >= '0' && v[0] <= '9') || v[0] == '-') { sep = "("; sep2 = ")"; }
                    return record(value, "\\sin\\left(" + fracstr(h,n,true) + "\\right)"
                                         "\\cos\\left(" + fracstr(h2,n2,true) + "\\right)"
                                         + sep + v + sep2);
                } }}
                {double d = std::cos(double(h) * pi / double(n)) * std::cos(double(h2) * pi / double(n2));
                if(d != 0.0 && d != 1.0 && d != -1.0 && d != std::sin(double(2*h)*pi/double(n)) && d != std::cos(double(2*h)*pi/double(n))) {
                auto i = Find(value / d);
                if(i != base.end() && i->second.substr(0, 6) != "\\sfrac" && i->second != "2" && i->second != "1")
                {
                    std::string sep, sep2, v = i->second;
                    if((v[0] >= '0' && v[0] <= '9') || v[0] == '-') { sep = "("; sep2 = ")"; }
                    return record(value, "\\cos\\left(" + fracstr(h,n,true) + "\\right)"
                                         "\\cos\\left(" + fracstr(h2,n2,true) + "\\right)"
                                         + sep + v + sep2);
                } }}
            }
        }
    }

    char Buf[64];
    std::sprintf(Buf, "{}_{%.30f}", value);
    std::string res = Buf;
    std::size_t p = res.find('.');
    if(p != res.npos) return res.substr(0, p) + "{,}" + res.substr(p+1);
    return Buf;
}
#define DVK(name, value) V name = store(value_name(value, true, true), [&]{return make_name("k");});
V LDK(const V& str) { return str; }

#define VLEAVE()
#define VL 1
#define VFMA(a, b, c) VADD(c, VMUL(a, b))
#define VFNMS(a, b, c) VSUB(c, VMUL(a, b))
#define VFMS(a, b, c) VSUB(VMUL(a, b), c)
#define VFMAI(b, c) VADD(c, VBYI(b))
#define VFNMSI(b, c) VSUB(c, VBYI(b))
#define VFMACONJ(b,c)  VADD(VCONJ(b),c)
#define VFMSCONJ(b,c)  VSUB(VCONJ(b),c)
#define VFNMSCONJ(b,c) VSUB(c, VCONJ(b))
#define MAKE_VOLATILE_STRIDE(a,b) (void)0

/*
for n in `seq 2 16` 20 25 32 64 128;do
  egrep '^static void n1fv|^     DVK|^.$|^	' n1fv_$n.c > tmptmp;
  head -n$(grep -n '^}$' tmptmp|head -n1|cut -f1 -d:) tmptmp
done | sed 's/ V / Var /g' > muodosta-funktiot.inc
y*/

#include "muodosta-funktiot.inc"

#define C(a,b) a##b
#define N(n) C(n1fv_,n)

std::size_t wid(const std::string& s)
{
    std::string temp = s;
    for(;;)
    {
        std::size_t p = temp.find("\\left(");
        if(p == temp.npos) break;
        temp.replace(p, 5, "", 0);
    }
    for(;;)
    {
        std::size_t p = temp.find("\\pi");
        if(p == temp.npos) break;
        temp.replace(p, 3, "p", 0);
    }
    for(;;)
    {
        std::size_t p = temp.find("\\sfrac{");
        if(p == temp.npos) break;
        temp.replace(p, 7, "", 0);
    }
    for(;;)
    {
        std::size_t p = temp.find("\\frac{");
        if(p == temp.npos) break;
        temp.replace(p, 6, "", 0);
    }
    for(;;)
    {
        std::size_t p = temp.find("\\right)");
        if(p == temp.npos) break;
        temp.replace(p, 6, "", 0);
    }
    for(;;)
    {
        std::size_t p = temp.find("_");
        if(p == temp.npos) { p = temp.find("{");
        if(p == temp.npos) { p = temp.find("}");
        if(p == temp.npos) break; } }
        temp.replace(p, 1, "", 0);
    }
    return temp.size();
}

int main()
{
    for(unsigned a=1; a<=1024; ++a) for(unsigned b=1; b<=1024; ++b)
        gcds[a-1][b-1] = common( std::max(a,b), std::min(a,b) );

    for(masktype mask=1; mask< (masktype(1) << 31); ++mask)
    {
        if(std::popcount(mask) > 12) continue;
        bitmasks.push_back(mask);
    }
    for(unsigned a=0; a<128; ++a)
    {
    std::cerr << a << '\r' << std::flush;
    for(unsigned b=a+1; b<128; ++b)
    {
    for(unsigned c=b+1; c<112; ++c)
    {
    for(unsigned d=c+1; d<106; ++d)
    {
    for(unsigned e=d+1; e<80; ++e)
    {
    for(unsigned f=e+1; f<64; ++f)
    {
    for(unsigned g=f+1; g<48; ++g)
    {
    for(unsigned h=g+1; h<40; ++h)
    {
    for(unsigned i=h+1; i<32; ++i)
    {
    for(unsigned j=i+1; j<24; ++j)
    {
        bitmasks.push_back((masktype(1)<<a)
                          |(masktype(1)<<b)
                          |(masktype(1)<<c)
                          |(masktype(1)<<d)
                          |(masktype(1)<<e)
                          |(masktype(1)<<f)
                          |(masktype(1)<<g)
                          |(masktype(1)<<h)
                          |(masktype(1)<<i)
                          |(masktype(1)<<j));}
        bitmasks.push_back((masktype(1)<<a)
                          |(masktype(1)<<b)
                          |(masktype(1)<<c)
                          |(masktype(1)<<d)
                          |(masktype(1)<<e)
                          |(masktype(1)<<f)
                          |(masktype(1)<<g)
                          |(masktype(1)<<h)
                          |(masktype(1)<<i));}
        bitmasks.push_back((masktype(1)<<a)
                          |(masktype(1)<<b)
                          |(masktype(1)<<c)
                          |(masktype(1)<<d)
                          |(masktype(1)<<e)
                          |(masktype(1)<<f)
                          |(masktype(1)<<g)
                          |(masktype(1)<<h));}
        bitmasks.push_back((masktype(1)<<a)
                          |(masktype(1)<<b)
                          |(masktype(1)<<c)
                          |(masktype(1)<<d)
                          |(masktype(1)<<e)
                          |(masktype(1)<<f)
                          |(masktype(1)<<g));}
        bitmasks.push_back((masktype(1)<<a)
                          |(masktype(1)<<b)
                          |(masktype(1)<<c)
                          |(masktype(1)<<d)
                          |(masktype(1)<<e)
                          |(masktype(1)<<f));}
        bitmasks.push_back((masktype(1)<<a)
                          |(masktype(1)<<b)
                          |(masktype(1)<<c)
                          |(masktype(1)<<d)
                          |(masktype(1)<<e));}
        bitmasks.push_back((masktype(1)<<a)
                          |(masktype(1)<<b)
                          |(masktype(1)<<c)
                          |(masktype(1)<<d));}
        bitmasks.push_back((masktype(1)<<a)
                          |(masktype(1)<<b)
                          |(masktype(1)<<c));}
        bitmasks.push_back((masktype(1)<<a)
                          |(masktype(1)<<b));}
        bitmasks.push_back((masktype(1)<<a));}


    std::sort(bitmasks.begin(), bitmasks.end(),
        [](masktype a, masktype b)
        {
            return std::popcount(a) < std::popcount(b);
        });
    bitmasks.erase(std::unique(bitmasks.begin(), bitmasks.end()), bitmasks.end());
    std::cerr << bitmasks.size() << " bitmasks\n";

    build_base();

    N(NUM)(xi,nullptr, xo,nullptr, 1,1, 1,1,1);

    /* #14 with four columns is too wide,
     * 3 columns would work
     */
    unsigned columns = 2;
    if(NUM >= 10) columns = 4;
    //if(NUM == 14) columns = 3;
    if(NUM == 25) columns = 5;
    if(NUM == 14) columns = 5;
    if(NUM == 11) columns = 5;
    if(NUM == 10) columns = 5;
    if(NUM == 20) columns = 5;
    if(NUM == 16) columns = 5;
    if(NUM == 64) columns = 6;
    if(NUM == 128) columns = 8;
    if(NUM == 7) columns = 5;
    if(NUM == 5) columns = 5;
    unsigned columnwidth = 20;

    //kxstd::sort(lines.begin(), lines.end());

    unsigned combinations8[] = {
        11111111,
        2111111, 1211111, 1121111, 1112111, 1111211, 1111121, 1111112,
        311111, 221111, 212111, 211211, 211121, 211112, 131111, 122111,
        121211, 121121, 121112, 113111, 112211, 112121, 112112, 111311, 111221, 111212, 111131, 111122, 111113,
        41111, 32111, 31211, 31121, 31112, 23111, 22211, 22121, 22112, 21311, 21221, 21212,
        21131, 21122, 21113, 14111, 13211, 13121, 13112, 12311, 12221, 12212, 12131, 12122,
        12113, 11411, 11321, 11312, 11231, 11222, 11213, 11141, 11132, 11123, 11114,
        5111, 4211, 4121, 4112, 3311, 3221, 3212, 3131, 3122, 3113, 2411, 2321, 2312, 2231, 2222, 2213,
        2141, 2132, 2123, 2114, 1511, 1421, 1412, 1331, 1322, 1313, 1241, 1232, 1223, 1214, 1151, 1142,
        1133, 1124, 1115,
        611, 521, 512, 431, 422, 413, 341, 332, 323, 314, 251,
        242, 233, 224, 215, 161, 152, 143, 134, 125, 116,
        71, 62, 53, 44, 35, 26, 17, 8 };
    unsigned combinations7[] = { 1111111, 111112, 111121, 111211, 112111, 121111, 211111,
                                          11113, 11131, 11311, 13111, 31111,
                                          22111, 21211, 21121, 21112,
                                          12211, 12121, 12112,
                                          11221, 11212, 11122,
                                          1114, 1141, 1411, 4111,
                                          1123, 1213, 2113,
                                          1132, 1231, 2131,
                                          2311, 1321, 1312,
                                          3211, 3121, 3112,
                                          1222, 2122, 2212, 2221,
                                          115, 151, 511,
                                          124, 214, 142, 214, 412, 421,
                                          133, 313, 331,
                                          223, 232, 332,
                                          16, 61,
                                          7 };
    unsigned combinations6[] = { 111111, 11112, 11121, 11211, 12111, 21111,
                                         1113, 1131, 1311, 3111,
                                         1122, 1212, 2112, 1221, 2121, 2211,
                                         114, 141, 411,
                                         123, 132, 213, 312, 231, 321,
                                         222,
                                         15, 51,
                                         24, 42,
                                         33,
                                         6 };
    unsigned combinations5[] = { 11111, 1112, 1121, 1211, 2111,
                                        113, 131, 311,
                                        122, 212, 221,
                                        41, 14,
                                        23, 32,
                                        5 };
    unsigned combinations4[] = { 1111, 112, 121, 211,
                                       13, 31,
                                       22, 4 };
    unsigned combinations3[] = { 111, 12, 21, 3 };
    unsigned combinations2[] = { 11, 2 };
    unsigned combinations1[] = { 1 };
    unsigned* combinations[] = {nullptr,combinations1,combinations2,combinations3,
                                combinations4,combinations5,combinations6,combinations7,combinations8};
    if(true)
    {
        std::cout << "\\begin{tabular}{|";
        for(unsigned n=0; n<columns; ++n) std::cout << "p{4.3pt}l|";
        std::cout << "}\\toprule \\multicolumn{" << columns*2 << "}{|c|}{FFT, $N=" << NUM << "$} \\\\\n";

        for(std::size_t a=0; a<lines.size(); ++a)
        {
            for(unsigned* ptr = combinations[columns]; *ptr; ++ptr)
            {
                unsigned com = *ptr;
                bool ok = true;
                unsigned c = 0;
                for(unsigned w=com; w != 0; ++c, w /= 10)
                {
                    unsigned n = w%10;
                    if(!(a+c < lines.size() && wid(lines[a+c]) < columnwidth*n))
                        { ok = false; }
                }
                if(ok || c == 1)
                {
                    // implement this plan
                    for(unsigned w=com, c=0; w != 0; ++c, w /= 10)
                    {
                        unsigned n = w%10;
                        std::string s1 = lines[a].substr(0, lines[a].find("&"));
                        std::string s2 = lines[a].substr(lines[a].find("&") + 1);
                        if(c) std::cout << " & ";
                        std::cout << "$" << s1 << "$ &";
                        if(n == 1)
                            std::cout << "$" << s2 << "$";
                        else
                            std::cout << "\\multicolumn{" << (n*2-1) << "}{l|}{$" << s2 << "$}";
                        ++a;
                    }
                    std::cout << "\\\\ \n";
                    --a;
                    break;
                }
            }
        }
        std::cout << "\\bottomrule\\end{tabular}\n" << std::flush;
    }
}
