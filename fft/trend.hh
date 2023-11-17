#include <cmath>
#include <cstdio>
#include <string>
#include <complex>
#include <vector>

#define PRECISION 4
#define PRECISION_EXP 3

namespace std {
    template<typename T> auto exp10(T&& v) { return std::pow(std::remove_cvref_t<T>(10), v); }
    template<typename T> auto cbrt(T&& v) { return std::pow(v, std::remove_cvref_t<T>(1)/std::remove_cvref_t<T>(3)); }
}

inline bool IsZero(const std::string& s, int sign=0)
{
    const char* p = s.c_str();
    if(*p == '+') { ++p; if(sign < 0) return false; }
    else if(*p == '-') { ++p; if(sign > 0) return false; }
    else { if(sign < 0) return false; }
    while(*p == '0') ++p;
    if(*p == '.')
    {
        ++p;
        while(*p == '0') ++p;
    }
    return !*p;
}
inline bool IsOne(const std::string& s, int sign=0)
{
    const char* p = s.c_str();
    if(*p == '+') { ++p; if(sign < 0) return false; }
    else if(*p == '-') { ++p; if(sign > 0) return false; }
    else { if(sign < 0) return false; }
    while(*p == '0') ++p;
    if(*p == '1') ++p; else return false;
    if(*p == '.')
    {
        ++p;
        while(*p == '0') ++p;
    }
    return !*p;
}

template<typename Float = double>
struct Trend
{
    virtual inline ~Trend();
    // Translate X and Y coordinates so that linear regression can be used
    virtual void estimate(const std::vector<std::pair<Float, Float>>& values) = 0;
    // Express the formula as a string
    virtual std::string translate(const std::string& var, int mulx = 0) const = 0;
    // Predict Y coordinate from X coordinate
    virtual Float predict(Float) const = 0;
    virtual const char* name() const = 0;

    std::string powstr(Float degree, const std::string& var) const
    {
        char Buf[32];
        if(degree == 0) return {};
        if(degree == 1.0) return var;
        std::sprintf(Buf, "^{%.*g}", PRECISION_EXP, degree);
        return var + Buf;
    }
};
template<typename Float>
inline Trend<Float>::~Trend()
{
}

template<unsigned Degree, typename Float = double>
struct TrendPoly: public Trend<Float>
{
    Float factors[Degree+1];

    void estimate(const std::vector<Float>& x, const std::vector<Float>& y, bool intercept = false, float interceptvalue = 0)
    {
        auto n = x.size();
        if(Degree == 1)
        {
            Float sumX = 0; for(auto v: x) sumX += v; Float avgX = sumX/n;
            Float sumY = 0; for(auto v: y) sumY += v - (intercept ? interceptvalue : 0); Float avgY = sumY/n;
            if(intercept) avgX = avgY = 0;
            Float sxx = 0, /*syy = 0, */sxy = 0;
            for(std::size_t a=0; a<n; ++a)
            {
                sxx += std::pow(x[a]-avgX, 2);
                //syy += std::pow(y[a]-avgY, 2);
                sxy += (x[a]-avgX) * (y[a]-avgY);
            }
            factors[1] = sxy / sxx;
            factors[0] = intercept ? interceptvalue : (avgY - factors[1]*avgX);
        }
        else
        {
            std::size_t NoValues = y.size(), NoPowers = intercept ? Degree : (Degree+1), minorsize = std::min(NoValues,NoPowers);
            std::vector<Float> yvector(y), qr(NoValues * NoPowers), diag(minorsize);
            if(intercept) for(auto& v: yvector) v -= interceptvalue;
            for(std::size_t j=0; j<NoPowers; ++j)
                for(std::size_t i=0; i<NoValues; ++i)
                    { qr[i + j*NoValues] = std::pow(x[i], int(j+(intercept?1:0))); }
            for(std::size_t minor=0; minor<minorsize; ++minor)
            {
                Float normsqr = 0;
                for(std::size_t x=minor; x<NoValues; ++x) { normsqr += std::pow(qr[x + minor*NoValues], 2); }
                Float a = diag[minor] = (qr[minor+minor*NoValues] > 0 ? -std::sqrt(normsqr) : std::sqrt(normsqr));
                if(a == 0.0) continue;
                qr[minor + minor*NoValues] -= a;
                for(std::size_t c = minor+1; c<NoPowers; ++c)
                {
                    Float alfa = 0;
                    for(std::size_t r = minor; r < NoValues; ++r) alfa -= qr[r + c*NoValues] * qr[r + minor*NoValues];
                    alfa /= a * qr[minor + minor*NoValues];
                    for(std::size_t r = minor; r < NoValues; ++r) qr[r + c*NoValues] -= alfa * qr[r + minor*NoValues];
                }
            }
            for(std::size_t minor=0; minor<minorsize; ++minor)
            {
                Float dot = 0;
                for(std::size_t r=minor; r<NoValues; ++r) dot += yvector[r] * qr[r + minor*NoValues];
                dot /= diag[minor]*qr[minor+minor*NoValues];
                for(std::size_t r=minor; r<NoValues; ++r) yvector[r] += dot * qr[r + minor*NoValues];
            }
            for(std::size_t r=diag.size(); r-- > 0; )
            {
                Float yr = factors[r + (intercept?1:0)] = (yvector[r] /= diag[r]);
                for(std::size_t i=0; i<r; ++i) yvector[i] -= yr * qr[i + r*NoValues];
            }
            if(intercept) factors[0] = interceptvalue;
        }
    }
    virtual void estimate(const std::vector<std::pair<Float, Float>>& values) override
    {
        std::vector<Float> x,y; x.reserve(values.size()); y.reserve(values.size());
        for(auto& v: values) { x.push_back(v.first); y.push_back(v.second); }
        estimate(x,y);
    }
    virtual std::string translate(const std::string& var = "x", int mulx=0) const override
    {
        std::string result;
        bool first = true;
        for(int degree=Degree; degree>=0; --degree)
        {
            int deg = degree + mulx;
            char Buf[32];
            if(!first)
                std::sprintf(Buf, "%+.*g", PRECISION, factors[degree]);
            else
                std::sprintf(Buf, "%.*g", PRECISION, factors[degree]);
            std::string k = Buf;
            if(IsZero(k)) continue;
            if(k[0] == '+' || k[0] == '-')
                result += " " + std::string(1, k[0]) + " " + k.substr(1);
            else
                result += k;
            first = false;
            result += Trend<Float>::powstr(deg, " "+var);
        }
        return result;
    }
    // Predict Y coordinate from X coordinate
    virtual Float predict(Float x) const override
    {
        Float result = 0, xfac = 1;
        for(std::size_t a=0; a<Degree+1; ++a, xfac *= x)
            result += factors[a] * xfac;
        return result;
    }
    virtual const char* name() const { return Degree==1?"Linear":(Degree==2?"Quadratic":"Polynomial"); }

    virtual TrendPoly<Degree - (Degree?1:0), Float> differentiate() const
    {
        TrendPoly<Degree - (Degree?1:0), Float> result;
        for(std::size_t i=0; i<Degree; ++i)
            result.factors[i] = factors[i+1] * (i+1);
        return result;
    }
    virtual Float solve_newton(Float x0, Float precision, std::size_t maxiter) const
    {
        auto diff = differentiate();
        for(std::size_t i = 0; i < maxiter; ++i)
        {
            auto f = predict(x0), g = diff.predict(x0);
            auto err = std::abs(g);
            if(err < Float(1e-40)) break;
            auto x1 = x0 - f/g;
            if(std::abs(x0-x1) < precision) { x0=x1; break; }
            x0 = x1;
        }
        return x0;
    }
    virtual std::vector<std::complex<Float>> solve() const
    {
        using cplx = std::complex<Float>;
        switch(Degree)
        {
            case 0: return {};
            case 1: return {-factors[0]/factors[1]}; // 5x+3 = 0
            case 2: { cplx c = factors[0], b = factors[1], a = factors[2];
                      return { -b+std::sqrt(b*b-a*c*Float(4))/(a+a),
                               -b-std::sqrt(b*b-a*c*Float(4))/(a+a) };
                    }
            case 3: { cplx d = factors[0], c = factors[1], b = factors[2], a = factors[3];
                      cplx d0 = b*b-Float(3)*a*c, d1 = Float(2)*b*b*b - Float(9)*a*b*c + Float(27)*a*a*d;
                      cplx C = std::cbrt((d1 + std::sqrt(d1*d1-Float(4)*d0*d0*d0)) / Float(2));
                      cplx k = (-1+std::cbrt(Float(-3)))/Float(2), f = -Float(1)/(Float(3)*a);
                      return { f*(b+std::pow(k,0)*C+d0/(std::pow(k,0)*C)),
                               f*(b+std::pow(k,1)*C+d0/(std::pow(k,1)*C)),
                               f*(b+std::pow(k,2)*C+d0/(std::pow(k,2)*C)) };
                    }
            case 4: { cplx E = factors[0], D = factors[1], C = factors[2], B = factors[3], A = factors[4];
                      cplx a = -Float(3)*B*B/(Float(8)*A*A)+C/A, b = B*B*B/(Float(8)*A*A*A)-B*C/(Float(2)*A*A)+D/A;
                      cplx g = -Float(3)*B*B*B*B/(Float(256)*A*A*A*A)+C*B*B/(Float(16)*A*A*A)-B*D/(Float(4)*A*A)+E/A;
                      cplx P = -a*a/Float(12)-g, Q = -a*a*a/Float(108)+a*g/Float(3) - b*b/Float(8);
                      cplx R = -Q/Float(2) + std::sqrt(Q*Q/Float(4)+P*P*P/Float(27)), U = std::cbrt(R);
                      cplx y = -Float(5)*a/Float(6) + (U!=Float{} ? U-P/(Float(3)*U) : -std::cbrt(Q));
                      cplx W = std::sqrt(a + Float(2)*y), k = -B/(Float(4)*A);
                      return { k + ((+W)+std::sqrt(-(Float(3)*a+Float(2)*y+Float(2)*b/W))) / Float(2),
                               k + ((+W)-std::sqrt(-(Float(3)*a+Float(2)*y+Float(2)*b/W))) / Float(2),
                               k + ((-W)+std::sqrt(-(Float(3)*a+Float(2)*y-Float(2)*b/W))) / Float(2),
                               k + ((-W)-std::sqrt(-(Float(3)*a+Float(2)*y-Float(2)*b/W))) / Float(2) };
                    }
        }
    }
};

template<typename Float=double>
struct TrendLinear: public TrendPoly<1, Float>
{
    virtual void estimate(const std::vector<std::pair<Float, Float>>& values) override
    {
        std::vector<Float> x,y; x.reserve(values.size()); y.reserve(values.size());
        for(auto& v: values)
        {
            x.push_back(v.first);
            y.push_back(v.second);
        }
        TrendPoly<1>::estimate(x, y, true, 0);
    }
};

template<typename Float=double>
struct TrendQuadratic: public TrendPoly<2, Float>
{
    virtual void estimate(const std::vector<std::pair<Float, Float>>& values) override
    {
        std::vector<Float> x,y; x.reserve(values.size()); y.reserve(values.size());
        for(auto& v: values)
        {
            x.push_back(v.first);
            y.push_back(v.second * v.first);
        }
        //TrendPoly<2>::estimate(x, y);
        TrendPoly<2>::estimate(x, y, true, 0);
    }
    virtual std::string translate(const std::string& var = "x", int mulx=0) const override
    {
        return TrendPoly<2>::translate(var, 0*mulx);
    }

    // Predict Y coordinate from X coordinate
    virtual Float predict(Float x) const override
    {
        return TrendPoly<2>::predict(x) / x;
    }
};

//using TrendQuadratic = TrendPoly<2>;

template<typename Float=double>
struct TrendLog: public Trend<Float>
{
    TrendPoly<1> slave;
    virtual void estimate(const std::vector<std::pair<Float, Float>>& values) override
    {
        std::vector<Float> x,y; x.reserve(values.size()); y.reserve(values.size());
        for(auto& v: values)
        {
            x.push_back(std::log(v.first));
            y.push_back(v.second);
        }
        slave.estimate(x, y);
        // Now, slave.predict(x) = y = a + b log(x)
    }
    virtual std::string translate(const std::string& var = "x", int mulx=0) const override
    {
        std::string s = slave.translate("\\ln "+var);
        if(mulx != 0)
            return Trend<Float>::powstr(mulx, var) + "\\left(" + s + "\\right)";
        return s;
    }
    virtual Float predict(Float x) const override
    {
        return slave.predict(std::log(x));
    }
    virtual const char* name() const { return "Logarithmic"; }
};

template<typename Float=double>
struct TrendPower: public Trend<Float>
{
    TrendPoly<1> slave;

    virtual void estimate(const std::vector<std::pair<Float, Float>>& values) override
    {
        std::vector<Float> x,y; x.reserve(values.size()); y.reserve(values.size());
        for(auto& v: values)
        {
            x.push_back(std::log(v.first));
            y.push_back(std::log(v.second));
        }
        slave.estimate(x, y);
        // Now, slave.predict(x) = log(y) = a + b log(x)
        //                                = a + log(x^b)
        //                             y  = exp(a + log(x^b))
        //                                = x^b * exp(a)
    }
    virtual std::string translate(const std::string& var = "x", int mulx = 0) const override
    {
        Float a = std::exp(slave.factors[0]);
        Float b = slave.factors[1] + mulx;

        std::string result;
        char Buf[32];
        std::sprintf(Buf, "%.*g", PRECISION, a);
        std::string k = Buf;
        if(!IsOne(k, 1))
        {
            if(IsOne(k, -1))
                result += '-';
            else
                result += Buf;
            result += " ";
        }
        result += Trend<Float>::powstr(b, var);
        return result;
    }
    virtual Float predict(Float x) const override
    {
        return std::exp(slave.predict(std::log(x)));
    }
    virtual const char* name() const { return "Power"; }
};
