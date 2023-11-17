#include <set>
#include <map>
#include <vector>
#include <string>
#include <algorithm>
#include <fstream>

#define DO_SVG 0

std::set<float> Rows; // Jaksopituus
std::map<std::string/*key*/, std::map<float/*jakso*/, float/*arvo*/>> Columns;
std::vector<std::string> Headers;

static constexpr char delimiter = ',';
static std::vector<std::string> split(const std::string& s)
{
    std::size_t begin = 0;
    std::vector<std::string> result;
    while(begin != s.size())
    {
        std::size_t p = s.find(delimiter, begin);
        if(p == s.npos) { result.emplace_back(s.begin() + begin, s.end()); break; }
        result.emplace_back(s.begin() + begin, s.begin() + p);
        begin = p+1;
    }
    return result;
}
static void ReadCSV(const std::string& filename)
{
    std::vector<std::string> headers;
    std::map<std::string, std::vector<std::string>> result;
    std::ifstream f(filename);
    std::string line;
    while(std::getline(f, line))
    {
        if(headers.empty())
        {
            headers = split(line);
            for(const auto& h: headers) { result[h]; }
        }
        else
        {
            unsigned h=0;
            for(auto& s: split(line))
            {
                if(h > headers.size()) break;
                if(s[0] == '"')
                {
                    s.erase(0,1);
                    s.erase(s.size()-1,1);
                }
                result[headers[h]].push_back(std::move(s));
                ++h;
            }
        }
    }
    auto& firstcol = result[headers[0]];
    headers.erase(headers.begin(), headers.begin()+1); // Jaksopituus

    for(auto h: firstcol) Rows.insert(std::stod(h));
    for(auto& name: headers)
    {
        auto& srccolumn = result[name];
        auto& tgtcolumn = Columns[name];
        for(std::size_t a=0; a<firstcol.size(); ++a)
        {
            if(!srccolumn[a].empty())
            {
                tgtcolumn[std::stod(firstcol[a])] = std::stod(srccolumn[a]);
            }
        }
    }
    Headers = std::move(headers);
}
static void ToUTF8(char* target, char32_t c)
{
    alignas(16) static constexpr unsigned S[4] = {0x7F,0x3F,0x3F,0x3F}, q[4] = {0xFF,0,0,0}, o[4] = {0,0x80,0x80,0x80};
    unsigned n = (c >= 0x80u) + (c >= 0x800u) + (c >= 0x10000u);
    unsigned val = 0xF0E0C000u >> (n*8u);
    alignas(16) unsigned w[4]={val,val,val,val};
    #pragma omp simd
    for(unsigned m=0; m<4; ++m)
        w[m] = ((w[m] & q[m]) + ((c >> ((n-m)*6u) & S[m]) + o[m])) << (m*8u);
    unsigned sum = 0;
    #pragma omp simd
    for(unsigned m=0; m<4; ++m) sum += w[m];
    for(unsigned m=0; m<4; ++m) target[m] = sum >> (m*8);
    target[n+1] = 0;
}


#include <cairo.h>
#include <cairo-svg.h>
#include <cmath>
#include <thread>
#include <gd.h>

#include "trend.hh"

struct Style
{
    double width, r,g,b, dashes[16];
    unsigned ndashes;

    double symsize; char32_t symbol; float symxdelta, symydelta;
};
static std::map<std::string, Style> styles
{
    { "DFT",           Style{ 2.0, .6, 0, .4,  {1,1}, 2, 15, U'◆',-8,4   } },
    { "FFTW-est",     Style{ 0.5, 0, .4, .8, {1,1}, 2, 10, U'■',-2,4   } },
    { "FFTW-exh",     Style{ 0.15,.8, .6, .2, {1,1}, 2, 10, U'W',-2,4   } },
    { "Radix2",        Style{ 3.0, 0, .6, .2, {1,1}, 2, 60, U'◉',-15,18 } },
    { "Radix2_fast",   Style{ 3.0, 0, .8, .2, {1,2}, 2, 60, U'◈',-35,18 } },
    { "TukeyR1_RecAlways", Style{ 0.8, .4, .2, .1, {1,1}, 2, 10, U'■',-2,4   } },
    { "TukeyR2_NoRec", Style{ 0.8, .2, .4, .1, {1,1}, 2, 10, U'■',-2,4   } },
    //
    { "Rader",         Style{ 0.6, .9, .4, .2, {1,1}, 2, 20, U'▲',-6,6 } },
    { "Bluestein_prime",     Style{ 1.2, 0, 0, .70,  {1,1}, 2, 20, U'▼',-6,6 } },
    { "Bluestein2_prime",    Style{ 1.9, 0, 0, .10,  {1,1}, 2, 20, U'◀',-12,6 } },
    { "Bluestein3_prime",    Style{ 1.9, 0, .6, .10,  {1,1}, 2, 20, U'▷',-12,6 } },
    { "Bluestein_any",       Style{ 1.2, 0, 0, .70,  {1,1}, 2, 20, U'▼',-6,6 } },
    { "Bluestein2_any",      Style{ 1.9, 0, 0, .10,  {1,1}, 2, 10, U'◀',-6,3 } },
    { "Bluestein3_any",      Style{ 1.9, 0, .6, .10,  {1,1}, 2, 10, U'▷',-6,3 } },
    //
    { "any",           Style{ 1.0, 1,1,1,      {1,1}, 2, 10, U'Z',-2,4   } },
    { "Any-est",       Style{ 1.0, .5,.2,.2,   {1,1}, 2, 10, U'X',-2,4   } },
};

#include <atomic>
static cairo_surface_t* latex_render(const std::string& formula)
{
    static std::atomic<unsigned> counter{};
    unsigned myid = ++counter;

    std::string prefix = "rendertemp_" + std::to_string(myid);
    std::string outfn  = "rendertemp_" + std::to_string(myid) + ".png";
    std::string infn = "rendertemp_" + std::to_string(myid) + ".tex";

    std::string tex = R"(\documentclass{standalone}
\begin{document}$\displaystyle )"
                   + formula + R"($\text{ ns}\end{document})";
    FILE* fp = std::fopen(infn.c_str(), "w");
    std::fwrite(tex.c_str(), tex.size(), 1, fp);
    std::fclose(fp);
    std::string command = "pdflatex --interaction=nonstopmode -jobname=" + prefix + " " + infn + " 2>/dev/null";
    std::system(command.c_str());
    std::remove(infn.c_str());
    std::remove((prefix+".aux").c_str());
    std::remove((prefix+".log").c_str());
    std::string command2 = "convert -background transparent -alpha remove -density 400 " + prefix + ".pdf " + outfn;
    std::system(command2.c_str());
    std::remove((prefix+".pdf").c_str());
    auto resvalue1 = cairo_image_surface_create_from_png(outfn.c_str());
    std::remove(outfn.c_str());
    return resvalue1;
}

template<typename F1>
static void Render(
    const std::string& outfn,
    F1 key_selector,
    std::vector<std::pair<std::string, std::vector<Trend<double>*>>> columns,
    unsigned xres, unsigned yres,
    float first_x, float last_x,
    float first_y, float last_y
)
{
#if !DO_SVG
    std::vector<unsigned> pixels(xres*yres);
    std::fill_n(&pixels[0], xres*yres, 0xFFFFFFFFu);
#endif
    int grid_start_x = 160, grid_end_x = xres-10;
    int grid_start_y = 10,  grid_end_y = yres-100, nub = 9;
    float grid_width = 1;

    //int legend_start_x = 2200; int legend_line_width  = 100;
    //int legend_start_y = 2000; int legend_item_height = -60;
    int legend_start_x =  190; int legend_line_width  = 100;
    int legend_start_y =   80; int legend_item_height =  60;

#if !DO_SVG
    cairo_surface_t* sfc = cairo_image_surface_create_for_data
        ((unsigned char*)&pixels[0], CAIRO_FORMAT_ARGB32, xres, yres, xres*sizeof(unsigned));
#else
    cairo_surface_t* sfc = cairo_svg_surface_create(
        outfn.c_str(),
        xres, yres);
#endif

    cairo_t* c   = cairo_create(sfc);
    cairo_set_antialias(c, CAIRO_ANTIALIAS_GRAY);

    auto ycoord = [=](float value) // Duration ns to Ycoord
    {
        float log_value = value > 0 ? std::log(value) : 0;
        float log_starty = std::log(first_y), log_endy = std::log(last_y);
        return (log_value - log_starty) * (grid_end_y - grid_start_y) / (log_endy - log_starty) + grid_start_y;
    };
    auto xcoord = [=](float value) // Sample size to Xcoord
    {
        float log_value = value > 0 ? std::log(value) : 0;
        float log_startx = std::log(first_x), log_endx = std::log(last_x);
        return (log_value - log_startx) * (grid_end_x - grid_start_x) / (log_endx - log_startx) + grid_start_x;
    };
    auto xreverse = [=](float value) // Xcoord to sample size
    {
        float log_startx = std::log(first_x), log_endx = std::log(last_x);
        float log_value = ((log_startx-log_endx)*value - grid_end_x*log_startx + grid_start_x*log_endx) / (grid_start_x-grid_end_x);
        return std::exp(log_value);
    };
    auto Line = [&](unsigned x1,unsigned y1, unsigned x2,unsigned y2, float width)
    {
        cairo_new_path(c);
        cairo_set_line_width(c, width);
        cairo_move_to(c, x1,y1);
        cairo_line_to(c, x2,y2);
        cairo_stroke(c);
        cairo_close_path(c);
    };
    auto Hline = [&](unsigned x1,unsigned x2,unsigned y, float width) { Line(x1,y,x2,y, width); };
    auto Vline = [&](unsigned y1,unsigned y2,unsigned x, float width) { Line(x,y1,x,y2, width); };
#if !DO_SVG
    auto SaveFrame = [xres,yres](const std::string& filename, const unsigned* pixels)
    {
        static std::mutex lock;
        static struct { std::vector<unsigned> pixels;
                        std::thread saver;
                      } data{ {}, {} };
        std::lock_guard<std::mutex> lk(lock);
        if(data.saver.joinable()) data.saver.join();
        if(!pixels) return;
        data.pixels.assign(pixels, pixels + xres*yres);
        data.saver = std::thread([filename,xres,yres]()
        {
            fprintf(stderr, "Saving %s...\n", filename.c_str());
            gdImagePtr im = gdImageCreateTrueColor(xres, yres);
            ////BgdImageSaveAlpha(im, 1);
            gdImageAlphaBlending(im, 0);
            for(unsigned p=0, y=0; y<yres; ++y)
                for(unsigned x=0; x<xres; ++x, ++p)
                    gdImageSetPixel(im, x,y, data.pixels[p]);// ^ 0xFF000000

            std::FILE* fp = std::fopen(filename.c_str(), "wb");
            if(!fp) std::perror(filename.c_str());
            if(fp)
            {
                gdImagePng(im, fp);
                std::fclose(fp);
            }
            gdImageDestroy(im);
        });
    };
#endif
    auto Save = [&]()
    {
        cairo_surface_flush(sfc);
#if DO_SVG
        cairo_surface_finish(sfc);
#else
        SaveFrame(outfn, &pixels[0]);
#endif
    };

    /* Render grid */
    cairo_set_source_rgb(c, .5, .5, .5); // grid line color
    cairo_select_font_face(c, "sans-serif", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_NORMAL);
    float grid_font_size = 45;
    Hline(grid_start_x,grid_end_x, grid_start_y, grid_width);
    Hline(grid_start_x,grid_end_x, grid_end_y,   grid_width);
    Vline(grid_start_y,grid_end_y, grid_start_x, grid_width);
    Vline(grid_start_y,grid_end_y, grid_end_x,   grid_width);

    cairo_set_font_size(c, grid_font_size*3/4);
    cairo_move_to(c, grid_start_x+(grid_end_x-grid_start_x)*5/12 -40, grid_end_y+grid_font_size*2.1);
    cairo_show_text(c, "Syötteen pituus [näytettä]");

    cairo_save(c);
    cairo_move_to(c, grid_start_x - 120, grid_end_y - 100);
    cairo_rotate(c, -3.141592653/2);
    cairo_show_text(c, "Kesto [ns / näyte], pienempi on parempi");
    cairo_restore(c);

    cairo_set_font_size(c, grid_font_size);
    for(int ly = std::round(std::log10(std::min(last_y, first_y)))-1,
            ey = std::round(std::log10(std::max(last_y, first_y)))+1;
            ly <= ey;
            ++ly)
        for(int n=1; n<=9; ++n)
        {
            float val = std::exp10(ly)*n;
            float y = ycoord(val);
            if(y > grid_start_y-1 && y < grid_end_y+1)
            {
                cairo_set_source_rgb(c, .5, .5, .5); // grid line color
                if(n == 1) Hline(grid_start_x-nub*2, grid_end_x,       y, grid_width);
                else       Hline(grid_start_x,       grid_start_x+nub, y, grid_width);
                if(n == 1)
                {
                    cairo_move_to(c, grid_start_x - ly*grid_font_size*2/5 - nub*6, y + grid_font_size*2/5);
                    cairo_set_source_rgb(c, 0,0,0); // text color
                    cairo_show_text(c, std::to_string(int(std::round(val))).c_str());
                }
            }
        }
    for(int lx = std::round(std::log10(std::min(last_x, first_x)))-1,
            ex = std::round(std::log10(std::max(last_x, first_x)))+1;
            lx <= ex;
            ++lx)
        for(int n=1; n<=9; ++n)
        {
            float val = std::exp10(lx) * n;
            float x   = xcoord(val);
            if(x > grid_start_x-1 && x < grid_end_x+1)
            {
                cairo_set_source_rgb(c, .5, .5, .5); // grid line color
                if(n == 1) Vline(grid_start_y, grid_end_y+nub*2, x, grid_width);
                else       Vline(grid_end_y,   grid_end_y+nub, x, grid_width);
                if(n == 1)
                {
                    cairo_move_to(c, x - lx*grid_font_size*2/5, grid_end_y + nub*2 + grid_font_size);
                    cairo_set_source_rgb(c, 0,0,0); // text color
                    cairo_show_text(c, std::to_string(int(std::round(val))).c_str());
                }
            }
        }

    cairo_rectangle(c, grid_start_x,            grid_start_y,
                       grid_end_x-grid_start_x, grid_end_y-grid_start_y);
    cairo_clip(c);

    auto LegendAdd = [&](const std::string& name,
                         const std::string& formula,
                         const char* symbol, float symxdelta=0, float symydelta=0,
                         const double* dashes,unsigned ndashes, float orig_width)
    {
        float font_size = 45;
        float line_width = 4;
        cairo_set_line_width(c, line_width);

        double dashes1[32];
        for(unsigned n=0; n<ndashes; ++n)
            dashes1[n] = dashes[n] * line_width/std::max(1.5f,orig_width);
        cairo_set_dash(c, dashes1, ndashes, 0);

        cairo_set_font_size(c, 40);
        cairo_new_path(c);
        cairo_move_to(c, legend_start_x, legend_start_y - font_size/3);
        cairo_line_to(c, legend_start_x + legend_line_width, legend_start_y - font_size/3);
        cairo_stroke(c);
        cairo_close_path(c);
        cairo_move_to(c, legend_start_x + legend_line_width / 2 + symxdelta, legend_start_y - font_size/3 + symydelta);
        cairo_show_text(c, symbol);
        cairo_select_font_face(c, "Liberation Sans Narrow", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_NORMAL);
        cairo_set_font_size(c, font_size);
        cairo_move_to(c, legend_start_x + legend_line_width, legend_start_y);
        cairo_show_text(c, name.c_str());
        int x,y;
        cairo_move_to(c, x=legend_start_x + legend_line_width + font_size*14, y=legend_start_y);
        if(!formula.empty())
        {
            std::string form = formula;
            for(;;)
            {
                auto s = form.find('.');
                if(s == form.npos) break;
                form.replace(s, 1, "{,}", 3);
            }
            for(;;)
            {
                auto s1 = form.find("e+"), s2 = form.find("e-");
                if(s1 == form.npos && s2 == form.npos) break;
                std::size_t s = s1, n = 2;
                if(s2 < s1 || s1 == form.npos) { n = 1; s = s2; }
                form.replace(s, n, "\\cdot 10^{", 10); s += 10;
                if(form[s] == '-') ++s;
                while(form[s] == '0') form.erase(s, 1);
                while(form[s] >= '0' && form[s] <= '9') ++s;
                form.replace(s, 0, "}", 1);
            }
            auto source = latex_render(form);
            cairo_set_source_surface(c, source, x - font_size*4.2,y - 40);
            cairo_paint(c);
            //cairo_show_text(c, formula.c_str());
        }
        legend_start_y += legend_item_height;
    };

    unsigned linecounter = 0;
    for(auto& [column_name, trends]: columns)
    {
        printf("Column %s selected\n", column_name.c_str());
        auto& column = Columns[column_name];

        cairo_new_path(c);
        //cairo_stroke(c);
        cairo_set_line_join(c, CAIRO_LINE_JOIN_ROUND);

        auto s = styles.find(column_name);
        Style style{ 0.5, 1,0,0, {1,1}, 2, 50, U'A',0,0 };
        if(s != styles.end()) style = s->second;

        cairo_set_source_rgb(c, style.r, style.g, style.b);
        cairo_set_line_width(c, style.width);
        cairo_set_dash(c, style.dashes, style.ndashes, 0);

        cairo_select_font_face(c, "NotoSansMono", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_NORMAL);
        cairo_set_font_size(c, style.symsize);
        char Buf[16] = {};
        ToUTF8(Buf, style.symbol);
        bool first = false;

        std::vector<std::pair<double,double>> values;

        for(auto i = column.begin(); i != column.end(); ++i)
            if(key_selector(i->first))
            {
                float size = i->first, time = i->second;
                // Time is microseconds (absolute). Translate to nanoseconds per sample.
                time = time * 1e6 / size;
                float x = xcoord(size);  // sample size
                float y = ycoord(time); // time
                if(!first)
                {
                    cairo_line_to(c, x,y);
                }
                cairo_move_to(c, x + style.symxdelta, y + style.symydelta);
                cairo_show_text(c, Buf);
                cairo_move_to(c, x,y);
                first = false;

                values.emplace_back(size, time);
            }
        cairo_stroke(c);
        cairo_close_path(c);

        LegendAdd(column_name, "", Buf, style.symxdelta*40/style.symsize, style.symydelta*40/style.symsize,
                  style.dashes,style.ndashes, style.width);

        static const char* const dashpatterns[16] =
        {
            "-",
            //
            "-.",
            "-..",
            "--.",
            //
            "-...",
            "---.",
            //
            "---..",
            "--...",
            "-....",
        };

        for(auto trend: trends)
        {
            const char* pattern = dashpatterns[linecounter++];

            std::vector<std::pair<double,double>> values2 = values;

            if(column_name == "FFTW-exh")
                values2.erase(std::remove_if(values2.begin(), values2.end(),
                    [](auto& a) { return !(a.first <= 30000); }), values2.end());

            //if(column_name == "FFTW-est")
            //    values2.erase(std::remove_if(values2.begin(), values2.end(),
            //        [](auto& a) { return !(a.first <= 30000); }), values2.end());

            trend->estimate(values2);

            std::string formula = trend->translate("N", 1);
            printf("t = %s (pattern=%s)\n", formula.c_str(), pattern);

            cairo_new_path(c);
            cairo_stroke(c);
            cairo_set_line_join(c, CAIRO_LINE_JOIN_ROUND);
            cairo_set_source_rgb(c, 0,0,0);
            float width = 6;
            if(column_name == "DFT") width = 2;
            cairo_set_line_width(c, width);
            double dashes[32];
            unsigned n = 0;
            for(; *pattern; ++pattern)
            {
                dashes[n++] = (*pattern == '-' ? 6 : 1) * width;
                dashes[n++] = 2*width;
            }
            double dist = 0;
            float prevy = 0;
            cairo_set_dash(c, dashes, n, 0);
            for(int x = grid_start_x; x <= grid_end_x; ++x)
            {
                float y = ycoord(trend->predict(xreverse(x)));
                if(x > grid_start_x)
                {
                    dist += std::sqrt(1*1 + (y-prevy)*(y-prevy));
                    cairo_set_dash(c, dashes, n, dist);
                    cairo_line_to(c, x, y);
                    cairo_stroke(c);
                }
                cairo_move_to(c, x, y);
                prevy = y;
            }
            cairo_close_path(c);

            Buf[0] = '\0';
            cairo_set_dash(c, dashes, n, 0);
            LegendAdd(std::string(trend->name()) + "(" + column_name + ")", "t \\approx " + formula, Buf,0,0,
                      dashes,n, width);
        }
    }

    Save();
    /* Terminate saving thread at end */
#if !DO_SVG
    SaveFrame("", nullptr);
#endif
}

#include "factor.hh"
int main()
{
    ReadCSV("/home/bisqwit/tmp.csv");

  #pragma omp parallel sections
  {
#if 1
    Render("img/dft-vs-radix2.png",
        [](float /*where*/) { return true; },
        {
            {"DFT",    {new TrendLinear}},
            {"Radix2", {new TrendLog}},
        },
           3840, 1000,
           1,1.35e8, 20e3,4
    );
#endif
   #pragma omp section
    Render("img/dft-vs-radix2-power.png",
        [](float /*where*/) { return true; },
        {
            {"DFT",    {new TrendPower}},
            {"Radix2", {new TrendPower}},
        },
           3840, 1000,
           1,1.35e8, 20e3,4
    );
   #pragma omp section
    Render("img/dft-vs-radix2-quadratic.png",
        [](float /*where*/) { return true; },
        {
            {"DFT",    {new TrendQuadratic}},
            {"Radix2", {new TrendQuadratic}},
        },
           3840, 1000,
           1,1.35e8, 20e3,4
    );
#if 1
   #pragma omp section
    Render("img/radix2-vs-radix2fast.png",
        [](float /*where*/) { return true; },
        {
            {"Radix2",      {new TrendLog}},
            {"Radix2_fast", {new TrendLog}},
        },
           3840, 800,
           1,1.35e8, 650,4
    );
   #pragma omp section
    Render("img/radix2-vs-radix2fast-power.png",
        [](float /*where*/) { return true; },
        {
            {"Radix2",      {new TrendPower}},
            {"Radix2_fast", {new TrendPower}},
        },
           3840, 800,
           1,1.35e8, 650,7
    );
   #pragma omp section
    Render("img/radix2fast-vs-tukey.png",
        [](float /*where*/)    { return true; },
        {
            {"Radix2_fast",    {new TrendPower}},
            //{"TukeyR2_NoRec",  {new TrendPower}},
            {"TukeyR1_RecAlways",{new TrendPower}},
            {"DFT",            {new TrendPower}},
        },
           3840, 1600,
           1,600e3, 200e3,7
    );
   #pragma omp section
    Render("img/radix2fast-vs-tukey-small.png",
        [](float where)        { return small_factors_only(int(where)); },
        {
            {"Radix2_fast",    {new TrendPower}},
            {"TukeyR2_NoRec",  {new TrendPower}},
            {"TukeyR1_RecAlways",{new TrendPower}},
            {"DFT",            {new TrendPower}},
        },
           3840, 1000,
           1,500e3, 4e3,7
    );
   #pragma omp section
    Render("img/radix2fast-vs-tukey-small-without-dft.png",
        [](float where)        { return small_factors_only(int(where)); },
        {
            {"Radix2_fast",    {new TrendLog}},
            //{"TukeyR2_NoRec",  {new TrendLog}},
            {"TukeyR1_RecAlways",{new TrendPower}},
        },
           3840, 1000,
           1,500e3, 4e3,7
    );
   #pragma omp section
    Render("img/radix2fast-dft-vs-rader-vs-bluestein.png",
        [](float w)             { return int(smallest_factor(int(w)))==int(w) || small_factors_only(int(w)); },
        {
            //{"DFT",            {new TrendPower}},
            {"Rader",           {new TrendPower}},
            {"Bluestein_prime", {new TrendPower}},
            {"Bluestein2_prime",{new TrendPower}},
            {"Bluestein3_prime",{new TrendPower}},
            {"Radix2_fast",     {new TrendPower}},
        },
           3840, 1600,
           1,1600e3, 6e3,7
    );
   #pragma omp section
    Render("img/any-vs-fftw.png",
        [](float /*where*/)   { return true; },
        {
            {"DFT",            {/*new TrendPower*/}},
            {"Bluestein2_any",  {new TrendPower}},
            //{"any",            {new TrendPower}},
            {"Any-est",        {new TrendPower}},
            //{"any",            {new TrendLog}},
            {"FFTW-est",        {new TrendPower}},
            //{"FFTW-exh",        {new TrendPower}},
            //{"FFTW3",          {new TrendLog}},
            {"Radix2_fast",    {new TrendPower}},
        },
           3840, 2000,
           1,1600e3, 6e3,1
    );
#endif
  }
}
