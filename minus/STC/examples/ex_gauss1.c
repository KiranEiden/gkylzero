#include <stdio.h>
#include <time.h>
#include <math.h>
#include <stc/crandom.h>
#include <stc/cstr.h>
#include <stc/cmap.h>
#include <stc/cvec.h>

// Declare int -> int hashmap. Uses typetag 'ii' for ints.
using_cmap(ii, int32_t, size_t);

// Declare int vector with map entries that can be sorted by map keys.
typedef struct {int first; size_t second;} mapval;
static int compare(mapval *a, mapval *b) {
    return c_default_compare(&a->first, &b->first);
}

using_cvec(pair, mapval, compare);

int main()
{
    enum {N = 10000000};
    const double Mean = -12.0, StdDev = 6.0, Scale = 74;

    printf("Demo of gaussian / normal distribution of %d random samples\n", N);

    // Setup random engine with normal distribution.
    uint64_t seed = time(NULL);
    stc64_t rng = stc64_init(seed);
    stc64_normalf_t dist = stc64_normalf_init(Mean, StdDev);

    // Create and init histogram vec and map with defered destructors:
    c_forauto (cvec_pair, histvec)
    c_forauto (cmap_ii, histmap)
    {
        c_forrange (N) {
            int index = (int) round( stc64_normalf(&rng, &dist) );
            cmap_ii_emplace(&histmap, index, 0).ref->second += 1;
        }

        // Transfer map to vec and sort it by map keys.
        c_foreach (i, cmap_ii, histmap)
            cvec_pair_push_back(&histvec, (mapval){i.ref->first, i.ref->second});

        cvec_pair_sort(&histvec);

        // Print the gaussian bar chart
        c_forauto (cstr, bar)
        c_foreach (i, cvec_pair, histvec) {
            size_t n = (size_t) (i.ref->second * StdDev * Scale * 2.5 / (float)N);
            if (n > 0) {
                cstr_resize(&bar, n, '*');
                printf("%4d %s\n", i.ref->first, bar.str);
            }
        }
    }
}