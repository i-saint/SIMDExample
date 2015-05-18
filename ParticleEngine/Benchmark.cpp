#include "ParticleEngine.h"
#include <cstdio>
#include <ctime>


const int ParticleCount = 1024 * 4;
const int LoopCount = 16;

template<class F>
void Test(const F &f, const char *name)
{
    float average = 0.0f;
    peContext *ctx = peCreateContext(ParticleCount);
    for (int i = 0; i < LoopCount; ++i) {
        clock_t t = clock();
        f(ctx, 1.0f / 60.0f);
        average += float(clock() - t) / CLOCKS_PER_SEC * 1000.0f;
    }
    peDestroyContext(ctx);

    printf("%s: average %.2fms\n", name, average / LoopCount);
}

int main(int argc, char *argv[])
{
    Test(peUpdate_Plain, "Test_Plain");
    Test(peUpdate_SIMD, "Test_SIMD");
    Test(peUpdate_SIMDSoA, "Test_SIMDSoA");
    Test(peUpdate_ISPC, "Test_ISPC");
}
