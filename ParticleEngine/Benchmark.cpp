#include "ParticleEngine.h"
#include <cstdio>
#include <ctime>


const int ParticleCount = 1024 * 4;
const int LoopCount = 16;

void Test(peUpdateRoutine r, bool mt, const char *name)
{
    float average = 0.0f;
    peContext *ctx = peCreateContext(ParticleCount);
    peEnableMultiThreading(ctx, mt);
    peSetUpdateRoutine(ctx, r);
    for (int i = 0; i < LoopCount; ++i) {
        clock_t t = clock();
        peUpdate(ctx, 1.0f / 60.0f);
        average += float(clock() - t) / CLOCKS_PER_SEC * 1000.0f;
    }
    peDestroyContext(ctx);

    printf("%s: average %.2fms\n", name, average / LoopCount);
}

int main(int argc, char *argv[])
{
    Test(peE_Plain,   false, "Test_Plain (ST)  " );
    Test(peE_SIMD,    false, "Test_SIMD (ST)   " );
    Test(peE_SIMDSoA, false, "Test_SIMDSoA (ST)" );
    Test(peE_ISPC,    false, "Test_ISPC (ST)   " );

    Test(peE_Plain,    true, "Test_Plain (MT)  " );
    Test(peE_SIMD,     true, "Test_SIMD (MT)   " );
    Test(peE_SIMDSoA,  true, "Test_SIMDSoA (MT)" );
    Test(peE_ISPC,     true, "Test_ISPC (MT)   " );
}
