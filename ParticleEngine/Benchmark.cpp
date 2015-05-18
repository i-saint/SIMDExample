#include "ParticleEngine.h"
#include <cstdio>
#include <ctime>


const int ParticleCount = 1024 * 4;
const int LoopCount = 16;

void Test(peUpdateRoutine r, const char *name)
{
    float average = 0.0f;
    peContext *ctx = peCreateContext(ParticleCount);
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
    Test(peE_Plain,     "Test_Plain"    );
    Test(peE_SIMD,      "Test_SIMD"     );
    Test(peE_SIMDSoA,   "Test_SIMDSoA"  );
    Test(peE_ISPC,      "Test_ISPC"     );
}
