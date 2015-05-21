#include "Vector.h"

#define peCLinkage extern "C"
#ifdef _MSC_VER
#define peExport __declspec(dllexport)
#else
#define peExport
#endif


const int peDefaultAlign = 32;
void* peAlignedAlloc(size_t size, size_t align = peDefaultAlign);
void peAlignedFree(void *p);

struct peParticle
{
    float4 position;
    float4 velocity;
};

struct peParticleSoA
{
    float *pos_x;
    float *pos_y;
    float *pos_z;
    float *vel_x;
    float *vel_y;
    float *vel_z;
};

class peContext;


class peCopyToTextureBase
{
public:
    virtual ~peCopyToTextureBase() {}
    virtual void copy(void *texptr, int width, int height, const void *data, int datasize) = 0;
};
peCopyToTextureBase* peCreateCopyToTextureD3D11(void *device);

typedef void(__stdcall *CSUpdateRoutine)(float dt, int begin, int end);

enum peUpdateRoutine
{
    peE_Plain,
    peE_SIMD,
    peE_SIMDSoA,
    peE_ISPC,
    peE_CSharp,
};

struct peParams
{
    peUpdateRoutine routine;
    bool multi_threading;
    float particle_size;
    float pressure_stiffness;
    float wall_stiffness;
    CSUpdateRoutine update_velocity;
    CSUpdateRoutine update_position;
};


peCLinkage peExport peContext*  peCreateContext(int n);
peCLinkage peExport void        peDestroyContext(peContext *ctx);

peCLinkage peExport void        peSetParams(peContext *ctx, peParams *v);
peCLinkage peExport void        peUpdate(peContext *ctx, float dt);
peCLinkage peExport void        peCopyDataToTexture(peContext *ctx, void *texture, int width, int height);
peCLinkage peExport peParticle* peGetParticles(peContext *ctx);
peCLinkage peExport void        peResetParticles(peContext *ctx);

peCLinkage peExport const char* peBenchmark(peContext *ctx, int loop_count);
