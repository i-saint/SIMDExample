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



peCLinkage peExport peContext*  peCreateContext(int n);
peCLinkage peExport void        peDestroyContext(peContext *ctx);
peCLinkage peExport void        peUpdate_Plain(peContext *ctx, float dt);
peCLinkage peExport void        peUpdate_SIMD(peContext *ctx, float dt);
peCLinkage peExport void        peUpdate_SIMDSoA(peContext *ctx, float dt);
peCLinkage peExport void        peUpdate_ISPC(peContext *ctx, float dt);
