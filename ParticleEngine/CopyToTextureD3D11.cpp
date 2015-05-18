#include <windows.h>
#include <d3d11.h>
#include "ParticleEngine.h"


inline int ceildiv(int v, int d) { return v / d + (v%d == 0 ? 0 : 1); }
#define peSafeRelease(obj) if(obj) { obj->Release(); obj=nullptr; }



class peCopyToTextureD3D11 : public peCopyToTextureBase
{
public:
    peCopyToTextureD3D11(void *dev);
    virtual ~peCopyToTextureD3D11();
    virtual void copy(void *texptr, int width, int height, const void *data, int datasize);

private:
    ID3D11Device        *m_pDevice;
    ID3D11DeviceContext *m_pImmediateContext;
};


peCopyToTextureBase* peCreateCopyToTextureD3D11(void *device) { return new peCopyToTextureD3D11(device); }

peCopyToTextureD3D11::peCopyToTextureD3D11(void *dev)
    : m_pDevice(nullptr)
    , m_pImmediateContext(nullptr)
{
    m_pDevice = (ID3D11Device*)dev;
    m_pDevice->GetImmediateContext(&m_pImmediateContext);
}

peCopyToTextureD3D11::~peCopyToTextureD3D11()
{
    peSafeRelease(m_pImmediateContext);
}


void peCopyToTextureD3D11::copy(void *texptr, int width, int height, const void *dataptr, int datasize)
{
    ID3D11Texture2D *tex = (ID3D11Texture2D*)texptr;

    D3D11_BOX box;
    box.left = 0;
    box.right = width;
    box.top = 0;
    box.bottom = ceildiv(datasize / sizeof(float4), width);
    box.front = 0;
    box.back = 1;
    m_pImmediateContext->UpdateSubresource(tex, 0, &box, dataptr, datasize, 0);
}

