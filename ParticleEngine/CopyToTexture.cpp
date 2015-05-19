#include <windows.h>
#include <d3d11.h>
#include "ParticleEngine.h"


peCopyToTextureBase* peCreateCopyToTextureD3D9(void *device);
peCopyToTextureBase* peCreateCopyToTextureD3D11(void *device);
peCopyToTextureBase* peCreateCopyToTextureOpenGL(void *device);

// Graphics device identifiers in Unity
enum GfxDeviceRenderer
{
    kGfxRendererOpenGL = 0,          // OpenGL
    kGfxRendererD3D9,                // Direct3D 9
    kGfxRendererD3D11,               // Direct3D 11
    kGfxRendererGCM,                 // Sony PlayStation 3 GCM
    kGfxRendererNull,                // "null" device (used in batch mode)
    kGfxRendererHollywood,           // Nintendo Wii
    kGfxRendererXenon,               // Xbox 360
    kGfxRendererOpenGLES,            // OpenGL ES 1.1
    kGfxRendererOpenGLES20Mobile,    // OpenGL ES 2.0 mobile variant
    kGfxRendererMolehill,            // Flash 11 Stage3D
    kGfxRendererOpenGLES20Desktop,   // OpenGL ES 2.0 desktop variant (i.e. NaCl)
    kGfxRendererCount
};

// Event types for UnitySetGraphicsDevice
enum GfxDeviceEventType {
    kGfxDeviceEventInitialize = 0,
    kGfxDeviceEventShutdown,
    kGfxDeviceEventBeforeReset,
    kGfxDeviceEventAfterReset,
};


extern peCopyToTextureBase *g_CopyToTexture;
void *g_graphicsDevice;
int g_deviceType;

peCLinkage peExport void* peGetGraphicsDevice() { return g_graphicsDevice; }
peCLinkage peExport int   peGetGraphicsDeviceType() { return g_deviceType; }
typedef void* (*peGetGraphicsDeviceT)();
typedef int (*peGetGraphicsDeviceTypeT)();


#define peSupportD3D11

peCLinkage peExport void UnitySetGraphicsDevice(void* device, int deviceType, int eventType)
{
    if (device == nullptr) { return; }

    g_graphicsDevice = device;
    g_deviceType = deviceType;
    if (eventType == kGfxDeviceEventInitialize) {
#ifdef peSupportD3D9
        if (deviceType == kGfxRendererD3D9)
        {
            // todo
        }
#endif // peSupportD3D9
#ifdef peSupportD3D11
        if (deviceType == kGfxRendererD3D11)
        {
            g_CopyToTexture = peCreateCopyToTextureD3D11(device);
        }
#endif // peSupportD3D11
#ifdef peSupportOpenGL
        if (deviceType == kGfxRendererOpenGL)
        {
            g_CopyToTexture = peCreateCopyToTextureOpenGL(device);
        }
#endif // peSupportOpenGL
    }

    if (eventType == kGfxDeviceEventShutdown) {
        delete g_CopyToTexture;
        g_CopyToTexture = nullptr;
    }
}

peCLinkage peExport void UnityRenderEvent(int eventID)
{
}


#ifndef peMaster

// PatchLibrary で突っ込まれたモジュールは UnitySetGraphicsDevice() が呼ばれないので、
// DLL_PROCESS_ATTACH のタイミングで先にロードされているモジュールからデバイスをもらって同等の処理を行う。
BOOL WINAPI DllMain(HINSTANCE module_handle, DWORD reason_for_call, LPVOID reserved)
{
    if (reason_for_call == DLL_PROCESS_ATTACH)
    {
        HMODULE m = GetModuleHandleA("ParticleEngine.dll");
        if (m) {
            auto p1 = (peGetGraphicsDeviceT)GetProcAddress(m, "peGetGraphicsDevice");
            auto p2 = (peGetGraphicsDeviceTypeT)GetProcAddress(m, "peGetGraphicsDeviceType");
            if (p1 && p2) {
                UnitySetGraphicsDevice(p1(), p2(), kGfxDeviceEventInitialize);
            }
        }
    }
    else if (reason_for_call == DLL_PROCESS_DETACH)
    {
    }
    return TRUE;
}

// "DllMain already defined in MSVCRT.lib" 対策
#ifdef _X86_
extern "C" { int _afxForceUSRDLL; }
#else
extern "C" { int __afxForceUSRDLL; }
#endif

#endif // peMaster
