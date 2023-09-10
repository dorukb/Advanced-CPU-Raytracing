#ifndef __DORKTRACER_PATHTRACINGPARAMS__
#define __DORKTRACER_PATHTRACINGPARAMS__

namespace DorkTracer
{
    class RendererParams
    {

    public:
        bool pathTracingEnabled = false;
        bool sampleUniform = false;
        bool sampleImportance = false;
        bool nextEventEstimationEnabled = false;
        bool russianRouletteEnabled = false;

        void SetParams(bool enablePathTracing, bool enableImportanceSampling, bool enableNextEventEstimation, bool enableRussianRoulette)
        {
            this->pathTracingEnabled = enablePathTracing;
            this->sampleUniform = !enableImportanceSampling;
            this->sampleImportance = enableImportanceSampling;
            this->nextEventEstimationEnabled = enableNextEventEstimation;
            this->russianRouletteEnabled =  enableRussianRoulette;
        }
    };
}

#endif