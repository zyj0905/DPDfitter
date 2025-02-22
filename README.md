# DPD fitter

The repository presents two versions of code:
- `gpuversion` is used to perform fit in the GPU, hence the part related to amplitude is coded in CUDA.
- `cpuversion` is used to sampling MC based on the fit result in CPU.

The Dalitz-plot decomposition related parts are:
1. `DPD.cxx`, used to calculate the part of Eq.7 in your article.
2. `Amplitude.cxx`, used to perform sum in Eq. 3 in your article.
3. `Event.cxx`, saving the event information and calculating the angles (alignment, scattering, wigner rotation...)

## Notebooks

Notebooks are run with [`Pluto`](https://github.com/fonsp/Pluto.jl).

Follow [simple instructions](https://plutojl.org/#install) to run:
1. Install julia
2. Add pluto
3. Open notebook
