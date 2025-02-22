# DPD fitter

The repository presents two versions of code:
- `gpuversion` is used to perform fit in the GPU, hence the part related to amplitude is coded in CUDA.
- `cpuversion` is used to sampling MC based on the fit result in CPU.

The Dalitz-plot decomposition related parts are:
1. `DPD.cxx`, used to calculate the part of Eq.7 in your article.
2. `Amplitude.cxx`, used to perform sum in Eq. 3 in your article.
3. `Event.cxx`, saving the event information and calculating the angles (alignment, scattering, wigner rotation...)

A example with single event for your check is attached in this email. In this sample,
- I am using the `gamma*->D*Dpi`, with `Zc3900->D*D` as a check.
- In the dir `DstpiD`, the order of daughter particles are `Dst(1) pi(2) D(3)`, and the isobar is set to be `(1,3)`,
- in the dir `piDDst`, the order of daughter particles are `pi(1) D(2) Dst(3)`, and the isobar is set to be `(2,3)`.

The resultant amplitudes are found to be not consistent.
