# DPD Fitter

**DPD Fitter** is a C++ implementation of three-body decay amplitudes based on the formalism developed in the [Dalitz-Plot Decomposition paper](https://inspirehep.net/literature/1758460).

The project provides:

- a CUDA/GPU implementation for amplitude fitting;
- a CPU implementation for evaluating fitted amplitudes and sampling Monte Carlo events;
- JSON-based configuration and parameter files;
- support for simultaneous fits to multiple decay channels;
- fit-fraction and interference-term calculations.

## Contents

- [DPD Fitter](#dpd-fitter)
  - [Contents](#contents)
  - [Validation](#validation)
  - [Examples](#examples)
  - [JSON configuration](#json-configuration)
    - [Spin-density matrix](#spin-density-matrix)
    - [Mother particle](#mother-particle)
    - [Daughter particles](#daughter-particles)
    - [Secondary decay](#secondary-decay)
    - [Likelihood model](#likelihood-model)
    - [Resonances and decay chains](#resonances-and-decay-chains)
    - [Input ROOT samples](#input-root-samples)
    - [Minuit strategy](#minuit-strategy)
    - [Parameter files](#parameter-files)
    - [Fit fractions](#fit-fractions)
    - [Output ROOT files](#output-root-files)
  - [Running a fit](#running-a-fit)
  - [Simultaneous fits and shared parameters](#simultaneous-fits-and-shared-parameters)
  - [Sampling Monte Carlo events](#sampling-monte-carlo-events)
  - [Citation](#citation)

## Validation

The implementation has been checked in several complementary ways:

- numerical comparisons with established amplitude-analysis tools such as **TF-PWA** for selected decay chains;
- validation under permutations of the final-state particle ordering;
- independent source-level review with ChatGPT 5.6 sol.

> The code was developed and validated by the authors before AI-assisted code review was added. Nevertheless, users are encouraged to perform process-specific checks before applying the fitter to a new decay channel or lineshape.

## Examples

Examples for the process

\[
e^+e^- \to D^*D\pi
\]

are organized as follows in the full project distribution:

```text
example/
├── runFitinGPU/       # Fit data with a GPU device
└── runSamplinginCPU/  # Sample MC events using fitted parameters
```

The GPU example demonstrates how to configure and run a fit. The CPU example demonstrates how to apply the fitted amplitude to reconstructed and truth-level phase-space MC samples.

## JSON configuration

Input configuration files and parameter lists use JSON. The examples below follow the structure of `Dst0DmPip.json`.

### Spin-density matrix

```json
"SDM": [1, 0, 0, 0, 0, 0, 0, 0, 1]
```

`SDM` is the spin-density matrix of the mother particle, stored as a flattened array. In this example it describes a virtual photon produced in an \(e^+e^-\) collision without longitudinal polarization.

For a mother particle with spin \(J\), the array must contain

\[
(2J+1)^2
\]

elements.

### Mother particle

```json
"Mom": {
  "name": "vpho",
  "spin": 1,
  "p-parity": -1
}
```

| Field | Meaning |
|---|---|
| `name` | Particle label used internally and in generated parameter names |
| `spin` | Particle spin; integer and half-integer values are supported |
| `p-parity` | Intrinsic parity, normally `1` or `-1` |

### Daughter particles

```json
"Daus": {
  "Dau1": {
    "name": "Dst0",
    "spin": 1,
    "p-parity": -1
  },
  "Dau2": {
    "name": "Dm",
    "spin": 0,
    "p-parity": -1
  },
  "Dau3": {
    "name": "Pip",
    "spin": 0,
    "p-parity": -1
  }
}
```

The order `Dau1`, `Dau2`, `Dau3` defines the particle indices 1, 2, and 3 throughout the configuration. The same order must be used for all four-momentum branches in `p4_final`.

### Secondary decay

```json
"set_sec_decay": {
  "is_set_sec": true,
  "type_sec": 0,
  "idx_sec": 1,
  "p4_sec": ["p4_sub_D", "p4_sub_Pi"]
}
```

This block describes a subsequent decay of one of the three daughter particles. In the example,

\[
D^* \to D\pi.
\]

| Field | Meaning |
|---|---|
| `is_set_sec` | Enable or disable the secondary-decay amplitude |
| `type_sec` | Secondary-decay model index implemented in `Amplitude.cu` |
| `idx_sec` | Index of the daughter that decays: `1`, `2`, or `3` |
| `p4_sec` | ROOT branches containing the two secondary-daughter four-momenta |

Currently, `type_sec = 0` implements a vector-to-two-pseudoscalar decay, \(V\to PP\). Additional decay types must be added to `Amplitude.cu` when needed.

Each four-momentum branch is read as

```cpp
SetPxPyPzE(p4[0], p4[1], p4[2], p4[3]);
```

and therefore follows the order

```text
[px, py, pz, E].
```

See `NLL_estimator.cu` for the ROOT branch handling.

### Likelihood model

```json
"set_fit_type": {
  "fit_type": "cFit",
  "PDF_bg": "PDF_bg",
  "bg_ratio": -1
}
```

Two likelihood constructions are available:

| Mode | Description |
|---|---|
| `nFit` | The background contribution is subtracted using the input background sample |
| `cFit` | Signal and background PDFs are combined in the event likelihood |

For `cFit`, the event probability is constructed schematically as

\[
P=(1-f_{\mathrm{bg}})P_{\mathrm{sig}}
  +f_{\mathrm{bg}}P_{\mathrm{bg}}.
\]

Additional requirements for `cFit`:

- `PDF_bg` gives the ROOT branch containing the background-PDF value for each data and MC event;
- `bg_ratio` gives the background fraction;
- if `bg_ratio` is `-1`, the code calculates it from the numbers of events in the input background and data samples.

### Resonances and decay chains

```json
"Res": {
  "A_Zc_3900": {
    "name": "A_Zc_3900",
    "res_par": [3.9015, 0.075, 0.135, 3.096, 2.0325, 1.865, 2.010],
    "spin": 1,
    "p-parity": 1,
    "dynamic_type": 3,
    "idx_isobar": 3
  }
}
```

Each member of `Res` defines one amplitude component or decay chain.

| Field | Meaning |
|---|---|
| `name` | Resonance/component name used in fit-parameter names |
| `res_par` | Lineshape parameters, such as mass, width, or Flatté-like couplings and thresholds |
| `spin` | Spin of the intermediate state |
| `p-parity` | Parity of the intermediate state |
| `dynamic_type` | Lineshape index implemented in `Dynamic.h` |
| `idx_isobar` | Index of the bachelor/spectator daughter |

The `idx_isobar` convention is:

| `idx_isobar` | Bachelor | Isobar pair |
|---:|---|---|
| 1 | daughter 1 | (23) |
| 2 | daughter 2 | (31) |
| 3 | daughter 3 | (12) |

For example, in

\[
\gamma^* \to \pi Z_c,
\]

if the pion is daughter 3, then `idx_isobar` is `3`.

The meaning and order of `res_par` depend on `dynamic_type`. Refer to `gpuversion/include/Dynamic.h` in the distributed project layout (or the corresponding `include/Dynamic.h` in the GPU source directory). New lineshapes can be added there.

### Input ROOT samples

```json
"data": {
  "filename": "./Dst0DmPip_data_withpdf.root",
  "chainname": "ana",
  "p4_final": ["p4_Dst0", "p4_Dm", "p4_Pip"]
}
```

| Field | Meaning |
|---|---|
| `filename` | Input ROOT file or `TChain`-compatible file pattern |
| `chainname` | Name of the input ROOT tree |
| `p4_final` | Four-momentum branches for daughters 1, 2, and 3 |

The order of `p4_final` must match the order in `Daus`. The same structure is used for `data`, `mc`, and `bg`.

### Minuit strategy

```json
"Fit_stragety": {
  "is_scan": false,
  "is_migrad": true,
  "is_hesse": true
}
```

This block selects the `TMinuit` operations to run:

- `is_scan`: run `SCAN`;
- `is_migrad`: run `MIGRAD`;
- `is_hesse`: run `HESSE`.

> The key is spelled `Fit_stragety` in the current implementation. The JSON file must use this exact spelling.

### Parameter files

```json
"para_list": {
  "random_number": -1,
  "listin": "./par_in/Dst0DmPip_par.json",
  "listout": "./par_out/Dst0DmPip_par.json"
}
```

| Field | Meaning |
|---|---|
| `random_number` | `-1` disables randomized initial values; any other value is used as the random seed |
| `listin` | Input parameter list; an empty or unavailable file leaves the code-generated defaults in place |
| `listout` | Output parameter list containing the fit result |

A free parameter is written as

```json
"parameter_name": [value, step, lower_bound, upper_bound]
```

Setting `step` to zero fixes the parameter.

### Fit fractions

```json
"Cal_FitFraction": "./out/output_FFs_Dst0DmPip.root"
```

When this optional field is present, the program calculates:

- diagonal fit fractions;
- interference terms between amplitude components;
- the corresponding statistical uncertainties.

The results are saved as matrices in the requested ROOT file and in an accompanying text file.

### Output ROOT files

```json
"save_root": {
  "data": "./out/output_data_Dst0DmPip.root",
  "mc": "./out/output_mc_Dst0DmPip.root",
  "bg": "./out/output_bg_Dst0DmPip.root"
}
```

The output trees contain calculated kinematic variables. The MC output additionally contains total and component event weights for plotting and projection studies.

## Running a fit

Run a single-channel fit with:

```bash
./Fit.exe Dst0DmPip.json
```

The executable loads the data, background, and normalization MC samples; constructs the configured amplitudes; runs the requested Minuit operations; and writes the fitted parameters and ROOT outputs.

## Simultaneous fits and shared parameters

Pass multiple configuration files to perform a simultaneous fit:

```bash
./Fit.exe Dst0DmPip.json DstmD0Pip.json
```

Parameter names start with the sample index: `s0_` for the first channel, `s1_` for the second channel, and so on.

A parameter can be constrained to a fixed multiple of another parameter:

```json
"s1_A_Zc_3900_LS1coff_L0S1.0_phi": [
  "s0_A_Zc_3900_LS1coff_L0S1.0_phi",
  1,
  0
]
```

This enforces

\[
\phi_{s1}=1\times\phi_{s0}.
\]

The array has the form

```text
["reference parameter", ratio, cached fixed value]
```

The third value is not used to establish the constraint when the parameter file is read. It is updated to the resolved fixed value in the output parameter file so that the result can subsequently be used by the CPU MC-sampling program.

## Sampling Monte Carlo events

After fitting, evaluate the fitted amplitude on phase-space MC with:

```bash
./SampleMC.exe 1 Dst0DmPip.json
```

The first argument controls component output:

| Value | Behavior |
|---:|---|
| `0` | Save only the total amplitude weight |
| `1` | Also save amplitude-component projections |

Before running the sampling job, update the configuration so that:

- `mc` points to the reconstructed phase-space MC sample;
- `Truth` points to the truth-level phase-space MC sample;
- `para_list.listin` points to the fitted parameter file;
- the output paths in `save_root` do not overwrite important inputs.

The CPU program then writes weighted reconstructed and truth-level samples that can be used for efficiency studies and fit projections.

## Citation

If this software is used in a physics analysis, please cite the Dalitz-plot decomposition formalism:

```text
M. Mikhasenko et al.,
Dalitz-plot decomposition for three-body decays,
Phys. Rev. D 101, 034033 (2020).
```

Paper record: <https://inspirehep.net/literature/1758460>
