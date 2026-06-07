# Dockerized dev environment for `solver_density_cuda`

## Overview

This directory is set up to use Docker as a disposable development environment for `solver_density_cuda`.

- The main target is CUDA development with `Dockerfile.cuda.dev`.
- Source code is bind-mounted into the container, so edits stay on the host.
- Build artifacts are generated in the mounted workspace, not inside an image layer.
- The intended workflow is `build image -> run container -> build inside container`.

In practice, the image is the reusable part. Containers are expected to be short-lived.

## Recommended host setup

If your goal is to get started quickly on Windows + WSL2, using Docker Desktop with WSL integration is acceptable.

If your goal is to stay close to AWS or a plain Linux server environment, prefer installing Docker Engine directly in Linux/WSL instead of depending on Docker Desktop long term.

Recommended approach:

1. Bring the environment up first with Docker Desktop if that is the fastest path.
2. Verify the CUDA container workflow for this repository.
3. Move to a direct Docker Engine setup later if you want parity with EC2 or Linux development hosts.

## Prerequisites

Before using the commands below, make sure these work on the host:

```bash
docker version
docker run --rm hello-world
```

For GPU usage, also confirm:

```bash
docker run --rm --gpus all nvidia/cuda:12.4.1-base-ubuntu22.04 nvidia-smi
```

If these commands fail, fix Docker or GPU runtime integration first.

## Build the image

From this directory:

```bash
cd /home/sano/work/forge/solver_density_cuda
docker build -f Dockerfile.cuda.dev -t forge-solver:cuda-dev .
```

This creates a local image named `forge-solver:cuda-dev`.

You only need to rebuild it when:

- `Dockerfile.cuda.dev` changes
- system packages in the image need refreshing
- you want to invalidate the previous build cache

## Start a development container

After the image exists, run:

```bash
docker run --rm -it --gpus all \
	-v "$(pwd)":/workspace \
	-w /workspace \
	forge-solver:cuda-dev bash
```

What this does:

- `--rm`: remove the container when you exit the shell
- `-it`: start an interactive shell
- `--gpus all`: expose available GPUs to the container
- `-v "$(pwd)":/workspace`: mount the current directory into the container
- `-w /workspace`: start in that mounted directory

This does not reopen an existing container. It starts a new container from the existing image each time.

That is the intended workflow here, because the source tree is mounted from the host.

## Build inside the container

Inside the container:

```bash
./tools/build.sh
```

This script configures CMake with the HDF5 include and library paths that are already exported by the image.

## Docker Compose option

There is also a minimal Compose definition in `compose.yml`.

```bash
docker compose run --rm dev-cuda bash
```

Then inside the container:

```bash
./tools/build.sh
```

Notes:

- `dev-cuda` matches the current repository state.
- `dev-cpu` is listed in `compose.yml`, but there is no `Dockerfile.dev` in the repository at the moment.
- Treat the CUDA path as the maintained path unless a CPU Dockerfile is added.

## Host shell helpers

The user shell configuration may define helper aliases such as:

- `forge-cuda`
- `forge-cuda2`
- `forge-gmsh`
- `forge-paraview`
- `forge-build`

These are convenience wrappers around `docker run`.

Be careful: those shortcuts are host-specific and may drift from the repository layout. For example, some helper functions assume directory names such as `cases`, while this workspace currently uses `case`.

Use the raw `docker build` and `docker run` commands in this document as the source of truth.

## GUI tools

The CUDA development image includes:

- `gmsh`
- `paraview`

If you want to launch them from the container, use the host-side wrapper scripts in `tools/`:

```bash
cd /home/sano/work/forge/solver_density_cuda
./tools/run_gmsh_gui.sh
./tools/run_paraview_gui.sh
```

You can also pass a repository-local file to open immediately:

```bash
./tools/run_gmsh_gui.sh /home/sano/work/forge/case/20.naca_ml/001.test/mesh/naca.geo
./tools/run_paraview_gui.sh /home/sano/work/forge/case/20.naca_ml/001.test/run_slau/res_0.xmf
```

These scripts are intended to be run on the host, not from inside the container.

On WSL2, they are set up for WSLg and do the following automatically:

- run the container as the host UID/GID instead of root
- create a private `XDG_RUNTIME_DIR` inside the container
- forward the required Wayland, Pulse, and X11 sockets

If toolbar icons are missing in ParaView, rebuild the image so the GUI dependencies from `Dockerfile.cuda.dev` are included:

```bash
cd /home/sano/work/forge/solver_density_cuda
docker build -f Dockerfile.cuda.dev -t forge-solver:cuda-dev .
```

On native Linux, you may still need to adjust environment variables or mounts depending on your desktop session.

## Profiling helpers

The repository now includes host-side wrapper scripts for NVIDIA profiling tools in `tools/`.

### Nsight Systems

Use this when you want a timeline view of:

- CUDA kernel launch order
- CPU wait time and synchronization
- memcpy activity
- overall GPU busy time

Command pattern:

```bash
cd /home/sano/work/forge/solver_density_cuda
./tools/run_nsys_profile.sh /home/sano/work/forge/case/08.bump/run_slau_mach1.65_4pctHt
```

Current repository state:

- the wrapper script is present
- the current `forge-solver:cuda-dev` image does not yet include `nsys`
- if you run the wrapper now, it will stop with a clear error instead of failing silently

That means the maintained path today is:

1. use the built-in `FORGE_PROFILE=1` runtime summary for coarse timings
2. use `run_ncu_profile.sh` for kernel-level GPU efficiency metrics
3. add Nsight Systems to the image later if timeline analysis becomes necessary

Outputs are written under the case directory:

```bash
case/08.bump/run_slau_mach1.65_4pctHt/profiles/nsys/
```

This wrapper enables `FORGE_PROFILE=1` by default, so the solver's runtime section summary is captured in the same run.

If the case directory already contains `res_*.h5` or `res_*.xmf` files owned by `root`, the wrapper will stop before launching and tell you to remove or `chown` them. This avoids waiting through startup only to fail at the first output write.

### Nsight Compute

Use this when you want kernel-level efficiency metrics such as:

- achieved occupancy
- SM utilization
- memory throughput
- scheduler stall reasons

Command pattern:

```bash
cd /home/sano/work/forge/solver_density_cuda
./tools/run_ncu_profile.sh /home/sano/work/forge/case/08.bump/run_slau_mach1.65_4pctHt
```

You can also narrow the profile to specific kernels:

```bash
./tools/run_ncu_profile.sh \
	/home/sano/work/forge/case/08.bump/run_slau_mach1.65_4pctHt \
	'SLAU_d|viscousFlux_d|calcGradient_2_d'
```

For a quick profiling run without touching the source case, use temporary overrides:

```bash
FORGE_NCU_NSTEP=1 \
FORGE_NCU_OUTSTEP_INTERVAL=1000000 \
NCU_LAUNCH_COUNT=1 \
./tools/run_ncu_profile.sh \
	/home/sano/work/forge/case/08.bump/run_slau_mach1.65_4pctHt \
	'SLAU_d'
```

Outputs are written under:

```bash
case/08.bump/run_slau_mach1.65_4pctHt/profiles/ncu/
```

Like the Nsight Systems wrapper, this script pre-checks whether existing result files in the case directory are writable by the current host user.

If `ncu` stops with `ERR_NVGPUCTRPERM`, the container setup is not the issue. That error means the host NVIDIA driver is blocking access to GPU performance counters. In that case, ask the host admin to allow profiling counters or configure the NVIDIA kernel module with `NVreg_RestrictProfilingToAdminUsers=0`.

For a practical workflow, use `run_nsys_profile.sh` first to identify the heaviest kernels, then use `run_ncu_profile.sh` with a narrow kernel regex to inspect GPU efficiency in detail.

Given the current image contents, replace that workflow with `FORGE_PROFILE=1` plus `run_ncu_profile.sh` unless you explicitly extend the image to add `nsys`.

## Windows Nsight Compute via WSL target

Docker is still the maintained execution environment for this repository, and it is fine for:

- building the solver
- running cases
- launching Gmsh and ParaView
- collecting the coarse runtime summary added through `FORGE_PROFILE=1`

The current issue is narrower than that: in this workspace's WSL2 + WDDM + Docker setup, `ncu` inside the container reaches the target process, but kernel profiling fails with `Failed to prepare kernel for profiling` / `Unknown Error on device 0`.

That does not prove Docker profiling is impossible in general. It only means this environment is not currently reliable for Nsight Compute kernel collection.

If you want Windows-side Nsight Compute to inspect GPU efficiency, the safer path is:

1. install Nsight Compute on Windows
2. build and run `forge` natively inside WSL instead of inside Docker
3. use Windows Nsight Compute to remote-launch the WSL target over SSH

Why use a native WSL target for Nsight Compute:

- it removes one moving part from the profiling path
- it matches NVIDIA's documented Windows host -> Linux target workflow
- it avoids depending on container-side tool behavior when the target kernel replay is already unstable

### WSL target prerequisites

Check the current WSL readiness with:

```bash
cd /home/sano/work/forge/solver_density_cuda
./tools/check_wsl_profile_target.sh
```

Typical packages needed in WSL:

```bash
sudo apt update
sudo apt install -y \
	build-essential cmake ninja-build gfortran pkg-config git \
	python3 python3-pip python3-h5py \
	libhdf5-dev libyaml-cpp-dev libmetis-dev \
	openssh-server
```

You also need a CUDA toolkit installed inside WSL so that `nvcc` exists. On WSL, use the WSL-Ubuntu CUDA toolkit path from NVIDIA and avoid installing Linux driver meta-packages such as `cuda`, `cuda-12-x`, or `cuda-drivers` inside WSL.

### Native WSL build

After the prerequisites are installed:

```bash
cd /home/sano/work/forge/solver_density_cuda
git submodule update --init --recursive
./tools/build_native_wsl.sh
```

This builds into:

```bash
solver_density_cuda/build/native-relwithdebinfo/
```

### Expose WSL as a profiling target

Windows Nsight Compute uses SSH for remote Linux targets, so enable `sshd` inside WSL:

```bash
sudo service ssh start
```

If you want it to come up automatically when WSL starts, configure your preferred WSL/systemd startup flow for `sshd`.

### Launch flow with Windows Nsight Compute

Once Windows Nsight Compute is installed:

1. add your WSL instance as a Linux remote target via SSH
2. set the application executable to the native WSL build output:

```bash
/home/sano/work/forge/solver_density_cuda/build/native-relwithdebinfo/forge
```

3. set the working directory to the case directory you want to run, for example:

```bash
/home/sano/work/forge/case/08.bump/run_slau_mach1.65_4pctHt
```

4. use a short run configuration first, then widen the metric set only after one `.ncu-rep` is collected successfully

Recommended interpretation:

- Docker remains the default path for normal development and case execution
- native WSL is the preferred path only for Nsight Compute collection if the container path stays unstable

## Implementation notes

- Source is bind-mounted from host to `/workspace`.
- HDF5 is provided by Ubuntu packages inside the image.
- `HDF5_INC=/usr/include/hdf5/serial`
- `HDF5_LIBDIR=/usr/lib/x86_64-linux-gnu/hdf5/serial`

These are passed to CMake because the current project still expects explicit HDF5 paths.

## Known gaps

- `README_docker.md` previously mentioned `Dockerfile.dev`, but that file is not present now.
- The current maintained path is the CUDA development image.
- If you want a CPU-only Docker workflow, add a matching `Dockerfile.dev` first.
