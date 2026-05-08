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
