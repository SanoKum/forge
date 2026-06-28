#include "periodicNode_d.cuh"
#include "cuda_forge/cudaWrapper.cuh"

// node-centered 周期境界 DOF 同一視 (median-dual M4, §4.5)。
//
// periodicRoot[ic] は ic の所属する周期 group の root (master) CV index。root/非周期では root==ic。
// 継ぎ目で割れた 1 つの CV の両側部分が同じ group に入る。継ぎ目には双対面を作らず (周期半割面は移流
// ループから除外)、両側の内部双対面が各 partial CV の res を組む。これらを group で足し合わせ、全員へ
// 同じ和を書き戻すことで「1 CV」の残差にする。合併体積 (buildPeriodicNodeGroups) と合わせ、両者が
// 同 res・同 vol で更新され同期する。free-stream 保存は closure と継ぎ目半割面の相殺で成立 (plans §4.5.3)。

// 第1パス: slave (root != self) の res を root へ atomicAdd で集約。
__global__ void periodicGatherToRoot_d
(
    geom_int nCells, geom_int* root,
    flow_float* res_ro, flow_float* res_roUx, flow_float* res_roUy, flow_float* res_roUz, flow_float* res_roe
)
{
    geom_int ic = blockDim.x*blockIdx.x + threadIdx.x;
    if (ic < nCells) {
        const geom_int r = root[ic];
        if (r != ic) {
            atomicAdd(&res_ro[r],   res_ro[ic]);
            atomicAdd(&res_roUx[r], res_roUx[ic]);
            atomicAdd(&res_roUy[r], res_roUy[ic]);
            atomicAdd(&res_roUz[r], res_roUz[ic]);
            atomicAdd(&res_roe[r],  res_roe[ic]);
        }
    }
}

// 第2パス: root の集約済み和を slave へ broadcast (group 全員が同じ和を持つ)。
__global__ void periodicBroadcastFromRoot_d
(
    geom_int nCells, geom_int* root,
    flow_float* res_ro, flow_float* res_roUx, flow_float* res_roUy, flow_float* res_roUz, flow_float* res_roe
)
{
    geom_int ic = blockDim.x*blockIdx.x + threadIdx.x;
    if (ic < nCells) {
        const geom_int r = root[ic];
        if (r != ic) {
            res_ro[ic]   = res_ro[r];
            res_roUx[ic] = res_roUx[r];
            res_roUy[ic] = res_roUy[r];
            res_roUz[ic] = res_roUz[r];
            res_roe[ic]  = res_roe[r];
        }
    }
}

void periodicNodeGather_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var)
{
    if (cfg.discretization != "node" || msh.periodicRoot_d == nullptr || msh.nPeriodicMembers == 0) return;

    periodicGatherToRoot_d<<<cuda_cfg.dimGrid_cell , cuda_cfg.dimBlock>>>(
        msh.nCells, msh.periodicRoot_d,
        var.c_d["res_ro"], var.c_d["res_roUx"], var.c_d["res_roUy"], var.c_d["res_roUz"], var.c_d["res_roe"]
    );
    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchkKernelSync();

    periodicBroadcastFromRoot_d<<<cuda_cfg.dimGrid_cell , cuda_cfg.dimBlock>>>(
        msh.nCells, msh.periodicRoot_d,
        var.c_d["res_ro"], var.c_d["res_roUx"], var.c_d["res_roUy"], var.c_d["res_roUz"], var.c_d["res_roe"]
    );
    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchkKernelSync();
}

// block-DPLUR sweep 後の dq ミラー (§4.5.7)。root の補正 dq を member へ broadcast し、
// 周期同一視ノードの dq drift (多値化→発散) を防ぐ。periodicBroadcastFromRoot_d を dq buffer に流用。
void periodicMirrorDq_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var)
{
    if (cfg.discretization != "node" || msh.periodicRoot_d == nullptr || msh.nPeriodicMembers == 0) return;

    if (cfg.blockDPLUR == 1) {
        periodicBroadcastFromRoot_d<<<cuda_cfg.dimGrid_cell , cuda_cfg.dimBlock>>>(
            msh.nCells, msh.periodicRoot_d,
            var.c_d["dq_block_old_0"], var.c_d["dq_block_old_1"], var.c_d["dq_block_old_2"],
            var.c_d["dq_block_old_3"], var.c_d["dq_block_old_4"]
        );
    } else {
        periodicBroadcastFromRoot_d<<<cuda_cfg.dimGrid_cell , cuda_cfg.dimBlock>>>(
            msh.nCells, msh.periodicRoot_d,
            var.c_d["dq_ro_old"], var.c_d["dq_roUx_old"], var.c_d["dq_roUy_old"],
            var.c_d["dq_roUz_old"], var.c_d["dq_roe_old"]
        );
    }
    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchkKernelSync();
}
