#include <iostream>
// ホスト側の関数
void myKernel(int* result1, int* result2) {
    // 処理して結果をresult1とresult2に格納
    *result1 = 10;
    *result2 = 20;
}

__device__ void calcDeltaIJ(geom_float pcx , geom_float pcy, geom_float pcz, 
                            flow_float dudx, flow_float dudy,flow_float dudz,
                            flow_float delu_max, flow_float delu_min ,
                            flow_float* delta, flow_float* deltap, flow_float* deltam
                            ) {
    *delta = res1;
    *deltap= res2;
    *deltam= res3;
}

int main() {
    int result1, result2;
    // ポインタを渡して関数を呼び出す
    myKernel(&result1, &result2);
    
    // 結果を出力
    std::cout << "Result 1: " << result1 << std::endl;
    std::cout << "Result 2: " << result2 << std::endl;

    return 0;
}
