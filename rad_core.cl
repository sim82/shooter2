__kernel void rad_kernel( int num_planes, __global float *fact, __global int *target, __global int *offsets, __global int *num, 
                          __global float4 *col_diff_, __global float4 *emit, __global float4 *rad, __global float4 *rad2, __global int *perm )
{
    unsigned int x = perm[get_global_id(0)];
    //unsigned int x = get_global_id(0);
   // rad[x] = emit[x];//(float4)( 0.5, 0.5, 0.5, 0.5);
    
    float4 rad_acc = (0.0, 0.0, 0.0, 0.0);
    float4 col_diff = col_diff_[x];
    int offs = offsets[x];
    for( int i = 0; i < num[x]; ++i ) {
        int t = target[offs+i];
        float4 trad = rad[t];
        float ff = fact[offs+i];
        rad_acc += (col_diff * trad) * ff;
        
    }
    
    rad2[x] = rad_acc + emit[x];
    //rad2[x] = (1.0,1.0,1.0,1.0) * 0.5;
    //rad2[x] = emit[x];
}