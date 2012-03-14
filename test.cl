

__kernel void convert_colors(__global uchar4* out, __global float *in, unsigned int size )
{
    
    unsigned int x = get_global_id(0);
    
    uchar r = clamp( in[x*3] * 255, (float)0, (float)255 );
    uchar g = clamp( in[x*3+1] * 255, (float)0, (float)255 );
    uchar b = clamp( in[x*3+2] * 255, (float)0, (float)255 );
    uchar4 outc = (uchar4)(r,g,b,0xff);
    
    out[x*4] = outc;
    out[x*4+1] = outc;
    out[x*4+2] = outc;
    out[x*4+3] = outc;
    
/*
    unsigned int r = in[x*3] * 255;
    unsigned int g = in[x*3+1] * 255;
    unsigned int b = in[x*3+2] * 255;
    
    
    unsigned int col = 0xff000000 | (b << 16) | (g<<8) | r;
    out[x*4] = col;
    out[x*4 + 1] = col;
    out[x*4 + 2] = col;
    out[x*4 + 3] = col;
  */  

/*
    unsigned int y = get_global_id(1);

    // calculate uv coordinates
    float u = x / (float) width;
    float v = y / (float) height;
    u = u*2.0f - 1.0f;
    v = v*2.0f - 1.0f;

    // calculate simple sine wave pattern
    float freq = 4.0f;
    float w = sin(u*freq + time) * cos(v*freq + time) * 0.5f;

    // write output vertex
    pos[y*width+x] = (float4)(u, w, v, 1.0f);
*/
}
