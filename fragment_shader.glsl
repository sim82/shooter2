varying vec4 Color;
varying vec2 TexCoord;

uniform sampler2D Texture;
uniform float Amount;
uniform float Timer;


uniform vec2 samp_offs[15];
uniform float samp_weight[15];

uniform float abber_amount;

float rand( float seed ) {
    return 0.0;
    
}
vec4 blur2() {
    vec4 fr_left = texture2D( Texture, TexCoord + vec2( abber_amount, 0 ) );
    vec4 fr_right = texture2D( Texture, TexCoord + vec2( -abber_amount, 0) );
    vec4 fr_up = texture2D( Texture, TexCoord + vec2( 0, -abber_amount ) );
    vec4 fr_down = texture2D( Texture, TexCoord + vec2( 0, abber_amount ) );
    vec4 fr_mid = texture2D( Texture, TexCoord );
    
    return fr_mid * 0.2 + fr_left * 0.2 + fr_right * 0.2 + fr_up * 0.2 + fr_down * 0.2;
    
    
    
}

vec4 blur() {
    vec4 c = vec4(0.0);
    for (int i = 0 ; i < 15; i++)
    {
        vec4 fr_up = texture2D( Texture, TexCoord + samp_offs[i] ) * samp_weight[i];
        vec4 fr_down = texture2D( Texture, TexCoord - samp_offs[i] ) * samp_weight[i];
        vec4 fr_mid = texture2D(Texture, (TexCoord)) * samp_weight[i];
        c += fr_up * vec4( 1.0, 0.0, 0.0, 0.0);
        c += fr_down * vec4( 0.0, 1.0, 0.0, 0.0);
        c += fr_mid * vec4( 0.0, 0.0, 1.0, 0.0);
        
    }
    c.a = 1.0;
    return c;
}

vec4 scanlines()
{
    //fragment.b = mix(fragment.b, ((fragment.b * sin(TexCoord.y * 3000.0)) * 3.0), Amount);
    
    vec4 fr_up = texture2D( Texture, TexCoord + vec2( 0, -abber_amount ) );
    vec4 fr_down = texture2D( Texture, TexCoord + vec2( 0, abber_amount ) );
    vec4 fr_mid = texture2D( Texture, TexCoord );
    
       
    return vec4( fr_up.r, fr_down.g, fr_mid.b, 1.0 );
}

vec4 noise()
{
    float x = TexCoord.x * TexCoord.y * 123456.0 * Timer;
    x = (mod(x, 13.0) * mod(x, 123.0));
    float dx = mod(x, 0.0015) * Amount;
    float dy = mod(x, 0.0005) * Amount;
    vec4 c = texture2D(Texture, (TexCoord + vec2(dx, dy)));
    return c;
}

vec4 grey(in vec4 fragment)
{
	vec3 grey = vec3(dot(fragment.rgb, vec3(0.3, 0.59, 0.11)));
	fragment.rgb = mix(fragment.rgb, grey, Amount);
	fragment.a = 1.0;
	return fragment;
}

void main() 
{
//     gl_FragColor = scanlines();
    gl_FragColor = blur();
	//gl_FragColor = grey(scanlines(noise()));
//	gl_FragColor = grey(scanlines(texture2D(Texture, TexCoord)));
}
