Shader "Hidden/RaymarchShader"
{
    Properties
    {
        _MainTex ("Texture", 2D) = "white" {}
    }
    SubShader
    {
        // No culling or depth
        Cull Off ZWrite Off ZTest Always

        Pass
        {
            CGPROGRAM
            #pragma vertex vert
            #pragma fragment frag
            #pragma target 3.0

            #include "UnityCG.cginc"
            #include "DistanceFunctions.cginc"
  
            sampler2D _MainTex;
            uniform sampler2D _CameraDepthTexture;
            uniform float4x4 _CamFrustum, _CamToWorld;
            uniform float _maxDistance, _box1round, _boxSphereSmooth, _sphereIntersectSmooth;
            uniform float4 _sphere1, _sphere2, _box1, _fractal;
            uniform float3 _modInterval;
            uniform float3 _LightDir;
            uniform fixed4 _mainColor, _secColor, _skyColor;
            uniform float _precision, _lightIntensity, _shadowIntensity, _aoIntensity, _forceFieldRad, _GlobalScale;
            uniform int _iterations, _functionNum;
            uniform int _useNormal;
            uniform int _useShadow;
            uniform float4 _forceFieldNormal, _forceFieldColor;
            uniform float3 _player;
            uniform float3 _modOffsetPos;
            uniform float3 _modOffsetRot;
            uniform float4x4 _iterationTransform;
            uniform float4x4 _globalTransform;
            uniform float _smoothRadius, _scaleFactor, _innerSphereRad, _power;

            struct appdata
            {
                float4 vertex : POSITION;
                float2 uv : TEXCOORD0;
                float3 ray : TEXCOORD1;
            };

            struct v2f
            {
                float2 uv : TEXCOORD0;
                float3 ray : TEXCOORD1;
                float4 vertex : SV_POSITION;
            };

            v2f vert (appdata v)
            {
                v2f o;
                half index = v.vertex.z;
                v.vertex.z = 0;
                o.vertex = UnityObjectToClipPos(v.vertex);
                o.uv = v.uv;

                o.ray = _CamFrustum[(int)index].xyz;

                o.ray /= abs(o.ray.z);

                o.ray = mul(_CamToWorld, o.ray);

                return o;
            }

            // the distancefunction for the forcefield around the player.
            float sdforceField(float3 p)
            {
                //simple sphere
                return sdSphere(p - _player , _forceFieldRad);
            }

            //Ray dicstance count
            float distanceField(float3 p){

                float2 dist;

                if(_functionNum == 1){
                    dist = sdMerger(p,_GlobalScale, _iterations,_modOffsetPos ,_iterationTransform, _globalTransform, _smoothRadius, _scaleFactor);
                }
                //merger cylinder
                else if(_functionNum == 2){
                    dist = sdMergerCyl(p,_GlobalScale, _iterations,_modOffsetPos ,_iterationTransform, _globalTransform, _smoothRadius, _scaleFactor);
                }
                //mergerPyr
                else if(_functionNum == 3){
                    //dist = sdtriangleCross(p, _GlobalScale);
                    dist = sdMergerPyr(p,_GlobalScale, _iterations,_modOffsetPos ,_iterationTransform, _globalTransform, _smoothRadius, _scaleFactor,
                                        float4x4(1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1));
                } 
                // neg sphere
                else if(_functionNum == 4){
                    dist = sdNegSphere(p,_GlobalScale, _iterations,_modOffsetPos ,_iterationTransform, _globalTransform, _innerSphereRad, _scaleFactor);
                }
                // Sierpinski
                else if(_functionNum == 5){
                    dist = sdSierpinski(p, _scaleFactor);
                }
                // Mandelbulb
                else if(_functionNum == 6){
                    dist = mandelbulb(p, _power,  _iterations, _smoothRadius);
                } 
                // Tower IFS
                else if(_functionNum == 7){
                    dist = towerIFS(p);
                }
                // Modern windows
                else if(_functionNum == 8){
                    dist = modernWindows(p);
                } 
                // Jungles 
                else if(_functionNum == 9){
                    dist = infinityJungles(p);
                }
                // Pseudo Kleinian
                else if(_functionNum == 10){
                    dist = pseudo_kleinian(p);
                }
                // Lampshade pattern
                else if(_functionNum == 11){
                    dist = terrain3SDF(p, float4(0,1,0,0));
                }
                return dist;
            }

            //Shadow func

            float shadowCalc( in float3 ro, in float3 rd, float mint, float maxt, float k )
            {
                float res = 1.0;
                float ph = 1e20;
                for( float t=mint; t<maxt; )
                {
                    float h = min(distanceField(ro + rd*t),sdforceField(ro + rd*t));
                    if( h<0.001 ) return 0.0;
                    float y = h*h/(2.0*ph);
                    float d = sqrt(h*h-y*y);
                    res = min( res, k*d/max(0.0,t-y) );
                    ph = h;
                    t += h;
                }
                return res;
            }

            // returns the normal in a single point of the fractal
            float3 getNormal(float3 p){
                float d = distanceField(p).x;
                const float2 e = float2(.01, 0);
                float3 n = d - float3(distanceField(p - e.xyy).x,distanceField(p - e.yxy).x,distanceField(p - e.yyx).x);
                return normalize(n);
            }

            // returns the normal of the forcefield
            float3 getNormalForceField(float3 p)
            {

              float d = sdforceField(p);
              const float2 e = float2(.01, 0);
              float3 n = d - float3(sdforceField(p - e.xyy),sdforceField(p - e.yxy),sdforceField(p - e.yyx));
              return normalize(n);

            }

            fixed4 raymarching(float3 ro, float3 rd, float depth){
                fixed4 result = fixed4(0, 0, 0, 0.5);
                const int max_iteration = 400;
                bool _forceFieldHit = false;
                float3 _forceFieldNormal;
                float t = 0; //distance travelled along the ray direction

                for (int i = 0; i < max_iteration; i++) {

                    //sends out ray from the camera
                    float3 p = ro + rd * t;

                    //return distance to forcefield
                    float _forceField = sdforceField(p);

                    //check for hit in distancefield
                    float2 d = distanceField(p);

                    if (d.x < _precision) { //We have hit smth
                        //shading

                        float3 colorDepth;
                        float light;
                        float shadow;

                        float3 color = float3(_mainColor.rgb*(_iterations-d.y)/_iterations + _secColor.rgb*(d.y)/_iterations);

                        if(_useNormal == 1){
                            float3 n = getNormal(p);
                             light = (dot(-_LightDir, n) * (1 - _lightIntensity) + _lightIntensity); //lambertian shading
                        }
                        else  light = 1;
                        
                        if(_useShadow == 1){
                             shadow = (shadowCalc(p, -_LightDir, 0.1, _maxDistance, 3) * (1 - _shadowIntensity) + _shadowIntensity); // soft shadows

                        }
                        else  shadow = 1;

                        float ao = (1 - 2 * i/float(max_iteration)) * (1 - _aoIntensity) + _aoIntensity; // ambient occlusion
                        float3 colorLight = float3 (color * light * shadow * ao); // multiplying all values between 0 and 1 to return final color
                        colorDepth = float3 (colorLight*(_maxDistance-t)/(_maxDistance) + _skyColor.rgb*(t)/(_maxDistance)); // Background color, multiplying with distance

                        if(_forceFieldHit == true)
                        {
                            colorDepth =dot(-rd, _forceFieldNormal)* colorDepth + (1-dot(-rd, _forceFieldNormal))*_forceFieldColor; // multiply by transparant forcefield
                            
                        }

                        result = fixed4(colorDepth ,1);
                        break;
                    }

                    // adds distance to the distance traveled and next point

                    if(_forceFieldHit == false)
                    {
                        // closer points get higher precicion to limit overstepping

                        if((d.x) < 10)
                        {
                            t+=  min(d.x * 0.75f, _forceField);
                        }
                        else if( abs(d.x) < 2)
                        {
                            t+= min(d.x * 0.5f, _forceField);
                        }
                        else t+= min(d.x, _forceField);
                        
                        
                    }
                    else t += d.x;
                    
                }
                return result;
            }


            fixed4 frag (v2f i) : SV_Target
            {
                float depth = LinearEyeDepth(tex2D(_CameraDepthTexture, i.uv).r);
                depth *= length(i.ray);
                fixed3 col = tex2D(_MainTex, i.uv);
                float3 rayDirection = normalize(i.ray.xyz);
                float3 rayOrigin = _WorldSpaceCameraPos;
                fixed4 result = raymarching(rayOrigin, rayDirection, depth);

                return fixed4(col * (1.0 - result.w) + result.xyz * result.w, 1.0);
            }
            ENDCG
        }
    }
}

