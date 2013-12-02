void main()
{
	vec3 posEye = vec3(gl_ModelViewMatrix * vec4(gl_Vertex.xyz, 1.0));
    float dist = length(posEye);
	if(gl_Vertex.w == 0)
	{
		gl_PointSize = 550.0/dist;
	}
	
	if(gl_Vertex.w == 1)
	{
		gl_PointSize = 350.0/dist;
	}
	gl_Vertex.w=1.0f;
    
	gl_TexCoord[0] = gl_MultiTexCoord0;
    gl_Position = gl_ModelViewProjectionMatrix * gl_Vertex;
    gl_FrontColor = gl_Color;
}