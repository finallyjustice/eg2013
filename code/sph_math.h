#ifndef __SPHMATH_H__
#define __SPHMATH_H__

#define EXP 2.71828f

float wavelet(float t)
{
	float p=2.0f;
	return (float)(2/(sqrt(3.0*p)*pow(PI, 0.25f)) * (1-pow(t, 2)/pow(p, 2)) * pow(EXP, (0.0f-pow(t, 2))/(2*pow(p, 2))));
}

float wavelet4(float t)
{
	float p=4.0f;
	float result= (float)(2/(sqrt(3.0*p)*pow(PI, 0.25f)) * (1-pow(t, 2)/pow(p, 2)) * pow(EXP, (0.0f-pow(t, 2))/(2*pow(p, 2))));

	if(result > 0.0f)
	{
		return result;
	}

	return 0.0f;
}

float polyk(float s)
{
	float result=pow((1.0f-s*s), 3);
	
	if(result > 0.0f)
	{
		return result;
	}

	return 0.0f;
}

#endif
