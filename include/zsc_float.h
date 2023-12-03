#ifndef FLOAT
#define FLOAT
const long double pi(3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679l), e(2.7182818284590452353602874713526624977572470936999595749669676277240766303535475945713821785251664274l);
float add(float x)noexcept//大于x最小的浮点数
{
	if ((*(unsigned*)(void*)&x & 8388607u) == 8388607u)
		++(*(unsigned*)(void*)&x >>= 23) <<= 23;
	else
		++*(unsigned*)(void*)&x;
	return x;
}
float& iadd(float& x)noexcept
{
	if ((*(unsigned*)(void*)&x & 8388607u) == 8388607u)
		++(*(unsigned*)(void*)&x >>= 23) <<= 23;
	else
		++*(unsigned*)(void*)&x;
	return x;
}
double add(double x)noexcept
{
	if ((*(unsigned long long*)(void*) & x & 4503599627370495ull) == 4503599627370495ull)
		++(*(unsigned long long*)(void*) & x >>= 23) <<= 23;
	else
		++* (unsigned long long*)(void*)& x;
	return x;
}
double& iadd(double& x)noexcept
{
	if ((*(unsigned long long*)(void*) & x & 4503599627370495ull) == 4503599627370495ull)
		++(*(unsigned long long*)(void*) & x >>= 23) <<= 23;
	else
		++* (unsigned long long*)(void*)& x;
	return x;
}
float fix(float x)noexcept//浮点意义下的取整,不会溢出
{
	unsigned t(*(unsigned*)(void*)&x);
	signed char e((t >> 23) - 127);
	if (e < 0)return 0.f;
	if (e > 22)return x;
	return *(float*)(void*)&(t &= ~((1 << 23 - e) - 1));
}
float& ifix(float& x)noexcept
{
	unsigned* t((unsigned*)(void*)&x);
	signed char e((*t >> 23) - 127);
	if (e < 0)return x = 0.f;
	if (e > 22)return x;
	*t &= ~((1 << 23 - e) - 1);
	return x;
}
double& ifix(double& x)noexcept
{
	unsigned long long* t((unsigned long long*)(void*) & x);
	short e((2047 & int(*t >> 52)) - 1023);
	if (e < 0)return x = 0.;
	if (e > 51)return x;
	*t &= ~((1 << 52 - e) - 1);
	return x;
}
double fix(double x)noexcept
{
	unsigned long long t(*(unsigned long long*)(void*) & x);
	short e((2047 & int(t >> 52)) - 1023);
	if (e < 0)return 0.;
	if (e > 51)return x;
	return *(double*)(void*)&(t &= ~((1 << 52 - e) - 1));
}
#endif