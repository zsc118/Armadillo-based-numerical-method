#ifndef ZSC_STR
#define ZSC_STR
namespace val
{
	long double ld(const char* s)
	{
		if (!*s)return 0.;
		long double x(0.), t(1.);
		bool f(true);
		if (*s == '-')
		{
			if (!*++s)return -1.;
			if (*s == '.')
				f = false;
			else
				if (f)
					x = '0' - *s;
				else
					x = ('0' - *s) * (t /= 10);
			while (*++s)
				if (*s == '.')
					f = false;
				else
					if (f)
						x = '0' - *s + x * 10;
					else
						x += ('0' - *s) * (t /= 10);
		}
		else
			do
				if (*s == '.')
					f = false;
				else
					if (f)
						x = *s - '0' + x * 10;
					else
						x += (*s - '0') * (t /= 10);
		while (*++s);
		return x;
	}
	double d(const char* s)
	{
		if (!*s)return 0.;
		double x(0.), t(1.);
		bool f(true);
		if (*s == '-')
		{
			if (!*++s)return -1.;
			if (*s == '.')
				f = false;
			else
				if (f)
					x = '0' - *s;
				else
					x = ('0' - *s) * (t /= 10);
			while (*++s)
				if (*s == '.')
					f = false;
				else
					if (f)
						x = '0' - *s + x * 10;
					else
						x += ('0' - *s) * (t /= 10);
		}
		else
			do
				if (*s == '.')
					f = false;
				else
					if (f)
						x = *s - '0' + x * 10;
					else
						x += (*s - '0') * (t /= 10);
			while (*++s);
		return x;
	}
	float f(const char* s)
	{
		if (!*s)return 0.f;
		float x(0.f), t(1.f);
		bool f(true);
		if (*s == '-')
		{
			if (!*++s)return -1.f;
			if (*s == '.')
				f = false;
			else
				if (f)
					x = '0' - *s;
				else
					x = ('0' - *s) * (t /= 10);
			while (*++s)
				if (*s == '.')
					f = false;
				else
					if (f)
						x = '0' - *s + x * 10;
					else
						x += ('0' - *s) * (t /= 10);
		}
		else
			do
				if (*s == '.')
					f = false;
				else
					if (f)
						x = *s - '0' + x * 10;
					else
						x += (*s - '0') * (t /= 10);
		while (*++s);
		return x;
	}
	int i(const char* s)
	{
		int x(0);
		if (*s == '-')
		{
			x = '0' - *++s;
			while (*++s)x = '0' - *s + x * 10;
		}
		else
			do x = *s - '0' + x * 10; while (*++s);
		return x;
	}
	long long ll(const char* s)
	{
		long long x(0);
		if (*s == '-')
		{
			x = '0' - *++s;
			while (*++s)x = '0' - *s + x * 10;
		}
		else
			do x = *s - '0' + x * 10; while (*++s);
		return x;
	}
	long l(const char* s)
	{
		long x(0);
		if (*s == '-')
		{
			x = '0' - *++s;
			while (*++s)x = '0' - *s + x * 10;
		}
		else
			do x = *s - '0' + x * 10; while (*++s);
		return x;
	}
	short s(const char* s)
	{
		short x(0);
		if (*s == '-')
		{
			x = '0' - *++s;
			while (*++s)x = '0' - *s + x * 10;
		}
		else
			do x = *s - '0' + x * 10; while (*++s);
		return x;
	}
	signed char sc(const char* s)
	{
		signed char x(0);
		if (*s == '-')
		{
			x = '0' - *++s;
			while (*++s)x = '0' - *s + x * 10;
		}
		else
			do x = *s - '0' + x * 10; while (*++s);
		return x;
	}
	unsigned u(const char* s)
	{
		unsigned x(0);
		do x = *s - '0' + x * 10; while (*++s);
		return x;
	}
	unsigned short us(const char* s)
	{
		unsigned short x(0);
		do x = *s - '0' + x * 10; while (*++s);
		return x;
	}
	unsigned char uc(const char* s)
	{
		unsigned char x(0);
		do x = *s - '0' + x * 10; while (*++s);
		return x;
	}
	unsigned long ul(const char* s)
	{
		unsigned long x(0);
		do x = *s - '0' + x * 10; while (*++s);
		return x;
	}
	unsigned long long ull(const char* s)
	{
		unsigned long long x(0);
		do x = *s - '0' + x * 10; while (*++s);
		return x;
	}
}
#endif