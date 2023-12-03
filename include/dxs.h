#ifndef ZSC_DXS
#define ZSC_DXS
#include <iostream>
#include "str.h"
#include "zsc_float.h"
//多项式类
//a[i]表示x^i前面的系数
template<typename T = float,typename N = unsigned char>
struct dxs
{
    T* a;
    N n;
    dxs()
    {
        a = new T[1];
        *a = 0;
        n = 0;
    }
    dxs(T x)
    {
        a = new T[1];
        *a = x;
        n = 0;
    }
    dxs(const T* b, const T* e)
    {
        n = e - b;
        a = new T[n];
        T* i = a + n;
        --n;
        do *--i = *--e; while (e > b);
    }
    dxs(T* p, N m)
    {
        a = p;
        n = m;
    }
    dxs(std::initializer_list<T> l)
    {
        const T* j(l.end());
        N k(n = --j - l.begin());
        T* i(a = new T[++k]);
        *i = *j;
        while (--k)*++i = *--j;
    }
    ~dxs()
    {
        delete[]a;
    }
    dxs(const dxs& x)
    {
        n = x.n;
        T* i((a = new T[n + 1]) + n);
        const T* j(x.a + n);
        *i = *j;
        do *--i = *--j; while (i > a);
    }
    operator bool()
    {
        return n || *a;
    }
    bool operator!()
    {
        return !n && !*a;
    }
    T operator()(T x)
    {
        T s(*a);
        if (!n)return s;
        T t(x), * p(a);
        N m(n);
        s += *++p * t;
        while (--m)s += *++p * (t *= x);
        return s;
    }
    dxs& operator=(const dxs& x)
    {
        delete[]a;
        N k;
        const T* j(x.a);
        T* i(a = new T[k = (n = x.n) + 1]);
        *i = *j;
        while (--k)*++i = *++j;
        return*this;
    }
    dxs& operator+=(const dxs& x)
    {
        if (n > x.n)
        {
            T* p(a + x.n);
            const T* q(x.a + x.n);
            *p += *q;
            while (--p >= a)*p += *--q;
            return*this;
        }
        if (x.n > n)
        {
            N k(n - x.n);
            T* p(new T[x.n + 1] + x.n);
            const T* i(a + n), * j(x.a + x.n);
        }
    }
};
template<class N>
dxs<float, N>& sim(dxs<float, N>& x, float e = 1e-4f)
{
    float* p(x.a), t(*p - fix(*p)), E(1.f - e);
    N n(x.n + 1);
    if (t<e && t>-e)ifix(*p);
    else if (t > E)++ifix(*p);
    else if (t < -E)--ifix(*p);
    while (--n)
    {
        ++p;
        t = *p - fix(*p);
        if (t<e && t>-e)ifix(*p);
        else if (t > E)++ifix(*p);
        else if (t < -E)--ifix(*p);
    }
    return x;
}
template<class N>
dxs<double, N>& sim(dxs<double, N>& x, double e = 1e-6)
{
    double* p(x.a), t(*p - fix(*p)), E(1. - e);
    N n(x.n + 1);
    if (t<e && t>-e)ifix(*p);
    else if (t > E)++ifix(*p);
    else if (t < -E)--ifix(*p);
    while (--n)
    {
        ++p;
        t = *p - fix(*p);
        if (t<e && t>-e)ifix(*p);
        else if (t > E)++ifix(*p);
        else if (t < -E)--ifix(*p);
    }
    return x;
}
template<class T, class N>
T dtgl(const dxs<T, N>& x, T a, T b)
{
    const T* p(x.a);
    T i(a), j(b), s(*p * b - *p * a);
    N n(1);
    while (n <= x.n)s += ((j *= b) - (i *= a)) * *++p / ++n;
    return s;
}
template<class T, class N>
dxs<T, N> tgl(const dxs<T, N>& x, T C)
{
    N n(x.n + 2);
    const T* p(x.a + x.n);
    T* a(new T[n--]);
    *a = C;
    *(a += n) = *p / n;
    while (--n)*--a = *--p / n;
    return dxs<T, N>(--a, x.n + 1);
}
template<class T, class N>
dxs<T, N>& itgl(dxs<T, N>& x, T C)
{
    const T* p(x.a + x.n);
    N n(++x.n);
    T* a(new T[n + 1]);
    *a = C;
    *(a += n) = *p / n;
    while (--n)*--a = *--p / n;
    delete[]x.a;
    x.a = --a;
    return x;
}
template<class T, class N>
dxs<T, N> der(const dxs<T, N>& x)
{
    if (!x.n)return dxs<T, N>();
    T* a(new T[x.n] + x.n);
    const T* p(x.a + x.n);
    N m(x.n);
    *--a = m * *p;
    while (--m)*--a = m * *--p;
    return dxs<T, N>(a, x.n - 1);
}
template<class T, class N>
dxs<T, N>& ider(dxs<T, N>& x)
{
    if (!x.n)
    {
        *x.a = 0;
        return x;
    }
    T* p(x.a + x.n);
    x.a = new T[x.n] + x.n;
    N n(x.n--);
    *--x.a = n * *p;
    while (--n)*--x.a = n * *--p;
    delete[]--p;
    return x;
}
template<class T, class N>
bool operator==(const dxs<T, N>& x, const dxs<T, N>& y)
{
    if (x.n != y.n)return false;
    const T* i(x.a + x.n), * j(y.a + y.n);
    do if (*--i != *--j)return false; while (i > x.a);
    return true;
}
template<class T, class N>
bool operator!=(const dxs<T, N>& x, const dxs<T, N>& y)
{
    if (x.n != y.n)return true;
    const T* i(x.a + x.n), * j(y.a + y.n);
    do if (*--i != *--j)return true; while (i > x.a);
    return false;
}
template<class T, class N>
dxs<T, N> operator%(const dxs<T, N>& x, const dxs<T, N>& y)
{
    if (!y.n)return 0;
    if (y.n > x.n)return y;
    if (y.n == x.n)
    {
        T* a(new T[x.n] + x.n), t(*x.a / *y.a);
        const T* i(x.a + x.n), * j(y.a + x.n);
        if ((*--a = *j - t * *i)<1e-4&&*a>-1e-4)
        {
            T v;
            do if (--j == y.a)return dxs<T, N>(); while ((v = *j - t * *--i) > 1e-4 || v < -1e-4);
            delete[]++(a -= x.n);
            N n(j - y.a);
            *(a = new T[n + 1] + n) = v;
            while (--j > y.a)*--a = *j - t * *--i;
            return dxs<T, N>(a, n);
        }
        while (--j > y.a)*--a = *j - t * *--i;
        return dxs<T, N>(a, x.n - 1);
    }
    N k(x.n - y.n), g(k);
    const T* i(x.a + x.n), * j(y.a + y.n);
    T* p(new T[k] + k), * q(new T[y.n] + y.n), t(*i / *j), * u(p), * v(q);
    do
    {
        *--p = *--i - t * *--j;
        if (j == y.a)
        {
            while (--g)*--p = *--i;
            do *--q = *--i; while (++g < y.n);
            goto F;
        }
    } while (--g);
    do
    {
        *--q = *--i - t * *--j;
        if (j == y.a)
        {
            while (++g < y.n)*--q = *--i;
            break;
        }
    } while (++g < y.n);
F:i = y.a + y.n;
    while (--k)
    {
        g = k;
        q = v;
        t = *--(p = --u) / *--(j = i);
        while (--g)
        {
            *--p -= t * *--j;
            if (j == y.a)goto G;
        }
        do *--q -= t * *--j; while (j > y.a);
    G:;
    }
    t = *--u / *i;
    delete[]u;
    g = y.n;
    do *--v -= t * *--i; while (i > y.a);
    return dxs<T, N>(v, y.n - 1);
}
template<class T, class N>
dxs<T, N> operator/(const dxs<T, N>& x, const dxs<T, N>& y)
{
    if (x.n < y.n)return dxs<T, N>();
    if (x.n == y.n)return x.a[x.n] / y.a[y.n];
    N n(x.n - y.n), k(n);
    T* a(new T[++k] + n), * q(new T[--k]), * p(q + n);
    const T* i(x.a + x.n), * j(y.a + y.n);
    *a = *i / *j;
    do *--p = *--i - *a * *--j; while (p > q);
    i = y.a + y.n;
    while (--k)
    {
        *--a = *(p = q + k) / *(j = i);
        do *--p -= *a * *--j; while (j > y.a && p > q);
    }
    *--a = *q / *i;
    delete[]q;
    return dxs<T, N>(a, n);
}
template<class T, class N>
dxs<T, N> gcd(const dxs<T, N>& x, const dxs<T, N>& y)
{

}
template<class T, class N>
dxs<T, N> operator*(const dxs<T, N>& x, const dxs<T, N>& y)
{
    N n(x.n + y.n), k(n), i;
    T* a(new T[n + 1] + n);
    *a = x.a[x.n] * y.a[y.n];
    while (--k)
    {
        i = k > x.n ? x.n : k;
        if (k > y.n + i)continue;
        *--a = x.a[i] * y.a[k - i];
        while (--i)
        {
            if (k > y.n + i)goto F;
            *a += x.a[i] * y.a[k - i];
        }
        if (k <= y.n + i)*a += *x.a * y.a[k];
        F:;
    }
    *--a = *x.a * *y.a;
    return dxs<T, N>(a, n);
}
template<class T, class N>
dxs<T, N> operator-(const dxs<T, N>& x)
{
    N n(x.n);
    T* a(new T[++n] + x.n);
    const T* i(x.a + x.n);
    *a = -*i;
    while (--n)*--a = -*--i;
    return dxs<T, N>(a, x.n);
}
template<class T, class N>
dxs<T, N> operator-(const dxs<T, N>& x, const dxs<T, N>& y)
{
    if (x.n > y.n)
    {
        N k(x.n - y.n);
        T* a(new T[x.n + 1] + x.n);
        const T* i(x.a + x.n), * j(y.a + y.n);
        *a = *i;
        while (--k)*--a = *--i;
        k = y.n + 1;
        *--a = *--i - *j;
        while (--k)*--a = *--i - *--j;
        return dxs<T, N>(a, x.n);
    }
    if (y.n > x.n)
    {
        N k(y.n - x.n);
        T* a(new T[y.n + 1] + y.n);
        const T* i(y.a + y.n), * j(x.a + x.n);
        *a = -*i;
        while (--k)*--a = -*--i;
        k = x.n + 1;
        *--a = *j - *--i;
        while (--k)*--a = *--j - *--i;
        return dxs<T, N>(a, y.n);
    }
    const T* i(x.a + x.n + 1), * j(y.a + y.n + 1);
    while (*--i - *--j < 1e-4 && *i - *j > -1e-4)if (i == x.a)return dxs<T, N>();
    N n(i - x.a), m(n);
    T* p(new T[++m] + n);
    *p = *i - *j;
    while (--m)*--p = *--i - *--j;
    return dxs<T, N>(p, n);
}
template<class T, class N>
dxs<T, N> operator+(const dxs<T, N>& x, const dxs<T, N>& y)
{
    if (x.n > y.n)
    {
        N k(x.n - y.n);
        T* a(new T[x.n + 1] + x.n);
        const T* i(x.a + x.n), * j(y.a + y.n);
        *a = *i;
        while (--k)*--a = *--i;
        k = y.n + 1;
        *--a = *--i + *j;
        while (--k)*--a = *--i + *--j;
        return dxs<T, N>(a, x.n);
    }
    if (y.n > x.n)
    {
        N k(y.n - x.n);
        T* a(new T[y.n + 1] + y.n);
        const T* i(y.a + y.n), * j(x.a + x.n);
        *a = *i;
        while (--k)*--a = *--i;
        k = x.n + 1;
        *--a = *j + *--i;
        while (--k)*--a = *--j + *--i;
        return dxs<T, N>(a, y.n);
    }
    const T* i(x.a + x.n + 1), * j(y.a + y.n + 1);
    while (*--i + *--j < 1e-4 && *i + *j > -1e-4)if (i == x.a)return dxs<T, N>();
    N n(i - x.a), m(n);
    T* p(new T[++m] + n);
    *p = *i + *j;
    while (--m)*--p = *--i + *--j;
    return dxs<T, N>(p, n);
}
template<class T,class N>
std::ostream& operator<<(std::ostream& c, const dxs<T, N>& x)
{
    const T* i(x.a + x.n);
    unsigned long long n(x.n);
    if (!x.n)
    {
        c << *i;
        return c;
    }
    if (*i != 1)
        if (*i == -1)
            c << '-';
        else
            c << *i;
    c << 'x';
    if (x.n == 1)
    {
        if (!*x.a)
            return c;
        else if (*x.a > 0)
            c << '+';
        c << *x.a;
        return c;
    }
    c << '^' << n;
    while(--n>1)
        if (*--i)
        {
            if (*i > 0)c << '+';
            if (*i != 1)
                if (*i == -1)
                    c << '-';
                else
                    c << *i;
            c << 'x' << '^' << n;
        }
    if (*--i)
    {
        if (*i > 0)c << '+';
        if (*i != 1)
            if (*i == -1)
                c << '-';
            else
                c << *i;
        c << 'x';
    }
    if (*x.a)
    {
        if (*x.a > 0)c << '+';
        c << *x.a;
    }
    return c;
}
template<typename N>
std::istream& operator>>(std::istream& c, dxs<float, N>& x)
{
    delete[]x.a;
    char s[41], * i(s);
    do c.get(*i); while (*i == ' ' || *i == '\t' || *i == '\n');
    if (*i == '+')c.get(*i);
    bool f(*i == 'x');
    while (*i != 'x')
    {
        c.get(*++i);
        if (*i == ' ' || *i == '\t' || *i == '\n')
        {
            *i = '\0';
            *(x.a = new float[1]) = val::f(s);
            x.n = 0;
            return c;
        }
    }
    *i = '\0';
    float* p, g(val::f(s));
    if (f)g = 1.f;
    c.get(*(i = s));
    if (*i == ' ' || *i == '\t' || *i == '\n')
    {
        *i = '\0';
        x.n = 1;
        *(x.a = new float[2]) = 0.f;
        x.a[1] = g;
        return c;
    }
    if (*i == '+')
    {
        x.n = 1;
        *(p = (x.a = new float[2]) + 1) = g;
        c.get(*i);
    }
    else if (*i == '-')
    {
        x.n = 1;
        *(p = (x.a = new float[2]) + 1) = g;
    }
    else
    {
        c.get(*i);
        do c.get(*++i); while (*i >= '0' && *i <= '9');
        if (*i == ' ' || *i == '\t' || *i == '\n')
        {
            *i = '\0';
            *(p = (x.a = new float[(x.n = val::ull(s)) + 1]) + x.n) = g;
            do *--p = 0.f; while (p > x.a);
            return c;
        }
        if (*i == '+')
        {
            *i = '\0';
            *(p = (x.a = new float[(x.n = val::ull(s)) + 1]) + x.n) = g;
            c.get(*(i = s));
        }
        else
        {
            *i = '\0';
            *(p = (x.a = new float[(x.n = val::ull(s)) + 1]) + x.n) = g;
            *(i = s) = '-';
        }
    }
    float* q(p);
    N w(x.n), v;
    do
    {
        c.get(*++i);
        if (*i == ' ' || *i == '\t' || *i == '\n')
        {
            *i = '\0';
            while (--p > x.a)*p = 0.f;
            *p = val::f(s);
            return c;
        }
        if (*i == 'x')
        {
            *i = '\0';
            g = val::f(s);
            F:c.get(*(i = s));
            if (*i == ' ' || *i == '\t' || *i == '\n')
            {
                while (--p > x.a)*p = 0.f;
                *p = g;
                *--p = 0.f;
                return c;
            }
            if (*i == '+')
            {
                *i = '\0';
                while (--p > x.a)*p = 0.f;
                *++p = g;
                c.get(*(i = s));
                do c.get(*++i); while (*i != ' ' && *i != '\t' && *i != '\n');
                *i = '\0';
                *x.a = val::f(s);
                return c;
            }
            if (*i == '-')
            {
                *i = '\0';
                while (--p > x.a)*p = 0.f;
                *++p = g;
                *(i = s) = '-';
                c.get(*++i);
                do c.get(*++i); while (*i != ' ' && *i != '\t' && *i != '\n');
                *i = '\0';
                *x.a = val::f(s);
                return c;
            }
            c.get(*i);
            do c.get(*++i); while (*i >= '0' && *i <= '9');
            if (*i == '+')
            {
                *i = '\0';
                q -= w - (v = val::ull(s));
                w = v;
                while (--p > q)*p = 0.f;
                *q = g;
                c.get(*(i = s));
            }
            else if(*i=='-')
            {
                *i = '\0';
                q -= w - (v = val::ull(s));
                w = v;
                while (--p > q)*p = 0.f;
                *q = g;
                *(i = s) = '-';
            }
            else
            {
                *i = '\0';
                q -= w - (v = val::ull(s));
                w = v;
                while (--p > q)*p = 0.f;
                *q = g;
                do *--q = 0.f; while (q > x.a);
                return c;
            }
            if (*i == 'x')
            {
                g = 1.f;
                goto F;
            }
        }
    } while (true);
}
template<typename N>
std::istream& operator>>(std::istream& c, dxs<double, N>& x)
{
    delete[]x.a;
    char s[41], * i(s);
    do c.get(*i); while (*i == ' ' || *i == '\t' || *i == '\n');
    if (*i == '+')c.get(*i);
    bool f(*i == 'x');
    while (*i != 'x')
    {
        c.get(*++i);
        if (*i == ' ' || *i == '\t' || *i == '\n')
        {
            *i = '\0';
            *(x.a = new double[1]) = val::d(s);
            x.n = 0;
            return c;
        }
    }
    *i = '\0';
    double* p, g(val::d(s));
    if (f)g = 1.;
    c.get(*(i = s));
    if (*i == ' ' || *i == '\t' || *i == '\n')
    {
        *i = '\0';
        x.n = 1;
        *(x.a = new double[2]) = 0.;
        x.a[1] = g;
        return c;
    }
    if (*i == '+')
    {
        x.n = 1;
        *(p = (x.a = new double[2]) + 1) = g;
        c.get(*i);
    }
    else if (*i == '-')
    {
        x.n = 1;
        *(p = (x.a = new double[2]) + 1) = g;
    }
    else
    {
        c.get(*i);
        do c.get(*++i); while (*i >= '0' && *i <= '9');
        if (*i == ' ' || *i == '\t' || *i == '\n')
        {
            *i = '\0';
            *(p = (x.a = new double[(x.n = val::ull(s)) + 1]) + x.n) = g;
            do *--p = 0.; while (p > x.a);
            return c;
        }
        if (*i == '+')
        {
            *i = '\0';
            *(p = (x.a = new double[(x.n = val::ull(s)) + 1]) + x.n) = g;
            c.get(*(i = s));
        }
        else
        {
            *i = '\0';
            *(p = (x.a = new double[(x.n = val::ull(s)) + 1]) + x.n) = g;
            *(i = s) = '-';
        }
    }
    double* q(p);
    N w(x.n), v;
    do
    {
        c.get(*++i);
        if (*i == ' ' || *i == '\t' || *i == '\n')
        {
            *i = '\0';
            while (--p > x.a)*p = 0.;
            *p = val::d(s);
            return c;
        }
        if (*i == 'x')
        {
            *i = '\0';
            g = val::d(s);
        F:c.get(*(i = s));
            if (*i == ' ' || *i == '\t' || *i == '\n')
            {
                while (--p > x.a)*p = 0.;
                *p = g;
                *--p = 0.;
                return c;
            }
            if (*i == '+')
            {
                *i = '\0';
                while (--p > x.a)*p = 0.;
                *++p = g;
                c.get(*(i = s));
                do c.get(*++i); while (*i != ' ' && *i != '\t' && *i != '\n');
                *i = '\0';
                *x.a = val::d(s);
                return c;
            }
            if (*i == '-')
            {
                *i = '\0';
                while (--p > x.a)*p = 0.;
                *++p = g;
                *(i = s) = '-';
                c.get(*++i);
                do c.get(*++i); while (*i != ' ' && *i != '\t' && *i != '\n');
                *i = '\0';
                *x.a = val::d(s);
                return c;
            }
            c.get(*i);
            do c.get(*++i); while (*i >= '0' && *i <= '9');
            if (*i == '+')
            {
                *i = '\0';
                q -= w - (v = val::ull(s));
                w = v;
                while (--p > q)*p = 0.;
                *q = g;
                c.get(*(i = s));
            }
            else if (*i == '-')
            {
                *i = '\0';
                q -= w - (v = val::ull(s));
                w = v;
                while (--p > q)*p = 0.;
                *q = g;
                *(i = s) = '-';
            }
            else
            {
                *i = '\0';
                q -= w - (v = val::ull(s));
                w = v;
                while (--p > q)*p = 0.;
                *q = g;
                do *--q = 0.; while (q > x.a);
                return c;
            }
            if (*i == 'x')
            {
                g = 1.;
                goto F;
            }
        }
    } while (true);
}
template<typename N>
std::istream& operator>>(std::istream& c, dxs<long double, N>& x)
{
    delete[]x.a;
    char s[41], * i(s);
    do c.get(*i); while (*i == ' ' || *i == '\t' || *i == '\n');
    if (*i == '+')c.get(*i);
    bool f(*i == 'x');
    while (*i != 'x')
    {
        c.get(*++i);
        if (*i == ' ' || *i == '\t' || *i == '\n')
        {
            *i = '\0';
            *(x.a = new long double[1]) = val::ld(s);
            x.n = 0;
            return c;
        }
    }
    *i = '\0';
    long double* p, g(val::ld(s));
    if (f)g = 1.;
    c.get(*(i = s));
    if (*i == ' ' || *i == '\t' || *i == '\n')
    {
        *i = '\0';
        x.n = 1;
        *(x.a = new long double[2]) = 0.;
        x.a[1] = g;
        return c;
    }
    if (*i == '+')
    {
        x.n = 1;
        *(p = (x.a = new long double[2]) + 1) = g;
        c.get(*i);
    }
    else if (*i == '-')
    {
        x.n = 1;
        *(p = (x.a = new long double[2]) + 1) = g;
    }
    else
    {
        c.get(*i);
        do c.get(*++i); while (*i >= '0' && *i <= '9');
        if (*i == ' ' || *i == '\t' || *i == '\n')
        {
            *i = '\0';
            *(p = (x.a = new long double[(x.n = val::ull(s)) + 1]) + x.n) = g;
            do *--p = 0.; while (p > x.a);
            return c;
        }
        if (*i == '+')
        {
            *i = '\0';
            *(p = (x.a = new long double[(x.n = val::ull(s)) + 1]) + x.n) = g;
            c.get(*(i = s));
        }
        else
        {
            *i = '\0';
            *(p = (x.a = new long double[(x.n = val::ull(s)) + 1]) + x.n) = g;
            *(i = s) = '-';
        }
    }
    long double* q(p);
    N w(x.n), v;
    do
    {
        c.get(*++i);
        if (*i == ' ' || *i == '\t' || *i == '\n')
        {
            *i = '\0';
            while (--p > x.a)*p = 0.;
            *p = val::ld(s);
            return c;
        }
        if (*i == 'x')
        {
            *i = '\0';
            g = val::ld(s);
        F:c.get(*(i = s));
            if (*i == ' ' || *i == '\t' || *i == '\n')
            {
                while (--p > x.a)*p = 0.;
                *p = g;
                *--p = 0.;
                return c;
            }
            if (*i == '+')
            {
                *i = '\0';
                while (--p > x.a)*p = 0.;
                *++p = g;
                c.get(*(i = s));
                do c.get(*++i); while (*i != ' ' && *i != '\t' && *i != '\n');
                *i = '\0';
                *x.a = val::ld(s);
                return c;
            }
            if (*i == '-')
            {
                *i = '\0';
                while (--p > x.a)*p = 0.;
                *++p = g;
                *(i = s) = '-';
                c.get(*++i);
                do c.get(*++i); while (*i != ' ' && *i != '\t' && *i != '\n');
                *i = '\0';
                *x.a = val::ld(s);
                return c;
            }
            c.get(*i);
            do c.get(*++i); while (*i >= '0' && *i <= '9');
            if (*i == '+')
            {
                *i = '\0';
                q -= w - (v = val::ull(s));
                w = v;
                while (--p > q)*p = 0.;
                *q = g;
                c.get(*(i = s));
            }
            else if (*i == '-')
            {
                *i = '\0';
                q -= w - (v = val::ull(s));
                w = v;
                while (--p > q)*p = 0.;
                *q = g;
                *(i = s) = '-';
            }
            else
            {
                *i = '\0';
                q -= w - (v = val::ull(s));
                w = v;
                while (--p > q)*p = 0.;
                *q = g;
                do *--q = 0.; while (q > x.a);
                return c;
            }
            if (*i == 'x')
            {
                g = 1.;
                goto F;
            }
        }
    } while (true);
}
template<typename N>
std::istream& operator>>(std::istream& c, dxs<int, N>& x)
{
    delete[]x.a;
    char s[41], * i(s);
    do c.get(*i); while (*i == ' ' || *i == '\t' || *i == '\n');
    if (*i == '+')c.get(*i);
    bool f(*i == 'x');
    while (*i != 'x')
    {
        c.get(*++i);
        if (*i == ' ' || *i == '\t' || *i == '\n')
        {
            *i = '\0';
            *(x.a = new int[1]) = val::i(s);
            x.n = 0;
            return c;
        }
    }
    *i = '\0';
    int* p, g(val::i(s));
    if (f)g = 1;
    c.get(*(i = s));
    if (*i == ' ' || *i == '\t' || *i == '\n')
    {
        *i = '\0';
        x.n = 1;
        *(x.a = new int[2]) = 0;
        x.a[1] = g;
        return c;
    }
    if (*i == '+')
    {
        x.n = 1;
        *(p = (x.a = new int[2]) + 1) = g;
        c.get(*i);
    }
    else if (*i == '-')
    {
        x.n = 1;
        *(p = (x.a = new int[2]) + 1) = g;
    }
    else
    {
        c.get(*i);
        do c.get(*++i); while (*i >= '0' && *i <= '9');
        if (*i == ' ' || *i == '\t' || *i == '\n')
        {
            *i = '\0';
            *(p = (x.a = new int[(x.n = val::ull(s)) + 1]) + x.n) = g;
            do *--p = 0; while (p > x.a);
            return c;
        }
        if (*i == '+')
        {
            *i = '\0';
            *(p = (x.a = new int[(x.n = val::ull(s)) + 1]) + x.n) = g;
            c.get(*(i = s));
        }
        else
        {
            *i = '\0';
            *(p = (x.a = new int[(x.n = val::ull(s)) + 1]) + x.n) = g;
            *(i = s) = '-';
        }
    }
    int* q(p);
    N w(x.n), v;
    do
    {
        c.get(*++i);
        if (*i == ' ' || *i == '\t' || *i == '\n')
        {
            *i = '\0';
            while (--p > x.a)*p = 0;
            *p = val::i(s);
            return c;
        }
        if (*i == 'x')
        {
            *i = '\0';
            g = val::i(s);
        F:c.get(*(i = s));
            if (*i == ' ' || *i == '\t' || *i == '\n')
            {
                while (--p > x.a)*p = 0;
                *p = g;
                *--p = 0;
                return c;
            }
            if (*i == '+')
            {
                *i = '\0';
                while (--p > x.a)*p = 0;
                *++p = g;
                c.get(*(i = s));
                do c.get(*++i); while (*i != ' ' && *i != '\t' && *i != '\n');
                *i = '\0';
                *x.a = val::i(s);
                return c;
            }
            if (*i == '-')
            {
                *i = '\0';
                while (--p > x.a)*p = 0;
                *++p = g;
                *(i = s) = '-';
                c.get(*++i);
                do c.get(*++i); while (*i != ' ' && *i != '\t' && *i != '\n');
                *i = '\0';
                *x.a = val::i(s);
                return c;
            }
            c.get(*i);
            do c.get(*++i); while (*i >= '0' && *i <= '9');
            if (*i == '+')
            {
                *i = '\0';
                q -= w - (v = val::ull(s));
                w = v;
                while (--p > q)*p = 0;
                *q = g;
                c.get(*(i = s));
            }
            else if (*i == '-')
            {
                *i = '\0';
                q -= w - (v = val::ull(s));
                w = v;
                while (--p > q)*p = 0;
                *q = g;
                *(i = s) = '-';
            }
            else
            {
                *i = '\0';
                q -= w - (v = val::ull(s));
                w = v;
                while (--p > q)*p = 0;
                *q = g;
                do *--q = 0; while (q > x.a);
                return c;
            }
            if (*i == 'x')
            {
                g = 1;
                goto F;
            }
        }
    } while (true);
}
template<typename N>
std::istream& operator>>(std::istream& c, dxs<long, N>& x)
{
    delete[]x.a;
    char s[41], * i(s);
    do c.get(*i); while (*i == ' ' || *i == '\t' || *i == '\n');
    if (*i == '+')c.get(*i);
    bool f(*i == 'x');
    while (*i != 'x')
    {
        c.get(*++i);
        if (*i == ' ' || *i == '\t' || *i == '\n')
        {
            *i = '\0';
            *(x.a = new long[1]) = val::l(s);
            x.n = 0;
            return c;
        }
    }
    *i = '\0';
    long* p, g(val::l(s));
    if (f)g = 1;
    c.get(*(i = s));
    if (*i == ' ' || *i == '\t' || *i == '\n')
    {
        *i = '\0';
        x.n = 1;
        *(x.a = new long[2]) = 0;
        x.a[1] = g;
        return c;
    }
    if (*i == '+')
    {
        x.n = 1;
        *(p = (x.a = new long[2]) + 1) = g;
        c.get(*i);
    }
    else if (*i == '-')
    {
        x.n = 1;
        *(p = (x.a = new long[2]) + 1) = g;
    }
    else
    {
        c.get(*i);
        do c.get(*++i); while (*i >= '0' && *i <= '9');
        if (*i == ' ' || *i == '\t' || *i == '\n')
        {
            *i = '\0';
            *(p = (x.a = new long[(x.n = val::ull(s)) + 1]) + x.n) = g;
            do *--p = 0; while (p > x.a);
            return c;
        }
        if (*i == '+')
        {
            *i = '\0';
            *(p = (x.a = new long[(x.n = val::ull(s)) + 1]) + x.n) = g;
            c.get(*(i = s));
        }
        else
        {
            *i = '\0';
            *(p = (x.a = new long[(x.n = val::ull(s)) + 1]) + x.n) = g;
            *(i = s) = '-';
        }
    }
    long* q(p);
    N w(x.n), v;
    do
    {
        c.get(*++i);
        if (*i == ' ' || *i == '\t' || *i == '\n')
        {
            *i = '\0';
            while (--p > x.a)*p = 0;
            *p = val::l(s);
            return c;
        }
        if (*i == 'x')
        {
            *i = '\0';
            g = val::l(s);
        F:c.get(*(i = s));
            if (*i == ' ' || *i == '\t' || *i == '\n')
            {
                while (--p > x.a)*p = 0;
                *p = g;
                *--p = 0;
                return c;
            }
            if (*i == '+')
            {
                *i = '\0';
                while (--p > x.a)*p = 0;
                *++p = g;
                c.get(*(i = s));
                do c.get(*++i); while (*i != ' ' && *i != '\t' && *i != '\n');
                *i = '\0';
                *x.a = val::l(s);
                return c;
            }
            if (*i == '-')
            {
                *i = '\0';
                while (--p > x.a)*p = 0;
                *++p = g;
                *(i = s) = '-';
                c.get(*++i);
                do c.get(*++i); while (*i != ' ' && *i != '\t' && *i != '\n');
                *i = '\0';
                *x.a = val::l(s);
                return c;
            }
            c.get(*i);
            do c.get(*++i); while (*i >= '0' && *i <= '9');
            if (*i == '+')
            {
                *i = '\0';
                q -= w - (v = val::ull(s));
                w = v;
                while (--p > q)*p = 0;
                *q = g;
                c.get(*(i = s));
            }
            else if (*i == '-')
            {
                *i = '\0';
                q -= w - (v = val::ull(s));
                w = v;
                while (--p > q)*p = 0;
                *q = g;
                *(i = s) = '-';
            }
            else
            {
                *i = '\0';
                q -= w - (v = val::ull(s));
                w = v;
                while (--p > q)*p = 0;
                *q = g;
                do *--q = 0; while (q > x.a);
                return c;
            }
            if (*i == 'x')
            {
                g = 1;
                goto F;
            }
        }
    } while (true);
}
template<typename N>
std::istream& operator>>(std::istream& c, dxs<short, N>& x)
{
    delete[]x.a;
    char s[41], * i(s);
    do c.get(*i); while (*i == ' ' || *i == '\t' || *i == '\n');
    if (*i == '+')c.get(*i);
    bool f(*i == 'x');
    while (*i != 'x')
    {
        c.get(*++i);
        if (*i == ' ' || *i == '\t' || *i == '\n')
        {
            *i = '\0';
            *(x.a = new short[1]) = val::s(s);
            x.n = 0;
            return c;
        }
    }
    *i = '\0';
    short* p, g(val::s(s));
    if (f)g = 1;
    c.get(*(i = s));
    if (*i == ' ' || *i == '\t' || *i == '\n')
    {
        *i = '\0';
        x.n = 1;
        *(x.a = new short[2]) = 0;
        x.a[1] = g;
        return c;
    }
    if (*i == '+')
    {
        x.n = 1;
        *(p = (x.a = new short[2]) + 1) = g;
        c.get(*i);
    }
    else if (*i == '-')
    {
        x.n = 1;
        *(p = (x.a = new short[2]) + 1) = g;
    }
    else
    {
        c.get(*i);
        do c.get(*++i); while (*i >= '0' && *i <= '9');
        if (*i == ' ' || *i == '\t' || *i == '\n')
        {
            *i = '\0';
            *(p = (x.a = new short[(x.n = val::ull(s)) + 1]) + x.n) = g;
            do *--p = 0; while (p > x.a);
            return c;
        }
        if (*i == '+')
        {
            *i = '\0';
            *(p = (x.a = new short[(x.n = val::ull(s)) + 1]) + x.n) = g;
            c.get(*(i = s));
        }
        else
        {
            *i = '\0';
            *(p = (x.a = new short[(x.n = val::ull(s)) + 1]) + x.n) = g;
            *(i = s) = '-';
        }
    }
    short* q(p);
    N w(x.n), v;
    do
    {
        c.get(*++i);
        if (*i == ' ' || *i == '\t' || *i == '\n')
        {
            *i = '\0';
            while (--p > x.a)*p = 0;
            *p = val::s(s);
            return c;
        }
        if (*i == 'x')
        {
            *i = '\0';
            g = val::s(s);
        F:c.get(*(i = s));
            if (*i == ' ' || *i == '\t' || *i == '\n')
            {
                while (--p > x.a)*p = 0;
                *p = g;
                *--p = 0;
                return c;
            }
            if (*i == '+')
            {
                *i = '\0';
                while (--p > x.a)*p = 0;
                *++p = g;
                c.get(*(i = s));
                do c.get(*++i); while (*i != ' ' && *i != '\t' && *i != '\n');
                *i = '\0';
                *x.a = val::s(s);
                return c;
            }
            if (*i == '-')
            {
                *i = '\0';
                while (--p > x.a)*p = 0;
                *++p = g;
                *(i = s) = '-';
                c.get(*++i);
                do c.get(*++i); while (*i != ' ' && *i != '\t' && *i != '\n');
                *i = '\0';
                *x.a = val::s(s);
                return c;
            }
            c.get(*i);
            do c.get(*++i); while (*i >= '0' && *i <= '9');
            if (*i == '+')
            {
                *i = '\0';
                q -= w - (v = val::ull(s));
                w = v;
                while (--p > q)*p = 0;
                *q = g;
                c.get(*(i = s));
            }
            else if (*i == '-')
            {
                *i = '\0';
                q -= w - (v = val::ull(s));
                w = v;
                while (--p > q)*p = 0;
                *q = g;
                *(i = s) = '-';
            }
            else
            {
                *i = '\0';
                q -= w - (v = val::ull(s));
                w = v;
                while (--p > q)*p = 0;
                *q = g;
                do *--q = 0; while (q > x.a);
                return c;
            }
            if (*i == 'x')
            {
                g = 1;
                goto F;
            }
        }
    } while (true);
}
template<typename N>
std::istream& operator>>(std::istream& c, dxs<signed char, N>& x)
{
    delete[]x.a;
    char s[41], * i(s);
    do c.get(*i); while (*i == ' ' || *i == '\t' || *i == '\n');
    if (*i == '+')c.get(*i);
    bool f(*i == 'x');
    while (*i != 'x')
    {
        c.get(*++i);
        if (*i == ' ' || *i == '\t' || *i == '\n')
        {
            *i = '\0';
            *(x.a = new signed char[1]) = val::sc(s);
            x.n = 0;
            return c;
        }
    }
    *i = '\0';
    signed char* p, g(val::sc(s));
    if (f)g = 1;
    c.get(*(i = s));
    if (*i == ' ' || *i == '\t' || *i == '\n')
    {
        *i = '\0';
        x.n = 1;
        *(x.a = new signed char[2]) = 0;
        x.a[1] = g;
        return c;
    }
    if (*i == '+')
    {
        x.n = 1;
        *(p = (x.a = new signed char[2]) + 1) = g;
        c.get(*i);
    }
    else if (*i == '-')
    {
        x.n = 1;
        *(p = (x.a = new signed char[2]) + 1) = g;
    }
    else
    {
        c.get(*i);
        do c.get(*++i); while (*i >= '0' && *i <= '9');
        if (*i == ' ' || *i == '\t' || *i == '\n')
        {
            *i = '\0';
            *(p = (x.a = new signed char[(x.n = val::ull(s)) + 1]) + x.n) = g;
            do *--p = 0; while (p > x.a);
            return c;
        }
        if (*i == '+')
        {
            *i = '\0';
            *(p = (x.a = new signed char[(x.n = val::ull(s)) + 1]) + x.n) = g;
            c.get(*(i = s));
        }
        else
        {
            *i = '\0';
            *(p = (x.a = new signed char[(x.n = val::ull(s)) + 1]) + x.n) = g;
            *(i = s) = '-';
        }
    }
    signed char* q(p);
    N w(x.n), v;
    do
    {
        c.get(*++i);
        if (*i == ' ' || *i == '\t' || *i == '\n')
        {
            *i = '\0';
            while (--p > x.a)*p = 0;
            *p = val::sc(s);
            return c;
        }
        if (*i == 'x')
        {
            *i = '\0';
            g = val::sc(s);
        F:c.get(*(i = s));
            if (*i == ' ' || *i == '\t' || *i == '\n')
            {
                while (--p > x.a)*p = 0;
                *p = g;
                *--p = 0;
                return c;
            }
            if (*i == '+')
            {
                *i = '\0';
                while (--p > x.a)*p = 0;
                *++p = g;
                c.get(*(i = s));
                do c.get(*++i); while (*i != ' ' && *i != '\t' && *i != '\n');
                *i = '\0';
                *x.a = val::sc(s);
                return c;
            }
            if (*i == '-')
            {
                *i = '\0';
                while (--p > x.a)*p = 0;
                *++p = g;
                *(i = s) = '-';
                c.get(*++i);
                do c.get(*++i); while (*i != ' ' && *i != '\t' && *i != '\n');
                *i = '\0';
                *x.a = val::sc(s);
                return c;
            }
            c.get(*i);
            do c.get(*++i); while (*i >= '0' && *i <= '9');
            if (*i == '+')
            {
                *i = '\0';
                q -= w - (v = val::ull(s));
                w = v;
                while (--p > q)*p = 0;
                *q = g;
                c.get(*(i = s));
            }
            else if (*i == '-')
            {
                *i = '\0';
                q -= w - (v = val::ull(s));
                w = v;
                while (--p > q)*p = 0;
                *q = g;
                *(i = s) = '-';
            }
            else
            {
                *i = '\0';
                q -= w - (v = val::ull(s));
                w = v;
                while (--p > q)*p = 0;
                *q = g;
                do *--q = 0; while (q > x.a);
                return c;
            }
            if (*i == 'x')
            {
                g = 1;
                goto F;
            }
        }
    } while (true);
}
template<typename N>
std::istream& operator>>(std::istream& c, dxs<unsigned char, N>& x)
{
    delete[]x.a;
    char s[41], * i(s);
    do c.get(*i); while (*i == ' ' || *i == '\t' || *i == '\n');
    if (*i == '+')c.get(*i);
    bool f(*i == 'x');
    while (*i != 'x')
    {
        c.get(*++i);
        if (*i == ' ' || *i == '\t' || *i == '\n')
        {
            *i = '\0';
            *(x.a = new unsigned char[1]) = val::uc(s);
            x.n = 0;
            return c;
        }
    }
    *i = '\0';
    unsigned char* p, g(val::uc(s));
    if (f)g = 1;
    c.get(*(i = s));
    if (*i == ' ' || *i == '\t' || *i == '\n')
    {
        *i = '\0';
        x.n = 1;
        *(x.a = new unsigned char[2]) = 0;
        x.a[1] = g;
        return c;
    }
    if (*i == '+')
    {
        x.n = 1;
        *(p = (x.a = new unsigned char[2]) + 1) = g;
        c.get(*i);
    }
    else if (*i == '-')
    {
        x.n = 1;
        *(p = (x.a = new unsigned char[2]) + 1) = g;
    }
    else
    {
        c.get(*i);
        do c.get(*++i); while (*i >= '0' && *i <= '9');
        if (*i == ' ' || *i == '\t' || *i == '\n')
        {
            *i = '\0';
            *(p = (x.a = new unsigned char[(x.n = val::ull(s)) + 1]) + x.n) = g;
            do *--p = 0; while (p > x.a);
            return c;
        }
        if (*i == '+')
        {
            *i = '\0';
            *(p = (x.a = new unsigned char[(x.n = val::ull(s)) + 1]) + x.n) = g;
            c.get(*(i = s));
        }
        else
        {
            *i = '\0';
            *(p = (x.a = new unsigned char[(x.n = val::ull(s)) + 1]) + x.n) = g;
            *(i = s) = '-';
        }
    }
    unsigned char* q(p);
    N w(x.n), v;
    do
    {
        c.get(*++i);
        if (*i == ' ' || *i == '\t' || *i == '\n')
        {
            *i = '\0';
            while (--p > x.a)*p = 0;
            *p = val::uc(s);
            return c;
        }
        if (*i == 'x')
        {
            *i = '\0';
            g = val::uc(s);
        F:c.get(*(i = s));
            if (*i == ' ' || *i == '\t' || *i == '\n')
            {
                while (--p > x.a)*p = 0;
                *p = g;
                *--p = 0;
                return c;
            }
            if (*i == '+')
            {
                *i = '\0';
                while (--p > x.a)*p = 0;
                *++p = g;
                c.get(*(i = s));
                do c.get(*++i); while (*i != ' ' && *i != '\t' && *i != '\n');
                *i = '\0';
                *x.a = val::uc(s);
                return c;
            }
            if (*i == '-')
            {
                *i = '\0';
                while (--p > x.a)*p = 0;
                *++p = g;
                *(i = s) = '-';
                c.get(*++i);
                do c.get(*++i); while (*i != ' ' && *i != '\t' && *i != '\n');
                *i = '\0';
                *x.a = val::uc(s);
                return c;
            }
            c.get(*i);
            do c.get(*++i); while (*i >= '0' && *i <= '9');
            if (*i == '+')
            {
                *i = '\0';
                q -= w - (v = val::ull(s));
                w = v;
                while (--p > q)*p = 0;
                *q = g;
                c.get(*(i = s));
            }
            else if (*i == '-')
            {
                *i = '\0';
                q -= w - (v = val::ull(s));
                w = v;
                while (--p > q)*p = 0;
                *q = g;
                *(i = s) = '-';
            }
            else
            {
                *i = '\0';
                q -= w - (v = val::ull(s));
                w = v;
                while (--p > q)*p = 0;
                *q = g;
                do *--q = 0; while (q > x.a);
                return c;
            }
            if (*i == 'x')
            {
                g = 1;
                goto F;
            }
        }
    } while (true);
}
template<typename N>
std::istream& operator>>(std::istream& c, dxs<long long, N>& x)
{
    delete[]x.a;
    char s[41], * i(s);
    do c.get(*i); while (*i == ' ' || *i == '\t' || *i == '\n');
    if (*i == '+')c.get(*i);
    bool f(*i == 'x');
    while (*i != 'x')
    {
        c.get(*++i);
        if (*i == ' ' || *i == '\t' || *i == '\n')
        {
            *i = '\0';
            *(x.a = new long long[1]) = val::ll(s);
            x.n = 0;
            return c;
        }
    }
    *i = '\0';
    long long* p, g(val::ll(s));
    if (f)g = 1;
    c.get(*(i = s));
    if (*i == ' ' || *i == '\t' || *i == '\n')
    {
        *i = '\0';
        x.n = 1;
        *(x.a = new long long[2]) = 0;
        x.a[1] = g;
        return c;
    }
    if (*i == '+')
    {
        x.n = 1;
        *(p = (x.a = new long long[2]) + 1) = g;
        c.get(*i);
    }
    else if (*i == '-')
    {
        x.n = 1;
        *(p = (x.a = new long long[2]) + 1) = g;
    }
    else
    {
        c.get(*i);
        do c.get(*++i); while (*i >= '0' && *i <= '9');
        if (*i == ' ' || *i == '\t' || *i == '\n')
        {
            *i = '\0';
            *(p = (x.a = new long long[(x.n = val::ull(s)) + 1]) + x.n) = g;
            do *--p = 0; while (p > x.a);
            return c;
        }
        if (*i == '+')
        {
            *i = '\0';
            *(p = (x.a = new long long[(x.n = val::ull(s)) + 1]) + x.n) = g;
            c.get(*(i = s));
        }
        else
        {
            *i = '\0';
            *(p = (x.a = new long long[(x.n = val::ull(s)) + 1]) + x.n) = g;
            *(i = s) = '-';
        }
    }
    long long* q(p);
    N w(x.n), v;
    do
    {
        c.get(*++i);
        if (*i == ' ' || *i == '\t' || *i == '\n')
        {
            *i = '\0';
            while (--p > x.a)*p = 0;
            *p = val::ll(s);
            return c;
        }
        if (*i == 'x')
        {
            *i = '\0';
            g = val::ll(s);
        F:c.get(*(i = s));
            if (*i == ' ' || *i == '\t' || *i == '\n')
            {
                while (--p > x.a)*p = 0;
                *p = g;
                *--p = 0;
                return c;
            }
            if (*i == '+')
            {
                *i = '\0';
                while (--p > x.a)*p = 0;
                *++p = g;
                c.get(*(i = s));
                do c.get(*++i); while (*i != ' ' && *i != '\t' && *i != '\n');
                *i = '\0';
                *x.a = val::ll(s);
                return c;
            }
            if (*i == '-')
            {
                *i = '\0';
                while (--p > x.a)*p = 0;
                *++p = g;
                *(i = s) = '-';
                c.get(*++i);
                do c.get(*++i); while (*i != ' ' && *i != '\t' && *i != '\n');
                *i = '\0';
                *x.a = val::ll(s);
                return c;
            }
            c.get(*i);
            do c.get(*++i); while (*i >= '0' && *i <= '9');
            if (*i == '+')
            {
                *i = '\0';
                q -= w - (v = val::ull(s));
                w = v;
                while (--p > q)*p = 0;
                *q = g;
                c.get(*(i = s));
            }
            else if (*i == '-')
            {
                *i = '\0';
                q -= w - (v = val::ull(s));
                w = v;
                while (--p > q)*p = 0;
                *q = g;
                *(i = s) = '-';
            }
            else
            {
                *i = '\0';
                q -= w - (v = val::ull(s));
                w = v;
                while (--p > q)*p = 0;
                *q = g;
                do *--q = 0; while (q > x.a);
                return c;
            }
            if (*i == 'x')
            {
                g = 1;
                goto F;
            }
        }
    } while (true);
}
template<typename N>
std::istream& operator>>(std::istream& c, dxs<unsigned short, N>& x)
{
    delete[]x.a;
    char s[41], * i(s);
    do c.get(*i); while (*i == ' ' || *i == '\t' || *i == '\n');
    if (*i == '+')c.get(*i);
    bool f(*i == 'x');
    while (*i != 'x')
    {
        c.get(*++i);
        if (*i == ' ' || *i == '\t' || *i == '\n')
        {
            *i = '\0';
            *(x.a = new unsigned short[1]) = val::us(s);
            x.n = 0;
            return c;
        }
    }
    *i = '\0';
    unsigned short* p, g(val::us(s));
    if (f)g = 1;
    c.get(*(i = s));
    if (*i == ' ' || *i == '\t' || *i == '\n')
    {
        *i = '\0';
        x.n = 1;
        *(x.a = new unsigned short[2]) = 0;
        x.a[1] = g;
        return c;
    }
    if (*i == '+')
    {
        x.n = 1;
        *(p = (x.a = new unsigned short[2]) + 1) = g;
        c.get(*i);
    }
    else if (*i == '-')
    {
        x.n = 1;
        *(p = (x.a = new unsigned short[2]) + 1) = g;
    }
    else
    {
        c.get(*i);
        do c.get(*++i); while (*i >= '0' && *i <= '9');
        if (*i == ' ' || *i == '\t' || *i == '\n')
        {
            *i = '\0';
            *(p = (x.a = new unsigned short[(x.n = val::ull(s)) + 1]) + x.n) = g;
            do *--p = 0; while (p > x.a);
            return c;
        }
        if (*i == '+')
        {
            *i = '\0';
            *(p = (x.a = new unsigned short[(x.n = val::ull(s)) + 1]) + x.n) = g;
            c.get(*(i = s));
        }
        else
        {
            *i = '\0';
            *(p = (x.a = new unsigned short[(x.n = val::ull(s)) + 1]) + x.n) = g;
            *(i = s) = '-';
        }
    }
    unsigned short* q(p);
    N w(x.n), v;
    do
    {
        c.get(*++i);
        if (*i == ' ' || *i == '\t' || *i == '\n')
        {
            *i = '\0';
            while (--p > x.a)*p = 0;
            *p = val::us(s);
            return c;
        }
        if (*i == 'x')
        {
            *i = '\0';
            g = val::us(s);
        F:c.get(*(i = s));
            if (*i == ' ' || *i == '\t' || *i == '\n')
            {
                while (--p > x.a)*p = 0;
                *p = g;
                *--p = 0;
                return c;
            }
            if (*i == '+')
            {
                *i = '\0';
                while (--p > x.a)*p = 0;
                *++p = g;
                c.get(*(i = s));
                do c.get(*++i); while (*i != ' ' && *i != '\t' && *i != '\n');
                *i = '\0';
                *x.a = val::us(s);
                return c;
            }
            if (*i == '-')
            {
                *i = '\0';
                while (--p > x.a)*p = 0;
                *++p = g;
                *(i = s) = '-';
                c.get(*++i);
                do c.get(*++i); while (*i != ' ' && *i != '\t' && *i != '\n');
                *i = '\0';
                *x.a = val::us(s);
                return c;
            }
            c.get(*i);
            do c.get(*++i); while (*i >= '0' && *i <= '9');
            if (*i == '+')
            {
                *i = '\0';
                q -= w - (v = val::ull(s));
                w = v;
                while (--p > q)*p = 0;
                *q = g;
                c.get(*(i = s));
            }
            else if (*i == '-')
            {
                *i = '\0';
                q -= w - (v = val::ull(s));
                w = v;
                while (--p > q)*p = 0;
                *q = g;
                *(i = s) = '-';
            }
            else
            {
                *i = '\0';
                q -= w - (v = val::ull(s));
                w = v;
                while (--p > q)*p = 0;
                *q = g;
                do *--q = 0; while (q > x.a);
                return c;
            }
            if (*i == 'x')
            {
                g = 1;
                goto F;
            }
        }
    } while (true);
}
template<typename N>
std::istream& operator>>(std::istream& c, dxs<unsigned, N>& x)
{
    delete[]x.a;
    char s[41], * i(s);
    do c.get(*i); while (*i == ' ' || *i == '\t' || *i == '\n');
    if (*i == '+')c.get(*i);
    bool f(*i == 'x');
    while (*i != 'x')
    {
        c.get(*++i);
        if (*i == ' ' || *i == '\t' || *i == '\n')
        {
            *i = '\0';
            *(x.a = new unsigned[1]) = val::u(s);
            x.n = 0;
            return c;
        }
    }
    *i = '\0';
    unsigned* p, g(val::u(s));
    if (f)g = 1;
    c.get(*(i = s));
    if (*i == ' ' || *i == '\t' || *i == '\n')
    {
        *i = '\0';
        x.n = 1;
        *(x.a = new unsigned[2]) = 0;
        x.a[1] = g;
        return c;
    }
    if (*i == '+')
    {
        x.n = 1;
        *(p = (x.a = new unsigned[2]) + 1) = g;
        c.get(*i);
    }
    else if (*i == '-')
    {
        x.n = 1;
        *(p = (x.a = new unsigned[2]) + 1) = g;
    }
    else
    {
        c.get(*i);
        do c.get(*++i); while (*i >= '0' && *i <= '9');
        if (*i == ' ' || *i == '\t' || *i == '\n')
        {
            *i = '\0';
            *(p = (x.a = new unsigned[(x.n = val::ull(s)) + 1]) + x.n) = g;
            do *--p = 0; while (p > x.a);
            return c;
        }
        if (*i == '+')
        {
            *i = '\0';
            *(p = (x.a = new unsigned[(x.n = val::ull(s)) + 1]) + x.n) = g;
            c.get(*(i = s));
        }
        else
        {
            *i = '\0';
            *(p = (x.a = new unsigned[(x.n = val::ull(s)) + 1]) + x.n) = g;
            *(i = s) = '-';
        }
    }
    unsigned* q(p);
    N w(x.n), v;
    do
    {
        c.get(*++i);
        if (*i == ' ' || *i == '\t' || *i == '\n')
        {
            *i = '\0';
            while (--p > x.a)*p = 0;
            *p = val::u(s);
            return c;
        }
        if (*i == 'x')
        {
            *i = '\0';
            g = val::u(s);
        F:c.get(*(i = s));
            if (*i == ' ' || *i == '\t' || *i == '\n')
            {
                while (--p > x.a)*p = 0;
                *p = g;
                *--p = 0;
                return c;
            }
            if (*i == '+')
            {
                *i = '\0';
                while (--p > x.a)*p = 0;
                *++p = g;
                c.get(*(i = s));
                do c.get(*++i); while (*i != ' ' && *i != '\t' && *i != '\n');
                *i = '\0';
                *x.a = val::u(s);
                return c;
            }
            if (*i == '-')
            {
                *i = '\0';
                while (--p > x.a)*p = 0;
                *++p = g;
                *(i = s) = '-';
                c.get(*++i);
                do c.get(*++i); while (*i != ' ' && *i != '\t' && *i != '\n');
                *i = '\0';
                *x.a = val::u(s);
                return c;
            }
            c.get(*i);
            do c.get(*++i); while (*i >= '0' && *i <= '9');
            if (*i == '+')
            {
                *i = '\0';
                q -= w - (v = val::ull(s));
                w = v;
                while (--p > q)*p = 0;
                *q = g;
                c.get(*(i = s));
            }
            else if (*i == '-')
            {
                *i = '\0';
                q -= w - (v = val::ull(s));
                w = v;
                while (--p > q)*p = 0;
                *q = g;
                *(i = s) = '-';
            }
            else
            {
                *i = '\0';
                q -= w - (v = val::ull(s));
                w = v;
                while (--p > q)*p = 0;
                *q = g;
                do *--q = 0; while (q > x.a);
                return c;
            }
            if (*i == 'x')
            {
                g = 1;
                goto F;
            }
        }
    } while (true);
}
template<typename N>
std::istream& operator>>(std::istream& c, dxs<unsigned long, N>& x)
{
    delete[]x.a;
    char s[41], * i(s);
    do c.get(*i); while (*i == ' ' || *i == '\t' || *i == '\n');
    if (*i == '+')c.get(*i);
    bool f(*i == 'x');
    while (*i != 'x')
    {
        c.get(*++i);
        if (*i == ' ' || *i == '\t' || *i == '\n')
        {
            *i = '\0';
            *(x.a = new unsigned long[1]) = val::ul(s);
            x.n = 0;
            return c;
        }
    }
    *i = '\0';
    unsigned long* p, g(val::ul(s));
    if (f)g = 1;
    c.get(*(i = s));
    if (*i == ' ' || *i == '\t' || *i == '\n')
    {
        *i = '\0';
        x.n = 1;
        *(x.a = new unsigned long[2]) = 0;
        x.a[1] = g;
        return c;
    }
    if (*i == '+')
    {
        x.n = 1;
        *(p = (x.a = new unsigned long[2]) + 1) = g;
        c.get(*i);
    }
    else if (*i == '-')
    {
        x.n = 1;
        *(p = (x.a = new unsigned long[2]) + 1) = g;
    }
    else
    {
        c.get(*i);
        do c.get(*++i); while (*i >= '0' && *i <= '9');
        if (*i == ' ' || *i == '\t' || *i == '\n')
        {
            *i = '\0';
            *(p = (x.a = new unsigned long[(x.n = val::ull(s)) + 1]) + x.n) = g;
            do *--p = 0; while (p > x.a);
            return c;
        }
        if (*i == '+')
        {
            *i = '\0';
            *(p = (x.a = new unsigned long[(x.n = val::ull(s)) + 1]) + x.n) = g;
            c.get(*(i = s));
        }
        else
        {
            *i = '\0';
            *(p = (x.a = new unsigned long[(x.n = val::ull(s)) + 1]) + x.n) = g;
            *(i = s) = '-';
        }
    }
    unsigned long* q(p);
    N w(x.n), v;
    do
    {
        c.get(*++i);
        if (*i == ' ' || *i == '\t' || *i == '\n')
        {
            *i = '\0';
            while (--p > x.a)*p = 0;
            *p = val::ul(s);
            return c;
        }
        if (*i == 'x')
        {
            *i = '\0';
            g = val::ul(s);
        F:c.get(*(i = s));
            if (*i == ' ' || *i == '\t' || *i == '\n')
            {
                while (--p > x.a)*p = 0;
                *p = g;
                *--p = 0;
                return c;
            }
            if (*i == '+')
            {
                *i = '\0';
                while (--p > x.a)*p = 0;
                *++p = g;
                c.get(*(i = s));
                do c.get(*++i); while (*i != ' ' && *i != '\t' && *i != '\n');
                *i = '\0';
                *x.a = val::ul(s);
                return c;
            }
            if (*i == '-')
            {
                *i = '\0';
                while (--p > x.a)*p = 0;
                *++p = g;
                *(i = s) = '-';
                c.get(*++i);
                do c.get(*++i); while (*i != ' ' && *i != '\t' && *i != '\n');
                *i = '\0';
                *x.a = val::ul(s);
                return c;
            }
            c.get(*i);
            do c.get(*++i); while (*i >= '0' && *i <= '9');
            if (*i == '+')
            {
                *i = '\0';
                q -= w - (v = val::ull(s));
                w = v;
                while (--p > q)*p = 0;
                *q = g;
                c.get(*(i = s));
            }
            else if (*i == '-')
            {
                *i = '\0';
                q -= w - (v = val::ull(s));
                w = v;
                while (--p > q)*p = 0;
                *q = g;
                *(i = s) = '-';
            }
            else
            {
                *i = '\0';
                q -= w - (v = val::ull(s));
                w = v;
                while (--p > q)*p = 0;
                *q = g;
                do *--q = 0; while (q > x.a);
                return c;
            }
            if (*i == 'x')
            {
                g = 1;
                goto F;
            }
        }
    } while (true);
}
template<typename N>
std::istream& operator>>(std::istream& c, dxs<unsigned long long, N>& x)
{
    delete[]x.a;
    char s[41], * i(s);
    do c.get(*i); while (*i == ' ' || *i == '\t' || *i == '\n');
    if (*i == '+')c.get(*i);
    bool f(*i == 'x');
    while (*i != 'x')
    {
        c.get(*++i);
        if (*i == ' ' || *i == '\t' || *i == '\n')
        {
            *i = '\0';
            *(x.a = new unsigned long long[1]) = val::ull(s);
            x.n = 0;
            return c;
        }
    }
    *i = '\0';
    unsigned long long* p, g(val::ull(s));
    if (f)g = 1;
    c.get(*(i = s));
    if (*i == ' ' || *i == '\t' || *i == '\n')
    {
        *i = '\0';
        x.n = 1;
        *(x.a = new unsigned long long[2]) = 0;
        x.a[1] = g;
        return c;
    }
    if (*i == '+')
    {
        x.n = 1;
        *(p = (x.a = new unsigned long long[2]) + 1) = g;
        c.get(*i);
    }
    else if (*i == '-')
    {
        x.n = 1;
        *(p = (x.a = new unsigned long long[2]) + 1) = g;
    }
    else
    {
        c.get(*i);
        do c.get(*++i); while (*i >= '0' && *i <= '9');
        if (*i == ' ' || *i == '\t' || *i == '\n')
        {
            *i = '\0';
            *(p = (x.a = new unsigned long long[(x.n = val::ull(s)) + 1]) + x.n) = g;
            do *--p = 0; while (p > x.a);
            return c;
        }
        if (*i == '+')
        {
            *i = '\0';
            *(p = (x.a = new unsigned long long[(x.n = val::ull(s)) + 1]) + x.n) = g;
            c.get(*(i = s));
        }
        else
        {
            *i = '\0';
            *(p = (x.a = new unsigned long long[(x.n = val::ull(s)) + 1]) + x.n) = g;
            *(i = s) = '-';
        }
    }
    unsigned long long* q(p);
    N w(x.n), v;
    do
    {
        c.get(*++i);
        if (*i == ' ' || *i == '\t' || *i == '\n')
        {
            *i = '\0';
            while (--p > x.a)*p = 0;
            *p = val::ull(s);
            return c;
        }
        if (*i == 'x')
        {
            *i = '\0';
            g = val::ull(s);
        F:c.get(*(i = s));
            if (*i == ' ' || *i == '\t' || *i == '\n')
            {
                while (--p > x.a)*p = 0;
                *p = g;
                *--p = 0;
                return c;
            }
            if (*i == '+')
            {
                *i = '\0';
                while (--p > x.a)*p = 0;
                *++p = g;
                c.get(*(i = s));
                do c.get(*++i); while (*i != ' ' && *i != '\t' && *i != '\n');
                *i = '\0';
                *x.a = val::ull(s);
                return c;
            }
            if (*i == '-')
            {
                *i = '\0';
                while (--p > x.a)*p = 0;
                *++p = g;
                *(i = s) = '-';
                c.get(*++i);
                do c.get(*++i); while (*i != ' ' && *i != '\t' && *i != '\n');
                *i = '\0';
                *x.a = val::ull(s);
                return c;
            }
            c.get(*i);
            do c.get(*++i); while (*i >= '0' && *i <= '9');
            if (*i == '+')
            {
                *i = '\0';
                q -= w - (v = val::ull(s));
                w = v;
                while (--p > q)*p = 0;
                *q = g;
                c.get(*(i = s));
            }
            else if (*i == '-')
            {
                *i = '\0';
                q -= w - (v = val::ull(s));
                w = v;
                while (--p > q)*p = 0;
                *q = g;
                *(i = s) = '-';
            }
            else
            {
                *i = '\0';
                q -= w - (v = val::ull(s));
                w = v;
                while (--p > q)*p = 0;
                *q = g;
                do *--q = 0; while (q > x.a);
                return c;
            }
            if (*i == 'x')
            {
                g = 1;
                goto F;
            }
        }
    } while (true);
}
template<class T, class N>
const dxs<T, N>& out(const dxs<T, N>& x, char ch = 'x')
{
    const T* i(x.a + x.n);
    unsigned long long n(x.n);
    if (!x.n)
    {
        std::cout << *i;
        return x;
    }
    if (*i != 1)
        if (*i == -1)
            std::cout << '-';
        else
            std::cout << *i;
    std::cout << ch;
    if (x.n == 1)
    {
        if (!*x.a)
            return x;
        else if (*x.a > 0)
            std::cout << '+';
        std::cout << *x.a;
        return x;
    }
    std::cout << '^' << n;
    while (--n > 1)
        if (*--i)
        {
            if (*i > 0)std::cout << '+';
            if (*i != 1)
                if (*i == -1)
                    std::cout << '-';
                else
                    std::cout << *i;
            std::cout << ch << '^' << n;
        }
    if (*--i)
    {
        if (*i > 0)std::cout << '+';
        if (*i != 1)
            if (*i == -1)
                std::cout << '-';
            else
                std::cout << *i;
        std::cout << ch;
    }
    if (*x.a)
    {
        if (*x.a > 0)std::cout << '+';
        std::cout << *x.a;
    }
    return x;
}
template<class T, class N>
dxs<T, N>& out(dxs<T, N>& x, char ch = 'x')
{
    T* i(x.a + x.n);
    unsigned long long n(x.n);
    if (!x.n)
    {
        std::cout << *i;
        return x;
    }
    if (*i != 1)
        if (*i == -1)
            std::cout << '-';
        else
            std::cout << *i;
    std::cout << ch;
    if (x.n == 1)
    {
        if (!*x.a)
            return x;
        else if (*x.a > 0)
            std::cout << '+';
        std::cout << *x.a;
        return x;
    }
    std::cout << '^' << n;
    while (--n > 1)
        if (*--i)
        {
            if (*i > 0)std::cout << '+';
            if (*i != 1)
                if (*i == -1)
                    std::cout << '-';
                else
                    std::cout << *i;
            std::cout << ch << '^' << n;
        }
    if (*--i)
    {
        if (*i > 0)std::cout << '+';
        if (*i != 1)
            if (*i == -1)
                std::cout << '-';
            else
                std::cout << *i;
        std::cout << ch;
    }
    if (*x.a)
    {
        if (*x.a > 0)std::cout << '+';
        std::cout << *x.a;
    }
    return x;
}
template<typename N>
dxs<float, N>& in(dxs<float, N>& x, char ch = 'x')
{
    delete[]x.a;
    char s[41], * i(s);
    do std::cin.get(*i); while (*i == ' ' || *i == '\t' || *i == '\n');
    if (*i == '+')std::cin.get(*i);
    bool f(*i == ch);
    while (*i != ch)
    {
        std::cin.get(*++i);
        if (*i == ' ' || *i == '\t' || *i == '\n')
        {
            *i = '\0';
            *(x.a = new float[1]) = val::f(s);
            x.n = 0;
            return x;
        }
    }
    *i = '\0';
    float* p, g(val::f(s));
    if (f)g = 1.f;
    std::cin.get(*(i = s));
    if (*i == ' ' || *i == '\t' || *i == '\n')
    {
        *i = '\0';
        x.n = 1;
        *(x.a = new float[2]) = 0.f;
        x.a[1] = g;
        return x;
    }
    if (*i == '+')
    {
        x.n = 1;
        *(p = (x.a = new float[2]) + 1) = g;
        std::cin.get(*i);
    }
    else if (*i == '-')
    {
        x.n = 1;
        *(p = (x.a = new float[2]) + 1) = g;
    }
    else
    {
        std::cin.get(*i);
        do std::cin.get(*++i); while (*i >= '0' && *i <= '9');
        if (*i == ' ' || *i == '\t' || *i == '\n')
        {
            *i = '\0';
            *(p = (x.a = new float[(x.n = val::ull(s)) + 1]) + x.n) = g;
            do *--p = 0.f; while (p > x.a);
            return x;
        }
        if (*i == '+')
        {
            *i = '\0';
            *(p = (x.a = new float[(x.n = val::ull(s)) + 1]) + x.n) = g;
            std::cin.get(*(i = s));
        }
        else
        {
            *i = '\0';
            *(p = (x.a = new float[(x.n = val::ull(s)) + 1]) + x.n) = g;
            *(i = s) = '-';
        }
    }
    float* q(p);
    N w(x.n), v;
    do
    {
        std::cin.get(*++i);
        if (*i == ' ' || *i == '\t' || *i == '\n')
        {
            *i = '\0';
            while (--p > x.a)*p = 0.f;
            *p = val::f(s);
            return x;
        }
        if (*i == ch)
        {
            *i = '\0';
            g = val::f(s);
        F:std::cin.get(*(i = s));
            if (*i == ' ' || *i == '\t' || *i == '\n')
            {
                while (--p > x.a)*p = 0.f;
                *p = g;
                *--p = 0.f;
                return x;
            }
            if (*i == '+')
            {
                *i = '\0';
                while (--p > x.a)*p = 0.f;
                *++p = g;
                std::cin.get(*(i = s));
                do std::cin.get(*++i); while (*i != ' ' && *i != '\t' && *i != '\n');
                *i = '\0';
                *x.a = val::f(s);
                return x;
            }
            if (*i == '-')
            {
                *i = '\0';
                while (--p > x.a)*p = 0.f;
                *++p = g;
                *(i = s) = '-';
                std::cin.get(*++i);
                do std::cin.get(*++i); while (*i != ' ' && *i != '\t' && *i != '\n');
                *i = '\0';
                *x.a = val::f(s);
                return x;
            }
            std::cin.get(*i);
            do std::cin.get(*++i); while (*i >= '0' && *i <= '9');
            if (*i == '+')
            {
                *i = '\0';
                q -= w - (v = val::ull(s));
                w = v;
                while (--p > q)*p = 0.f;
                *q = g;
                std::cin.get(*(i = s));
            }
            else if (*i == '-')
            {
                *i = '\0';
                q -= w - (v = val::ull(s));
                w = v;
                while (--p > q)*p = 0.f;
                *q = g;
                *(i = s) = '-';
            }
            else
            {
                *i = '\0';
                q -= w - (v = val::ull(s));
                w = v;
                while (--p > q)*p = 0.f;
                *q = g;
                do *--q = 0.f; while (q > x.a);
                return x;
            }
            if (*i == ch)
            {
                g = 1.f;
                goto F;
            }
        }
    } while (true);
}
template<typename N>
dxs<double, N>& in(dxs<double, N>& x, char ch = 'x')

{
    delete[]x.a;
    char s[41], * i(s);
    do std::cin.get(*i); while (*i == ' ' || *i == '\t' || *i == '\n');
    if (*i == '+')std::cin.get(*i);
    bool f(*i == ch);
    while (*i != ch)
    {
        std::cin.get(*++i);
        if (*i == ' ' || *i == '\t' || *i == '\n')
        {
            *i = '\0';
            *(x.a = new double[1]) = val::d(s);
            x.n = 0;
            return x;
        }
    }
    *i = '\0';
    double* p, g(val::d(s));
    if (f)g = 1.;
    std::cin.get(*(i = s));
    if (*i == ' ' || *i == '\t' || *i == '\n')
    {
        *i = '\0';
        x.n = 1;
        *(x.a = new double[2]) = 0.;
        x.a[1] = g;
        return x;
    }
    if (*i == '+')
    {
        x.n = 1;
        *(p = (x.a = new double[2]) + 1) = g;
        std::cin.get(*i);
    }
    else if (*i == '-')
    {
        x.n = 1;
        *(p = (x.a = new double[2]) + 1) = g;
    }
    else
    {
        std::cin.get(*i);
        do std::cin.get(*++i); while (*i >= '0' && *i <= '9');
        if (*i == ' ' || *i == '\t' || *i == '\n')
        {
            *i = '\0';
            *(p = (x.a = new double[(x.n = val::ull(s)) + 1]) + x.n) = g;
            do *--p = 0.; while (p > x.a);
            return x;
        }
        if (*i == '+')
        {
            *i = '\0';
            *(p = (x.a = new double[(x.n = val::ull(s)) + 1]) + x.n) = g;
            std::cin.get(*(i = s));
        }
        else
        {
            *i = '\0';
            *(p = (x.a = new double[(x.n = val::ull(s)) + 1]) + x.n) = g;
            *(i = s) = '-';
        }
    }
    double* q(p);
    N w(x.n), v;
    do
    {
        std::cin.get(*++i);
        if (*i == ' ' || *i == '\t' || *i == '\n')
        {
            *i = '\0';
            while (--p > x.a)*p = 0.;
            *p = val::d(s);
            return x;
        }
        if (*i == ch)
        {
            *i = '\0';
            g = val::d(s);
        F:std::cin.get(*(i = s));
            if (*i == ' ' || *i == '\t' || *i == '\n')
            {
                while (--p > x.a)*p = 0.;
                *p = g;
                *--p = 0.;
                return x;
            }
            if (*i == '+')
            {
                *i = '\0';
                while (--p > x.a)*p = 0.;
                *++p = g;
                std::cin.get(*(i = s));
                do std::cin.get(*++i); while (*i != ' ' && *i != '\t' && *i != '\n');
                *i = '\0';
                *x.a = val::d(s);
                return x;
            }
            if (*i == '-')
            {
                *i = '\0';
                while (--p > x.a)*p = 0.;
                *++p = g;
                *(i = s) = '-';
                std::cin.get(*++i);
                do std::cin.get(*++i); while (*i != ' ' && *i != '\t' && *i != '\n');
                *i = '\0';
                *x.a = val::d(s);
                return x;
            }
            std::cin.get(*i);
            do std::cin.get(*++i); while (*i >= '0' && *i <= '9');
            if (*i == '+')
            {
                *i = '\0';
                q -= w - (v = val::ull(s));
                w = v;
                while (--p > q)*p = 0.;
                *q = g;
                std::cin.get(*(i = s));
            }
            else if (*i == '-')
            {
                *i = '\0';
                q -= w - (v = val::ull(s));
                w = v;
                while (--p > q)*p = 0.;
                *q = g;
                *(i = s) = '-';
            }
            else
            {
                *i = '\0';
                q -= w - (v = val::ull(s));
                w = v;
                while (--p > q)*p = 0.;
                *q = g;
                do *--q = 0.; while (q > x.a);
                return x;
            }
            if (*i == ch)
            {
                g = 1.;
                goto F;
            }
        }
    } while (true);
}
template<typename N>
dxs<long double, N>& in(dxs<long double, N>& x, char ch = 'x')
{
    delete[]x.a;
    char s[41], * i(s);
    do std::cin.get(*i); while (*i == ' ' || *i == '\t' || *i == '\n');
    if (*i == '+')std::cin.get(*i);
    bool f(*i == ch);
    while (*i != ch)
    {
        std::cin.get(*++i);
        if (*i == ' ' || *i == '\t' || *i == '\n')
        {
            *i = '\0';
            *(x.a = new long double[1]) = val::ld(s);
            x.n = 0;
            return x;
        }
    }
    *i = '\0';
    long double* p, g(val::ld(s));
    if (f)g = 1.;
    std::cin.get(*(i = s));
    if (*i == ' ' || *i == '\t' || *i == '\n')
    {
        *i = '\0';
        x.n = 1;
        *(x.a = new long double[2]) = 0.;
        x.a[1] = g;
        return x;
    }
    if (*i == '+')
    {
        x.n = 1;
        *(p = (x.a = new long double[2]) + 1) = g;
        std::cin.get(*i);
    }
    else if (*i == '-')
    {
        x.n = 1;
        *(p = (x.a = new long double[2]) + 1) = g;
    }
    else
    {
        std::cin.get(*i);
        do std::cin.get(*++i); while (*i >= '0' && *i <= '9');
        if (*i == ' ' || *i == '\t' || *i == '\n')
        {
            *i = '\0';
            *(p = (x.a = new long double[(x.n = val::ull(s)) + 1]) + x.n) = g;
            do *--p = 0.; while (p > x.a);
            return x;
        }
        if (*i == '+')
        {
            *i = '\0';
            *(p = (x.a = new long double[(x.n = val::ull(s)) + 1]) + x.n) = g;
            std::cin.get(*(i = s));
        }
        else
        {
            *i = '\0';
            *(p = (x.a = new long double[(x.n = val::ull(s)) + 1]) + x.n) = g;
            *(i = s) = '-';
        }
    }
    long double* q(p);
    N w(x.n), v;
    do
    {
        std::cin.get(*++i);
        if (*i == ' ' || *i == '\t' || *i == '\n')
        {
            *i = '\0';
            while (--p > x.a)*p = 0.;
            *p = val::ld(s);
            return x;
        }
        if (*i == ch)
        {
            *i = '\0';
            g = val::ld(s);
        F:std::cin.get(*(i = s));
            if (*i == ' ' || *i == '\t' || *i == '\n')
            {
                while (--p > x.a)*p = 0.;
                *p = g;
                *--p = 0.;
                return x;
            }
            if (*i == '+')
            {
                *i = '\0';
                while (--p > x.a)*p = 0.;
                *++p = g;
                std::cin.get(*(i = s));
                do std::cin.get(*++i); while (*i != ' ' && *i != '\t' && *i != '\n');
                *i = '\0';
                *x.a = val::ld(s);
                return x;
            }
            if (*i == '-')
            {
                *i = '\0';
                while (--p > x.a)*p = 0.;
                *++p = g;
                *(i = s) = '-';
                std::cin.get(*++i);
                do std::cin.get(*++i); while (*i != ' ' && *i != '\t' && *i != '\n');
                *i = '\0';
                *x.a = val::ld(s);
                return x;
            }
            std::cin.get(*i);
            do std::cin.get(*++i); while (*i >= '0' && *i <= '9');
            if (*i == '+')
            {
                *i = '\0';
                q -= w - (v = val::ull(s));
                w = v;
                while (--p > q)*p = 0.;
                *q = g;
                std::cin.get(*(i = s));
            }
            else if (*i == '-')
            {
                *i = '\0';
                q -= w - (v = val::ull(s));
                w = v;
                while (--p > q)*p = 0.;
                *q = g;
                *(i = s) = '-';
            }
            else
            {
                *i = '\0';
                q -= w - (v = val::ull(s));
                w = v;
                while (--p > q)*p = 0.;
                *q = g;
                do *--q = 0.; while (q > x.a);
                return x;
            }
            if (*i == ch)
            {
                g = 1.;
                goto F;
            }
        }
    } while (true);
}
template<typename N>
dxs<int, N>& in(dxs<int, N>& x, char ch = 'x')
{
    delete[]x.a;
    char s[41], * i(s);
    do std::cin.get(*i); while (*i == ' ' || *i == '\t' || *i == '\n');
    if (*i == '+')std::cin.get(*i);
    bool f(*i == ch);
    while (*i != ch)
    {
        std::cin.get(*++i);
        if (*i == ' ' || *i == '\t' || *i == '\n')
        {
            *i = '\0';
            *(x.a = new int[1]) = val::i(s);
            x.n = 0;
            return x;
        }
    }
    *i = '\0';
    int* p, g(val::i(s));
    if (f)g = 1;
    std::cin.get(*(i = s));
    if (*i == ' ' || *i == '\t' || *i == '\n')
    {
        *i = '\0';
        x.n = 1;
        *(x.a = new int[2]) = 0;
        x.a[1] = g;
        return x;
    }
    if (*i == '+')
    {
        x.n = 1;
        *(p = (x.a = new int[2]) + 1) = g;
        std::cin.get(*i);
    }
    else if (*i == '-')
    {
        x.n = 1;
        *(p = (x.a = new int[2]) + 1) = g;
    }
    else
    {
        std::cin.get(*i);
        do std::cin.get(*++i); while (*i >= '0' && *i <= '9');
        if (*i == ' ' || *i == '\t' || *i == '\n')
        {
            *i = '\0';
            *(p = (x.a = new int[(x.n = val::ull(s)) + 1]) + x.n) = g;
            do *--p = 0; while (p > x.a);
            return x;
        }
        if (*i == '+')
        {
            *i = '\0';
            *(p = (x.a = new int[(x.n = val::ull(s)) + 1]) + x.n) = g;
            std::cin.get(*(i = s));
        }
        else
        {
            *i = '\0';
            *(p = (x.a = new int[(x.n = val::ull(s)) + 1]) + x.n) = g;
            *(i = s) = '-';
        }
    }
    int* q(p);
    N w(x.n), v;
    do
    {
        std::cin.get(*++i);
        if (*i == ' ' || *i == '\t' || *i == '\n')
        {
            *i = '\0';
            while (--p > x.a)*p = 0;
            *p = val::i(s);
            return x;
        }
        if (*i == ch)
        {
            *i = '\0';
            g = val::i(s);
        F:std::cin.get(*(i = s));
            if (*i == ' ' || *i == '\t' || *i == '\n')
            {
                while (--p > x.a)*p = 0;
                *p = g;
                *--p = 0;
                return x;
            }
            if (*i == '+')
            {
                *i = '\0';
                while (--p > x.a)*p = 0;
                *++p = g;
                std::cin.get(*(i = s));
                do std::cin.get(*++i); while (*i != ' ' && *i != '\t' && *i != '\n');
                *i = '\0';
                *x.a = val::i(s);
                return x;
            }
            if (*i == '-')
            {
                *i = '\0';
                while (--p > x.a)*p = 0;
                *++p = g;
                *(i = s) = '-';
                std::cin.get(*++i);
                do std::cin.get(*++i); while (*i != ' ' && *i != '\t' && *i != '\n');
                *i = '\0';
                *x.a = val::i(s);
                return x;
            }
            std::cin.get(*i);
            do std::cin.get(*++i); while (*i >= '0' && *i <= '9');
            if (*i == '+')
            {
                *i = '\0';
                q -= w - (v = val::ull(s));
                w = v;
                while (--p > q)*p = 0;
                *q = g;
                std::cin.get(*(i = s));
            }
            else if (*i == '-')
            {
                *i = '\0';
                q -= w - (v = val::ull(s));
                w = v;
                while (--p > q)*p = 0;
                *q = g;
                *(i = s) = '-';
            }
            else
            {
                *i = '\0';
                q -= w - (v = val::ull(s));
                w = v;
                while (--p > q)*p = 0;
                *q = g;
                do *--q = 0; while (q > x.a);
                return x;
            }
            if (*i == ch)
            {
                g = 1;
                goto F;
            }
        }
    } while (true);
}
template<typename N>
dxs<long, N>& in(dxs<long, N>& x, char ch = 'x')
{
    delete[]x.a;
    char s[41], * i(s);
    do std::cin.get(*i); while (*i == ' ' || *i == '\t' || *i == '\n');
    if (*i == '+')std::cin.get(*i);
    bool f(*i == ch);
    while (*i != ch)
    {
        std::cin.get(*++i);
        if (*i == ' ' || *i == '\t' || *i == '\n')
        {
            *i = '\0';
            *(x.a = new long[1]) = val::l(s);
            x.n = 0;
            return x;
        }
    }
    *i = '\0';
    long* p, g(val::l(s));
    if (f)g = 1;
    std::cin.get(*(i = s));
    if (*i == ' ' || *i == '\t' || *i == '\n')
    {
        *i = '\0';
        x.n = 1;
        *(x.a = new long[2]) = 0;
        x.a[1] = g;
        return x;
    }
    if (*i == '+')
    {
        x.n = 1;
        *(p = (x.a = new long[2]) + 1) = g;
        std::cin.get(*i);
    }
    else if (*i == '-')
    {
        x.n = 1;
        *(p = (x.a = new long[2]) + 1) = g;
    }
    else
    {
        std::cin.get(*i);
        do std::cin.get(*++i); while (*i >= '0' && *i <= '9');
        if (*i == ' ' || *i == '\t' || *i == '\n')
        {
            *i = '\0';
            *(p = (x.a = new long[(x.n = val::ull(s)) + 1]) + x.n) = g;
            do *--p = 0; while (p > x.a);
            return x;
        }
        if (*i == '+')
        {
            *i = '\0';
            *(p = (x.a = new long[(x.n = val::ull(s)) + 1]) + x.n) = g;
            std::cin.get(*(i = s));
        }
        else
        {
            *i = '\0';
            *(p = (x.a = new long[(x.n = val::ull(s)) + 1]) + x.n) = g;
            *(i = s) = '-';
        }
    }
    long* q(p);
    N w(x.n), v;
    do
    {
        std::cin.get(*++i);
        if (*i == ' ' || *i == '\t' || *i == '\n')
        {
            *i = '\0';
            while (--p > x.a)*p = 0;
            *p = val::l(s);
            return x;
        }
        if (*i == ch)
        {
            *i = '\0';
            g = val::l(s);
        F:std::cin.get(*(i = s));
            if (*i == ' ' || *i == '\t' || *i == '\n')
            {
                while (--p > x.a)*p = 0;
                *p = g;
                *--p = 0;
                return x;
            }
            if (*i == '+')
            {
                *i = '\0';
                while (--p > x.a)*p = 0;
                *++p = g;
                std::cin.get(*(i = s));
                do std::cin.get(*++i); while (*i != ' ' && *i != '\t' && *i != '\n');
                *i = '\0';
                *x.a = val::l(s);
                return x;
            }
            if (*i == '-')
            {
                *i = '\0';
                while (--p > x.a)*p = 0;
                *++p = g;
                *(i = s) = '-';
                std::cin.get(*++i);
                do std::cin.get(*++i); while (*i != ' ' && *i != '\t' && *i != '\n');
                *i = '\0';
                *x.a = val::l(s);
                return x;
            }
            std::cin.get(*i);
            do std::cin.get(*++i); while (*i >= '0' && *i <= '9');
            if (*i == '+')
            {
                *i = '\0';
                q -= w - (v = val::ull(s));
                w = v;
                while (--p > q)*p = 0;
                *q = g;
                std::cin.get(*(i = s));
            }
            else if (*i == '-')
            {
                *i = '\0';
                q -= w - (v = val::ull(s));
                w = v;
                while (--p > q)*p = 0;
                *q = g;
                *(i = s) = '-';
            }
            else
            {
                *i = '\0';
                q -= w - (v = val::ull(s));
                w = v;
                while (--p > q)*p = 0;
                *q = g;
                do *--q = 0; while (q > x.a);
                return x;
            }
            if (*i == ch)
            {
                g = 1;
                goto F;
            }
        }
    } while (true);
}
template<typename N>
dxs<short, N>& in(dxs<short, N>& x, char ch = 'x')
{
    delete[]x.a;
    char s[41], * i(s);
    do std::cin.get(*i); while (*i == ' ' || *i == '\t' || *i == '\n');
    if (*i == '+')std::cin.get(*i);
    bool f(*i == ch);
    while (*i != ch)
    {
        std::cin.get(*++i);
        if (*i == ' ' || *i == '\t' || *i == '\n')
        {
            *i = '\0';
            *(x.a = new short[1]) = val::s(s);
            x.n = 0;
            return x;
        }
    }
    *i = '\0';
    short* p, g(val::s(s));
    if (f)g = 1;
    std::cin.get(*(i = s));
    if (*i == ' ' || *i == '\t' || *i == '\n')
    {
        *i = '\0';
        x.n = 1;
        *(x.a = new short[2]) = 0;
        x.a[1] = g;
        return x;
    }
    if (*i == '+')
    {
        x.n = 1;
        *(p = (x.a = new short[2]) + 1) = g;
        std::cin.get(*i);
    }
    else if (*i == '-')
    {
        x.n = 1;
        *(p = (x.a = new short[2]) + 1) = g;
    }
    else
    {
        std::cin.get(*i);
        do std::cin.get(*++i); while (*i >= '0' && *i <= '9');
        if (*i == ' ' || *i == '\t' || *i == '\n')
        {
            *i = '\0';
            *(p = (x.a = new short[(x.n = val::ull(s)) + 1]) + x.n) = g;
            do *--p = 0; while (p > x.a);
            return x;
        }
        if (*i == '+')
        {
            *i = '\0';
            *(p = (x.a = new short[(x.n = val::ull(s)) + 1]) + x.n) = g;
            std::cin.get(*(i = s));
        }
        else
        {
            *i = '\0';
            *(p = (x.a = new short[(x.n = val::ull(s)) + 1]) + x.n) = g;
            *(i = s) = '-';
        }
    }
    short* q(p);
    N w(x.n), v;
    do
    {
        std::cin.get(*++i);
        if (*i == ' ' || *i == '\t' || *i == '\n')
        {
            *i = '\0';
            while (--p > x.a)*p = 0;
            *p = val::s(s);
            return x;
        }
        if (*i == ch)
        {
            *i = '\0';
            g = val::s(s);
        F:std::cin.get(*(i = s));
            if (*i == ' ' || *i == '\t' || *i == '\n')
            {
                while (--p > x.a)*p = 0;
                *p = g;
                *--p = 0;
                return x;
            }
            if (*i == '+')
            {
                *i = '\0';
                while (--p > x.a)*p = 0;
                *++p = g;
                std::cin.get(*(i = s));
                do std::cin.get(*++i); while (*i != ' ' && *i != '\t' && *i != '\n');
                *i = '\0';
                *x.a = val::s(s);
                return x;
            }
            if (*i == '-')
            {
                *i = '\0';
                while (--p > x.a)*p = 0;
                *++p = g;
                *(i = s) = '-';
                std::cin.get(*++i);
                do std::cin.get(*++i); while (*i != ' ' && *i != '\t' && *i != '\n');
                *i = '\0';
                *x.a = val::s(s);
                return x;
            }
            std::cin.get(*i);
            do std::cin.get(*++i); while (*i >= '0' && *i <= '9');
            if (*i == '+')
            {
                *i = '\0';
                q -= w - (v = val::ull(s));
                w = v;
                while (--p > q)*p = 0;
                *q = g;
                std::cin.get(*(i = s));
            }
            else if (*i == '-')
            {
                *i = '\0';
                q -= w - (v = val::ull(s));
                w = v;
                while (--p > q)*p = 0;
                *q = g;
                *(i = s) = '-';
            }
            else
            {
                *i = '\0';
                q -= w - (v = val::ull(s));
                w = v;
                while (--p > q)*p = 0;
                *q = g;
                do *--q = 0; while (q > x.a);
                return x;
            }
            if (*i == ch)
            {
                g = 1;
                goto F;
            }
        }
    } while (true);
}
template<typename N>
dxs<signed char, N>& in(dxs<signed char, N>& x, char ch = 'x')
{
    delete[]x.a;
    char s[41], * i(s);
    do std::cin.get(*i); while (*i == ' ' || *i == '\t' || *i == '\n');
    if (*i == '+')std::cin.get(*i);
    bool f(*i == ch);
    while (*i != ch)
    {
        std::cin.get(*++i);
        if (*i == ' ' || *i == '\t' || *i == '\n')
        {
            *i = '\0';
            *(x.a = new signed char[1]) = val::sc(s);
            x.n = 0;
            return x;
        }
    }
    *i = '\0';
    signed char* p, g(val::sc(s));
    if (f)g = 1;
    std::cin.get(*(i = s));
    if (*i == ' ' || *i == '\t' || *i == '\n')
    {
        *i = '\0';
        x.n = 1;
        *(x.a = new signed char[2]) = 0;
        x.a[1] = g;
        return x;
    }
    if (*i == '+')
    {
        x.n = 1;
        *(p = (x.a = new signed char[2]) + 1) = g;
        std::cin.get(*i);
    }
    else if (*i == '-')
    {
        x.n = 1;
        *(p = (x.a = new signed char[2]) + 1) = g;
    }
    else
    {
        std::cin.get(*i);
        do std::cin.get(*++i); while (*i >= '0' && *i <= '9');
        if (*i == ' ' || *i == '\t' || *i == '\n')
        {
            *i = '\0';
            *(p = (x.a = new signed char[(x.n = val::ull(s)) + 1]) + x.n) = g;
            do *--p = 0; while (p > x.a);
            return x;
        }
        if (*i == '+')
        {
            *i = '\0';
            *(p = (x.a = new signed char[(x.n = val::ull(s)) + 1]) + x.n) = g;
            std::cin.get(*(i = s));
        }
        else
        {
            *i = '\0';
            *(p = (x.a = new signed char[(x.n = val::ull(s)) + 1]) + x.n) = g;
            *(i = s) = '-';
        }
    }
    signed char* q(p);
    N w(x.n), v;
    do
    {
        std::cin.get(*++i);
        if (*i == ' ' || *i == '\t' || *i == '\n')
        {
            *i = '\0';
            while (--p > x.a)*p = 0;
            *p = val::sc(s);
            return x;
        }
        if (*i == ch)
        {
            *i = '\0';
            g = val::sc(s);
        F:std::cin.get(*(i = s));
            if (*i == ' ' || *i == '\t' || *i == '\n')
            {
                while (--p > x.a)*p = 0;
                *p = g;
                *--p = 0;
                return x;
            }
            if (*i == '+')
            {
                *i = '\0';
                while (--p > x.a)*p = 0;
                *++p = g;
                std::cin.get(*(i = s));
                do std::cin.get(*++i); while (*i != ' ' && *i != '\t' && *i != '\n');
                *i = '\0';
                *x.a = val::sc(s);
                return x;
            }
            if (*i == '-')
            {
                *i = '\0';
                while (--p > x.a)*p = 0;
                *++p = g;
                *(i = s) = '-';
                std::cin.get(*++i);
                do std::cin.get(*++i); while (*i != ' ' && *i != '\t' && *i != '\n');
                *i = '\0';
                *x.a = val::sc(s);
                return x;
            }
            std::cin.get(*i);
            do std::cin.get(*++i); while (*i >= '0' && *i <= '9');
            if (*i == '+')
            {
                *i = '\0';
                q -= w - (v = val::ull(s));
                w = v;
                while (--p > q)*p = 0;
                *q = g;
                std::cin.get(*(i = s));
            }
            else if (*i == '-')
            {
                *i = '\0';
                q -= w - (v = val::ull(s));
                w = v;
                while (--p > q)*p = 0;
                *q = g;
                *(i = s) = '-';
            }
            else
            {
                *i = '\0';
                q -= w - (v = val::ull(s));
                w = v;
                while (--p > q)*p = 0;
                *q = g;
                do *--q = 0; while (q > x.a);
                return x;
            }
            if (*i == ch)
            {
                g = 1;
                goto F;
            }
        }
    } while (true);
}
template<typename N>
dxs<unsigned char, N>& in(dxs<unsigned char, N>& x, char ch = 'x')
{
    delete[]x.a;
    char s[41], * i(s);
    do std::cin.get(*i); while (*i == ' ' || *i == '\t' || *i == '\n');
    if (*i == '+')std::cin.get(*i);
    bool f(*i == ch);
    while (*i != ch)
    {
        std::cin.get(*++i);
        if (*i == ' ' || *i == '\t' || *i == '\n')
        {
            *i = '\0';
            *(x.a = new unsigned char[1]) = val::uc(s);
            x.n = 0;
            return x;
        }
    }
    *i = '\0';
    unsigned char* p, g(val::uc(s));
    if (f)g = 1;
    std::cin.get(*(i = s));
    if (*i == ' ' || *i == '\t' || *i == '\n')
    {
        *i = '\0';
        x.n = 1;
        *(x.a = new unsigned char[2]) = 0;
        x.a[1] = g;
        return x;
    }
    if (*i == '+')
    {
        x.n = 1;
        *(p = (x.a = new unsigned char[2]) + 1) = g;
        std::cin.get(*i);
    }
    else if (*i == '-')
    {
        x.n = 1;
        *(p = (x.a = new unsigned char[2]) + 1) = g;
    }
    else
    {
        std::cin.get(*i);
        do std::cin.get(*++i); while (*i >= '0' && *i <= '9');
        if (*i == ' ' || *i == '\t' || *i == '\n')
        {
            *i = '\0';
            *(p = (x.a = new unsigned char[(x.n = val::ull(s)) + 1]) + x.n) = g;
            do *--p = 0; while (p > x.a);
            return x;
        }
        if (*i == '+')
        {
            *i = '\0';
            *(p = (x.a = new unsigned char[(x.n = val::ull(s)) + 1]) + x.n) = g;
            std::cin.get(*(i = s));
        }
        else
        {
            *i = '\0';
            *(p = (x.a = new unsigned char[(x.n = val::ull(s)) + 1]) + x.n) = g;
            *(i = s) = '-';
        }
    }
    unsigned char* q(p);
    N w(x.n), v;
    do
    {
        std::cin.get(*++i);
        if (*i == ' ' || *i == '\t' || *i == '\n')
        {
            *i = '\0';
            while (--p > x.a)*p = 0;
            *p = val::uc(s);
            return x;
        }
        if (*i == ch)
        {
            *i = '\0';
            g = val::uc(s);
        F:std::cin.get(*(i = s));
            if (*i == ' ' || *i == '\t' || *i == '\n')
            {
                while (--p > x.a)*p = 0;
                *p = g;
                *--p = 0;
                return x;
            }
            if (*i == '+')
            {
                *i = '\0';
                while (--p > x.a)*p = 0;
                *++p = g;
                std::cin.get(*(i = s));
                do std::cin.get(*++i); while (*i != ' ' && *i != '\t' && *i != '\n');
                *i = '\0';
                *x.a = val::uc(s);
                return x;
            }
            if (*i == '-')
            {
                *i = '\0';
                while (--p > x.a)*p = 0;
                *++p = g;
                *(i = s) = '-';
                std::cin.get(*++i);
                do std::cin.get(*++i); while (*i != ' ' && *i != '\t' && *i != '\n');
                *i = '\0';
                *x.a = val::uc(s);
                return x;
            }
            std::cin.get(*i);
            do std::cin.get(*++i); while (*i >= '0' && *i <= '9');
            if (*i == '+')
            {
                *i = '\0';
                q -= w - (v = val::ull(s));
                w = v;
                while (--p > q)*p = 0;
                *q = g;
                std::cin.get(*(i = s));
            }
            else if (*i == '-')
            {
                *i = '\0';
                q -= w - (v = val::ull(s));
                w = v;
                while (--p > q)*p = 0;
                *q = g;
                *(i = s) = '-';
            }
            else
            {
                *i = '\0';
                q -= w - (v = val::ull(s));
                w = v;
                while (--p > q)*p = 0;
                *q = g;
                do *--q = 0; while (q > x.a);
                return x;
            }
            if (*i == ch)
            {
                g = 1;
                goto F;
            }
        }
    } while (true);
}
template<typename N>
dxs<long long, N>& in(dxs<long long, N>& x, char ch = 'x')
{
    delete[]x.a;
    char s[41], * i(s);
    do std::cin.get(*i); while (*i == ' ' || *i == '\t' || *i == '\n');
    if (*i == '+')std::cin.get(*i);
    bool f(*i == ch);
    while (*i != ch)
    {
        std::cin.get(*++i);
        if (*i == ' ' || *i == '\t' || *i == '\n')
        {
            *i = '\0';
            *(x.a = new long long[1]) = val::ll(s);
            x.n = 0;
            return x;
        }
    }
    *i = '\0';
    long long* p, g(val::ll(s));
    if (f)g = 1;
    std::cin.get(*(i = s));
    if (*i == ' ' || *i == '\t' || *i == '\n')
    {
        *i = '\0';
        x.n = 1;
        *(x.a = new long long[2]) = 0;
        x.a[1] = g;
        return x;
    }
    if (*i == '+')
    {
        x.n = 1;
        *(p = (x.a = new long long[2]) + 1) = g;
        std::cin.get(*i);
    }
    else if (*i == '-')
    {
        x.n = 1;
        *(p = (x.a = new long long[2]) + 1) = g;
    }
    else
    {
        std::cin.get(*i);
        do std::cin.get(*++i); while (*i >= '0' && *i <= '9');
        if (*i == ' ' || *i == '\t' || *i == '\n')
        {
            *i = '\0';
            *(p = (x.a = new long long[(x.n = val::ull(s)) + 1]) + x.n) = g;
            do *--p = 0; while (p > x.a);
            return x;
        }
        if (*i == '+')
        {
            *i = '\0';
            *(p = (x.a = new long long[(x.n = val::ull(s)) + 1]) + x.n) = g;
            std::cin.get(*(i = s));
        }
        else
        {
            *i = '\0';
            *(p = (x.a = new long long[(x.n = val::ull(s)) + 1]) + x.n) = g;
            *(i = s) = '-';
        }
    }
    long long* q(p);
    N w(x.n), v;
    do
    {
        std::cin.get(*++i);
        if (*i == ' ' || *i == '\t' || *i == '\n')
        {
            *i = '\0';
            while (--p > x.a)*p = 0;
            *p = val::ll(s);
            return x;
        }
        if (*i == ch)
        {
            *i = '\0';
            g = val::ll(s);
        F:std::cin.get(*(i = s));
            if (*i == ' ' || *i == '\t' || *i == '\n')
            {
                while (--p > x.a)*p = 0;
                *p = g;
                *--p = 0;
                return x;
            }
            if (*i == '+')
            {
                *i = '\0';
                while (--p > x.a)*p = 0;
                *++p = g;
                std::cin.get(*(i = s));
                do std::cin.get(*++i); while (*i != ' ' && *i != '\t' && *i != '\n');
                *i = '\0';
                *x.a = val::ll(s);
                return x;
            }
            if (*i == '-')
            {
                *i = '\0';
                while (--p > x.a)*p = 0;
                *++p = g;
                *(i = s) = '-';
                std::cin.get(*++i);
                do std::cin.get(*++i); while (*i != ' ' && *i != '\t' && *i != '\n');
                *i = '\0';
                *x.a = val::ll(s);
                return x;
            }
            std::cin.get(*i);
            do std::cin.get(*++i); while (*i >= '0' && *i <= '9');
            if (*i == '+')
            {
                *i = '\0';
                q -= w - (v = val::ull(s));
                w = v;
                while (--p > q)*p = 0;
                *q = g;
                std::cin.get(*(i = s));
            }
            else if (*i == '-')
            {
                *i = '\0';
                q -= w - (v = val::ull(s));
                w = v;
                while (--p > q)*p = 0;
                *q = g;
                *(i = s) = '-';
            }
            else
            {
                *i = '\0';
                q -= w - (v = val::ull(s));
                w = v;
                while (--p > q)*p = 0;
                *q = g;
                do *--q = 0; while (q > x.a);
                return x;
            }
            if (*i == ch)
            {
                g = 1;
                goto F;
            }
        }
    } while (true);
}
template<typename N>
dxs<unsigned short, N>& in(dxs<unsigned short, N>& x, char ch = 'x')
{
    delete[]x.a;
    char s[41], * i(s);
    do std::cin.get(*i); while (*i == ' ' || *i == '\t' || *i == '\n');
    if (*i == '+')std::cin.get(*i);
    bool f(*i == ch);
    while (*i != ch)
    {
        std::cin.get(*++i);
        if (*i == ' ' || *i == '\t' || *i == '\n')
        {
            *i = '\0';
            *(x.a = new unsigned short[1]) = val::us(s);
            x.n = 0;
            return x;
        }
    }
    *i = '\0';
    unsigned short* p, g(val::us(s));
    if (f)g = 1;
    std::cin.get(*(i = s));
    if (*i == ' ' || *i == '\t' || *i == '\n')
    {
        *i = '\0';
        x.n = 1;
        *(x.a = new unsigned short[2]) = 0;
        x.a[1] = g;
        return x;
    }
    if (*i == '+')
    {
        x.n = 1;
        *(p = (x.a = new unsigned short[2]) + 1) = g;
        std::cin.get(*i);
    }
    else if (*i == '-')
    {
        x.n = 1;
        *(p = (x.a = new unsigned short[2]) + 1) = g;
    }
    else
    {
        std::cin.get(*i);
        do std::cin.get(*++i); while (*i >= '0' && *i <= '9');
        if (*i == ' ' || *i == '\t' || *i == '\n')
        {
            *i = '\0';
            *(p = (x.a = new unsigned short[(x.n = val::ull(s)) + 1]) + x.n) = g;
            do *--p = 0; while (p > x.a);
            return x;
        }
        if (*i == '+')
        {
            *i = '\0';
            *(p = (x.a = new unsigned short[(x.n = val::ull(s)) + 1]) + x.n) = g;
            std::cin.get(*(i = s));
        }
        else
        {
            *i = '\0';
            *(p = (x.a = new unsigned short[(x.n = val::ull(s)) + 1]) + x.n) = g;
            *(i = s) = '-';
        }
    }
    unsigned short* q(p);
    N w(x.n), v;
    do
    {
        std::cin.get(*++i);
        if (*i == ' ' || *i == '\t' || *i == '\n')
        {
            *i = '\0';
            while (--p > x.a)*p = 0;
            *p = val::us(s);
            return x;
        }
        if (*i == ch)
        {
            *i = '\0';
            g = val::us(s);
        F:std::cin.get(*(i = s));
            if (*i == ' ' || *i == '\t' || *i == '\n')
            {
                while (--p > x.a)*p = 0;
                *p = g;
                *--p = 0;
                return x;
            }
            if (*i == '+')
            {
                *i = '\0';
                while (--p > x.a)*p = 0;
                *++p = g;
                std::cin.get(*(i = s));
                do std::cin.get(*++i); while (*i != ' ' && *i != '\t' && *i != '\n');
                *i = '\0';
                *x.a = val::us(s);
                return x;
            }
            if (*i == '-')
            {
                *i = '\0';
                while (--p > x.a)*p = 0;
                *++p = g;
                *(i = s) = '-';
                std::cin.get(*++i);
                do std::cin.get(*++i); while (*i != ' ' && *i != '\t' && *i != '\n');
                *i = '\0';
                *x.a = val::us(s);
                return x;
            }
            std::cin.get(*i);
            do std::cin.get(*++i); while (*i >= '0' && *i <= '9');
            if (*i == '+')
            {
                *i = '\0';
                q -= w - (v = val::ull(s));
                w = v;
                while (--p > q)*p = 0;
                *q = g;
                std::cin.get(*(i = s));
            }
            else if (*i == '-')
            {
                *i = '\0';
                q -= w - (v = val::ull(s));
                w = v;
                while (--p > q)*p = 0;
                *q = g;
                *(i = s) = '-';
            }
            else
            {
                *i = '\0';
                q -= w - (v = val::ull(s));
                w = v;
                while (--p > q)*p = 0;
                *q = g;
                do *--q = 0; while (q > x.a);
                return x;
            }
            if (*i == ch)
            {
                g = 1;
                goto F;
            }
        }
    } while (true);
}
template<typename N>
dxs<unsigned, N>& in(dxs<unsigned, N>& x, char ch = 'x')
{
    delete[]x.a;
    char s[41], * i(s);
    do std::cin.get(*i); while (*i == ' ' || *i == '\t' || *i == '\n');
    if (*i == '+')std::cin.get(*i);
    bool f(*i == ch);
    while (*i != ch)
    {
        std::cin.get(*++i);
        if (*i == ' ' || *i == '\t' || *i == '\n')
        {
            *i = '\0';
            *(x.a = new unsigned[1]) = val::u(s);
            x.n = 0;
            return x;
        }
    }
    *i = '\0';
    unsigned* p, g(val::u(s));
    if (f)g = 1;
    std::cin.get(*(i = s));
    if (*i == ' ' || *i == '\t' || *i == '\n')
    {
        *i = '\0';
        x.n = 1;
        *(x.a = new unsigned[2]) = 0;
        x.a[1] = g;
        return x;
    }
    if (*i == '+')
    {
        x.n = 1;
        *(p = (x.a = new unsigned[2]) + 1) = g;
        std::cin.get(*i);
    }
    else if (*i == '-')
    {
        x.n = 1;
        *(p = (x.a = new unsigned[2]) + 1) = g;
    }
    else
    {
        std::cin.get(*i);
        do std::cin.get(*++i); while (*i >= '0' && *i <= '9');
        if (*i == ' ' || *i == '\t' || *i == '\n')
        {
            *i = '\0';
            *(p = (x.a = new unsigned[(x.n = val::ull(s)) + 1]) + x.n) = g;
            do *--p = 0; while (p > x.a);
            return x;
        }
        if (*i == '+')
        {
            *i = '\0';
            *(p = (x.a = new unsigned[(x.n = val::ull(s)) + 1]) + x.n) = g;
            std::cin.get(*(i = s));
        }
        else
        {
            *i = '\0';
            *(p = (x.a = new unsigned[(x.n = val::ull(s)) + 1]) + x.n) = g;
            *(i = s) = '-';
        }
    }
    unsigned* q(p);
    N w(x.n), v;
    do
    {
        std::cin.get(*++i);
        if (*i == ' ' || *i == '\t' || *i == '\n')
        {
            *i = '\0';
            while (--p > x.a)*p = 0;
            *p = val::u(s);
            return x;
        }
        if (*i == ch)
        {
            *i = '\0';
            g = val::u(s);
        F:std::cin.get(*(i = s));
            if (*i == ' ' || *i == '\t' || *i == '\n')
            {
                while (--p > x.a)*p = 0;
                *p = g;
                *--p = 0;
                return x;
            }
            if (*i == '+')
            {
                *i = '\0';
                while (--p > x.a)*p = 0;
                *++p = g;
                std::cin.get(*(i = s));
                do std::cin.get(*++i); while (*i != ' ' && *i != '\t' && *i != '\n');
                *i = '\0';
                *x.a = val::u(s);
                return x;
            }
            if (*i == '-')
            {
                *i = '\0';
                while (--p > x.a)*p = 0;
                *++p = g;
                *(i = s) = '-';
                std::cin.get(*++i);
                do std::cin.get(*++i); while (*i != ' ' && *i != '\t' && *i != '\n');
                *i = '\0';
                *x.a = val::u(s);
                return x;
            }
            std::cin.get(*i);
            do std::cin.get(*++i); while (*i >= '0' && *i <= '9');
            if (*i == '+')
            {
                *i = '\0';
                q -= w - (v = val::ull(s));
                w = v;
                while (--p > q)*p = 0;
                *q = g;
                std::cin.get(*(i = s));
            }
            else if (*i == '-')
            {
                *i = '\0';
                q -= w - (v = val::ull(s));
                w = v;
                while (--p > q)*p = 0;
                *q = g;
                *(i = s) = '-';
            }
            else
            {
                *i = '\0';
                q -= w - (v = val::ull(s));
                w = v;
                while (--p > q)*p = 0;
                *q = g;
                do *--q = 0; while (q > x.a);
                return x;
            }
            if (*i == ch)
            {
                g = 1;
                goto F;
            }
        }
    } while (true);
}
template<typename N>
dxs<unsigned long, N>& in(dxs<unsigned long, N>& x, char ch = 'x')
{
    delete[]x.a;
    char s[41], * i(s);
    do std::cin.get(*i); while (*i == ' ' || *i == '\t' || *i == '\n');
    if (*i == '+')std::cin.get(*i);
    bool f(*i == ch);
    while (*i != ch)
    {
        std::cin.get(*++i);
        if (*i == ' ' || *i == '\t' || *i == '\n')
        {
            *i = '\0';
            *(x.a = new unsigned long[1]) = val::ul(s);
            x.n = 0;
            return x;
        }
    }
    *i = '\0';
    unsigned long* p, g(val::ul(s));
    if (f)g = 1;
    std::cin.get(*(i = s));
    if (*i == ' ' || *i == '\t' || *i == '\n')
    {
        *i = '\0';
        x.n = 1;
        *(x.a = new unsigned long[2]) = 0;
        x.a[1] = g;
        return x;
    }
    if (*i == '+')
    {
        x.n = 1;
        *(p = (x.a = new unsigned long[2]) + 1) = g;
        std::cin.get(*i);
    }
    else if (*i == '-')
    {
        x.n = 1;
        *(p = (x.a = new unsigned long[2]) + 1) = g;
    }
    else
    {
        std::cin.get(*i);
        do std::cin.get(*++i); while (*i >= '0' && *i <= '9');
        if (*i == ' ' || *i == '\t' || *i == '\n')
        {
            *i = '\0';
            *(p = (x.a = new unsigned long[(x.n = val::ull(s)) + 1]) + x.n) = g;
            do *--p = 0; while (p > x.a);
            return x;
        }
        if (*i == '+')
        {
            *i = '\0';
            *(p = (x.a = new unsigned long[(x.n = val::ull(s)) + 1]) + x.n) = g;
            std::cin.get(*(i = s));
        }
        else
        {
            *i = '\0';
            *(p = (x.a = new unsigned long[(x.n = val::ull(s)) + 1]) + x.n) = g;
            *(i = s) = '-';
        }
    }
    unsigned long* q(p);
    N w(x.n), v;
    do
    {
        std::cin.get(*++i);
        if (*i == ' ' || *i == '\t' || *i == '\n')
        {
            *i = '\0';
            while (--p > x.a)*p = 0;
            *p = val::ul(s);
            return x;
        }
        if (*i == ch)
        {
            *i = '\0';
            g = val::ul(s);
        F:std::cin.get(*(i = s));
            if (*i == ' ' || *i == '\t' || *i == '\n')
            {
                while (--p > x.a)*p = 0;
                *p = g;
                *--p = 0;
                return x;
            }
            if (*i == '+')
            {
                *i = '\0';
                while (--p > x.a)*p = 0;
                *++p = g;
                std::cin.get(*(i = s));
                do std::cin.get(*++i); while (*i != ' ' && *i != '\t' && *i != '\n');
                *i = '\0';
                *x.a = val::ul(s);
                return x;
            }
            if (*i == '-')
            {
                *i = '\0';
                while (--p > x.a)*p = 0;
                *++p = g;
                *(i = s) = '-';
                std::cin.get(*++i);
                do std::cin.get(*++i); while (*i != ' ' && *i != '\t' && *i != '\n');
                *i = '\0';
                *x.a = val::ul(s);
                return x;
            }
            std::cin.get(*i);
            do std::cin.get(*++i); while (*i >= '0' && *i <= '9');
            if (*i == '+')
            {
                *i = '\0';
                q -= w - (v = val::ull(s));
                w = v;
                while (--p > q)*p = 0;
                *q = g;
                std::cin.get(*(i = s));
            }
            else if (*i == '-')
            {
                *i = '\0';
                q -= w - (v = val::ull(s));
                w = v;
                while (--p > q)*p = 0;
                *q = g;
                *(i = s) = '-';
            }
            else
            {
                *i = '\0';
                q -= w - (v = val::ull(s));
                w = v;
                while (--p > q)*p = 0;
                *q = g;
                do *--q = 0; while (q > x.a);
                return x;
            }
            if (*i == ch)
            {
                g = 1;
                goto F;
            }
        }
    } while (true);
}
template<typename N>
dxs<unsigned long long, N>& in(dxs<unsigned long long, N>& x, char ch = 'x')
{
    delete[]x.a;
    char s[41], * i(s);
    do std::cin.get(*i); while (*i == ' ' || *i == '\t' || *i == '\n');
    if (*i == '+')std::cin.get(*i);
    bool f(*i == ch);
    while (*i != ch)
    {
        std::cin.get(*++i);
        if (*i == ' ' || *i == '\t' || *i == '\n')
        {
            *i = '\0';
            *(x.a = new unsigned long long[1]) = val::ull(s);
            x.n = 0;
            return x;
        }
    }
    *i = '\0';
    unsigned long long* p, g(val::ull(s));
    if (f)g = 1;
    std::cin.get(*(i = s));
    if (*i == ' ' || *i == '\t' || *i == '\n')
    {
        *i = '\0';
        x.n = 1;
        *(x.a = new unsigned long long[2]) = 0;
        x.a[1] = g;
        return x;
    }
    if (*i == '+')
    {
        x.n = 1;
        *(p = (x.a = new unsigned long long[2]) + 1) = g;
        std::cin.get(*i);
    }
    else if (*i == '-')
    {
        x.n = 1;
        *(p = (x.a = new unsigned long long[2]) + 1) = g;
    }
    else
    {
        std::cin.get(*i);
        do std::cin.get(*++i); while (*i >= '0' && *i <= '9');
        if (*i == ' ' || *i == '\t' || *i == '\n')
        {
            *i = '\0';
            *(p = (x.a = new unsigned long long[(x.n = val::ull(s)) + 1]) + x.n) = g;
            do *--p = 0; while (p > x.a);
            return x;
        }
        if (*i == '+')
        {
            *i = '\0';
            *(p = (x.a = new unsigned long long[(x.n = val::ull(s)) + 1]) + x.n) = g;
            std::cin.get(*(i = s));
        }
        else
        {
            *i = '\0';
            *(p = (x.a = new unsigned long long[(x.n = val::ull(s)) + 1]) + x.n) = g;
            *(i = s) = '-';
        }
    }
    unsigned long long* q(p);
    N w(x.n), v;
    do
    {
        std::cin.get(*++i);
        if (*i == ' ' || *i == '\t' || *i == '\n')
        {
            *i = '\0';
            while (--p > x.a)*p = 0;
            *p = val::ull(s);
            return x;
        }
        if (*i == ch)
        {
            *i = '\0';
            g = val::ull(s);
        F:std::cin.get(*(i = s));
            if (*i == ' ' || *i == '\t' || *i == '\n')
            {
                while (--p > x.a)*p = 0;
                *p = g;
                *--p = 0;
                return x;
            }
            if (*i == '+')
            {
                *i = '\0';
                while (--p > x.a)*p = 0;
                *++p = g;
                std::cin.get(*(i = s));
                do std::cin.get(*++i); while (*i != ' ' && *i != '\t' && *i != '\n');
                *i = '\0';
                *x.a = val::ull(s);
                return x;
            }
            if (*i == '-')
            {
                *i = '\0';
                while (--p > x.a)*p = 0;
                *++p = g;
                *(i = s) = '-';
                std::cin.get(*++i);
                do std::cin.get(*++i); while (*i != ' ' && *i != '\t' && *i != '\n');
                *i = '\0';
                *x.a = val::ull(s);
                return x;
            }
            std::cin.get(*i);
            do std::cin.get(*++i); while (*i >= '0' && *i <= '9');
            if (*i == '+')
            {
                *i = '\0';
                q -= w - (v = val::ull(s));
                w = v;
                while (--p > q)*p = 0;
                *q = g;
                std::cin.get(*(i = s));
            }
            else if (*i == '-')
            {
                *i = '\0';
                q -= w - (v = val::ull(s));
                w = v;
                while (--p > q)*p = 0;
                *q = g;
                *(i = s) = '-';
            }
            else
            {
                *i = '\0';
                q -= w - (v = val::ull(s));
                w = v;
                while (--p > q)*p = 0;
                *q = g;
                do *--q = 0; while (q > x.a);
                return x;
            }
            if (*i == ch)
            {
                g = 1;
                goto F;
            }
        }
    } while (true);
}
#endif
