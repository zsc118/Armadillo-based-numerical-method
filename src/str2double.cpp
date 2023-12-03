#include <string>
#include <vector>
#include <math.h>
#include <armadillo>
#include <QStringList>
using namespace std;
static char str_data[24];
using arma::vec;
namespace CAL {
vector<string> func{"sqrt","sin","log","exp","pow","abs","cos","tan","asin","acos","atan"};
}
static inline string dou_to_str(const double& x)
{
    sprintf(str_data,"%.15E",x);
    return str_data;
}
namespace _MYFUNCTION {
extern string &substring_replace(string&,const vector<string>&,const vector<string>&,const vector<string>&);
extern vector<string> vec_to_string(const vector<double>&);
vector<string> vec_to_string(const vec& x)
{
    vector<string> y(x.size());
    auto i=x.cbegin();
    auto k=y.begin();
    while(i<x.cend())*k++=*i++;
    return y;
}
// 仅支持四则运算的计算器
double _cal(string s)
{
    string s1;
    string::iterator a, b, i; // a记录'(',b记录')'
    string::size_type n = 0;
    int l;
    if (false)
    F:
        s.replace(a, i + 1, dou_to_str(_cal(string(a + 1, i))).c_str()); // 递归计算括号
R:
    i = s.begin();
    while (i < s.end())
    {
        if (*i == ' ' || *i == '\t')
        {
            s.erase(i);
            goto R;
        }
        if (*i == '(')
        {
            if (!n++)
                a = i;
        }
        else if (*i == ')')
        {
            if (!--n)
                goto F;
        }
        else if (*i != '\0' && *i != '.' && *i != 'e' && *i != '*' && *i != '-' && *i != '+' && *i != '/' && (*i < '0' || *i > '9') || (*i == '/' && *(i + 1) == '0' && *(i + 2) != '.'))
            return NAN; // 不是数学表达式
        ++i;
    }
    if (n)
        return NAN; // 左括号个数不等于右括号个数
    vector<string::size_type> u;
    /* k用作bool数组
     * 第一位表示第一个数字的正负性
     * 第二位表示第二个数字的正负性
     * 第三位表示第一个数字是否为字符串首个数字
     */
    unsigned char k('\0');
    // 计算乘除
    for (b = a = s.begin(); a < s.end(); ++a)
        if (*a == '*' || *a == '/')
            u.push_back(a - b);
    while (u.size())
    {
        i = a = b = s.begin() + u[0];
        while ((*--b >= '0' && *b <= '9') || *b == '.' || *b == 'e')
            if (b == s.begin())
            {
                k = 4;
                goto G;
            }
        // n记录第一个数字的首位的位置
        k += *b == '-';
        if (*b == '-')
            n = b - s.begin();
        else
            n = ++b - s.begin();
    G:
        double t = stod(string(b, a++));
        do
            if (++a == s.end())
                break;
        while ((*a >= '0' && *a <= '9') || *a == '.' || *a == 'e');
        double w = stod(string(i + 1, a));
        k += (w < 0) * 2;
        if (*i == '*')
        {
            l = (s1 = dou_to_str(t * w)).size() - (a - b);
            s.replace(b, a, s1.c_str());
        }
        else if (w)
        {
            l = (s1 = dou_to_str(t / w)).size() - (a - b);
            s.replace(b, a, s1.c_str());
        }
        else
            return NAN; // 除数为0
        switch (k)
        {
        case 2:
            s.erase(s.begin() + n - 1);
            break; // 删除-号前的+
        case 3:
            s.insert(s.begin() + n, '+'); // 在计算结果之前加一个+
        }
        u.erase(u.begin());
        for (auto &q : u)
            q += l;
    }
    // 计算加减
H:
    for (a = (b = s.begin()) + 1; a < s.end(); ++a)
        if (*a == '+')
        {
            if (*(a + 1) == '-' || *(a + 1) == '+')
            {
                s.erase(a);
                goto H;
            }
            u.push_back(a - b);
        }
        else if (*a == '-')
        {
            if (*(a + 1) == '-')
            {
                s.replace(a, a + 2, "+");
                goto H;
            }
            else if (*(a + 1) == '+')
            {
                s.erase(a + 1);
                goto H;
            }
            u.push_back(a - b);
        }
        else if(*a=='e'||*a=='E')
            if(*(a+1)=='-')
                ++a;
    while (u.size())
    {
        i = a = (b = s.begin()) + u[0];
        double t = stod(string(b, a++));
        do
            if (++a == s.end())
                break;
        while ((*a >= '0' && *a <= '9') || *a == '.' || *a == 'e');
        if (*i++ == '+')
            l = (s1 = dou_to_str(t + stod(string(i, a)))).size() - (a - b);
        else
            l = (s1 = dou_to_str(t - stod(string(i, a)))).size() - (a - b);
        s.replace(b, a, s1.c_str());
        u.erase(u.begin());
        for (auto &q : u)
            q += l;
    }
    return stod(s);
}
}
/*
 * 计算器
 * 计算字符串s代表的值
 * 支持的函数：sqrt,sin,log,exp,pow,abs,cos,tan,asin,acos,atan
 */
double calStr(string s)
{
    if(s.empty())return 0;
    string::iterator i(s.begin());
    do
        switch (*i)
        {
        case 's':
        {
            string::iterator j(i);
            if (*++j == 'i')
            {
                if (*++j != 'n')
                    return NAN;
                if (*++j != '(')
                    return NAN;
                string::iterator k(++j);
                string::size_type n = 1;
                if (*k == '(')
                    ++n;
                else if (*k == ')')
                    return NAN;
                while (++k != s.end())
                    switch (*k)
                    {
                    case '(':
                        ++n;
                        break;
                    case ')':
                        if (!--n)
                            return calStr(s.replace(i, k + 1, dou_to_str(sin(calStr(string(j, k))))));
                    }
                return NAN;
            }
            if (*j != 'q')
                return NAN;
            if (*++j != 'r')
                return NAN;
            if (*++j != 't')
                return NAN;
            if (*++j != '(')
                return NAN;
            string::iterator k(++j);
            string::size_type n = 1;
            if (*k == '(')
                ++n;
            else if (*k == ')')
                return NAN;
            while (++k != s.end())
                switch (*k)
                {
                case '(':
                    ++n;
                    break;
                case ')':
                    if (!--n)
                    {
                        double res(calStr(string(j, k)));
                        if(isnan(res)||res<0)return NAN;
                        return calStr(s.replace(i, k + 1, dou_to_str(sqrt(res))));
                    }
                }
            return NAN;
        }
        case 'l':
        {
            string::iterator j(i);
            if (*++j != 'o')
                return NAN;
            if (*++j != 'g')
                return NAN;
            if (*++j != '(')
                return NAN;
            string::iterator k(++j);
            string::size_type n = 1;
            if (*k == '(')
                ++n;
            else if (*k == ')')
                return NAN;
            while (++k != s.end())
                switch (*k)
                {
                case '(':
                    ++n;
                    break;
                case ')':
                    if (!--n)
                    {
                        double res(calStr(string(j, k)));
                        if(isnan(res)||res<=0)return NAN;
                        return calStr(s.replace(i, k + 1, dou_to_str(log(res))));
                    }
                }
            return NAN;
        }
        case 'e':
        {
            string::iterator j(i);
            if (*++j != 'x')
                continue;
            if (*++j != 'p')
                return NAN;
            if (*++j != '(')
                return NAN;
            string::iterator k(++j);
            string::size_type n = 1;
            if (*k == '(')
                ++n;
            else if (*k == ')')
                return NAN;
            while (++k != s.end())
                switch (*k)
                {
                case '(':
                    ++n;
                    break;
                case ')':
                    if (!--n)
                        return calStr(s.replace(i, k + 1, dou_to_str(exp(calStr(string(j, k))))));
                }
            return NAN;
        }
        case 'p':
        {
            string::iterator j(i);
            if (*++j != 'o')
            {
                if (*j == 'i')
                    return calStr(s.replace(i, ++j, "3.1415926535897932384626433832795"));
                return NAN;
            }
            if (*++j != 'w')
                return NAN;
            if (*++j != '(')
                return NAN;
            string::iterator k(++j), x;
            string::size_type n = 1;
            if (*k == '(')
                ++n;
            else if (*k == ')')
                return NAN;
            while (++k != s.end())
                switch (*k)
                {
                case ',':
                    x = k;
                    break;
                case '(':
                    ++n;
                    break;
                case ')':
                    if (!--n)
                        return calStr(s.replace(i, k + 1, dou_to_str(pow(calStr(string(j, x)), calStr(string(x + 1, k))))));
                }
            return NAN;
        }
        case 'c':
        {
            string::iterator j(i);
            if (*++j != 'o')
                return NAN;
            if (*++j != 's')
                return NAN;
            if (*++j != '(')
                return NAN;
            string::iterator k(++j);
            string::size_type n = 1;
            if (*k == '(')
                ++n;
            else if (*k == ')')
                return NAN;
            while (++k != s.end())
                switch (*k)
                {
                case '(':
                    ++n;
                    break;
                case ')':
                    if (!--n)
                        return calStr(s.replace(i, k + 1, dou_to_str(cos(calStr(string(j, k))))));
                }
            return NAN;
        }
        case 't':
        {
            string::iterator j(i);
            if (*++j != 'a')
                return NAN;
            if (*++j != 'n')
                return NAN;
            if (*++j != '(')
                return NAN;
            string::iterator k(++j);
            string::size_type n = 1;
            if (*k == '(')
                ++n;
            else if (*k == ')')
                return NAN;
            while (++k != s.end())
                switch (*k)
                {
                case '(':
                    ++n;
                    break;
                case ')':
                    if (!--n)
                        return calStr(s.replace(i, k + 1, dou_to_str(tan(calStr(string(j, k))))));
                }
            return NAN;
        }
        case 'a':
        {
            string::iterator j(i);
            switch (*++j)
            {
            case 's':
            {
                if (*++j != 'i')
                    return NAN;
                if (*++j != 'n')
                    return NAN;
                if (*++j != '(')
                    return NAN;
                string::iterator k(++j);
                string::size_type n = 1;
                if (*k == '(')
                    ++n;
                else if (*k == ')')
                    return NAN;
                while (++k != s.end())
                    switch (*k)
                    {
                    case '(':
                        ++n;
                        break;
                    case ')':
                        if (!--n)
                            return calStr(s.replace(i, k + 1, dou_to_str(asin(calStr(string(j, k))))));
                    }
                return NAN;
            }
            case 'c':
            {
                if (*++j != 'o')
                    return NAN;
                if (*++j != 's')
                    return NAN;
                if (*++j != '(')
                    return NAN;
                string::iterator k(++j);
                string::size_type n = 1;
                if (*k == '(')
                    ++n;
                else if (*k == ')')
                    return NAN;
                while (++k != s.end())
                    switch (*k)
                    {
                    case '(':
                        ++n;
                        break;
                    case ')':
                        if (!--n)
                            return calStr(s.replace(i, k + 1, dou_to_str(acos(calStr(string(j, k))))));
                    }
                return NAN;
            }
            case 't':
            {
                if (*++j != 'a')
                    return NAN;
                if (*++j != 'n')
                    return NAN;
                if (*++j != '(')
                    return NAN;
                string::iterator k(++j);
                string::size_type n = 1;
                if (*k == '(')
                    ++n;
                else if (*k == ')')
                    return NAN;
                while (++k != s.end())
                    switch (*k)
                    {
                    case '(':
                        ++n;
                        break;
                    case ')':
                        if (!--n)
                            return calStr(s.replace(i, k + 1, dou_to_str(atan(calStr(string(j, k))))));
                    }
                return NAN;
            }
            }
            if (*j != 'b')
                return NAN;
            if (*++j != 's')
                return NAN;
            if (*++j != '(')
                return NAN;
            string::iterator k(++j);
            string::size_type n = 1;
            if (*k == '(')
                ++n;
            else if (*k == ')')
                return NAN;
            while (++k != s.end())
                switch (*k)
                {
                case '(':
                    ++n;
                    break;
                case ')':
                    if (!--n)
                        return calStr(s.replace(i, k + 1, dou_to_str(abs(calStr(string(j, k))))));
                }
            return NAN;
        }
        }
    while (++i != s.end());
    return _MYFUNCTION::_cal(s);
}
//在表达式中带有x
double calStr_x(string s,const double& x)
{
    for(auto p(s.begin());p!=s.end();++p)
        if(*p=='x')
            if(p==s.begin()||p+1==s.end()||*(p-1)!='e'||*(p+1)!='p')
                return calStr_x(s.replace(p,p+1,dou_to_str(x)),x);
    return calStr(s);
}
//在表达式中带有x,y
double calStr_xy(string s,const double& x,const double& y)
{
    for(auto p(s.begin());p!=s.end();++p)
        if(*p=='y')
            return calStr_xy(s.replace(p,p+1,dou_to_str(y)),x,y);
        else if(*p=='x')
            if(p==s.begin()||p+1==s.end()||*(p-1)!='e'||*(p+1)!='p')
                return calStr_xy(s.replace(p,p+1,dou_to_str(x)),x,y);
    return calStr(s);
}
//在表达式中带有多个未知数
double calStr_f(string s,const vector<string>& x,const vector<string>& name)
{
    return calStr(_MYFUNCTION::substring_replace(s,name,x,CAL::func));
}
//在表达式中带有多个未知数
double calStr_f(string s,const vector<double>& x,const vector<string>& name)
{
    return calStr(_MYFUNCTION::substring_replace(s,name,_MYFUNCTION::vec_to_string(x),CAL::func));
}
//在表达式中带有多个未知数
double calStr_f(string s,const vec& x,const vector<string>& name)
{
    return calStr(_MYFUNCTION::substring_replace(s,name,_MYFUNCTION::vec_to_string(x),CAL::func));
}
//计算向量函数
vec calStr_F(const QStringList& S,const vec& x,const vector<string>& name)
{
    unsigned n(S.size());
    if(!n)return vec();
    auto p(S.cbegin());
    vec res(n,arma::fill::none);
    double* p_res(&res.at(0));
    *p_res=calStr_f(p->toStdString(),x,name);
    while(--n)*++p_res=calStr_f((++p)->toStdString(),x,name);
    return res;
}
