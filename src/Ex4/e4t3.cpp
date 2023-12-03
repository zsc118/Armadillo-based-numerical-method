#include "e4t3.h"
#include "ui_e4t3.h"
#include <stdio.h>
#include <QStandardPaths>
#include <string>
#include <algorithm>
#include <QMessageBox>
#include <QFileDialog>
#include "e4t3_1.h"
using namespace arma;
using namespace std;
#ifdef DEBUG43
#include <QDebug>
#endif
namespace _MYFUNCTION {
extern bool check_unique(const vec&);
//const unsigned max_unsigned(-1);
}
namespace Ex4 {
#ifdef DEBUG43E1
extern void Newton(const vec& x0,const vec& y0,const vec& x1,vec& y1);
#endif
/*
 * 分段线性插值
 * x0:观测点
 * y0:观测数据
 * x1:插值点
 * y1:预测值
 */
void piecewise_linear_interpolation(const vec& x0,const mat& y0,const vec&x1,vec& y1)
{
    if(x1.empty())
    {
        y1.clear();
        return;
    }
    if(x0.empty()||y0.empty())throw "观测集为空！";
    if(y0.n_cols!=1)throw "分段线性插值不需要导数信息！";
    if(x0.n_elem!=y0.n_elem)throw "观测点与观测数据向量大小不一致！";
    if(x0.n_elem==1)throw "观测点数量小于2！";
    if(_MYFUNCTION::check_unique(x0))throw "观测集包含重复元素！";
    const unsigned &n0(x0.n_elem),&n1(x1.n_elem);
    vec vx0(x0);
    uvec index(sort_index(x0));
    vx0.elem(index);
    mat vy0(y0.rows(index));
    y1.set_size(n1);
    unsigned u(n1);
    double *begin(&vx0.at(0)),*end(&vx0.at(n0-1)),*ybegin(&vy0.at(0)),*yend(&vy0.at(n0-1));
    do
    {
        const double& x(x1.at(--u));
//        unsigned k(_MYFUNCTION::max_unsigned);//顺序查找
//        while(++k!=x0.n_elem)
//            if(vx0.at(k)>x)
//                break;
        double*i(begin),*j(end);
        if(x<*i)
            y1.at(u)=*ybegin+(*(ybegin+1)-*ybegin)*(x-*begin)/(*(begin+1)-*begin);
        else if(x>*j)
            y1.at(u)=*yend+(*(yend-1)-*yend)*(x-*end)/(*(end-1)-*end);
        else
        {
            do//二分查找
            {
                double*m(i+(j-i)/2);
                if(*m>x)
                    --(j=m);
                else
                    ++(i=m);
            }while(i<=j);
            unsigned k(i-begin);
            double*yi(ybegin+k),*yj(yi-1);
            y1.at(u)=*yj+(*yi-*yj)*(x-*j)/(*i-*j);
        }
    }while(u);
}
/*
 * 分段Hermite插值
 * x0:观测点
 * y0:观测数据(第i列为i阶导数)
 * x1:插值点
 * y1:预测值
 */
void piecewise_Hermite_interpolation(const vec& x0,const mat& y0,const vec&x1,vec& y1)
{
    if(x1.empty())
    {
        y1.clear();
        return;
    }
    if(x0.empty()||y0.empty())throw "观测集为空！";
    if(y0.n_cols!=2)throw "导数信息有误！";
    if(x0.n_elem!=y0.n_rows)throw "观测点与观测数据向量大小不一致！";
    if(x0.n_elem==1)throw "观测点数量小于2！";
    if(_MYFUNCTION::check_unique(x0))throw "观测集包含重复元素！";
    const unsigned &n0(x0.n_elem),&n1(x1.n_elem);
    vec vx0(x0);
    uvec index(sort_index(x0));
    vx0.elem(index);
    mat vy0(y0.rows(index));
    y1.set_size(n1);
    unsigned u(n1);
    double *begin(&vx0.at(0)),*end(&vx0.at(n0-1)),*ybegin(&vy0.at(0)),*yend(&vy0.at(n0-1)),*dybegin(&vy0.at(0,1)),*dyend(&vy0.at(n0-1,1)),&begin1(*(begin+1)),&ybegin1(*(ybegin+1)),&dybegin1(*(dybegin+1)),&end1(*(end-1)),&yend1(*(yend-1)),&dyend1(*(dyend-1));
    do
    {
        const double& x(x1.at(--u));
        double*i(begin),*j(end);
        if(x<*i)
        {
            double d0(x-*begin),d1(x-begin1),t1(d0/(begin1-*begin)),t0(d1/(*begin-begin1)),t12(t1),t02(t0);
            t12*=t1;
            t02*=t0;
            y1.at(u)=ybegin1*(1+2*t0)*t12+*ybegin*(1+2*t1)*t02+*dybegin*d0*t02+dybegin1*d1*t12;
        }
        else if(x>*j)
        {
            double d0(x-*end),d1(x-end1),t0(d1/(*end-end1)),t1(d0/(end1-*end)),t02(t0),t12(t1);
            t02*=t0;
            t12*=t1;
            y1.at(u)=*yend*(1+2*t1)*t02+yend1*(1+2*t0)*t12+*dyend*d0*t02+dyend1*d1*t12;
        }
        else
        {
            do//二分查找
            {
                double*m(i+(j-i)/2);
                if(*m>x)
                    --(j=m);
                else
                    ++(i=m);
            }while(i<=j);
            unsigned k(i-begin),l(j-begin);
            double d0(x-*j),d1(x-*i),t0(d1/(*j-*i)),t1(d0/(*i-*j)),t02(t0),t12(t1);
            t02*=t0,t12*=t1;
            y1.at(u)=ybegin[l]*(1+2*t1)*t02+ybegin[k]*(1+2*t0)*t12+dybegin[l]*d0*t02+dybegin[k]*d1*t12;
        }
    }while(u);
}
/*
 * 三弯矩算法
 * x0:观测点
 * y0:观测数据
 * x1:插值点
 * y1:预测值
 * B :边界条件
 */
void Three_moment_algorithm(const vec& x0,const mat& y0,const vec& x1,vec& y1,const boundary_conditions& B)
{
    if(x1.empty())
    {
        y1.clear();
        return;
    }
    if(x0.empty()||y0.empty())throw "观测集为空！";
    if(y0.n_cols!=1)throw "三次样条插值不需要导数信息！";
    if(x0.n_elem!=y0.n_elem)throw "观测点与观测数据向量大小不一致！";
    if(x0.n_elem<4u)throw "观测点数量小于4！";
    if(_MYFUNCTION::check_unique(x0))throw "观测集包含重复元素！";
    vec vx0(x0);
    uvec index(sort_index(x0));
    vx0.elem(index);
    mat vy0(y0.rows(index));
    const unsigned n0(x0.n_elem),n1(x1.n_elem),n(n0-1);
    double *lambda(new double[n]),*mu(new double[n]),*M(new double[n0]),*d(new double[n0]),*l_p(lambda+n),*m_p(mu+n),*M_p(M),*d_p(d+n),*diff(new double[n]+n),*diffp(diff),*h(new double[n]+n),*h_p(h);
    double *xb(&vx0.at(0)),*yb(&vy0.at(0)),*xe(xb+n0),*ye(yb+n0),*xp(xe),*yp(ye);
    while(--yp,--xp!=xb)
        *--diff=(*yp-*(yp-1))/(*--h=*xp-*(xp-1));
    if(B.t&'\1')
    {
        double h_0n(*h+*--h_p);
        *--m_p=1-(*--l_p=*h/h_0n);
        *--d_p=6*(*diff-*--diffp)/h_0n;
        do
        {
            const double &th(*h_p),&tdiff(*diffp);
            double t(*--h_p+th);
            *--m_p=1-(*--l_p=th/t);
            *--d_p=6*(tdiff-*--diffp)/t;
        }while(m_p!=mu);
        if(!*l_p)
        {
            delete[]lambda;
            delete[]mu;
            delete[]M;
            delete[]d;
            delete[]diff;
            delete[]h;
            throw "分母为零！";
        }
        if(n==3u)
        {//直接使用Cramer法则求解
            double &mu0(*mu),&mu1(*++m_p),&mu2(*++m_p),&lambda0(*lambda),&lambda1(*++l_p),&lambda2(*++l_p),&d0(*d),&d1(*++d_p),&d2(*++d_p),D(8+mu0*mu1*mu2+lambda0*lambda1*lambda2-2*lambda2*mu0-mu1*lambda0*2-2*mu2*lambda1);
            if(!D)
            {
                delete[]lambda;
                delete[]mu;
                delete[]M;
                delete[]d;
                delete[]diff;
                delete[]h;
                throw "分母为零！";
            }
            *M=*(M_p+=3)=(4*d2+mu1*mu2*d0+lambda2*lambda0*d1-2*lambda2*d0-mu1*lambda0*d2-2*mu2*d1)/D;
            *--M_p=(4*d1+mu1*d2*mu0+mu2*d0*lambda1-mu2*d1*mu0-mu1*d0*2-2*d2*lambda1)/D;
            *--M_p=(4*d0+d1*mu2*mu0+d2*lambda0*lambda1-d2*2*mu0-d1*lambda0*2-d0*mu2*lambda1)/D;
        }
        else
        {
            //Sherman-Morrison方法求解循环三对角方程组
            //接下来diff没用了, 分别将diff和t看作向量y和q
            double *t(new double[n]),*t_p(t),&t_2(*++t_p=2/ *l_p);
            *t=-2*t_2;
            *++d_p-=2*(*++diffp=*d/ *l_p++);
            *diff=*++m_p;
            *++d_p-=*++m_p**diffp;
            *++t_p=-*m_p*t_2;
            *++diffp=2;
            unsigned i(n),k(i-=3);
            while(--i)
            {
                if(!*diffp)
                {
                    delete[]lambda;
                    delete[]mu;
                    delete[]M;
                    delete[]d;
                    delete[]diff;
                    delete[]h;
                    delete[]t;
                    throw "分母为零！";
                }
                double m(*++m_p/ *diffp),&d_k(*d_p);
                *++diffp=2-m**++l_p;
                *++d_p-=m*d_k;
                *++t_p=-2*m;
            }
            //注意到x_k的表达式只与d_k,b_k有关, 与d_{k+1},b_{k+1}无关
            if(!*diffp)
            {
                delete[]lambda;
                delete[]mu;
                delete[]M;
                delete[]d;
                delete[]diff;
                delete[]h;
                delete[]t;
                throw "分母为零！";
            }
            double m(*++m_p/ *diffp);
            double &d_k(*d_p),&l_k(*++l_p),vn(*mu/2),Ann(2-vn**++l_p);
            if(!(Ann-=m*l_k))
            {
                delete[]lambda;
                delete[]mu;
                delete[]M;
                delete[]d;
                delete[]diff;
                delete[]h;
                delete[]t;
                throw "分母为零！";
            }
            double &xk(*++diffp=(*++d_p-=m*d_k)/Ann),&bk(*--diffp),&qdk(*t_p),&qxk(*(++t_p)--=(*l_p-m*qdk)/Ann);
            (qdk-=*l_p*qxk)/=bk;
            bk=(*--d_p-*l_p*xk)/bk;
            while(--k)
            {
                double &xk(*diffp),&bk(*--diffp),&qxk(*t_p),&qdk(*--t_p);
                (qdk-=*--l_p*qxk)/=bk;
                bk=(*--d_p-*l_p*xk)/bk;
            }
            (*t-=*--l_p**t_p)/=*diff;
            *diff=(*--d_p-*l_p**diffp)/ *diff;
            --(i=n);
            double dot_c((*diff+vn**(diffp=diff+i))/(1+*t+vn**(t_p=t+i)));
            *M=*(M_p+=n)=*(diffp=diff+i)-dot_c**(t_p=t+i);
            do *--M_p=*--diffp-dot_c**--t_p;while(--i);
            delete[]t;
        }
    }
    else
    {
        --h_p,--diffp;
        if(B.t&'\4')
        {
            *d_p=2*B.y2;
            *--m_p=0;
        }
        else
        {
            *d_p=6*(B.y2-*diffp)/ *h_p;
            *--m_p=1;
        }
        if(B.t&'\2')
        {
            *d=2*B.y1;
            *lambda=0;
        }
        else
        {
            *d=6*(B.y1-*diff)/ *h;
            *lambda=1;
        }
        do
        {
            double &hi(*h_p),&diffi(*diffp),t(*--h_p+hi);
            *--m_p=1-(*--l_p=hi/t);
            *--d_p=6*(diffi-*--diffp)/t;
        }while(mu!=m_p);
        //追赶法求解
        unsigned i(n);
        double m(*m_p/(*M=2));
        *++M_p=2-m**--l_p;
        *d_p-=m**d;
        while(--i)
        {
            double m(*++m_p/ *M_p),&td(*d_p);
            *++M_p=2-m**++l_p;
            *++d_p-=m*td;
        }
        *M_p=*d_p/ *M_p;
        ++l_p;
        do
        {
            double &preM(*M_p--);
            *M_p=(*--d_p-*--l_p*preM)/ *M_p;
        }while(M_p!=M);
    }
    delete[]lambda;
    delete[]mu;
    delete[]d;
    delete[]diff;
    y1.set_size(n1);
    unsigned u(n1);
    double &M0(*M),&M1(M[1]),&X0(*xb),&X1(xb[1]),&Y0(*yb),&Y1(yb[1]),&XE(*(xe-1)),&XE1(*(xe-2)),&YE(*(ye-1)),&YE1(*(ye-2)),&ME(M[n]),&ME1(M[n-1]);
    do
    {
        const double& x(x1.at(--u));
        if(*xb>=x)
        {
            double Dx1(X1-x),Dx13(Dx1),Dx0(x-X0),Dx03(Dx0),&hi(*h),hi2(hi*hi/6);
            (y1.at(u)=(M0*((Dx13*=Dx1)*=Dx1)+M1*((Dx03*=Dx0)*=Dx0))/6+(Y0-M0*hi2)*Dx1+(Y1-M1*hi2)*Dx0)/=hi;
            continue;
        }
        if(XE<x)
        {
            double Dx1(XE-x),Dx13(Dx1),Dx0(x-XE1),Dx03(Dx0),&hi(h[n]),hi2(hi*hi/6);
            (y1.at(u)=(ME1*((Dx13*=Dx1)*=Dx1)+ME*((Dx03*=Dx0)*=Dx0))/6+(YE1-ME1*hi2)*Dx1+(YE-ME*hi2)*Dx0)/=hi;
            continue;
        }
        double *j(lower_bound(xb,xe,x)),*i(j),Dx1(*j-x),Dx13(Dx1),Dx0(x-*--i),Dx03(Dx0);
        unsigned k(i-xb),l(k);
        double &hi(h[k]),hi2(hi*hi/6),&M0(M[k]),&M1(M[++k]);
        (y1.at(u)=(M0*((Dx13*=Dx1)*=Dx1)+M1*((Dx03*=Dx0)*=Dx0))/6+(yb[l]-M0*hi2)*Dx1+(yb[k]-M1*hi2)*Dx0)/=hi;
    }while(u);
    delete[]M;
    delete[]h;
}
/*
 * 三转角算法
 * x0:观测点
 * y0:观测数据
 * x1:插值点
 * y1:预测值
 * B :边界条件
 */
void Three_angle_algorithm(const vec& x0,const mat& y0,const vec& x1,vec& y1,const boundary_conditions& B)//还有bug
/*错误信息:
Thread 1 received signal SIGTRAP, Trace/breakpoint trap.
0x00007ff94acef633 in ntdll!RtlIsZeroMemory () from C:\Windows\SYSTEM32\ntdll.dll
Execute debugger commands using "-exec <command>", for example "-exec info registers" will list registers in use (when GDB is the debugger)*/
//如果删去所有delete语句可正常运行, 但运行结果也不正确
{
    if(x1.empty())
    {
        y1.clear();
        return;
    }
    if(x0.empty()||y0.empty())throw "观测集为空！";
    if(y0.n_cols!=1)throw "三次样条插值不需要导数信息！";
    if(x0.n_elem!=y0.n_elem)throw "观测点与观测数据向量大小不一致！";
    if(x0.n_elem<4u)throw "观测点数量小于4！";
    if(_MYFUNCTION::check_unique(x0))throw "观测集包含重复元素！";
    vec vx0(x0);
    uvec index(sort_index(x0));
    vx0.elem(index);
    mat vy0(y0.rows(index));
    const unsigned n0(x0.n_elem),n1(x1.n_elem),n(n0-1);
    double *lambda(new double[n]),*mu(new double[n]),*M(new double[n0]),*d(new double[n0]),*l_p(lambda+n),*m_p(mu+n),*M_p(M),*d_p(d+n),*diff(new double[n]+n),*diffp(diff),*h(new double[n]+n),*h_p(h);
    double *xb(&vx0.at(0)),*yb(&vy0.at(0)),*xe(xb+n0),*ye(yb+n0),*xp(xe),*yp(ye);
    while(--yp,--xp!=xb)
        *--diff=(*yp-*(yp-1))/(*--h=*xp-*(xp-1));
    if(B.t&'\1')
    {
        //为了适应后面的方程组求解, 这里把lambda和mu互换一下
        *--l_p=1-(*--m_p=*h/(*h+*--h_p));
        *--d_p=3*(*l_p**diff+*m_p**--diffp);
        do
        {
            const double &th(*h_p),&tdiff(*diffp);
            *--l_p=1-(*--m_p=th/(*--h_p+th));
            *--d_p=3*(*m_p**--diffp+*l_p*tdiff);
        }while(m_p!=mu);
        if(!*l_p)
        {
            delete[]lambda;
            delete[]mu;
            delete[]M;
            delete[]d;
            delete[]diff;
            delete[]h;
            throw "分母为零！";
        }
        if(n==3u)
        {//直接使用Cramer法则求解
            double &mu0(*mu),&mu1(*++m_p),&mu2(*++m_p),&lambda0(*lambda),&lambda1(*++l_p),&lambda2(*++l_p),&d0(*d),&d1(*++d_p),&d2(*++d_p),D(8+mu0*mu1*mu2+lambda0*lambda1*lambda2-2*lambda2*mu0-mu1*lambda0*2-2*mu2*lambda1);
            if(!D)
            {
                delete[]lambda;
                delete[]mu;
                delete[]M;
                delete[]d;
                delete[]diff;
                delete[]h;
                throw "分母为零！";
            }
            *M=*(M_p+=3)=(4*d2+mu1*mu2*d0+lambda2*lambda0*d1-2*lambda2*d0-mu1*lambda0*d2-2*mu2*d1)/D;
            *--M_p=(4*d1+mu1*d2*mu0+mu2*d0*lambda1-mu2*d1*mu0-mu1*d0*2-2*d2*lambda1)/D;
            *--M_p=(4*d0+d1*mu2*mu0+d2*lambda0*lambda1-d2*2*mu0-d1*lambda0*2-d0*mu2*lambda1)/D;
        }
        else
        {
            //Sherman-Morrison方法求解循环三对角方程组
            //接下来diff没用了, 分别将diff和t看作向量y和q
            double *t(new double[n]),*t_p(t),&t_2(*++t_p=2/ *l_p);
            *t=-2*t_2;
            *++d_p-=2*(*++diffp=*d/ *l_p++);
            *diff=*++m_p;
            *++d_p-=*++m_p**diffp;
            *++t_p=-*m_p*t_2;
            *++diffp=2;
            unsigned i(n),k(i-=3);
            while(--i)
            {
                if(!*diffp)
                {
                    delete[]lambda;
                    delete[]mu;
                    delete[]M;
                    delete[]d;
                    delete[]diff;
                    delete[]h;
                    delete[]t;
                    throw "分母为零！";
                }
                double m(*++m_p/ *diffp),&d_k(*d_p);
                *++diffp=2-m**++l_p;
                *++d_p-=m*d_k;
                *++t_p=-2*m;
            }
            //注意到x_k的表达式只与d_k,b_k有关, 与d_{k+1},b_{k+1}无关
            if(!*diffp)
            {
                delete[]lambda;
                delete[]mu;
                delete[]M;
                delete[]d;
                delete[]diff;
                delete[]h;
                delete[]t;
                throw "分母为零！";
            }
            double m(*++m_p/ *diffp);
            double &d_k(*d_p),&l_k(*++l_p),vn(*mu/2),Ann(2-vn**++l_p);
            if(!(Ann-=m*l_k))
            {
                delete[]lambda;
                delete[]mu;
                delete[]M;
                delete[]d;
                delete[]diff;
                delete[]h;
                delete[]t;
                throw "分母为零！";
            }
            double &xk(*++diffp=(*++d_p-=m*d_k)/Ann),&bk(*--diffp),&qdk(*t_p),&qxk(*(++t_p)--=(*l_p-m*qdk)/Ann);
            (qdk-=*l_p*qxk)/=bk;
            bk=(*--d_p-*l_p*xk)/bk;
            while(--k)
            {
                double &xk(*diffp),&bk(*--diffp),&qxk(*t_p),&qdk(*--t_p);
                (qdk-=*--l_p*qxk)/=bk;
                bk=(*--d_p-*l_p*xk)/bk;
            }
            (*t-=*--l_p**t_p)/=*diff;
            *diff=(*--d_p-*l_p**diffp)/ *diff;
            --(i=n);
            double dot_c((*diff+vn**(diffp=diff+i))/(1+*t+vn**(t_p=t+i)));
            *M=*(M_p+=n)=*(diffp=diff+i)-dot_c**(t_p=t+i);
            do *--M_p=*--diffp-dot_c**--t_p;while(--i);
            delete[]t;
        }
    }
    else
    {
        --h_p,--diffp;
        if(B.t&'\4')
        {
            *d_p=3**diffp+*h_p*B.y2/2;
            *--l_p=1;
        }
        else
        {
            *d_p=2*B.y2;
            *--l_p=0;
        }
        if(B.t&'\2')
        {
            *d=3**diff-*h*B.y1/2;
            *mu=1;
        }
        else
        {
            *d=2*B.y1;
            *mu=0;
        }
        do
        {
            double &hi(*h_p),&diffi(*diffp);
            *--l_p=1-(*--m_p=hi/(*--h_p+hi));
            *--d_p=3*(*m_p**--diffp+*l_p*diffi);
        }while(mu!=m_p);
        unsigned i(n);
        double m(*m_p/(*M=2));
        *++M_p=2-m**--l_p;
        *d_p-=m**d;
        while(--i)
        {
            double m(*++m_p/ *M_p),&td(*d_p);
            *++M_p=2-m**++l_p;
            *++d_p-=m*td;
        }
        *M_p=*d_p/ *M_p;
        ++l_p;
        do
        {
            double &preM(*M_p--);
            *M_p=(*--d_p-*--l_p*preM)/ *M_p;
        }while(M_p!=M);
    }
    delete[]lambda;
    delete[]mu;
    delete[]d;
    delete[]diff;
    y1.set_size(n1);
    unsigned u(n1);
    double &M0(*M),&M1(M[1]),&X0(*xb),&X1(xb[1]),&Y0(*yb),&Y1(yb[1]),&XE(*(xe-1)),&XE1(*(xe-2)),&YE(*(ye-1)),&YE1(*(ye-2)),&ME(M[n]),&ME1(M[n-1]);
    do
    {
        const double& x(x1.at(--u));
        if(*xb>=x)
        {
            double Dx0(x-X0),Dx02(Dx0*Dx0),Dx1(x-X1),Dx12(Dx1*Dx1),&hi(*h),hi2(hi*hi);
            y1.at(u)=((Y0*(hi+2*Dx0)*Dx12+Y1*(hi-2*Dx1)*Dx02)/hi+Dx0*Dx12*M0+Dx1*Dx02*M1)/hi2;
            continue;
        }
        if(XE<x)
        {
            double Dx0(x-XE1),Dx02(Dx0*Dx0),Dx1(x-XE),Dx12(Dx1*Dx1),&hi(h[n]),hi2(hi*hi);
            y1.at(u)=((YE1*(hi+2*Dx0)*Dx12+YE*(hi-2*Dx1)*Dx02)/hi+Dx0*Dx12*ME1+Dx1*Dx02*ME)/hi2;
            continue;
        }
        double *j(lower_bound(xb,xe,x)),*i(j),Dx1(*j-x),Dx12(Dx1*Dx1),Dx0(x-*--i),Dx02(Dx0*Dx0);
        unsigned k(i-xb),l(k);
        double &hi(h[k]),hi2(hi*hi),&M0(M[k]),&M1(M[++k]);
        y1.at(u)=((yb[l]*(hi+2*Dx0)*Dx12+yb[k]*(hi-2*Dx1)*Dx02)/hi+Dx0*Dx12*M0+Dx1*Dx02*M1)/hi2;
    }while(u);
    delete[]h;
    delete[]M;
}
}
E4t3::E4t3(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::E4t3)
{
    ui->setupUi(this);
#ifdef DEBUG43
    auto f=[](double& x){x=1./(x*x*25.+1.);};
    auto f1=[](double& x){
        double t(25*x*x+1);
        x=-50*x/(t*=t);
    };//f的导数
    vec x0(linspace(-1,1,11)),x1(linspace(-1,1,1000)),y1,ans(x1);
    mat y00(x0),y01(x0);
    mat y0(x0);
    y0.for_each(f);
    y00.for_each(f);
    y01.for_each(f1);
    ans.for_each(f);
    Ex4::boundary_conditions B1,B2,B3('\6');
    vec y2,y3;
    try {
//        Ex4::piecewise_linear_interpolation(x0,y0,x1,y1);
        Ex4::piecewise_Hermite_interpolation(x0,join_rows(y00,y01),x1,y1);
//        Ex4::Three_moment_algorithm(x0,y0,x1,y1,B1);
//        Ex4::Three_moment_algorithm(x0,y0,x1,y2,Ex4::Lagrange_boundary(x0,y0,B2));
//        Ex4::Three_moment_algorithm(x0,y0,x1,y3,B3);
//        Ex4::Three_angle_algorithm(x0,y0,x1,y1,B1);
//        Ex4::Three_angle_algorithm(x0,y0,x1,y2,Ex4::Lagrange_boundary(x0,y0,B2));
//        Ex4::Three_angle_algorithm(x0,y0,x1,y3,B3);
        QFile file("D:/1.csv");
        file.open(QIODevice::WriteOnly);
        for(unsigned i(0);i!=1000;++i)
            file.write((to_string(x1.at(i))+','+to_string(y1.at(i))+','+to_string(ans.at(i))+'\n').c_str());
//        for(unsigned i(0);i!=1001;++i)
//            file.write((to_string(x1.at(i))+','+to_string(y1.at(i))+','+to_string(y2.at(i))+','+to_string(y3.at(i))+','+to_string(ans.at(i))+'\n').c_str());
    } catch (const char* s) {
        qDebug()<<s;
    }
#else
#ifdef DEBUG43E1
    vec x0(linspace(-5,5,11)),y0(x0),x1(linspace(-5,5,1000)),y1,y2,y3,y4,y5,ans(x1);
    auto f=[](double& x){
        x=1/(1+x*x);
    };
    y0.for_each(f);
    ans.for_each(f);
    Ex4::Newton(x0,y0,x1,y1);
    Ex4::boundary_conditions B1,B2,B3('\6');
    Ex4::piecewise_linear_interpolation(x0,y0,x1,y2);
    Ex4::Three_moment_algorithm(x0,y0,x1,y3,B1);
    Ex4::Three_moment_algorithm(x0,y0,x1,y4,Ex4::Lagrange_boundary(x0,y0,B2));
    Ex4::Three_moment_algorithm(x0,y0,x1,y5,B3);
    FILE* file;
    fopen_s(&file,"D:/1.csv","w");
    fprintf(file,"x,f(x),多项式插值,分段线性插值,周期边界条件,Lagrange型边界条件,自然边界条件");
    for(unsigned i(0);i!=1001;++i)
        fprintf(file,"\n%.14f,%.14f,%.14f,%.14f,%.14f,%.14f,%.14f",x1.at(i),ans.at(i),y1.at(i),y2.at(i),y3.at(i),y4.at(i),y5.at(i));
#else
    ui->widget->setTitle("观测点");
    ui->widget_2->setTitle("观测数据");
    ui->widget_3->setTitle("插值点");
    connect(ui->pushButton,&QPushButton::clicked,[=](){
        if(ui->widget->x.empty())
        {
            QMessageBox::warning(this,"分段低次插值","插值点为空！");
            return;
        }
        vec t;
        try {
            FILE* file;
            if(fopen_s(&file,QFileDialog::getSaveFileName(this,"请选择csv文件保存路径",QStandardPaths::writableLocation(QStandardPaths::DesktopLocation),"CSV Files (.csv);;All Files ()").toLatin1().data(),"w"))
                throw "文件保存失败！";
            if(ui->radioButton->isChecked())
                Ex4::piecewise_linear_interpolation(ui->widget->x,ui->widget_2->x,ui->widget_3->x,t);
            else if(ui->radioButton_2->isChecked())
                Ex4::piecewise_Hermite_interpolation(ui->widget->x,ui->widget_2->x,ui->widget_3->x,t);
            else
            {
                const vec& x0(ui->widget->x);
                const mat& y0(ui->widget_2->x);
                if(x0.empty()||y0.empty())throw "观测集为空！";
                if(y0.n_cols!=1)throw "三次样条插值不需要导数信息！";
                if(x0.n_elem!=y0.n_elem)throw "观测点与观测数据向量大小不一致！";
                if(x0.n_elem<4u)throw "观测点数量小于4！";
                Ex4::boundary_conditions B;
                E4t3_1* D=new E4t3_1(this,B,x0,y0);
                D->setModal(true);
                D->setAttribute(Qt::WA_DeleteOnClose);
                D->exec();
                Ex4::Three_moment_algorithm(x0,y0,ui->widget_3->x,t,B);
            }
            fprintf(file,"x,y");
            vec& x1(ui->widget_3->x);
            for(unsigned i(0);i!=x1.n_elem;++i)
                fprintf(file,"\n%.14lf,%.14lf",x1.at(i),t.at(i));
            fclose(file);
        } catch (const char* s) {
            QMessageBox::critical(this,"分段低次插值",s);
            return;
        }
        QMessageBox::information(this,"分段低次插值","导出成功！");
    });
#endif
#endif
}

E4t3::~E4t3()
{
    delete ui;
}
