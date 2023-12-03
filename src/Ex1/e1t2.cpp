#include "e1t2.h"
#include "ui_e1t2.h"
#include <QInputDialog>
#include <functional>
#include <math.h>
#include <QMessageBox>
#include <QFileDialog>
#include <QStandardPaths>
#include <QFile>
#include <QBitArray>
#include <vector>
#include <string>
extern double calStr(std::string s);
namespace CAL {
extern std::vector<std::string> func;
}
namespace _MYFUNCTION {
// double向量转string向量
std::vector<std::string> vec_to_string(const std::vector<double>& x)
{
    std::vector<std::string> y(x.size());
    auto i=x.cbegin();
    auto k=y.begin();
    while(i<x.cend())*k++=*i++;
    return y;
}
/*
 * 替换子串
 * str   : 要替换的字符串
 * oldStr: 旧子串
 * newStr: 新子串
 */
std::string &substring_replace(std::string &str,const std::string& oldStr,const std::string& newStr)
{
    if(oldStr.empty())
    {
        str.clear();
        return str;
    }
    size_t pos(0);
    while ((pos = str.find(oldStr, pos)) != std::string::npos)
    {
        str.replace(pos, oldStr.length(), newStr);
        pos += newStr.length();
    }
    return str;
}
/*
 * 替换子串(1个)
 * str   : 要替换的字符串
 * oldStr: 旧子串
 * newStr: 新子串
 * except: 排除的子串
 */
std::string &substring_replace(std::string &str,const std::string& oldStr,const std::string& newStr,const std::vector<std::string>& except)
{
    if(oldStr.empty())
    {
        str.clear();
        return str;
    }
    QBitArray a(str.length());
    for(const auto& i:except)
    {
        size_t pos(0);
        while ((pos = str.find(i, pos)) != std::string::npos)
        {
            a.fill(true,pos,pos+i.length());
            pos += i.length();
        }
    }
    size_t pos(0);
F:
    while ((pos = str.find(oldStr, pos)) != std::string::npos)
    {
        size_t i(oldStr.length()),p(pos);
        do
            if(a.at(p++))
            {
                pos+=oldStr.length();
                goto F;
            }
        while(--i);
        QBitArray t(a.size()+newStr.size()-oldStr.size());
        do t.setBit(i,a.at(i));while(++i<pos);
        p=pos+oldStr.length();
        str.replace(pos, oldStr.length(), newStr);
        if(i!=(pos += newStr.length()))do t.setBit(i,false);while(++i!=pos);
        while(p!=a.size())t.setBit(i++,a.at(p++));
        a=t;
    }
    return str;
}
/*
 * 替换子串(多个)
 * str   : 要替换的字符串
 * oldStr: 旧子串
 * newStr: 新子串
 * except: 排除的子串
 */
std::string &substring_replace(std::string &str,const std::vector<std::string>& oldStr,const std::vector<std::string>& newStr,const std::vector<std::string>& except)
{
    if(oldStr.empty())
    {
        str.clear();
        return str;
    }
//    if(oldStr.size()!=newStr.size())
//        if(oldStr.size()>newStr.size())
//            oldStr.erase(oldStr.begin()+newStr.size(),oldStr.end());
//        else
//            newStr.erase(newStr.begin()+oldStr.size(),newStr.end());
    size_t n(oldStr.size()>newStr.size()?newStr.size():oldStr.size());
    QBitArray a(str.length());
    for(const auto& i:except)
    {
        size_t pos(0);
        while ((pos = str.find(i, pos)) != std::string::npos)
        {
            a.fill(true,pos,pos+i.length());
            pos += i.length();
        }
    }
    size_t pos(0);
    for(size_t k(0);k<n;++k)
    {
        const std::string& pre=oldStr[k],&next=newStr[k];
F:
        while ((pos = str.find(pre, pos)) != std::string::npos)
        {
            size_t i(pre.length()),p(pos);
            do
                if(a.at(p++))
                {
                    pos+=pre.length();
                    goto F;
                }
            while(--i);
            QBitArray t(a.size()+next.size()-pre.size());
            do t.setBit(i,a.at(i));while(++i<pos);
            p=pos+pre.length();
            str.replace(pos, pre.length(), next);
            if(i!=(pos += next.length()))do t.setBit(i,false);while(++i!=pos);
            while(p!=a.size())t.setBit(i++,a.at(p++));
            a=t;
        }
    }
    return str;
}
// 判断是否为英文字母
inline bool is_not_alpha(char ch)noexcept
{
    return ch<'A'||ch>'z'||ch>'Z'&&ch<'a';
}
// 判断是否为数字
inline bool is_not_number(char ch)noexcept
{
    return ch<'0'||ch>'9';
}
// 判断是否为未知数字母
inline bool is_not_xch(char ch)noexcept
{
    return is_not_alpha(ch)&&is_not_number(ch)&&ch!='_';
}
}
namespace Ex1 {
/*
 * 迭代法(模板)
 * x   : 迭代自变量
 * func: 迭代函数
 * k   : 迭代次数
 */
template <class T>
T &iterative_method(T& x,std::function<T&(T&)> func, unsigned k)
{
    if(k)
        do func(x);while(--k);
    return x;
}
/*
 * 迭代法(模板)
 * x   : 迭代自变量
 * func: 迭代函数
 * k   : 迭代次数
 */
template <class T>
T &iterative_method(T& x,std::function<void(T&)> func, unsigned k)
{
    if(k)
        do func(x);while(--k);
    return x;
}
/*
 * 迭代法(模板)
 * x   : 迭代自变量
 * func: 迭代函数
 * k   : 迭代次数
 */
template <class T>
T &iterative_method(T& x,std::function<T(T)> func, unsigned k)
{
    if(k)
        do x=func(x);while(--k);
    return x;
}
/*
 * 迭代法(模板)
 * x   : 迭代自变量
 * func: 迭代函数
 * exit: 终止标准
 *  true  - 继续迭代
 *  false - 输出结果
 */
template <class T>
T &iterative_method(T& x,std::function<T(T)> func, std::function<bool(T)> exit)
{
    while(exit(x))x=func(x);
    return x;
}
/*
 * 迭代法(模板)
 * x   : 迭代自变量
 * func: 迭代函数
 * exit: 终止标准
 *  true  - 继续迭代
 *  false - 输出结果
 */
template <class T>
T &iterative_method(T& x,std::function<void(T&)> func, std::function<bool(T)> exit)
{
    while(exit(x))func(x);
    return x;
}
/*
 * 迭代法(模板)
 * x   : 迭代自变量
 * func: 迭代函数
 * exit: 终止标准
 *  true  - 继续迭代
 *  false - 输出结果
 */
template <class T>
T &iterative_method(T& x,std::function<T&(T&)> func, std::function<bool(T)> exit)
{
    while(exit(x))func(x);
    return x;
}
/*
 * 迭代法(模板)
 * x   : 迭代自变量
 * func: 迭代函数
 * exit: 终止标准
 *  true  - 继续迭代
 *  false - 输出结果
 *  第1个参数: x_n
 *  第2个参数: x_{n-1}
 */
template <class T>
T &iterative_method(T& x,std::function<T&(T&)> func, std::function<bool(T,T)> exit)
{
    T t(x);
    while(exit(func(x),t))t=x;
    return x;
}
/*
 * 迭代法(模板)
 * x   : 迭代自变量
 * func: 迭代函数
 * exit: 终止标准
 *  true  - 继续迭代
 *  false - 输出结果
 *  第1个参数: x_n
 *  第2个参数: x_{n-1}
 */
template <class T>
T &iterative_method(T& x,std::function<void(T&)> func, std::function<bool(T,T)> exit)
{
    T t(x);
    func(x);
    while(exit(x,t))
    {
        t=x;
        func(x);
    }
    return x;
}
/*
 * 迭代法(模板)
 * x   : 迭代自变量
 * func: 迭代函数
 * exit: 终止标准
 *  true  - 继续迭代
 *  false - 输出结果
 *  第1个参数: x_n
 *  第2个参数: x_{n-1}
 */
template <class T>
T &iterative_method(T& x,std::function<T(T)> func, std::function<bool(T,T)> exit)
{
    T t(x);
    while(exit(x=func(x),t))t=x;
    return x;
}
/*
 * 迭代法(lambda表达式)
 * x   : 迭代向量
 * func: 迭代函数
 */
inline std::vector<double> &iterative_method(std::vector<double>& x,std::function<std::vector<double>(std::vector<double>)> func)
{
    return x=func(x);
}
/*
 * 迭代法(模板)
 * x   : 迭代自变量
 * func: 迭代函数
 * exit: 终止标准
 *  true  - 继续迭代
 *  false - 输出结果
 * out :输出
 */
template <class T>
T &iterative_method(T& x,std::function<void(T&)> func,std::function<bool(const T&)> exit,std::function<void(const T&)> out)
{
    out(x);
    while(exit(x))
    {
        func(x);
        out(x);
    }
    return x;
}
/*
 * 迭代法(字符串向量)
 * current_vec: 迭代向量
 * func       : 迭代函数
 * x          : 未知数
 *
 * 返回(bool):
 *  true : 迭代失败
 *  false: 迭代成功
 */
bool iterative_method(std::vector<double>& current_vec,const std::vector<std::string> &f, const std::vector<std::string> &x)
{
    std::vector<double> x0(current_vec);
    unsigned n(current_vec.size());
    do
    {
        std::string toCal(f[--n]);
        _MYFUNCTION::substring_replace(toCal,x,_MYFUNCTION::vec_to_string(current_vec),CAL::func);
        if(isnan(current_vec[n]=calStr(toCal)))
        {
            current_vec=x0;
            return true;
        }
    }while(n);
    return false;
}
/*
 * 迭代法(多次迭代)
 * current_vec: 迭代向量
 * func       : 迭代函数
 * x          : 未知数
 * k          : 迭代次数
 * mid        : 迭代中间结果导出
 *
 * 返回(bool):
 *  true : 迭代失败
 *  false: 迭代成功
 */
bool iterative_method(std::vector<double>& current_vec,const std::vector<std::string> &func, const std::vector<std::string> &x,unsigned k,QString& mid)
{
    for(auto& i:x)mid.append(QString::fromStdString(i)).append(',');
    *(mid.end()-1)='\n';
    for(auto& i:current_vec)mid.append(QString::number(i,'g',14)).append(',');
    std::vector<double> x0(current_vec);
    if(k)
    {
        do
        {
            if(iterative_method(current_vec,func,x))
            {
                current_vec=x0;
                return true;
            }
            *(mid.end()-1)='\n';
            for(auto& i:current_vec)mid.append(QString::number(i,'g',14)).append(',');
        }while(--k);
    }
    mid.chop(1);
    return false;
}
}
e1t2::e1t2(QWidget *parent) : QWidget(parent),ui(new Ui::e1t2),func({"x","y","z"}),x({"x","y","z"}),current_vec(3,0.)
{
    ui->setupUi(this);
    connect(ui->pushButton,&QPushButton::clicked,[=](){
        QString s("迭代函数：");
        for(unsigned i(0);i<x.size();++i)
            s.append("\n   ").append(QString::fromStdString(x[i])).append("\'=").append(QString::fromStdString(func[i]=QInputDialog::getText(this,"设置迭代函数",QString::fromStdString(x[i]+"\'=")).toStdString()));
        ui->label->setText(s);
    });
    connect(ui->pushButton_2,&QPushButton::clicked,[=](){
        QString s1("当前向量值：   ("),s2("=(");
        double t;
        std::vector<double> x0(current_vec);
        for(unsigned i(0);i<x.size();++i)
        {
            if(isnan(t=calStr(QInputDialog::getText(this,"设置当前向量",QString::fromStdString(x[i]+'=')).toStdString())))
            {
                current_vec=x0;
                QMessageBox::critical(this,"错误","输入有误！");
                return;
            }
            s1.append(QString::fromStdString(x[i])).append(',');
            s2.append(QString::number(current_vec[i]=t)).append(',');
        }
        *(s2.end()-1)=*(s1.end()-1)=')';
        ui->label_2->setText(s1+s2);
    });
    connect(ui->pushButton_3,&QPushButton::clicked,[=](){
        int n(QInputDialog::getInt(this,"设置未知数","请输入向量维数:"));
        if(n<=0)
        {
            QMessageBox::critical(this,"错误","输入有误！");
            return;
        }
        std::vector<std::string> t;
        int i(0);
        do
        {
            std::string s(QInputDialog::getText(this,"设置未知数","请设置第"+QString::number(++i)+"个未知数：").toStdString());
            if(_MYFUNCTION::is_not_alpha(s[0]))
            {
                QMessageBox::critical(this,"错误","输入有误！");
                return;
            }
            for(auto p=s.begin()+1;p<s.end();++p)
                if(_MYFUNCTION::is_not_xch(*p))
                {
                    QMessageBox::critical(this,"错误","输入有误！");
                    return;
                }
            t.push_back(s);
        }while(i<n);
        func=x=t;
        current_vec=std::vector<double>(n,0.);
        ui->label_3->setText("向量维数：\n   "+QString::number(n));
        i=0;
        std::string s1("迭代函数："),s21("当前向量值：   ("),s22("=("),s4("未知数记号：");
        auto p=x.begin();
        do
        {
            s1+="\n   "+*p+"\'="+*p;
            s21+=*p+',';
            s22+="0,";
            s4+="\n  "+*p;
        }while(++p!=x.end());
        *(s21.end()-1)=*(s22.end()-1)=')';
        ui->label->setText(QString::fromStdString(s1));
        ui->label_2->setText(QString::fromStdString(s21+s22));
        ui->label_4->setText(QString::fromStdString(s4));
    });
    connect(ui->pushButton_4,&QPushButton::clicked,[=](){
        if(Ex1::iterative_method(current_vec,func,x))
        {
            QMessageBox::critical(this,"错误","输入有误！");
            return;
        }
        std::string s1("当前向量值：   ("),s2("=(");
        for(unsigned i(0);i<x.size();++i)
        {
            s1+=x[i]+',';
            s2+=std::to_string(current_vec[i])+',';
        }
        *(s1.end()-1)=*(s2.end()-1)=')';
        ui->label_2->setText(QString::fromStdString(s1+s2));
    });
    connect(ui->pushButton_5,&QPushButton::clicked,[=](){
        int k=QInputDialog::getInt(this,"迭代多次","请输入迭代次数：",-1);
        if(k<0)
        {
            QMessageBox::critical(this,"错误","输入有误！");
            return;
        }
        QString s0;
        QFile file(QFileDialog::getSaveFileName(this,"请选择csv文件保存路径",QStandardPaths::writableLocation(QStandardPaths::DesktopLocation),"CSV Files (.csv);;All Files ()"));
        if(Ex1::iterative_method(current_vec,func,x,k,s0))
        {
            QMessageBox::critical(this,"错误","输入有误！");
            return;
        }
        if(!file.open(QIODevice::WriteOnly))
        {
            QMessageBox::critical(this,"错误","文件保存失败！");
            return;
        }
        file.write(s0.toUtf8());
        file.close();
        std::string s1("当前向量值：   ("),s2("=(");
        for(unsigned i(0);i<x.size();++i)
        {
            s1+=x[i]+',';
            s2+=std::to_string(current_vec[i])+',';
        }
        *(s1.end()-1)=*(s2.end()-1)=')';
        ui->label_2->setText(QString::fromStdString(s1+s2));
        QMessageBox::information(this,"迭代法","导出成功！");
    });
}

e1t2::~e1t2()
{
    delete ui;
}
