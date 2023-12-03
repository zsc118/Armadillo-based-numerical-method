#include "ex1.h"
#include "e1t1.h"
#include "e1t2.h"
ex1::ex1(QWidget *parent) : QTabWidget(parent)
{
    E1t1*T1=new E1t1(this);
    addTab(T1,"有效数字与误差限");
    e1t2*T2=new e1t2(this);
    addTab(T2,"迭代法");
}
