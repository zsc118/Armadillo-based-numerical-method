#include "mainwindow.h"
#include<QMenu>
#include<QMenuBar>
#include<QAction>
#include<QStatusBar>
#include<QLabel>
#include<QStackedWidget>
#include "e1t1.h"
#include "e1t2.h"
#include "e1t3.h"
#include "e1t4.h"
#include "e2t1.h"
#include "e2t2.h"
#include "e2t3.h"
#include "e2t8.h"
#include "e2t5.h"
#include "e2t6.h"
#include "e2t7.h"
#include "e3t1.h"
#include "e3t2.h"
#include "e3t3.h"
#include "e3t4.h"
#include "e4t1.h"
#include "e4t2.h"
#include "e4t3.h"
#include "e4t4.h"
#include "e5t1.h"
#include "e5t2.h"
#include "e5t3.h"
#include "e6t1.h"
#include "e6t2.h"
MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
{
    setFixedSize(700,500);
    setWindowTitle("数学软件实验");
    QMenuBar*menu=menuBar();
    setMenuBar(menu);
    QMenu*menu1=menu->addMenu("实验一");
    QAction*A11=menu1->addAction("有效数字与误差限");
    menu1->addSeparator();
    QAction*A12=menu1->addAction("迭代法");
    QAction*A13=menu1->addAction("线性随机IFS迭代法");
    QAction*A14=menu1->addAction("复平面上的迭代");
    QMenu*menu2=menu->addMenu("实验二");
    QAction*A21=menu2->addAction("图解法");
    QAction*A22=menu2->addAction("二分法与试位法");
    menu2->addSeparator();
    QAction*A23=menu2->addAction("迭代法");
    QAction*A28=menu2->addAction("迭代法(选择)");
    menu2->addSeparator();
    QAction*A25=menu2->addAction("Sturm序列法");
    QAction*A26=menu2->addAction("劈因子法");
    QAction*A27=menu2->addAction("矩阵特征值法");
    QMenu*menu3=menu->addMenu("实验三");
    QAction*A31=menu3->addAction("Gauss消去法");
    QAction*A32=menu3->addAction("LU分解");
    QAction*A34=menu3->addAction("追赶法");
    menu3->addSeparator();
    QAction*A33=menu3->addAction("迭代法");
    QMenu*menu4=menu->addMenu("实验四");
    QAction*A41=menu4->addAction("多项式插值");
    QAction*A42=menu4->addAction("Hermite插值");
    QAction*A43=menu4->addAction("分段低次插值");
    QAction*A44=menu4->addAction("高维插值");
    QMenu*menu5=menu->addMenu("实验五");
    QAction*A51=menu5->addAction("插值型求积法");
    QAction*A52=menu5->addAction("复化求积法");
    QAction*A53=menu5->addAction("Romberg求积法");
    QMenu*menu6=menu->addMenu("实验六");
    QAction*A61=menu6->addAction("常微分方程数值解法");
    QAction*A62=menu6->addAction("常微分方程组数值解法");
    QMenu*menu7=menu->addMenu("实验七");
    QMenu*menu8=menu->addMenu("实验八");
    QMenu*menu9=menu->addMenu("实验九");
    QStatusBar*stBar=statusBar();
    setStatusBar(stBar);
    QLabel*leftStatus=new QLabel("实验一：数值计算的基本概念",this);
    QLabel*rightStatus=new QLabel("有效数字与误差限",this);
    stBar->addWidget(leftStatus);
    stBar->addPermanentWidget(rightStatus);
    QStackedWidget*mainWidget=new QStackedWidget(this);
    setCentralWidget(mainWidget);
    E1t1*W11=new E1t1(this);
    int I11=mainWidget->addWidget(W11);
    e1t2*W12=new e1t2(this);
    int I12=mainWidget->addWidget(W12);
    E1t3*W13=new E1t3(this);
    int I13=mainWidget->addWidget(W13);
    E1t4*W14=new E1t4(this);
    int I14=mainWidget->addWidget(W14);
    E2t1*W21=new E2t1(this);
    int I21=mainWidget->addWidget(W21);
    E2t2*W22=new E2t2(this);
    int I22=mainWidget->addWidget(W22);
    E2t3*W23=new E2t3(this);
    int I23=mainWidget->addWidget(W23);
    E2t8*W28=new E2t8(this);
    int I28=mainWidget->addWidget(W28);
    E2t5*W25=new E2t5(this);
    int I25=mainWidget->addWidget(W25);
    E2t6*W26=new E2t6(this);
    int I26=mainWidget->addWidget(W26);
    E2t7*W27=new E2t7(this);
    int I27=mainWidget->addWidget(W27);
    E3t1*W31=new E3t1(this);
    int I31=mainWidget->addWidget(W31);
    E3t2*W32=new E3t2(this);
    int I32=mainWidget->addWidget(W32);
    E3t3*W33=new E3t3(this);
    int I33=mainWidget->addWidget(W33);
    E4t1*W41=new E4t1(this);
    int I41=mainWidget->addWidget(W41);
    E4t2*W42=new E4t2(this);
    int I42=mainWidget->addWidget(W42);
    E4t3*W43=new E4t3(this);
    int I43=mainWidget->addWidget(W43);
    E4t4*W44=new E4t4(this);
    int I44=mainWidget->addWidget(W44);
    E5t1*W51=new E5t1(this);
    int I51=mainWidget->addWidget(W51);
    E5t2*W52=new E5t2(this);
    int I52=mainWidget->addWidget(W52);
    E5t3*W53=new E5t3(this);
    int I53=mainWidget->addWidget(W53);
    E6t1*W61=new E6t1(this);
    int I61=mainWidget->addWidget(W61);
    E6t2*W62=new E6t2(this);
    int I62=mainWidget->addWidget(W62);
    E3t4*W34=new E3t4(this);
    int I34=mainWidget->addWidget(W34);
    connect(A11,&QAction::triggered,[=](){
        mainWidget->setCurrentIndex(I11);
        leftStatus->setText("实验一：数值计算的基本概念");
        rightStatus->setText("有效数字与误差限");
    });
    connect(A12,&QAction::triggered,[=](){
        mainWidget->setCurrentIndex(I12);
        leftStatus->setText("实验一：数值计算的基本概念");
        rightStatus->setText("迭代法");
    });
    connect(A13,&QAction::triggered,[=](){
        mainWidget->setCurrentIndex(I13);
        leftStatus->setText("实验一：数值计算的基本概念");
        rightStatus->setText("线性随机IFS迭代法");
    });
    connect(A14,&QAction::triggered,[=](){
        mainWidget->setCurrentIndex(I14);
        leftStatus->setText("实验一：数值计算的基本概念");
        rightStatus->setText("复平面上的迭代");
    });
    connect(A21,&QAction::triggered,[=](){
        mainWidget->setCurrentIndex(I21);
        leftStatus->setText("实验二：非线性方程求根");
        rightStatus->setText("图解法");
    });
    connect(A22,&QAction::triggered,[=](){
        mainWidget->setCurrentIndex(I22);
        leftStatus->setText("实验二：非线性方程求根");
        rightStatus->setText("二分法与试位法");
    });
    connect(A23,&QAction::triggered,[=](){
        mainWidget->setCurrentIndex(I23);
        leftStatus->setText("实验二：非线性方程求根");
        rightStatus->setText("迭代法");
    });
    connect(A28,&QAction::triggered,[=](){
        mainWidget->setCurrentIndex(I28);
        leftStatus->setText("实验二：非线性方程求根");
        rightStatus->setText("迭代法(选择)");
    });
    connect(A25,&QAction::triggered,[=](){
        mainWidget->setCurrentIndex(I25);
        leftStatus->setText("实验二：非线性方程求根");
        rightStatus->setText("Sturm序列法");
    });
    connect(A26,&QAction::triggered,[=](){
        mainWidget->setCurrentIndex(I26);
        leftStatus->setText("实验二：非线性方程求根");
        rightStatus->setText("劈因子法");
    });
    connect(A27,&QAction::triggered,[=](){
        mainWidget->setCurrentIndex(I27);
        leftStatus->setText("实验二：非线性方程求根");
        rightStatus->setText("矩阵特征值法");
    });
    connect(A31,&QAction::triggered,[=](){
        mainWidget->setCurrentIndex(I31);
        leftStatus->setText("实验三：线性方程组的数值解法");
        rightStatus->setText("Gauss消去法");
    });
    connect(A32,&QAction::triggered,[=](){
        mainWidget->setCurrentIndex(I32);
        leftStatus->setText("实验三：线性方程组的数值解法");
        rightStatus->setText("LU分解");
    });
    connect(A33,&QAction::triggered,[=](){
        mainWidget->setCurrentIndex(I33);
        leftStatus->setText("实验三：线性方程组的数值解法");
        rightStatus->setText("迭代法");
    });
    connect(A41,&QAction::triggered,[=](){
        mainWidget->setCurrentIndex(I41);
        leftStatus->setText("实验四：函数逼近的插值、拟合");
        rightStatus->setText("多项式插值");
    });
    connect(A42,&QAction::triggered,[=](){
        mainWidget->setCurrentIndex(I42);
        leftStatus->setText("实验四：函数逼近的插值、拟合");
        rightStatus->setText("Hermite插值");
    });
    connect(A43,&QAction::triggered,[=](){
        mainWidget->setCurrentIndex(I43);
        leftStatus->setText("实验四：函数逼近的插值、拟合");
        rightStatus->setText("分段低次插值");
    });
    connect(A44,&QAction::triggered,[=](){
        mainWidget->setCurrentIndex(I44);
        leftStatus->setText("实验四：函数逼近的插值、拟合");
        rightStatus->setText("高维插值");
    });
    connect(A51,&QAction::triggered,[=](){
        mainWidget->setCurrentIndex(I51);
        leftStatus->setText("实验五：数值积分");
        rightStatus->setText("插值型求积法");
    });
    connect(A52,&QAction::triggered,[=](){
        mainWidget->setCurrentIndex(I52);
        leftStatus->setText("实验五：数值积分");
        rightStatus->setText("复化求积法");
    });
    connect(A53,&QAction::triggered,[=](){
        mainWidget->setCurrentIndex(I53);
        leftStatus->setText("实验五：数值积分");
        rightStatus->setText("Romberg求积法");
    });
    connect(A61,&QAction::triggered,[=](){
        mainWidget->setCurrentIndex(I61);
        leftStatus->setText("实验六：常微分方程数值解");
        rightStatus->setText("常微分方程数值解法");
    });
    connect(A62,&QAction::triggered,[=](){
        mainWidget->setCurrentIndex(I62);
        leftStatus->setText("实验六：常微分方程数值解");
        rightStatus->setText("常微分方程组数值解法");
    });
    connect(A34,&QAction::triggered,[=](){
        mainWidget->setCurrentIndex(I34);
        leftStatus->setText("实验三：线性方程组的数值解法");
        rightStatus->setText("追赶法");
    });
}

MainWindow::~MainWindow()
{
}

