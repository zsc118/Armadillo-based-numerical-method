#include "vecwidget.h"
#include "ui_vecwidget.h"
#include "showwidget.h"
#include "inputmat.h"
#include <QPushButton>
#include <QInputDialog>
#include <QMessageBox>

void vecWidget::setTitle(const char *s)
{
    ui->groupBox->setTitle(s);
}

vecWidget::vecWidget(QWidget *parent) :
    QWidget(parent),ui(new Ui::vecWidget)
{
    ui->setupUi(this);
    connect(ui->pushButton,&QPushButton::clicked,[=](){
        showWidget*D=new showWidget(parent,x);
        D->setModal(true);
        D->setAttribute(Qt::WA_DeleteOnClose);
        D->exec();
    });
    connect(ui->pushButton_2,&QPushButton::clicked,[=](){
        int n(QInputDialog::getInt(parent,"重定义向量","请输入向量维数:"));
        if(n<=0)
        {
            QMessageBox::critical(parent,"错误","输入向量维数有误！");
            return;
        }
        inputMat*D=new inputMat(parent,n,1,x);
        D->setModal(true);
        D->setAttribute(Qt::WA_DeleteOnClose);
        D->exec();
    });
}

vecWidget::~vecWidget()
{
    delete ui;
}
