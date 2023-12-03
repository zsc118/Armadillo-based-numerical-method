#include "matwidget.h"
#include "ui_matwidget.h"
#include "showwidget.h"
#include "inputmatrowcol.h"
#include <QPushButton>

void matWidget::setTitle(const char* s)
{
    ui->groupBox->setTitle(s);
}

matWidget::matWidget(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::matWidget)
{
    ui->setupUi(this);
    connect(ui->pushButton,&QPushButton::clicked,[=](){
        showWidget*D=new showWidget(parent,x);
        D->setModal(true);
        D->setAttribute(Qt::WA_DeleteOnClose);
        D->exec();
    });
    connect(ui->pushButton_2,&QPushButton::clicked,[=](){
        inputMatRowCol*D=new inputMatRowCol(parent,x);
        D->setModal(true);
        D->setAttribute(Qt::WA_DeleteOnClose);
        D->exec();
    });
}

matWidget::~matWidget()
{
    delete ui;
}
