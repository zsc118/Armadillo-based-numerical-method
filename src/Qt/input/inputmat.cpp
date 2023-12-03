#include "inputmat.h"
#include "ui_inputmat.h"
#include <QPushButton>
#include <QMessageBox>
#include <math.h>
extern double calStr(std::string s);
using namespace arma;
inputMat::inputMat(QWidget *parent, unsigned m, unsigned n, mat& x) :
    QDialog(parent), ui(new Ui::inputMat)
{
    ui->setupUi(this);
    ui->tableWidget->setRowCount(m);
    ui->tableWidget->setColumnCount(n);
    ui->buttonBox->button(QDialogButtonBox::Ok)->setText("确定");
    ui->buttonBox->button(QDialogButtonBox::Cancel)->setText("取消");
    connect(ui->buttonBox->button(QDialogButtonBox::Ok),&QPushButton::clicked,[=,&x](){
        x.zeros(m,n);
        for(unsigned i(0);i<m;++i)
            for(unsigned j(0);j<n;++j)
                if(isnan(x.at(i,j)=calStr(ui->tableWidget->item(i,j)->text().toStdString())))
                {
                    QMessageBox::critical(parent,"错误","输入错误！");
                    return;
                }
        delete this;
    });
    connect(ui->buttonBox->button(QDialogButtonBox::Cancel),&QPushButton::clicked,[=](){
        delete this;
    });
}

inputMat::~inputMat()
{
    delete ui;
}
