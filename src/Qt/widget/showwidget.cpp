#include "showwidget.h"
#include "ui_showwidget.h"
#include <QString>
#include <QPushButton>
#include <QMessageBox>
#include <math.h>
using arma::mat;
extern double calStr(std::string s);
showWidget::showWidget(QWidget *parent,mat& x) :
    QDialog(parent),ui(new Ui::showWidget)
{
    ui->setupUi(this);
    ui->tableWidget->setRowCount(x.n_rows);
    ui->tableWidget->setColumnCount(x.n_cols);
    for(unsigned i(0);i<x.n_rows;++i)
        for(unsigned j(0);j<x.n_cols;++j)
            ui->tableWidget->setItem(i,j,new QTableWidgetItem(QString::number(x.at(i,j))));
    ui->buttonBox->button(QDialogButtonBox::Ok)->setText("确定");
    ui->buttonBox->button(QDialogButtonBox::Cancel)->setText("取消");
    connect(ui->buttonBox->button(QDialogButtonBox::Cancel),&QPushButton::clicked,[=](){
        delete this;
    });
    connect(ui->buttonBox->button(QDialogButtonBox::Ok),&QPushButton::clicked,[=,&x](){
        for(unsigned i(0);i<x.n_rows;++i)
            for(unsigned j(0);j<x.n_cols;++j)
                if(isnan(x.at(i,j)=calStr(ui->tableWidget->item(i,j)->text().toStdString())))
                {
                    QMessageBox::critical(parent,"错误","输入错误！");
                    return;
                }
        delete this;
    });
}

showWidget::~showWidget()
{
    delete ui;
}
