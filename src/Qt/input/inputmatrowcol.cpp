#include "inputmatrowcol.h"
#include "ui_inputmatrowcol.h"
#include "inputmat.h"
#include <QPushButton>
#include <QMessageBox>
using namespace arma;
inputMatRowCol::inputMatRowCol(QWidget *parent, mat& x) :
    QDialog(parent), ui(new Ui::inputMatRowCol)
{
    ui->setupUi(this);
    ui->buttonBox->button(QDialogButtonBox::Ok)->setText("确定");
    ui->buttonBox->button(QDialogButtonBox::Cancel)->setText("取消");
    connect(ui->buttonBox->button(QDialogButtonBox::Cancel),&QPushButton::clicked,[=](){
        delete this;
    });
    connect(ui->buttonBox->button(QDialogButtonBox::Ok),&QPushButton::clicked,[=,&x](){
        unsigned m(ui->row->text().toUInt()),n(ui->col->text().toUInt());
        if(m&&n)
        {
            inputMat*D=new inputMat(parent,m,n,x);
            D->setModal(true);
            D->setAttribute(Qt::WA_DeleteOnClose);
            D->exec();
            delete this;
        }
        else
            QMessageBox::critical(this,"输入错误","矩阵行列数不能为零或非数字");
    });
}

inputMatRowCol::~inputMatRowCol()
{
    delete ui;
}
