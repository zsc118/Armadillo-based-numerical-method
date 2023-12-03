#include "e4t3_1.h"
#include "ui_e4t3_1.h"
#include <QPushButton>
#include <string>
#include <math.h>
#include <QMessageBox>
using namespace arma;
using namespace std;
extern double calStr(string);
namespace Ex4 {
boundary_conditions& Lagrange_boundary(const vec& x,const mat& y,boundary_conditions& B)noexcept
{
    B.t='\0';
    const double *xb(&x.at(0)),*xe(xb+x.n_elem),&x0(*xb),&x1(*++xb),&x2(*++xb),&x3(*++xb),*yb(&y.at(0)),*ye(yb+x.n_elem),&y0(*yb),&y1(*++yb),&y2(*++yb),&y3(*++yb),&u3(*--xe),&u2(*--xe),&u1(*--xe),&u0(*--xe),&v3(*--ye),&v2(*--ye),&v1(*--ye),&v0(*--ye);
    double x02(x0-x2),x03(x0-x3),x01(x0-x1),x12(x1-x2),x13(x1-x3),x23(x2-x3);
    B.y1=y0/x02+y0/x03+y0/x01-y1*x02*x03/x12/x13/x01+y2*x01*x03/x23/x02/x12-y3*x01*x02/x03/x13/x23;
    double u01(u0-u1),u02(u0-u2),u03(u0-u3),u12(u1-u2),u13(u1-u3),u23(u2-u3);
    B.y2=v0*u13*u23/u01/u02/u03-v3/u13-v3/u23-v1*u03*u23/u12/u13/u01+v2*u03*u13/u23/u02/u12-v3/u03;
    return B;
}
}
E4t3_1::E4t3_1(QWidget *parent,Ex4::boundary_conditions& B,const vec& x0,const mat& y0) :
    QDialog(parent),
    ui(new Ui::E4t3_1)
{
    ui->setupUi(this);
    ui->buttonBox->button(QDialogButtonBox::Cancel)->setText("取消");
    ui->buttonBox->button(QDialogButtonBox::Ok)->setText("确定");
    connect(ui->radioButton,&QRadioButton::clicked,[=](){
        ui->groupBox->setCheckable(false);
        ui->plainTextEdit->setEnabled(false);
        ui->groupBox_2->setCheckable(false);
        ui->plainTextEdit_2->setEnabled(false);
        ui->checkBox->setEnabled(false);
    });
    connect(ui->radioButton_2,&QRadioButton::clicked,[=](){
        ui->groupBox->setCheckable(true);
        ui->plainTextEdit->setEnabled(true);
        ui->groupBox_2->setCheckable(true);
        ui->plainTextEdit_2->setEnabled(true);
        ui->checkBox->setEnabled(true);
    });
    connect(ui->buttonBox->button(QDialogButtonBox::Cancel),&QPushButton::clicked,[=](){
        delete this;
    });
    connect(ui->buttonBox->button(QDialogButtonBox::Ok),&QPushButton::clicked,[=,&B](){
        if(ui->radioButton->isChecked())
        {
            B.t='\1';
            delete this;
            return;
        }
        if(ui->checkBox->isChecked())
        {
            Ex4::Lagrange_boundary(x0,y0,B);
            delete this;
            return;
        }
        if(isnan(B.y1=calStr(ui->plainTextEdit->toPlainText().toStdString())))
        {
            QMessageBox::warning(this,"边界条件","y_0输入有误！");
            return;
        }
        if(isnan(B.y2=calStr(ui->plainTextEdit_2->toPlainText().toStdString())))
        {
            QMessageBox::warning(this,"边界条件","y_n输入有误！");
            return;
        }
        B.t='\0';
        if(ui->radioButton_4->isChecked())B.t|='\2';
        if(ui->radioButton_6->isChecked())B.t|='\4';
        delete this;
    });
}

E4t3_1::~E4t3_1()
{
    delete ui;
}
