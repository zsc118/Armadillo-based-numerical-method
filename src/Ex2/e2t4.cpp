#include "e2t4.h"
#include "ui_e2t4.h"

E2t4::E2t4(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::E2t4)
{
    ui->setupUi(this);
}

E2t4::~E2t4()
{
    delete ui;
}
