#include "e2t7.h"
#include "ui_e2t7.h"

E2t7::E2t7(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::E2t7)
{
    ui->setupUi(this);
}

E2t7::~E2t7()
{
    delete ui;
}
