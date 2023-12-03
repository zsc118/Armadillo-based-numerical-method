#include "e2t6.h"
#include "ui_e2t6.h"

E2t6::E2t6(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::E2t6)
{
    ui->setupUi(this);
}

E2t6::~E2t6()
{
    delete ui;
}
