#ifndef E3T2_H
#define E3T2_H

#include <QWidget>

namespace Ui {
class E3t2;
}

class E3t2 : public QWidget
{
    Q_OBJECT

public:
    explicit E3t2(QWidget *parent = nullptr);
    ~E3t2();

private:
    Ui::E3t2 *ui;
};

#endif // E3T2_H
