#ifndef E3T4_H
#define E3T4_H

#include <QWidget>

namespace Ui {
class E3t4;
}

class E3t4 : public QWidget
{
    Q_OBJECT

public:
    explicit E3t4(QWidget *parent = nullptr);
    ~E3t4();

private:
    Ui::E3t4 *ui;
};

#endif // E3T4_H
