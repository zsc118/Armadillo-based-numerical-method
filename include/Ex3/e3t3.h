#ifndef E3T3_H
#define E3T3_H

#include <QWidget>

namespace Ui {
class E3t3;
}

class E3t3 : public QWidget
{
    Q_OBJECT

public:
    explicit E3t3(QWidget *parent = nullptr);
    ~E3t3();

private:
    Ui::E3t3 *ui;
};

#endif // E3T3_H
