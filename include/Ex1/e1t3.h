#ifndef E1T3_H
#define E1T3_H

#include <QWidget>

namespace Ui {
class E1t3;
}

class E1t3 : public QWidget
{
    Q_OBJECT

public:
    explicit E1t3(QWidget *parent = nullptr);
    ~E1t3();

private:
    Ui::E1t3 *ui;
};

#endif // E1T3_H
