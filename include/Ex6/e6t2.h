#ifndef E6T2_H
#define E6T2_H

#include <QWidget>

namespace Ui {
class E6t2;
}

class E6t2 : public QWidget
{
    Q_OBJECT

public:
    explicit E6t2(QWidget *parent = nullptr);
    ~E6t2();

private:
    Ui::E6t2 *ui;
};

#endif // E6T2_H
