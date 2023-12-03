#ifndef E2T1_H
#define E2T1_H

#include <QWidget>

namespace Ui {
class E2t1;
}

class E2t1 : public QWidget
{
    Q_OBJECT

public:
    explicit E2t1(QWidget *parent = nullptr);
    ~E2t1();

private:
    Ui::E2t1 *ui;
};

#endif // E2T1_H
