#ifndef E2T2_H
#define E2T2_H

#include <QWidget>

namespace Ui {
class E2t2;
}

class E2t2 : public QWidget
{
    Q_OBJECT

public:
    explicit E2t2(QWidget *parent = nullptr);
    ~E2t2();

private:
    Ui::E2t2 *ui;
};

#endif // E2T2_H
