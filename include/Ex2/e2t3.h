#ifndef E2T3_H
#define E2T3_H

#include <QWidget>

namespace Ui {
class E2t3;
}

class E2t3 : public QWidget
{
    Q_OBJECT

public:
    explicit E2t3(QWidget *parent);
    ~E2t3();

private:
    Ui::E2t3 *ui;
    unsigned n;
};

#endif // E2T3_H
