#ifndef E2T6_H
#define E2T6_H

#include <QWidget>

namespace Ui {
class E2t6;
}

class E2t6 : public QWidget
{
    Q_OBJECT

public:
    explicit E2t6(QWidget *parent = nullptr);
    ~E2t6();

private:
    Ui::E2t6 *ui;
};

#endif // E2T6_H
