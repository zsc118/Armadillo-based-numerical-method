#ifndef E2T7_H
#define E2T7_H

#include <QWidget>

namespace Ui {
class E2t7;
}

class E2t7 : public QWidget
{
    Q_OBJECT

public:
    explicit E2t7(QWidget *parent = nullptr);
    ~E2t7();

private:
    Ui::E2t7 *ui;
};

#endif // E2T7_H
