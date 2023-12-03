#ifndef E2T5_H
#define E2T5_H

#include <QWidget>

namespace Ui {
class E2t5;
}

class E2t5 : public QWidget
{
    Q_OBJECT

public:
    explicit E2t5(QWidget *parent = nullptr);
    ~E2t5();

private:
    Ui::E2t5 *ui;
};
namespace Ex2 {
class Interval
{
    double a,b;
public:
    constexpr Interval()noexcept:a(0.),b(1.){}
    Interval(double t)noexcept:a(t),b(t+1){}
    Interval(double a0,double b0)noexcept:a(a0),b(b0){}
    double& get_a()noexcept
    {
        return a;
    }
    double& get_b()noexcept
    {
        return b;
    }
};
}
#endif // E2T5_H
