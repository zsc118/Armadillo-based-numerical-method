#include "veclabel.h"
#include <QString>
using namespace arma;
vecLabel::vecLabel(QWidget* parent):QLabel(parent)
{
    setScaledContents(true);
    setAlignment(Qt::AlignCenter);
}

void vecLabel::showVec(const vec &x,const char* first)
{
    if(x.empty())
    {
        QString s(first);
        setText(s.append("："));
        return;
    }
    QString s(first);
    s.append("：(");
    unsigned n(x.n_elem-1);
    for(unsigned i(0);i<n;++i)
        s.append(QString::number(x.at(i))).append(',');
    setText(s.append(QString::number(x.at(n))).append(")^T"));
}
