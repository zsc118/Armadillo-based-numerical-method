#include "matlabel.h"
#include "showwidget.h"
#include <QString>
#include <QMouseEvent>
namespace _MYFUNCTION {
void mat_append(QString& s,const double& x)
{
    QString t(QString::number(x,'g',6));
    s.append(t).append(QString(' ',7-t.length()));
}
}
matLabel::matLabel(QWidget*parent):QLabel(parent)
{
    setScaledContents(true);
    setAlignment(Qt::AlignCenter);
}
void matLabel::showMat(const arma::mat &x)
{
    if(x.empty())
    {
        setText("");
        return;
    }
    QString s;
    unsigned i(-1),j,m(x.n_rows-1),n(x.n_cols-1);
    while(++i<m)
    {
        j=-1;
        while(++j<n)
            _MYFUNCTION::mat_append(s,x.at(i,j));
        s.append(QString::number(x.at(i,n))).append('\n');
    }
    j=-1;
    while(++j<n)
        _MYFUNCTION::mat_append(s,x.at(m,j));
    setText(s.append(QString::number((a=x).at(m,n),'g',6)));
}
void matLabel::mousePressEvent(QMouseEvent* ev)
{
    if(ev->type()==QEvent::MouseButtonDblClick)
    {
        showWidget*D=new showWidget(this,a);
        D->setModal(true);
        D->setAttribute(Qt::WA_DeleteOnClose);
        D->exec();
    }
}
