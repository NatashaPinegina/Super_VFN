#ifndef KOL_IRI_SPUT_H
#define KOL_IRI_SPUT_H

#include <QDialog>
#include "mainwindow.h"

namespace Ui {
class Kol_IRI_Sput;
}

class Table_KOL: public QAbstractTableModel
{
    Q_OBJECT

public:
    Table_KOL(QObject *parent = 0);

    void populateData(const QList<QString> &Name,const QList<QString> &Znach);

    int rowCount(const QModelIndex &parent = QModelIndex()) const Q_DECL_OVERRIDE;
    int columnCount(const QModelIndex &parent = QModelIndex()) const Q_DECL_OVERRIDE;

    QVariant data(const QModelIndex &index, int role = Qt::DisplayRole) const Q_DECL_OVERRIDE;
    QVariant headerData(int section, Qt::Orientation orientation, int role = Qt::DisplayRole) const Q_DECL_OVERRIDE;

private:
    QList<QString> tm_name;
    QList<QString> tm_znach;
};

class Kol_IRI_Sput : public QDialog
{
    Q_OBJECT

public:
    explicit Kol_IRI_Sput(QWidget *parent = nullptr);
    ~Kol_IRI_Sput();

private:
    Ui::Kol_IRI_Sput *ui;
};

#endif // KOL_IRI_SPUT_H
