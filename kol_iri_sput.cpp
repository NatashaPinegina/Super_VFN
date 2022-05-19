#include "kol_iri_sput.h"
#include "ui_kol_iri_sput.h"
#include "GenerationSignal.h"

GenerationSignal kol_IRI;

Kol_IRI_Sput::Kol_IRI_Sput(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::Kol_IRI_Sput)
{
    ui->setupUi(this);

    QList<QString> Name;
    QList<QString> Znach;

    Name.push_back(QString::fromStdString("Кол-во ИРИ"));
    Name.push_back(QString::fromStdString("Кол-во спутников"));


    Znach.push_back(QString::number(kol_IRI.NumberOfSources));
    Znach.push_back(QString::number(kol_IRI.NumberOfSatellites));



        // Create model:
        Table_KOL *Param = new Table_KOL(this);

        // Populate model with data:
        Param->populateData(Name,Znach);

        // Connect model to table view:
        ui->tableView->setModel(Param);

        // Make table header visible and display table:
        ui->tableView->horizontalHeader()->setVisible(true);
        ui->tableView->horizontalHeader()->setSectionResizeMode(0, QHeaderView::Stretch);
        ui->tableView->show();
}

Kol_IRI_Sput::~Kol_IRI_Sput()
{
    delete ui;
}

Table_KOL::Table_KOL(QObject *parent) : QAbstractTableModel(parent)
{
}

// Create a method to populate the model with data:
void Table_KOL::populateData(const QList<QString> &Name,const QList<QString> &Znach)
{
    tm_name.clear();
    tm_name = Name;
    tm_znach.clear();
    tm_znach = Znach;
    return;
}

int Table_KOL::rowCount(const QModelIndex &parent) const
{
    Q_UNUSED(parent);
    return tm_name.length();
}

int Table_KOL::columnCount(const QModelIndex &parent) const
{
    Q_UNUSED(parent);
    return 2;
}

QVariant Table_KOL::data(const QModelIndex &index, int role) const
{
    if (!index.isValid() || role != Qt::DisplayRole) {
        return QVariant();
    }
    if (index.column() == 0)
    {
        return tm_name[index.row()];
    } else if (index.column() == 1)
    {
        return tm_znach[index.row()];
    }
    return QVariant();
}

QVariant Table_KOL::headerData(int section, Qt::Orientation orientation, int role) const
{
    if (role == Qt::DisplayRole && orientation == Qt::Horizontal)
    {
        if (section == 0) return QString("Название");
        else if (section == 1) return QString("Значение");
    }
    return QVariant();
}

