QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets printsupport

CONFIG += c++11

# You can make your code fail to compile if it uses deprecated APIs.
# In order to do so, uncomment the following line.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

SOURCES += \
    Generation.cpp \
    GenerationSignal.cpp \
    kol_iri_sput.cpp \
    main.cpp \
    mainwindow.cpp \
    paramSignal.cpp \
    qcustomplot.cpp

HEADERS += \
    Generation.h \
    GenerationSignal.h \
    Methods.h \
    Signal.h \
    kol_iri_sput.h \
    mainwindow.h \
    paramSignal.h \
    qcustomplot.h

FORMS += \
    kol_iri_sput.ui \
    mainwindow.ui \
    paramSignal.ui

# Default rules for deployment.
qnx: target.path = /tmp/$${TARGET}/bin
else: unix:!android: target.path = /opt/$${TARGET}/bin
!isEmpty(target.path): INSTALLS += target

RESOURCES += \
    Res.qrc


