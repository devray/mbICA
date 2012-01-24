# -------------------------------------------------
# Project created by QtCreator 2012-01-04T21:23:34
# -------------------------------------------------
QT -= core \
    gui

TARGET = unit_tests

TEMPLATE = app
CONFIG += console

SOURCES += \
    pca_test.cpp
HEADERS += \
    whitening_test.h

LIBS += -lmbica -larmadillo -lboost_unit_test_framework

INCLUDEPATH += ../include
QMAKE_LIBDIR += ../lib

DEFINES += ARMA_USE_LAPACK \
    BOOST_PARAMETER_MAX_ARITY=7
