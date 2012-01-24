# -------------------------------------------------
# Project created by QtCreator 2012-01-04T21:23:34
# -------------------------------------------------
QT -= core \
    gui

TARGET = ../lib/mbica

TEMPLATE = lib
CONFIG += staticlib

SOURCES += \
    icaseparator.cpp \
    nonlinearities.cpp \
    utils.cpp \
    policies.cpp
HEADERS += ../include/mbica.h \
    ../include/icaseparator.h \
    ../include/nonlinearities.h \
    ../include/utils.h \
    ../include/policies.h

INCLUDEPATH += ../include
LIBS += -larmadillo -lboost_unit_test_framework

DEFINES += ARMA_USE_LAPACK \
    BOOST_PARAMETER_MAX_ARITY=7
