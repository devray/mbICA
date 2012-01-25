# -------------------------------------------------
# Project created by QtCreator 2012-01-04T21:23:34
# -------------------------------------------------
QT -= core \
    gui

TARGET = mbica
DESTDIR = ../lib

TEMPLATE = lib
CONFIG += staticlib

SOURCES += \
    utils.cpp \
    policies.cpp \

HEADERS += ../include/mbica.h \
    ../include/icaseparator.h \
    ../include/nonlinearities.h \
    ../include/utils.h \
    ../include/policies.h

INCLUDEPATH += ../include
LIBS += -lboost_unit_test_framework
win32:LIBS += blas_win32_MT.lib lapack_win32_MT.lib
unix:LIBS += -larmadillo

DEFINES += ARMA_USE_LAPACK \
    BOOST_PARAMETER_MAX_ARITY=7
