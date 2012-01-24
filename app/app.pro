# -------------------------------------------------
# Project created by QtCreator 2012-01-04T21:23:34
# -------------------------------------------------
QT -= core \
    gui

TARGET = mbica
DESTDIR = ../bin

TEMPLATE = app 
CONFIG += console

SOURCES += \
    main.cpp \

LIBS += -lmbica
linux:LIBS += -larmadillo
win32:LIBS += ../lib/blas_win32_MT.lib ../lib/lapack_win32_MT.lib

INCLUDEPATH += ../include
QMAKE_LIBDIR += ../lib

DEFINES += ARMA_USE_LAPACK \
    BOOST_PARAMETER_MAX_ARITY=7
