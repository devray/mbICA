# -------------------------------------------------
# Project created by QtCreator 2012-01-04T21:23:34
# -------------------------------------------------
QT -= core \
    gui
TARGET = mbica

# TEMPLATE = lib
# CONFIG += staticlib
SOURCES += mbica.cpp \
    icaseparator.cpp \
    nonlinearities.cpp \
    main.cpp \
    utils.cpp
HEADERS += mbica.h \
    icaseparator.h \
    nonlinearities.h \
    utils.h
LIBS += -larmadillo
