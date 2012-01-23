# -------------------------------------------------
# Project created by QtCreator 2012-01-04T21:23:34
# -------------------------------------------------
QT -= core \
    gui
TARGET = mbica

# TEMPLATE = lib
# CONFIG += staticlib
SOURCES += \
    icaseparator.cpp \
    nonlinearities.cpp \
    main.cpp \
    utils.cpp \
    policies.cpp
HEADERS += mbica.h \
    icaseparator.h \
    nonlinearities.h \
    utils.h \
    policies.h
LIBS += -larmadillo

DEFINES += ARMA_USE_LAPACK \
    BOOST_PARAMETER_MAX_ARITY=7
