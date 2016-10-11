TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    initialize_electrons.cpp \
    potential.cpp

LIBS += -llapack -lblas -larmadillo

HEADERS += \
    catch.h

