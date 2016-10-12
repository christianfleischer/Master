TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    system.cpp

LIBS += -llapack -lblas -larmadillo

HEADERS += \
    catch.h \
    system.h

