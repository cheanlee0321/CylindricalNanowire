TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    cylindricalcoordinate.cpp \
    cylindricalquantumdd.cpp

QMAKE_CXXFLAGS += -fopenmp -lgsl -lgslcblas
LIBS += -fopenmp -lgsl -lgslcblas

HEADERS += \
    Parameter.h \
    cylindricalcoordinate.h \
    cylindricalquantumdd.h
