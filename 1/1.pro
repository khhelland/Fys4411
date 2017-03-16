TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    harmonicoscillator2d.cpp \
    hermite.cpp \
    integration.cpp \
    integrator.cpp \
    coulomb_functions.cpp \
    hartreefock.cpp

HEADERS += \
    harmonicoscillator2d.h \
    hermite.h \
    integrator.h \
    coulomb_functions.h \
    hartreefock.h

LIBS += -llapack -lblas -larmadillo

QMAKE_CXXFLAGS+= -fopenmp
QMAKE_LFLAGS +=  -fopenmp
