TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    harmonicoscillator2d.cpp \
    hermite.cpp \
    gausshermitequadrature.cpp

HEADERS += \
    harmonicoscillator2d.h \
    hermite.h

QMAKE_CXXFLAGS+= -fopenmp
QMAKE_LFLAGS +=  -fopenmp
