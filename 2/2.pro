TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    metropolis.cpp \
    montecarlointegrator.cpp

HEADERS += \
    montecarlointegrator.h

