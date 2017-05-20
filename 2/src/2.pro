TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
    vmc.cpp \
    slatervmc.cpp \
    ho2d.cpp \
    main.cpp \
    test.cpp

HEADERS += \
    vmc.h \
    slatervmc.h \
    ho2d.h

LIBS += -llapack -lblas -larmadillo

QMAKE_CXXFLAGS+= -fopenmp
QMAKE_LFLAGS +=  -fopenmp
