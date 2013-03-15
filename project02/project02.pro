TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    src/mainapplication.cpp \
    src/CAtom.cpp \
    src/CState.cpp \
    src/lib.cpp \
    src/CBox.cpp \
    src/CStatisticsSampler.cpp

HEADERS += \
    src/mainapplication.h \
    src/CAtom.h \
    src/CState.h \
    src/lib.h \
    src/CBox.h \
    src/inlines.h \
    src/linkedList.h \
    src/stack.h \
    src/linkedList2.h \
    src/CStatisticsSampler.h

LIBS += -larmadillo

release {
    DEFINES += ARMA_NO_DEBUG
    QMAKE_LFLAGS -= -O1
    QMAKE_LFLAGS += -O3
    QMAKE_LFLAGS_RELEASE -= -O1
    QMAKE_LFLAGS_RELEASE += -O3
    QMAKE_CXXFLAGS -= -O2
    QMAKE_CXXFLAGS += -O3
    QMAKE_CXXFLAGS_RELEASE -= -O2
    QMAKE_CXXFLAGS_RELEASE += -O3
}


QMAKE_LFLAGS -= -O1
QMAKE_CXXFLAGS -= -O2
QMAKE_LFLAGS_RELEASE -= -O1
QMAKE_CXXFLAGS_RELEASE -= -O2
