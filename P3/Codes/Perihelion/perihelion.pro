TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
    main.cpp \
    planet.cpp \
    solarsystem.cpp \
    verlet.cpp

HEADERS += \
    planet.h \
    solarsystem.h \
    verlet.h
