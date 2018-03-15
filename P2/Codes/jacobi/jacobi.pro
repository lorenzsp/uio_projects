TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
    define_matrices.cpp \
    jacobi.cpp \
    unit_test.cpp \
    unit_test_main.cpp

HEADERS += \
    catch.hpp \
    catch.hpp \
    jacobi.h \
    define_matrices.h
