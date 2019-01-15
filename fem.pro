QT -= gui

CONFIG += c++17 console
CONFIG -= app_bundle

QMAKE_CXXFLAGS += -std=c++1z
# The following define makes your compiler emit warnings if you use
# any Qt feature that has been marked deprecated (the exact warnings
# depend on your compiler). Please consult the documentation of the
# deprecated API in order to know how to port your code away from it.
DEFINES += QT_DEPRECATED_WARNINGS

# You can also make your code fail to compile if it uses deprecated APIs.
# In order to do so, uncomment the following line.
# You can also select to disable deprecated APIs only up to a certain version of Qt.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

SOURCES += \
        main.cpp \
    mesh_generation.cpp \
    mesh.cpp \
    ndimmatrix/matrix.cpp

# Default rules for deployment.
qnx: target.path = /tmp/$${TARGET}/bin
else: unix:!android: target.path = /opt/$${TARGET}/bin
!isEmpty(target.path): INSTALLS += target

HEADERS += \
    node.h \
    mesh_generation.h \
    mesh.h \
    ndimmatrix/matrix.h \
    ndimmatrix/matrix_impl.h \
    element.h

unix {
    target.path = /usr/lib
    #INSTALLS += target

    INCLUDEPATH += /opt/intel/parallel_studio_xe_2019.1.053/compilers_and_libraries_2019/linux/mkl/include/
    INCLUDEPATH += /usr/include/x86_64-linux-gnu/c++/8/
    LIBS += -L/opt/intel/parallel_studio_xe_2019.1.053/compilers_and_libraries_2019/linux/mkl/lib/intel64/ \
    -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core \
    -L/opt/intel/parallel_studio_xe_2019.1.053/compilers_and_libraries_2019/linux/compiler/lib/intel64/  \
    -liomp5 -lpthread -lm #-dl
}
