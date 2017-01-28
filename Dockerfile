FROM ubuntu:latest
MAINTAINER Yi Gao

RUN apt-get update
RUN apt-get install -y \
  build-essential \
  curl \
  git \
  cmake \
  libgsl-dev \
  libx11-dev \
  libxt-dev \
  libgl1-mesa-dev \
  libglu1-mesa-dev \
  libxrender-dev



# Instal VTK
WORKDIR /tmp/
RUN git clone http://vtk.org/VTK.git
WORKDIR /tmp/VTK/
RUN git pull
WORKDIR /tmp/VTK/build/
RUN cmake -DCMAKE_BUILD_TYPE:STRING=Release -DBUILD_EXAMPLES:BOOL=OFF -DBUILD_TESTING:BOOL=OFF ..
RUN make -j$(nproc)


# Instal ITK
WORKDIR /tmp/
RUN git clone http://itk.org/ITK.git
WORKDIR /tmp/ITK/
RUN git pull
WORKDIR /tmp/ITK/build/
RUN cmake -DCMAKE_BUILD_TYPE:STRING=Release -DVTK_DIR:STRING=/tmp/VTK/build -DBUILD_EXAMPLES:BOOL=OFF -DBUILD_TESTING:BOOL=OFF -DBUILD_TESTS:BOOL=OFF -DModule_ITKReview:BOOL=ON -DModule_ITKVtkGlue:BOOL=ON ..
RUN make -j$(nproc)


# Download alglib
WORKDIR /tmp/
RUN curl -O http://www.alglib.net/translator/re/alglib-3.10.0.cpp.gpl.tgz
RUN tar xf alglib-3.10.0.cpp.gpl.tgz
RUN mv cpp alglib


# compile SPoM
WORKDIR /tmp
RUN git clone https://github.com/gaoyi/SPoM.git
WORKDIR /tmp/SPoM/
RUN git pull
WORKDIR /tmp/SPoM/build
RUN cmake -DCMAKE_BUILD_TYPE:STRING=Release -DVTK_DIR:STRING=/tmp/VTK/build -DITK_DIR:STRING=/tmp/ITK/build ../src
RUN make -j$(nproc)




WORKDIR /tmp/SPoM/build/sumOfTwoPoisson


CMD ["/bin/bash"]
