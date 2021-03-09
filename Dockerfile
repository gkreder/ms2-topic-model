FROM bioconductor/bioconductor_docker
MAINTAINER Gabriel Reder gkreder@gmail.com

RUN R -e "install.packages('rcdk')"
RUN R -e "install.packages('docopt')"
RUN R -e "install.packages('stringr')"

ENV PATH /root/miniconda3/bin:$PATH
ENV PATH /usr/local/bin/R:$PATH
ENV LD_LIBRARY_PATH /usr/local/lib/R/lib:$LD_LIBRARY_PATH
ENV PKG_CONFIG_PATH /usr/local/lib/R/lib/pkgconfig/:$PKG_CONFIG_PATH


RUN apt-get update --fix-missing && \
    apt-get install -y wget bzip2 ca-certificates curl git

RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /root/miniconda3 && \
    rm ~/miniconda.sh && \
    /root/miniconda3/bin/conda clean -tipsy && \
    ln -s /root/miniconda3/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /root/miniconda3/etc/profile.d/conda.sh" >> ~/.bashrc && \
    echo "conda activate base" >> ~/.bashrc


RUN pip install tomotopy
RUN conda install -c conda-forge rdkit
RUN pip install scikit-learn scipy jupyterlab ipywidgets matplotlib seaborn lxml pyteomics xlrd tqdm molmass
RUN jupyter nbextension enable --py widgetsnbextension
RUN conda install -c conda-forge nodejs
RUN jupyter labextension install @jupyter-widgets/jupyterlab-manager
RUN pip install rpy2

CMD bash


