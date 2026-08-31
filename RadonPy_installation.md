

1.  Create conda env (如果遇到网络问题，把 -c conda-forge删掉)

```
conda create -n radonpy python=3.13 -c conda-forge
conda activate radonpy
# conda deactivate
```

2.  Installation of requirement packages

```

conda install -c conda-forge rdkit psi4 dftd3-python resp mdtraj psutil scipy pandas matplotlib pip
conda install -c conda-forge ambertools intermol
git clone https://github.com/RadonPy/RadonPy.git
# git clone https://github.com/AmonB/RadonPy.git
pip install --no-index --find-links=./RadonPy/dist/ radonpy-pypi
# force reinstall if you update your own fork
# pip install --no-index --find-links=./RadonPy/dist/ --force-reinstall radonpy-pypi
conda install -c conda-forge ambertools intermol
```

3.  Installation of LAMMPS by conda (有网络问题，而且好像不支持gpu加速)

```
conda install -c conda-forge lammps
```



3.1.  手动编译LAMMPS GPU、 PLUMED package，这里的GPU_ARCH是根据显卡架构来的，我的是2060对应sm_75，4090对应sm_89，PKG是安装额外可选的包**有一个隐形的坑，需要GPU MOLECULE  KSPACE  EXTRA-MOLECULE EXTRA-DUMP EXTRA-FIX package**，但安装指南没有给出。Ubuntu 22.04 內建的 OpenMPI 4.x 預設逐步淘汰了舊的 MPI C++ 綁定，因而導致在連結靜態庫 libplumed.a 時出現 undefined reference，暂时不 编译plumed。

```
wget https://download.lammps.org/tars/lammps-stable.tar.gz
tar -xzvf lammps*.tar.gz
cd lammps-22Jul2025
mkdir build
cd build
cmake ../cmake/ -D CMAKE_INSTALL_PREFIX=/home/amon/LAMMPS -D Kokkos_ARCH_HOSTARCH=yes -D Kokkos_ARCH_GPUARCH=yes -D Kokkos_ENABLE_CUDA=yes -D Kokkos_ENABLE_OPENMP=yes -D GPU_API=cuda -D GPU_ARCH=sm_75 -D PKG_GPU=yes -D PKG_KSPACSE=yes -D PKG_RIGID=yes -D PKG_TALLY=yes -D PKG_MOLECULE=yes -D PKG_EXTRA-MOLECULE=yes -D PKG_EXTRA-DUMP=yes -D PKG_KSPACE=yes -D PKG_EXTRA-FIX=yes
make -j 8
make install
echo 'export PATH=$PATH:/home/amon/LAMMPS/bin/' >> ~/.bashrc
source ~/.bashrc
lmp -h
```



编译plumed的版本

```
wget https://download.lammps.org/tars/lammps-stable.tar.gz
tar -xzvf lammps*.tar.gz
cd lammps-22Jul2025
mkdir build
cd build
cmake ../cmake/ -D CMAKE_INSTALL_PREFIX=/home/amon/LAMMPS -D Kokkos_ARCH_HOSTARCH=yes -D Kokkos_ARCH_GPUARCH=yes -D Kokkos_ENABLE_CUDA=yes -D Kokkos_ENABLE_OPENMP=yes -D GPU_API=cuda -D GPU_ARCH=sm_75 -D PKG_GPU=yes -D PKG_KSPACSE=yes -D PKG_RIGID=yes -D PKG_TALLY=yes -D PKG_MOLECULE=yes -D PKG_EXTRA-MOLECULE=yes -D PKG_EXTRA-DUMP=yes -D PKG_KSPACE=yes -D PKG_EXTRA-FIX=yes -D PKG_PLUMED=yes -D DOWNLOAD_PLUMED=yes
make -j 8
make install
echo 'export PATH=$PATH:/home/amon/LAMMPS/bin/' >> ~/.bashrc
source ~/.bashrc
lmp -h
```



3.2 another pre-installation

```
sudo apt install ffmpeg
sudo apt-get install libopenblas-dev
sudo apt-get install libgsl-dev
```

build Plumed (maybe don't need it)

```
wget https://github.com/plumed/plumed2/releases/download/v2.10.0/plumed-src-2.10.0.tgz
tar -xzvf plumed*.tar.gz
 ./configure --prefix=/home/amon/plumed
 make -j 8
 make install
 echo '# PLUMED' >> ~/.bashrc
 echo 'export PATH=$PATH:/home/amon/plumed/bin' >> ~/.bashrc
 echo 'export PATH=$PATH:/home/amon/plumed/include' >> ~/.bashrc
 echo 'export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/amon/plumed/lib' >> ~/.bashrc
 source ~/.bashrc
```

4.  手动编译LAMMPS，需要在当前的conda env 导出LAMMPS路径，实际安装完之后，不一定需要

```
export LAMMPS_EXEC=/home/amon/LAMMPS/bin/
```

   
