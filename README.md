### EIC Neutral Current Deep Inelastic Scattering (NC-DIS) using PYTHIA-8

Very simple setup to generate NC-DIS events from [PYTHIA 8](https://pythia.org/latest-manual/Welcome.html) with somewhat realistic beam effects.  Includes jet clustering and matching with deep-inelastic scattered parton. Beam crossing angle effects are impemented using Pythia8::BeamShape derived class eicBeamShape from [eicSimuBeamEffects/Pythia8](https://github.com/eic/eicSimuBeamEffects/tree/master/Pythia8) repository. 

## Installation and Compilation

Assumes PYTHIA 8, [FastJet](https://www.fastjet.fr/) and [ROOT](https://root.cern.ch/) already installed. If not please install them before proceeding. I am providing one way to install them here tailored to work with the code here. But feel free to install them anyway you like, but take care to set the paths in `Makefile` properly. 

Before you start make a folder
```
mkdir -p $HOME/heplibs
cd $HOME/heplibs
```
All the subsequent installations will have their outputs localized in this folder.
# PYTHIA8:
```
wget https://pythia.org/download/pythia83/pythia8312.tgz
tar xvfz pythia8312.tgz
cd pythia8312
make -j4
```
The xmldoc folder to initialize new `Pythia8::Pythia` objects will be in `$HOME/heplibs/pythia8310/share/Pythia8/xmldoc`

# FastJet:
```
wget https://fastjet.fr/repo/fastjet-3.4.3.tar.gz  
tar zxvf fastjet-3.4.3.tar.gz  
cd fastjet-3.4.3/
./configure --prefix=$PWD/../fastjet-install
make  
make check  
make install  
cd ..
```

# ROOT:
Follow steps laid out here: [Dependencies - ROOT](https://root.cern/install/dependencies/) to install all dependencies. Then select a release compiled for your OS and system here: [Releases - ROOT](https://root.cern/install/all_releases/). I want the latest release as of writing (6.34) this on my Ubuntu 22.04 computer so I would do:
```
wget https://root.cern/download/root_v6.34.00-rc1.Linux-ubuntu22.04-x86_64-gcc11.4.tar.gz
tar -xzvf root_v6.34.00-rc1.Linux-ubuntu22.04-x86_64-gcc11.4.tar.gz
```
This should add a folder `root` 

If you follow these instructions correctly then you should be able to compile the program without any issues. Just do:
```
mkdir -p $HOME/eicPythia8NCDIS
cd $HOME/eicPythia8NCDIS
git clone https://github.com/TanmayPani/EIC-NCDIS-setup.git
cd EIC-NCDIS-setup
make
./eicJetMaker.exe
```