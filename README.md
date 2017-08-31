# Falcon
A very fast hadronization and detector simulator based on ideas pioneered by Bruce Knuteson, in particular, the use of a lookup table to map events at the parton shower level to events at the reconstruction or analysis levels.

---

## [](#header-2)Acknowledgement

The codes for the Falcon package were originally developed by:
* Sergei Gleyzer, University of Florida and CERN;
* Harrison Prosper, Florida State University; and 
* Omar Zapata, UdeA and ITM.

The official Falcon package can be found [here](http://oproject.org/falcon).  The simulator package is mentioned in an academic publication: [Les Houches 2015: Physics at TeV colliders - new physics working group report](http://inspirehep.net/record/1456803#).

## [](#header-2)Background

The current version of Falcon started out as a simulator program that implements a lookup table to record exact matches from parton-level showers to particle-level jets based on geometric proximity, essentially doing what is idiomatically called jet reconstruction.  However, there are two favorable elements missing:
* Compatibility with TMVA
* Learning jet reconstruction rules

## [](#header-2)Project Objectives

In order to improve Falcon on these two areas, this Google Summer of Code 2017 project aims to deploy existing TMVA tools to empower Falcon with the ability to learn the mapping from parton-level showers to particle-level jets.  With this learning capability, TMVA engineers can then scale up the Falcon package to take care of more jet reconstruction tasks with different detector configurations.
 
## [](#header-2)Implementations

Among a basket of plausible methods, implementing a neural network with the MLP (Multi-Layer Perceptron) method of TMVA seems to stand out with its flexibility and fast learning ability.  Users can customize their own settings to train their neural networks for optimal jet reconstruction.  The program exports an Excel (xml) file that stores the learnt weights from the neural network and offers an option for users to launch the TMVA GUI to display regression metrics with a complete set of generated histograms and scatter plots.

## [](#header-2)Future Plans

While the core of the project has been finished, near the end of the project I tried something a bit more adventurous: I tried to take one step further to integrate the histogram customization feature in [TMVARegressionApplication.C](https://root.cern.ch/doc/v610/TMVARegressionApplication_8C.html) into Falcon.  However, since a critical TMVA class (TMVA::Reader) has not yet supported certain data types used by the data files for jet reconstruction, useful codes are stored in a branch separate from the master branch, namely TMVARegressionApplication, for future completion.  The master franch is kept free from the above attempt.

Right now, Falcon's implementation is still tailored towards a specific data input format.  Going forward, engineers may find making Falcon serve a more general purpose useful.


## [](#header-2)Building

1. Clone the repository of ROOT and build ROOT in the directory `~/build/` following the instructions [here](https://github.com/root-project/root).

2. Set up the environment variables for ROOT

$ cd build
$ source bin/thisroot.sh

3. Clone the repositories of Falcon and its data.

$ cd ..
$ mkdir falcon
$ cd falcon
$ git clone falcon@oproject.org:falcon
$ git clone falcon@oproject.org:falcondata data

4. Set up the environment variables for Falcon

$ source setup.sh

5. Make the Delphes package inside Falcon

$ cd falcon/delphes
$ ./configure
$ make

6. Make the Falcon package

$ cd ..
$ make

## [](#header-2)Running

1. Set up the environment variables for ROOT

$ cd build
$ source bin/thisroot.sh

4. Set up the environment variables for Falcon

$ cd ../falcon
$ source setup.sh

5. Build the lookup table from the data first

$ cd delphes
$ python ../test/testFalcon.py

6. Learn the mapping

$ python ../test/testFalcon.py learn [cut-options] [prepare-options] [train_options]

7. Run ROOT to launch TMVA GUI

$ root
$ Type in the ROOT interface: TMVA::TMVARegGui("FalconTMVAReg.root")

8. Click buttons on the launched TMVA GUI for histograms

## [](#header-2)Questions

Should you have any questions, please do not hesitate to contact Brian Lui (the author of this page) at brianlui@stanford.edu
