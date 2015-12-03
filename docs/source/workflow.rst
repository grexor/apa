**********************************
Quick start
**********************************

A quick start example to show you how to use the platform to process your data.

Setup the annotation.tab file and experiment folders
----------------------------------------------------
Experiments are organized in libraries. For this example we will call our library example_lib, containing 3 experiments (with id 1, 2 and 3).
Upload one FASTQ file per experiment into the library folder structure

| data_folder / example_lib / e1 / example_lib_e1.fastq.gz
| data_folder / example_lib / e2 / example_lib_e2.fastq.gz
| data_folder / example_lib / e3 / example_lib_e3.fastq.gz

The data_folder path is defined in the :doc:`path module </apa.path>`.
