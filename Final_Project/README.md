# Determining the Hubble Constant!

This python library contains the necessary data and analysis tools to calculate the Hubble Constant.

## Usage

data_extraction.py contains tools to read in the 3 provided datasets and extract from them distance and recession velocity measurements.

curve_fitting.py contains code to perform either a least-squares fit or a Markov Chain Monte Carlo (Metropolis-Hastings Algorithm) to solve for the value of the Hubble Constant.

Determining_Hubble_Const_Final_Project_Drachler.ipynb is a notebook that walks you through the analysis. This analysis can be easily be extended by adding more data or using a different fitting technique. 

## Contributing
If you have some data you'd like to fit, adding that data to the pipeline can be done by editing data_extraction.py to extract the distance and recessional velocity data. If you look at the code, you'll notice a template that is used by my 3 chosen datasets. Please follow this template.

If you'd like to add a new fitting technique, this can easily be added into curve_fitting.py. Your implementation of this is a little more open-ended. 

Please let me know if you find bugs or have suggestions! 

Author - Brendan Drachler, PhD Candidate, Rochester Institute of Technology - bcd3735@rit.edu