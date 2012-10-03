This Python script calculates the Lyman-alpha transmission function along a single line of sight. 
Required packages are Python, numpy and matplotlib

To use it, you must first create a data file containing the density, ionization and temperature data for your sightline. Look at the supplied example_data.dat to see how this is done.

You must also make a parameters file. Look at the supplied parameters.txt to see the format of this file.

To calculate the transmission, run:
'./Simpletransfer1D.py -p myparameters.txt -d mydata.dat'

The transmission will be saved to the file specified in your parameters file.

To plot the results, run:
'./plot_transmission.py -p myparameters.txt -t mytransmission.txt'

you can also run
'./plot_transmission.py -p myparameters.txt -t mytransmission.txt -d mydata.dat'
to get a plot of the density and temperature as well.

To view help, run:
'./Simpletransfer1D.py -h'
or 
'./plot_transmission.py -h'
