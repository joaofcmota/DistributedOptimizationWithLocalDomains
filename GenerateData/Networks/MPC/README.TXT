
This folder contains files for generating networks for the MPC problem in generic networks.

RawNetwork:
    Contains the US power grid network with about 5000 nodes. 


Scripts:
    Contains scripts for either creating networks (in case of random networks) or processing networks
    (in case there are networks in the RawNetwork folder)

        READ_FROM_RAW.m:
            Main file; it is a Matlab script that executes the sage scripts in the Scripts folder
            and creates a Networks struct in Matlab. These structs and the network statistics are 
            stored in the ProcessedNetwork folder

        processPowerGrid.sage:
            Reads the .gml file in the RawNetwork folder for the power grid and stores its information
            in a tmp.txt file (to be processed and deleted by READ_FROM_RAW.m)

        Barabasi.sage:
            Creates a Barabasi network with a given parameter and number of nodes


ProcessedNetwork:
    Contains .mat files with Networks structs and respective statistics.

            
