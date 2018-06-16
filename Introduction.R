library("Plum")

Folder="Directory where data is located"
Filename="Name data file in CVS format"
Data=Data_sim() #This creates a simulated core in home directory
Data
#Depth(bottom)  Density     210Pb     sd(210Pb)   Sample thickness    226Ra       sd(226Ra)
# 1               0.1000913  92.16821 10.000000    1                 12.370891    3
# 2               0.1006381 154.00736  9.689655    1                  18.359602    3
# 3               0.1017258 190.35590  9.379310    1                  17.964667    3
# ...              ...        ...       ...         ...                 ...         ...
#first row should have the described labels
#If 226Ra was no measured leave those column blank


runPlum() #Will create a simulation in ~/Plum/Sumulation-(a random number) and it will run plum in this simulation.





runPlum(folder = Folder,#  "Directory where data is located",
        Data =Filename , # "Name data file in CVS format",
        iterations = "number of iterations ",
        by =  "Chronology section  length", 
        number_supported = "Use when alpha spectrometry is use and you know how many bottom samples should be use to infer supported activity" ,
        detection_limit = "Deafult .05",
        memory_shape = "prior shape of memory parameter",
        memory_mean = "prior mean of memory parameter",
        acc_shape = "prior shape of accumulation parameter",
        acc_mean = "prior mean of accumulation parameter",
        fi_mean = "prior mean of supply",
        fi_acc = "prior shape of supply",
        As_mean = "prior mean of supported 210Pb",
        As_acc = "prior shape of supported 210Pb",
        seeds = "Simulations seed "
        )

#Most of the rumPlum variables have default
#folder and Data are the only necesary information 
#other variables should be change according to prior information
#Data is the name of the file as "Corecode.csv"




