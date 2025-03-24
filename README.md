***CAUTION***
Using this tool will result in the creation of hundreds of gigabytes of data. Ensure that the location you have selected has ample space. We recommend at least 1.5TB of free space.

Instructions:
1) 	Place this folder in it's entirety on the (external) location that you which to create the files in
2) 	Open the Julia code and ensure the current path is set to the location this folder is placed in
	pwd() can be used to check the current location 
	cd() can be used to change the directory  
	

	for example:
	If you want the data to be created at "D:/datafiles/benchmark_data"
	the downloaded folder needs to be placed at "D:/datafiles/benchmark_data/Milani_et_al_benchemark_data"
	to ensure proper functioning check the path in julia using "pwd()":

		julia> pwd()
		"C:\\Users\\analyst\\Documents\\Julia"

	This is not the location we need so change it using "cd()":

		julia> cd("D:\\datafiles\\benchmark_data")

	And recheck:

		julia> pwd()
		"D:\\datafiles\\benchmark_data"

	when correct, proceed. 

3)	Run the Julia file in it's entirety. 	
	Depending on the speed of your pc and your selected drive, this step may take several hours. We recommend (temporary) preventing the computer from going to sleep, and if possible to run it over night so it does not interfere with other activities. 
4)	After the script is finished, the files will be located in chosen path. In the example above this would be "D:/datafiles/benchmark_data". You may now remove the Milani_et_al_benchemark_data from this location if you please. 
